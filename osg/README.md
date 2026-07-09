# OSG/OSPool pipeline for the Gtau tilt sensitivity study

Distributes `main_simu_gtau_tilt_sensitivity.R` over HTCondor: one job = one
shard (a contiguous block of replication indices within one `(setup, n)` cell).
The per-replication logic is `src/gtau_tilt_core.R`, **shared with the local
serial driver**, and shard configs carry per-rep seed slices pre-drawn with the
driver's canonical scheme (`set.seed(123)` at `total_R`) — so a distributed run
computes exactly what the serial driver would, replication by replication.

Adapted from the mixC-era pipeline (kept locally as `osg_old/`, untracked).
Runs inside the Singularity image `osdf:///ospool/ap41/data/yuyao.wang/dCP1.sif`
(needs only R + survival + xgboost). Workers pin `OMP_NUM_THREADS=1` (run.sh)
to match `request_cpus = 1`.

## Files

| file | role |
|---|---|
| `prepare_gtau_jobs.R` | enumerate cells, draw + slice seeds, write shard configs, `manifest.tsv`, `queue.txt`, `job_ids.txt`, `plan.rds`; pre-create result/log dirs |
| `prepare_probe.sh` / `prepare_full_grid.sh` | wrappers with the agreed probe / production settings (env-var overridable) |
| `gtau_job.R` | worker: one shard via `run_one_gtau_rep()`; validates seeds and registry values; records R/package versions |
| `run.sh` | Condor executable: locate config, pin OMP, invoke `gtau_job.R` |
| `submit.sub` | one Condor job per `queue.txt` line, container-only R, `max_retries = 2` |
| `collect_gtau_results.R` | merge shards; for complete cells rebuild the CANONICAL result rds (identical schema to the local driver's) |

## Workflow (on the OSG access point, from the repo root)

1. Get the latest code:

```
git pull
```

2. Probe — replication 1 of every cell (6 one-rep jobs; seeds identical to the
   production study, so outputs double as a cross-machine consistency check):

```
bash osg/prepare_probe.sh
cd osg
condor_submit experiment_name=gtau_probe submit.sub
condor_watch_q
```

   When done, look at runtime and memory per cell (`condor_history -limit 6
   -af RequestMemory MemoryUsage RemoteWallClockTime` or the `.log` files under
   `logs/osg/gtau_probe/condor/`), then set `MEMORY_MAP` accordingly.

3. Production — full grid (linWB1+linWB2 × n ∈ {300,500,1000}, R=200,
   n_test=3000; 76 jobs with the default shard map):

```
bash osg/prepare_full_grid.sh
cd osg
condor_submit experiment_name=gtau_grid submit.sub
```

4. Collect (tolerates missing shards; reports them):

```
Rscript osg/collect_gtau_results.R --experiment-name gtau_grid
```

   Complete cells produce canonical files under
   `results/osg/gtau_grid/collected/<setup>/gtau_tilt/` — the exact structure
   the plot script reads. The collector NEVER writes into the top-level
   `results/<setup>/` tree; promoting a file there is a deliberate manual `cp`.

5. Resubmit failed jobs only (ids from the collector's report):

```
grep -E '^(0007|0012) ' jobs/gtau_grid/queue.txt > jobs/gtau_grid/queue_retry.txt
condor_submit experiment_name=gtau_grid queue_file=queue_retry.txt submit.sub
```

6. Bring results home (from the laptop):

```
rsync -av <ap-user>@<access-point>:<repo-path>/results/osg/gtau_grid/collected/ results/osg/gtau_grid/collected/
```

7. Plot without moving files (optional 6th argument = results root):

```
Rscript main_plot_gtau_tilt_sensitivity.R 200 1000 3000 0.1 linWB2 results/osg/gtau_grid/collected
```

## Reproducibility notes

- Same seeds as the local driver ⇒ same replications. Bitwise identity across
  machines additionally requires matching R + xgboost + survival versions; the
  worker records its versions in `meta` (`pkg_versions`, `r_version`) and the
  probe comparison against a local run measures the actual discrepancy.
- `--total-r` is the canonical study size and must stay 200 even for partial
  materializations (`--rep-subset`); the collector cross-checks every shard's
  seed slices against the canonical vectors and refuses on mismatch.
- New cells (other n, other setups) = new `--experiment-name`; never reuse a
  jobs dir (the prepare script refuses to overwrite).
