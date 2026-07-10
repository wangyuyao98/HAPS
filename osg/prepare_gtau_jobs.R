#!/usr/bin/env Rscript

## =====================================================================
## Prepare OSG/HTCondor jobs for the Gtau tilt sensitivity study.
##
## Enumerates (setup x n) cells, pre-draws the three per-replication seed
## vectors for each cell with the LOCAL DRIVER'S canonical scheme —
##     set.seed(123); seeds <- sample(1e7, total_R); seeds_test <- ...;
##     seeds_split <- ...
## (fixed by design; no master-seed flag, because changing it would break
## equivalence with main_simu_gtau_tilt_sensitivity.R) — slices them into
## contiguous shards, and writes one config .rds per shard plus
## manifest.tsv / queue.txt / job_ids.txt / plan.rds.
##
## IMPORTANT: seed vectors depend on total_R (the three sample() calls share
## one RNG stream), so --total-r must be the CANONICAL study size (200) even
## when only a subset of replications is materialized via --rep-subset
## (probe runs, targeted re-runs of failed reps).
##
## Usage (from anywhere; paths resolve relative to this script):
##   Rscript osg/prepare_gtau_jobs.R --experiment-name gtau_grid \
##     --setups linWB1,linWB2 --sample-sizes 300,500,1000 \
##     --n-test 3000 --total-r 200 --shard-map 300:25,500:20,1000:10 \
##     --memory-map 300:2GB,500:2GB,1000:3GB --disk-map all:2GB
##   Rscript osg/prepare_gtau_jobs.R --experiment-name gtau_probe \
##     --total-r 200 --rep-subset 1 --memory-map all:4GB
## =====================================================================

## ------------------------- Defaults + CLI -----------------------------
opts <- list(
    experiment_name = "gtau_grid",
    setups          = "linWB1,linWB2",
    sample_sizes    = "300,500,1000",
    n_test          = 3000L,
    total_r         = 200L,
    alpha           = 0.1,
    models          = "cox",       # a case from gtau_model_cases(); non-default
                                   # cases are tagged into result paths/filenames
    shard_map       = "300:25,500:20,1000:10",
    memory_map      = "300:2GB,500:2GB,1000:3GB",
    disk_map        = "all:2GB",
    rep_subset      = "",          # e.g. "1" or "3,17,42"; empty = all reps
    root_dir        = ""
)
args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
        stop("Arguments must be --key value pairs; got: ", key)
    }
    field <- gsub("-", "_", sub("^--", "", key))
    if (!field %in% names(opts)) stop("Unknown option: ", key)
    opts[[field]] <- args[[i + 1L]]
    i <- i + 2L
}
opts$n_test  <- as.integer(opts$n_test)
opts$total_r <- as.integer(opts$total_r)
opts$alpha   <- as.numeric(opts$alpha)
stopifnot(opts$total_r >= 1L, opts$n_test >= 10L,
          is.finite(opts$alpha), opts$alpha > 0, opts$alpha < 1)

script_dir <- tryCatch({
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else "."
}, error = function(e) ".")
root_dir <- if (nzchar(opts$root_dir)) normalizePath(opts$root_dir) else normalizePath(file.path(script_dir, ".."))
osg_dir  <- file.path(root_dir, "osg")
if (!file.exists(file.path(root_dir, "src", "linWB_dgm_registry.R"))) {
    stop("root_dir does not look like the HAPS repo root: ", root_dir)
}
source(file.path(root_dir, "src", "linWB_dgm_registry.R"))
source(file.path(root_dir, "src", "gtau_tilt_core.R"))

parse_csv <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
parse_map <- function(x, what, sizes) {
    entries <- parse_csv(x)
    kv <- strsplit(entries, ":", fixed = TRUE)
    bad <- vapply(kv, function(p) length(p) != 2L, logical(1))
    if (any(bad)) stop("Malformed ", what, " entry: ", paste(entries[bad], collapse = ", "))
    keys <- vapply(kv, `[[`, "", 1L); vals <- vapply(kv, `[[`, "", 2L)
    out <- setNames(vals, keys)
    resolve <- function(n) {
        if (as.character(n) %in% names(out)) out[[as.character(n)]]
        else if ("all" %in% names(out)) out[["all"]]
        else stop(what, " has no entry for n=", n, " and no 'all' fallback")
    }
    setNames(vapply(sizes, resolve, ""), as.character(sizes))
}

setups <- parse_csv(opts$setups)
for (s in setups) invisible(get_linWB_dgm_par(s))    # validates each setup
if (!opts$models %in% names(gtau_model_cases())) {
    stop("Unknown --models case '", opts$models, "'. Available: ",
         paste(names(gtau_model_cases()), collapse = ", "))
}
sizes  <- as.integer(parse_csv(opts$sample_sizes))
stopifnot(all(is.finite(sizes)), all(sizes >= 10L))
shard_size <- vapply(parse_map(opts$shard_map, "--shard-map", sizes),
                     function(x) as.integer(x), integer(1))
stopifnot(all(shard_size >= 1L))
memory_by_n <- parse_map(opts$memory_map, "--memory-map", sizes)
disk_by_n   <- parse_map(opts$disk_map,   "--disk-map",   sizes)

rep_subset <- if (nzchar(opts$rep_subset)) sort(unique(as.integer(parse_csv(opts$rep_subset)))) else seq_len(opts$total_r)
if (any(!is.finite(rep_subset)) || any(rep_subset < 1L) || any(rep_subset > opts$total_r)) {
    stop("--rep-subset entries must be replication indices in 1..total_r")
}

## ---------------------- Job directory (refuse overwrite) --------------
jobs_dir    <- file.path(osg_dir, "jobs", opts$experiment_name)
configs_dir <- file.path(jobs_dir, "configs")
if (dir.exists(jobs_dir) && length(list.files(jobs_dir, recursive = TRUE))) {
    stop("Jobs directory already exists and is not empty: ", jobs_dir,
         "\nChoose a different --experiment-name or remove it yourself first.")
}
dir.create(configs_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------- Shards: contiguous runs, chunked -----------------
## rep_subset is split into runs of consecutive indices; each run is chunked
## by the per-n shard size. With the default full subset this is just
## contiguous blocks 1..shard, shard+1..2*shard, ...
chunk_runs <- function(subset, size) {
    runs <- split(subset, cumsum(c(1L, diff(subset) != 1L)))
    out <- list()
    for (run in runs) {
        starts <- seq(1L, length(run), by = size)
        for (s in starts) out[[length(out) + 1L]] <- run[s:min(s + size - 1L, length(run))]
    }
    out
}

## ------------------------- Build all jobs -----------------------------
manifest <- list()
job_counter <- 0L
for (setup in setups) {
    dgm <- get_linWB_dgm_par(setup)
    registry_check <- list(tau_grid = c(0, dgm$change_times), delta_grid = dgm$delta_grid,
                           change_times = dgm$change_times, rho = dgm$rho, T_max = dgm$T_max)
    for (n in sizes) {
        ## Canonical driver seeding for this cell (identical for every cell,
        ## exactly as the serial driver would draw them at R = total_r).
        set.seed(123)
        seeds       <- sample(1e7, opts$total_r, replace = FALSE)
        seeds_test  <- sample(1e7, opts$total_r, replace = FALSE)
        seeds_split <- sample(1e7, opts$total_r, replace = FALSE)

        shards <- chunk_runs(rep_subset, shard_size[[as.character(n)]])
        for (reps in shards) {
            job_counter <- job_counter + 1L
            job_id <- sprintf("%04d", job_counter)
            cell_dir <- if (identical(opts$models, "cox")) {
                file.path(setup, sprintf("n%d", n))       # legacy layout
            } else {
                file.path(setup, opts$models, sprintf("n%d", n))
            }
            relpath <- file.path(cell_dir,
                                 sprintf("job_%s_reps_%03d_%03d.rds",
                                         job_id, reps[[1L]], reps[[length(reps)]]))
            config <- list(
                job_id = job_id,
                experiment_name = opts$experiment_name,
                driver_script = "gtau_job.R",
                setup = setup, models = opts$models,
                n = n, n_test = opts$n_test, alpha = opts$alpha,
                total_R = opts$total_r,
                rep_start = reps[[1L]], rep_end = reps[[length(reps)]],
                rep_count = length(reps), rep_indices = reps,
                result_relpath = relpath,
                seeds = seeds[reps], seeds_test = seeds_test[reps],
                seeds_split = seeds_split[reps],
                registry_check = registry_check
            )
            saveRDS(config, file.path(configs_dir, sprintf("job_%s.rds", job_id)))
            manifest[[job_counter]] <- data.frame(
                job_id = job_id, setup = setup, models = opts$models,
                n = n, n_test = opts$n_test,
                alpha = opts$alpha, total_R = opts$total_r,
                rep_start = reps[[1L]], rep_end = reps[[length(reps)]],
                rep_count = length(reps), result_relpath = relpath,
                request_memory = memory_by_n[[as.character(n)]],
                request_disk = disk_by_n[[as.character(n)]],
                stringsAsFactors = FALSE)
        }
    }
}
manifest <- do.call(rbind, manifest)

## --------------------------- Write outputs ----------------------------
write.table(manifest, file.path(jobs_dir, "manifest.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(sprintf("%s %s %s %s", manifest$job_id, manifest$result_relpath,
                   manifest$request_memory, manifest$request_disk),
           file.path(jobs_dir, "queue.txt"))
writeLines(manifest$job_id, file.path(jobs_dir, "job_ids.txt"))
saveRDS(list(options = opts, setups = setups, sizes = sizes,
             rep_subset = rep_subset, manifest = manifest,
             created = format(Sys.time(), tz = "UTC", usetz = TRUE)),
        file.path(jobs_dir, "plan.rds"))

## Pre-create the result/log trees Condor writes into (submit files do not
## create directories; a missing dir makes every job go on hold).
for (rel in unique(dirname(manifest$result_relpath))) {
    dir.create(file.path(root_dir, "results", "osg", opts$experiment_name, "raw", rel),
               recursive = TRUE, showWarnings = FALSE)
}
dir.create(file.path(root_dir, "results", "osg", opts$experiment_name, "collected"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "results", "osg", opts$experiment_name, "combined"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "logs", "osg", opts$experiment_name, "condor"),
           recursive = TRUE, showWarnings = FALSE)

## ----------------------------- Summary --------------------------------
cells <- unique(manifest[, c("setup", "n")])
cat(sprintf("Prepared experiment '%s' (models=%s): %d jobs over %d cells (%s), reps %s of total_R=%d.\n",
            opts$experiment_name, opts$models, nrow(manifest), nrow(cells),
            paste(sprintf("%s/n%d", cells$setup, cells$n), collapse = ", "),
            if (length(rep_subset) == opts$total_r) "1..all"
            else paste(range(rep_subset), collapse = ".."),
            opts$total_r))
cat("Wrote:", file.path(jobs_dir, c("manifest.tsv", "queue.txt", "job_ids.txt", "plan.rds")),
    sep = "\n  ")
cat("\nNext (from the OSG access point):\n")
cat(sprintf("  cd %s\n", osg_dir))
cat(sprintf("  condor_submit experiment_name=%s submit.sub\n", opts$experiment_name))
cat(sprintf("Then collect with:\n  Rscript %s --experiment-name %s\n",
            file.path(osg_dir, "collect_gtau_results.R"), opts$experiment_name))
