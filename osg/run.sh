#!/bin/bash

set -euo pipefail

JOB_ID="${1:-}"

if [[ -z "${JOB_ID}" ]]; then
    echo "Usage: run.sh <job_id>" >&2
    exit 1
fi

export TMPDIR="${_CONDOR_SCRATCH_DIR:-$PWD}"

# One core is requested (request_cpus = 1); pin OpenMP so xgboost does not
# oversubscribe the slot (also makes worker threading uniform across sites).
export OMP_NUM_THREADS=1

CONFIG_PATH="$(find . -name "job_${JOB_ID}.rds" -print -quit)"

if [[ -z "${CONFIG_PATH}" ]]; then
    echo "Could not find job config for job_id=${JOB_ID}" >&2
    find . -maxdepth 3 -type f >&2 || true
    exit 1
fi

Rscript gtau_job.R --config "${CONFIG_PATH}" --output result.rds
