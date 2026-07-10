#!/bin/bash

# Production run: the full grid linWB1+linWB2 x n in {300, 500, 1000},
# R=200 replications, n_test=3000 (the study configuration agreed for the
# general-Gtau paper). Adjust MEMORY_MAP/DISK_MAP from the probe's observed
# peaks before submitting.
#
# Override any of these via environment variables, e.g.:
#   MEMORY_MAP=300:2GB,500:2GB,1000:4GB bash osg/prepare_full_grid.sh
#
# OSG access points have no system R; run the R step through the container:
#   RSCRIPT="apptainer exec /ospool/ap41/data/yuyao.wang/dCP1.sif Rscript" bash osg/prepare_full_grid.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${RSCRIPT:-Rscript}"

EXPERIMENT_NAME="${EXPERIMENT_NAME:-gtau_grid}"
SETUPS="${SETUPS:-linWB1,linWB2}"
SAMPLE_SIZES="${SAMPLE_SIZES:-300,500,1000}"
N_TEST="${N_TEST:-3000}"
TOTAL_R="${TOTAL_R:-200}"
SHARD_MAP="${SHARD_MAP:-300:25,500:20,1000:10}"
MEMORY_MAP="${MEMORY_MAP:-300:2GB,500:2GB,1000:3GB}"
DISK_MAP="${DISK_MAP:-all:2GB}"

${RSCRIPT} "${SCRIPT_DIR}/prepare_gtau_jobs.R" \
  --experiment-name "${EXPERIMENT_NAME}" \
  --setups "${SETUPS}" \
  --sample-sizes "${SAMPLE_SIZES}" \
  --n-test "${N_TEST}" \
  --total-r "${TOTAL_R}" \
  --shard-map "${SHARD_MAP}" \
  --memory-map "${MEMORY_MAP}" \
  --disk-map "${DISK_MAP}"
