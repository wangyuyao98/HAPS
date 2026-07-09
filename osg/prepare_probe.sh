#!/bin/bash

# Probe run: replication 1 of every (setup, n) cell, one rep per job, with
# conservative resource requests. Seeds are drawn at the CANONICAL total_R
# (200), so probe outputs are the exact rep-1 results of the production
# study — usable both for resource calibration and for a cross-machine
# consistency check against local runs.
#
# Override any of these via environment variables, e.g.:
#   MEMORY_MAP=all:6GB bash osg/prepare_probe.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

EXPERIMENT_NAME="${EXPERIMENT_NAME:-gtau_probe}"
SETUPS="${SETUPS:-linWB1,linWB2}"
SAMPLE_SIZES="${SAMPLE_SIZES:-300,500,1000}"
N_TEST="${N_TEST:-3000}"
TOTAL_R="${TOTAL_R:-200}"
REP_SUBSET="${REP_SUBSET:-1}"
MEMORY_MAP="${MEMORY_MAP:-all:4GB}"
DISK_MAP="${DISK_MAP:-all:2GB}"

Rscript "${SCRIPT_DIR}/prepare_gtau_jobs.R" \
  --experiment-name "${EXPERIMENT_NAME}" \
  --setups "${SETUPS}" \
  --sample-sizes "${SAMPLE_SIZES}" \
  --n-test "${N_TEST}" \
  --total-r "${TOTAL_R}" \
  --rep-subset "${REP_SUBSET}" \
  --shard-map all:1 \
  --memory-map "${MEMORY_MAP}" \
  --disk-map "${DISK_MAP}"
