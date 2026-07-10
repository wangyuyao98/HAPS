#!/usr/bin/env Rscript

## =====================================================================
## OSG worker for the Gtau tilt sensitivity study (one shard = a block of
## replication indices within one (setup, n) cell).
##
## Contract (same as the old osg_old/mixC_job.R):
##   Rscript gtau_job.R --config <job_XXXX.rds> --output <result.rds>
##
## The per-replication logic is src/gtau_tilt_core.R::run_one_gtau_rep(),
## SHARED with the local serial driver main_simu_gtau_tilt_sensitivity.R.
## All study constants come from gtau_study_cfg(); the shard config only
## contributes (setup, n, n_test, alpha), the replication indices, and the
## three per-rep seed slices pre-drawn at prepare time with the driver's
## canonical scheme (set.seed(123) at total_R), so a shard reproduces the
## corresponding replications of the serial driver exactly (same R and
## package versions).
## =====================================================================

## ------------------------- CLI parsing --------------------------------
args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
output_path <- "result.rds"
i <- 1L
while (i <= length(args)) {
    if (args[[i]] == "--config" && i < length(args)) {
        config_path <- args[[i + 1L]]; i <- i + 2L
    } else if (args[[i]] == "--output" && i < length(args)) {
        output_path <- args[[i + 1L]]; i <- i + 2L
    } else {
        stop("Unknown argument: ", args[[i]], " (usage: gtau_job.R --config <rds> [--output <rds>])")
    }
}
if (is.null(config_path)) stop("--config is required")
if (!file.exists(config_path)) stop("Config file not found: ", config_path)

## --------------------- Locate and source the repo code ----------------
## On an OSG worker, Condor drops the transferred `src` directory (and this
## script) at the top of the scratch dir; locally the script is run from the
## repository root (or via its path inside osg/).
script_dir <- tryCatch({
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else "."
}, error = function(e) ".")
src_candidates <- c("src", file.path(script_dir, "src"), file.path(script_dir, "..", "src"))
src_dir <- NULL
for (cand in src_candidates) {
    if (file.exists(file.path(cand, "gen_ICML_simu.R"))) { src_dir <- cand; break }
}
if (is.null(src_dir)) {
    stop("Could not locate the src/ directory (tried: ",
         paste(src_candidates, collapse = ", "), ")")
}
source(file.path(src_dir, "gen_ICML_simu.R"))
source(file.path(src_dir, "helpers.R"))
source(file.path(src_dir, "helpers_AIPCW.R"))
source(file.path(src_dir, "dynamicCP.R"))
source(file.path(src_dir, "dynamicCP_AIPCW.R"))
source(file.path(src_dir, "gtau_eval_helpers.R"))
source(file.path(src_dir, "linWB_dgm_registry.R"))
source(file.path(src_dir, "gtau_tilt_core.R"))

suppressMessages(library(survival))
if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required (model.Xi = 'xgb_reg').")
}
## NOTE: unlike the local driver, the worker does NOT attach dplyr (never
## exercised on this path); randomForestSRC is attached below only when the
## models case needs it.

## ------------------------- Config validation --------------------------
config <- readRDS(config_path)
required <- c("job_id", "experiment_name", "setup", "n", "n_test", "alpha",
              "total_R", "rep_start", "rep_end", "rep_count", "rep_indices",
              "result_relpath", "seeds", "seeds_test", "seeds_split",
              "registry_check")
missing <- setdiff(required, names(config))
if (length(missing)) stop("Job config is missing fields: ", paste(missing, collapse = ", "))
stopifnot(
    length(config$rep_indices) == config$rep_count,
    length(config$seeds)       == config$rep_count,
    length(config$seeds_test)  == config$rep_count,
    length(config$seeds_split) == config$rep_count,
    config$rep_indices[[1L]] == config$rep_start,
    config$rep_indices[[config$rep_count]] == config$rep_end
)

## models case: configs prepared before the model-case feature carry no field
## and are the default all-cox case.
models <- if (!is.null(config$models)) config$models else "cox"
cfg <- gtau_study_cfg(config$setup, config$n, config$n_test, config$alpha,
                      models = models)
if (isTRUE(cfg$needs_rsf)) {
    if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
        stop("Models case '", models, "' needs package 'randomForestSRC'.")
    }
    suppressMessages(library(randomForestSRC))
}

## Drift guard: the registry values embedded at prepare time must match what
## this worker's copy of src/ resolves — otherwise the shard would silently
## compute a different study than the one prepared/submitted.
rc <- config$registry_check
for (fld in c("tau_grid", "delta_grid", "change_times", "rho", "T_max")) {
    if (!isTRUE(all.equal(rc[[fld]], cfg[[fld]]))) {
        stop("Registry drift for '", fld, "': prepare-time value (",
             paste(rc[[fld]], collapse = ","), ") != worker value (",
             paste(cfg[[fld]], collapse = ","), "). Re-run prepare_gtau_jobs.R ",
             "against the same code version as the worker.")
    }
}

cat(sprintf("gtau shard | exp=%s job=%s | setup=%s models=%s n=%d n_test=%d alpha=%s | reps %d-%d (%d of total_R=%d)\n",
            config$experiment_name, config$job_id, config$setup, models, config$n,
            config$n_test, format(config$alpha), config$rep_start, config$rep_end,
            config$rep_count, config$total_R))

## --------------------------- Run the shard ----------------------------
rows <- list()
rep_seconds <- numeric(config$rep_count)
t0 <- proc.time()[["elapsed"]]
for (k in seq_len(config$rep_count)) {
    r <- config$rep_indices[[k]]
    t1 <- proc.time()[["elapsed"]]
    rows <- c(rows, run_one_gtau_rep(r, config$seeds[[k]], config$seeds_test[[k]],
                                     config$seeds_split[[k]], cfg))
    rep_seconds[[k]] <- proc.time()[["elapsed"]] - t1
    cat(sprintf("  rep %d done (%.1f s)\n", r, rep_seconds[[k]]))
}
elapsed <- proc.time()[["elapsed"]] - t0

## ------------------------------ Save ----------------------------------
## `rows` is kept as the LIST of per-(tau, arm) data.frames — the collector
## concatenates shard lists in replication order and rbinds ONCE, exactly as
## the serial driver does, so the combined `results` table is identical.
save_obj <- list(
    job = list(job_id = config$job_id, experiment_name = config$experiment_name,
               setup = config$setup, models = models,
               n = config$n, n_test = config$n_test,
               alpha = config$alpha, total_R = config$total_R,
               rep_start = config$rep_start, rep_end = config$rep_end,
               rep_indices = config$rep_indices,
               result_relpath = config$result_relpath),
    seeds = list(seeds = config$seeds, seeds_test = config$seeds_test,
                 seeds_split = config$seeds_split),
    results = list(rows = rows, res = do.call(rbind, rows)),
    meta = list(script = "osg/gtau_job.R",
                elapsed_sec = elapsed,
                rep_seconds = setNames(rep_seconds, config$rep_indices),
                hostname = tryCatch(Sys.info()[["nodename"]], error = function(e) NA_character_),
                r_version = as.character(getRversion()),
                pkg_versions = c(
                    survival = as.character(utils::packageVersion("survival")),
                    xgboost  = as.character(utils::packageVersion("xgboost"))))
)
saveRDS(save_obj, output_path)
cat(sprintf("shard done in %.1f min | saved: %s\n", elapsed / 60, output_path))
