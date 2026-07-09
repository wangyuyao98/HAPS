#!/usr/bin/env Rscript

## =====================================================================
## Collect OSG shard outputs of the Gtau tilt sensitivity study.
##
## For every (setup, n) cell whose shards are ALL present and jointly cover
## replications 1..total_R, rebuilds the CANONICAL result object (identical
## structure to what main_simu_gtau_tilt_sensitivity.R saves — the plot
## script reads it unchanged) and writes it under
##   results/osg/<experiment>/collected/<setup>/gtau_tilt/<canonical name>
## NEVER under the top-level results/<setup>/ tree: promoting a collected
## file into the canonical location is a deliberate manual step (protects
## locally produced files from being overwritten).
##
## Incomplete cells (probe runs, failed shards) are reported with their
## missing job_ids and skipped; a combined archive with everything found so
## far is always written to results/osg/<experiment>/combined/.
##
## Usage: Rscript osg/collect_gtau_results.R --experiment-name <name>
## =====================================================================

opts <- list(experiment_name = "", root_dir = "")
args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) stop("Arguments must be --key value pairs; got: ", key)
    field <- gsub("-", "_", sub("^--", "", key))
    if (!field %in% names(opts)) stop("Unknown option: ", key)
    opts[[field]] <- args[[i + 1L]]
    i <- i + 2L
}
if (!nzchar(opts$experiment_name)) stop("--experiment-name is required")

script_dir <- tryCatch({
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else "."
}, error = function(e) ".")
root_dir <- if (nzchar(opts$root_dir)) normalizePath(opts$root_dir) else normalizePath(file.path(script_dir, ".."))
source(file.path(root_dir, "src", "linWB_dgm_registry.R"))
source(file.path(root_dir, "src", "gtau_tilt_core.R"))

exp <- opts$experiment_name
jobs_dir <- file.path(root_dir, "osg", "jobs", exp)
raw_root <- file.path(root_dir, "results", "osg", exp, "raw")
manifest_path <- file.path(jobs_dir, "manifest.tsv")
if (!file.exists(manifest_path)) stop("No manifest found at ", manifest_path)
manifest <- read.delim(manifest_path, colClasses = "character")
for (col in c("n", "n_test", "total_R", "rep_start", "rep_end", "rep_count")) {
    manifest[[col]] <- as.integer(manifest[[col]])
}
manifest$alpha <- as.numeric(manifest$alpha)
manifest$raw_path <- file.path(raw_root, manifest$result_relpath)
manifest$found <- file.exists(manifest$raw_path)

cat(sprintf("Experiment '%s': %d/%d shard outputs present.\n",
            exp, sum(manifest$found), nrow(manifest)))
if (any(!manifest$found)) {
    cat("Missing job_ids:", paste(manifest$job_id[!manifest$found], collapse = " "), "\n")
}

## ------------------------- Per-cell collection ------------------------
cells <- unique(manifest[, c("setup", "n", "n_test", "alpha", "total_R")])
collected_files <- character(0)
cell_status <- list()
for (ci in seq_len(nrow(cells))) {
    cell <- cells[ci, ]
    sub <- manifest[manifest$setup == cell$setup & manifest$n == cell$n &
                    manifest$n_test == cell$n_test & manifest$total_R == cell$total_R, ]
    sub <- sub[order(sub$rep_start), ]
    covered <- sort(unlist(mapply(seq, sub$rep_start[sub$found], sub$rep_end[sub$found],
                                  SIMPLIFY = FALSE)))
    complete <- all(sub$found) && identical(as.integer(covered), seq_len(cell$total_R))
    label <- sprintf("%s/n%d", cell$setup, cell$n)
    cell_status[[label]] <- list(cell = cell, jobs = nrow(sub), found = sum(sub$found),
                                 complete = complete,
                                 missing_job_ids = sub$job_id[!sub$found])
    if (!complete) {
        cat(sprintf("  %-14s %d/%d shards, reps covered %d/%d -> skipped (incomplete)\n",
                    label, sum(sub$found), nrow(sub), length(covered), cell$total_R))
        next
    }

    ## Recompute the canonical seed vectors and cross-check every shard slice.
    set.seed(123)
    seeds       <- sample(1e7, cell$total_R, replace = FALSE)
    seeds_test  <- sample(1e7, cell$total_R, replace = FALSE)
    seeds_split <- sample(1e7, cell$total_R, replace = FALSE)

    rows_all <- list()
    elapsed_total <- 0
    worker_meta <- list()
    for (j in seq_len(nrow(sub))) {
        shard <- readRDS(sub$raw_path[[j]])
        idx <- shard$job$rep_indices
        if (!identical(as.numeric(shard$seeds$seeds), as.numeric(seeds[idx])) ||
            !identical(as.numeric(shard$seeds$seeds_test), as.numeric(seeds_test[idx])) ||
            !identical(as.numeric(shard$seeds$seeds_split), as.numeric(seeds_split[idx]))) {
            stop("Seed mismatch in shard job_", sub$job_id[[j]], " (", label,
                 "): shard slices do not match the canonical set.seed(123) vectors ",
                 "at total_R=", cell$total_R, ". Was the shard prepared with a ",
                 "different total_R?")
        }
        rows_all <- c(rows_all, shard$results$rows)
        elapsed_total <- elapsed_total + shard$meta$elapsed_sec
        worker_meta[[j]] <- c(r = shard$meta$r_version, shard$meta$pkg_versions)
    }

    cfg <- gtau_study_cfg(cell$setup, cell$n, cell$n_test, cell$alpha)
    obj <- build_gtau_result_object(
        rows_all, cfg, cell$total_R, seeds, seeds_test, seeds_split,
        elapsed_sec = elapsed_total, completed_reps = cell$total_R,
        script = "main_simu_gtau_tilt_sensitivity.R",
        meta_extra = list(via_osg = TRUE, experiment_name = exp,
                          n_jobs = nrow(sub),
                          worker_versions = unique(worker_meta)))
    outfile <- file.path(root_dir, gtau_sensitivity_filename(
        cell$setup, cell$total_R, cell$n, cell$n_test, cell$alpha,
        results_root = file.path("results", "osg", exp, "collected")))
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    saveRDS(obj, outfile)
    collected_files <- c(collected_files, outfile)
    cat(sprintf("  %-14s complete (%d shards, %.1f CPU-h) -> %s\n",
                label, nrow(sub), elapsed_total / 3600, outfile))
}

## --------------------------- Combined archive -------------------------
combined_file <- file.path(root_dir, "results", "osg", exp, "combined",
                           sprintf("%s_combined.rds", exp))
dir.create(dirname(combined_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(manifest = manifest, cell_status = cell_status,
             collected_files = collected_files,
             collected_at = format(Sys.time(), tz = "UTC", usetz = TRUE)),
        combined_file)
cat(sprintf("Combined archive: %s\n", combined_file))
if (length(collected_files)) {
    cat("Collected canonical files (promote into results/<setup>/gtau_tilt/ manually if desired):\n")
    cat(sprintf("  %s\n", collected_files))
}
