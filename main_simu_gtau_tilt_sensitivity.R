rm(list = ls())

## =====================================================================
## Sensitivity analysis: robustness of HAPS-A to a censoring-distribution
## shift at the prediction time tau, using the tilted G_tau option.
##
## Setup: linWB1 / DGM1, correctly specified Cox prediction + Cox G + Cox S
## (so the ONLY thing varying is the Gtau handling, isolating its effect).
##
## For each prediction time tau we calibrate HAPS-A under:
##   - Gtau_mode = "one"      (the paper method; ignores any population shift)
##   - Gtau_mode = "tilted"   at each delta_cal in a grid (delta_cal = 0 is the
##                             "estimated" true-law anchor)
## and then EVALUATE every calibrated interval on the at-risk population
##   { T > tau, C_tilde > tau },   C_tilde | path ~ tilted TRUE censoring law,
## for each delta_eval in the grid, via the oracle weighting identity
##   coverage = sum_i w_i 1{T_i in C_hat_i} / sum_i w_i,   w_i = tilde G_{delta_eval}(tau | path_i)
## computed on UNCENSORED test data (full covariate path retained).
##
## The (delta_cal x delta_eval) sweep gives, in one pass:
##   * matched tilt   (delta_cal == delta_eval, the diagonal): should stay ~nominal,
##   * mismatched tilt (off-diagonal): robustness to choosing the wrong delta,
##   * "one" and "estimated" (delta_cal = 0) rows for reference.
##
## The per-replication logic lives in src/gtau_tilt_core.R and is SHARED with
## the OSG worker (osg/gtau_job.R); this script is the serial local driver.
## =====================================================================

if (!file.exists("src/gen_ICML_simu.R")) {
    stop("Please run main_simu_gtau_tilt_sensitivity.R from the repository root directory.")
}

source("src/gen_ICML_simu.R")
source("src/helpers.R")
source("src/helpers_AIPCW.R")
source("src/dynamicCP.R")
source("src/dynamicCP_AIPCW.R")
source("src/gtau_eval_helpers.R")
source("src/linWB_dgm_registry.R")
source("src/gtau_tilt_core.R")

suppressMessages({
    library(survival)
    library(randomForestSRC)
    library(dplyr)
})
if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required (model.Xi = 'xgb_reg').")
}

## ----------------------------- Config --------------------------------
## CLI: R n n_test [setup]   with setup in {linWB1 (default), linWB2}.
## linWB1 = original paper DGM1 (unbounded support; widest set from the training
## grid; positivity via trim.C -- see docs/dgm_positivity_notes.md).
## linWB2 = positivity-respecting DGM (T truncated at T_max = 20; retuned
## censoring; widest candidate set = (tau, T_max], so Algorithm 1 is feasible
## by construction).
setup  <- "linWB1"

R      <- 200
n      <- 1000
n_test <- 1000

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && nzchar(args[[1L]])) R      <- as.integer(args[[1L]])
if (length(args) >= 2L && nzchar(args[[2L]])) n      <- as.integer(args[[2L]])
if (length(args) >= 3L && nzchar(args[[3L]])) n_test <- as.integer(args[[3L]])
if (length(args) >= 4L && nzchar(args[[4L]])) setup  <- args[[4L]]
stopifnot(is.finite(R), R >= 1L, is.finite(n), n >= 10L, is.finite(n_test), n_test >= 10L)

## All study constants (models, grids, DGM parameters, variable names) come from
## the shared cfg builder -- the single source of truth used by the OSG worker
## too (src/gtau_tilt_core.R).
cfg        <- gtau_study_cfg(setup, n, n_test)   # validates setup; alpha = 0.1
alpha      <- cfg$alpha
tau_grid   <- cfg$tau_grid
delta_grid <- cfg$delta_grid

folder <- file.path("results", setup, "gtau_tilt")
if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)

cat(sprintf("Gtau tilt sensitivity | setup=%s | R=%d n=%d n_test=%d | tau={%s} | delta={%s}\n",
            setup, R, n, n_test, paste(tau_grid, collapse=","), paste(delta_grid, collapse=",")))

## --------------------------- Seeds -----------------------------------
set.seed(123)
seeds       <- sample(1e7, R, replace = FALSE)
seeds_test  <- sample(1e7, R, replace = FALSE)
seeds_split <- sample(1e7, R, replace = FALSE)

## --------------------- Output file + checkpointing -------------------
## The directory is re-created defensively at EVERY save (an external cleanup
## of results/ during a long run must not be able to lose finished work), and
## partial results are checkpointed every `checkpoint_every` replications so a
## crash costs minutes, not hours. meta$complete marks the final save.
outfile <- gtau_sensitivity_filename(setup, R, n, n_test, alpha)
checkpoint_every <- 10L
save_results <- function(rows_list, completed_reps, elapsed_sec) {
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    out <- build_gtau_result_object(
        rows_list, cfg, R, seeds, seeds_test, seeds_split,
        elapsed_sec = elapsed_sec, completed_reps = completed_reps)
    saveRDS(out, outfile)
    invisible(out)
}

## --------------------------- Main loop -------------------------------
rows <- list()
t0 <- proc.time()[["elapsed"]]
for (r in seq_len(R)) {
    if (r %% 10 == 0 || r == 1) cat("rep", r, "of", R, "\n")

    rows <- c(rows, run_one_gtau_rep(r, seeds[r], seeds_test[r], seeds_split[r], cfg))

    if (r %% checkpoint_every == 0L && r < R) {
        save_results(rows, r, proc.time()[["elapsed"]] - t0)
        cat(sprintf("  [checkpoint saved through rep %d]\n", r))
    }
}
res <- do.call(rbind, rows)
elapsed <- proc.time()[["elapsed"]] - t0
cat(sprintf("done in %.1f min\n", elapsed / 60))

## --------------------------- Save ------------------------------------
save_results(rows, R, elapsed)
cat("Saved:", outfile, "\n")

## ----------------- Quick console summary (matched vs one) ------------
# Mean over reps of the matched-tilt diagonal (delta_cal == delta_eval) and of
# the 'one' method, per (tau, delta_eval). Computed with base ave() to avoid
# aggregate() dropping the NA delta_cal rows of the 'one' method.
ok <- res[res$ok, ]
mean_cov <- function(sub) tapply(sub$coverage, sub$delta_eval,
                                 function(x) mean(x, na.rm = TRUE))
matched_sub <- ok[ok$method %in% c("estimated", "tilted") &
                  abs(ok$delta_cal - ok$delta_eval) < 1e-9, ]
one_sub <- ok[ok$method == "one", ]
cat("\nMean coverage by delta_eval={", paste(delta_grid, collapse=","),
    "} (nominal ", 1 - alpha, "):\n", sep = "")
for (tt in tau_grid) {
    md <- mean_cov(matched_sub[matched_sub$tau == tt, ])
    oa <- mean_cov(one_sub[one_sub$tau == tt, ])
    cat(sprintf(" tau=%g  matched-tilt: %s | 'one': %s\n", tt,
        paste(sprintf("%.3f", md[as.character(delta_grid)]), collapse = " "),
        paste(sprintf("%.3f", oa[as.character(delta_grid)]), collapse = " ")))
}

## Infeasibility rate: share of calibrations per (tau, arm) where Algorithm 1
## had no feasible theta (ok = FALSE; these reps are EXCLUDED from the means
## above, which makes affected cells slightly optimistic — report alongside).
cal_rows <- res[abs(res$delta_eval) < 1e-9, ]          # one row per calibration
arm_lab <- ifelse(cal_rows$method == "tilted",
                  sprintf("tilted(dc=%+.2f)", cal_rows$delta_cal), cal_rows$method)
infeas <- tapply(!cal_rows$ok, list(cal_rows$tau, arm_lab), mean)
cat("\nInfeasible-calibration rate by tau x arm (excluded from coverage means):\n")
print(round(infeas, 4))
