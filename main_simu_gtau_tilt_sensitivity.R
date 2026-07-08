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

suppressMessages({
    library(survival)
    library(randomForestSRC)
    library(dplyr)
})
if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required (model.Xi = 'xgb_reg').")
}

## ----------------------------- Config --------------------------------
setup    <- "linWB1"
dgm_name <- "linear_weibull"
folder   <- file.path("results", setup, "gtau_tilt")
if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)

rho    <- 0.3
n      <- 1000
n_test <- 1000
R      <- 200

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && nzchar(args[[1L]])) R      <- as.integer(args[[1L]])
if (length(args) >= 2L && nzchar(args[[2L]])) n      <- as.integer(args[[2L]])
if (length(args) >= 3L && nzchar(args[[3L]])) n_test <- as.integer(args[[3L]])
stopifnot(is.finite(R), R >= 1L, is.finite(n), n >= 10L, is.finite(n_test), n_test >= 10L)

# variable names
start.name <- "start"; stop.name <- "stop"; event.name <- "event"
event.C.name <- "event.C"; id.name <- "id"; TT.name <- "TT"

# covariates
covname.pred.timevarying <- c("L")
covname.C.timevarying    <- c("L")

# models: correctly-specified Cox everywhere; regression Xi (multiclass is
# invalid once g_tau != 1, so use xgb_reg uniformly across ALL arms so that the
# only difference between arms is the Gtau handling, not the Xi learner).
model.pred <- "cox"; model.C <- "cox"; model.S <- "cox"; model.Xi <- "xgb_reg"

alpha             <- 0.1
theta_grid_AIPCW  <- seq(0, 0.5, by = 0.01)
train_ratio       <- 0.5

# DGM1 parameters (linWB1)
change_times <- c(3, 6)
dgm_parK <- list(
    shape_T = c(4, 4, 5), shape_C = c(3, 3, 3),
    beta_T0 = c(-8, -8, -5), beta_TL = c(1, 2, 3),
    beta_C0 = c(-6, -6, -5), beta_CL = c(2, 2, 2)
)
tau_grid   <- c(0, change_times)                 # 0, 3, 6
delta_grid <- c(-0.30, -0.15, 0, 0.15, 0.30)     # tilt grid (cal and eval share it)

cat(sprintf("Gtau tilt sensitivity | R=%d n=%d n_test=%d | tau={%s} | delta={%s}\n",
            R, n, n_test, paste(tau_grid, collapse=","), paste(delta_grid, collapse=",")))

## --------------------------- Seeds -----------------------------------
set.seed(123)
seeds       <- sample(1e7, R, replace = FALSE)
seeds_test  <- sample(1e7, R, replace = FALSE)
seeds_split <- sample(1e7, R, replace = FALSE)

## ------------------------- Calibration wrapper -----------------------
calibrate <- function(dat_long, dat_test, tau, gtau_mode, gtau_delta, seed_split) {
    tryCatch(
        dynamicCP_AIPCW_split(
            dat = dat_long, dat_test = dat_test, tau = tau,
            start.name = start.name, stop.name = stop.name,
            event.name = event.name, event.C.name = event.C.name, id.name = id.name,
            covname.pred.timevarying = covname.pred.timevarying,
            model.pred = model.pred,
            model.C = model.C, covname.C.timevarying = covname.C.timevarying,
            model.S = model.S, model.Xi = model.Xi,
            theta_grid = theta_grid_AIPCW, alpha = alpha, trim.C = 0.05,
            TT.name = TT.name, visit_times = tau_grid,
            seed = seed_split, train_ratio = train_ratio,
            localization = FALSE, censoring_method = "piecewise",
            Gtau_mode = gtau_mode, Gtau_delta = gtau_delta,
            return_on_aipcw_fail = TRUE
        ),
        error = function(e) list(error = TRUE, msg = conditionMessage(e))
    )
}

# Turn one calibrated fit into per-delta_eval evaluation rows (weighted coverage
# on {T>tau, C_tilde>tau} for each delta_eval). Bounds/covered_T are computed once;
# only the weight vector changes with delta_eval.
eval_rows <- function(fit, method, delta_cal, tau, r, L_full, TT_by_id) {
    if (is.null(fit) || isTRUE(fit$error)) {
        return(lapply(delta_grid, function(de) data.frame(
            rep = r, tau = tau, method = method, delta_cal = delta_cal, delta_eval = de,
            coverage = NA_real_, med_len = NA_real_, ess = NA_real_,
            theta_hat = NA_real_, ok = FALSE, stringsAsFactors = FALSE)))
    }
    ids   <- as.character(fit$id_test_tau)
    lower <- fit$lower_AIPCW_test
    upper <- fit$upper_AIPCW_test
    TT    <- TT_by_id[ids]
    keep  <- is.finite(TT) & TT > tau              # survivors at tau (uncensored test)
    ids <- ids[keep]; lower <- lower[keep]; upper <- upper[keep]; TT <- TT[keep]
    covered <- as.integer(lower <= TT & TT <= upper)
    len     <- upper - lower
    Lsub    <- L_full[ids, , drop = FALSE]
    lapply(delta_grid, function(de) {
        w <- compute_Gtau_true_tilted_linWB1(Lsub, tau, de, change_times, dgm_parK)
        s <- weighted_coverage_summary(covered, len, w)
        data.frame(rep = r, tau = tau, method = method, delta_cal = delta_cal,
                   delta_eval = de, coverage = s$coverage, med_len = s$med_len,
                   ess = s$ess, theta_hat = fit$theta_hat_AIPCW, ok = TRUE,
                   stringsAsFactors = FALSE)
    })
}

## --------------------------- Main loop -------------------------------
rows <- list()
t0 <- proc.time()[["elapsed"]]
for (r in seq_len(R)) {
    if (r %% 10 == 0 || r == 1) cat("rep", r, "of", R, "\n")

    dat_tr <- simulate_dataset_long(
        n = n, seed = seeds[r], change_times = change_times, tau_max = Inf,
        no_censoring = FALSE, par = dgm_parK, rho = rho, dgm_name = dgm_name)
    dat_te <- simulate_dataset_long(
        n = n_test, seed = seeds_test[r], change_times = change_times, tau_max = Inf,
        no_censoring = TRUE, par = dgm_parK, rho = rho, dgm_name = dgm_name,
        return_L_full = TRUE)
    L_full   <- attr(dat_te, "L_full")
    rownames(L_full) <- as.character(seq_len(n_test))
    TT_by_id <- setNames(dat_te[[TT.name]][match(seq_len(n_test), dat_te[[id.name]])],
                         as.character(seq_len(n_test)))

    for (tau in tau_grid) {
        # calibrations: "one" + tilt family (delta_cal = 0 is the "estimated" anchor)
        fit_one  <- calibrate(dat_tr, dat_te, tau, "one", 0, seeds_split[r])
        rows[[length(rows) + 1L]] <- do.call(rbind,
            eval_rows(fit_one, "one", NA_real_, tau, r, L_full, TT_by_id))

        for (dcal in delta_grid) {
            method <- if (dcal == 0) "estimated" else "tilted"
            fit <- calibrate(dat_tr, dat_te, tau, "tilted", dcal, seeds_split[r])
            rows[[length(rows) + 1L]] <- do.call(rbind,
                eval_rows(fit, method, dcal, tau, r, L_full, TT_by_id))
        }
    }
}
res <- do.call(rbind, rows)
elapsed <- proc.time()[["elapsed"]] - t0
cat(sprintf("done in %.1f min\n", elapsed / 60))

## --------------------------- Save ------------------------------------
out <- list(
    results = res,
    config = list(setup = setup, dgm_name = dgm_name, rho = rho, n = n, n_test = n_test,
                  R = R, alpha = alpha, tau_grid = tau_grid, delta_grid = delta_grid,
                  model.pred = model.pred, model.C = model.C, model.S = model.S,
                  model.Xi = model.Xi, change_times = change_times, dgm_parK = dgm_parK),
    seeds = list(seeds = seeds, seeds_test = seeds_test, seeds_split = seeds_split),
    meta = list(script = "main_simu_gtau_tilt_sensitivity.R", elapsed_sec = elapsed)
)
outfile <- file.path(folder, sprintf("gtau_tilt_sensitivity_R%d_n%d_ntest%d_alpha%s.rds",
                                     R, n, n_test, format(alpha)))
saveRDS(out, outfile)
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
