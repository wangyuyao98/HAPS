rm(list = ls())

# Run this script from the repo root.
if (!file.exists("src/gen_ICML_simu.R")) {
    stop("Please run main_simu_dynamic_linWB1_cox_multiclass_loc0.R from the repository root directory.")
}

source("src/gen_ICML_simu.R")
source("src/helpers.R")
source("src/helpers_AIPCW.R")
source("src/dynamicCP.R")   # IPCW (previous implementation)
# Preserve IPCW-only implementation before loading the AIPCW implementation.
dynamicCP_IPCW_split <- dynamicCP_split
source("src/dynamicCP_AIPCW.R")      # IPCW + AIPCW implementation
source("src/gtau_eval_helpers.R")    # oracle Gtau weights for Gtau_mode="tilted" evaluation

library(survival)
library(randomForestSRC)
library(dplyr)

if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required for model.Xi='xgb_multiclass'.")
}


## ----------------------------- Setup ---------------------------------
setup <- "linWB1"
dgm_name <- "linear_weibull"   # "linear_weibull" or "mixC_uniform_cox"
method.list <- c("dCP_IPCW", "dCP_IPCW_new", "dCP_AIPCW")

folder <- file.path("results", setup)
if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
}

rho <- 0.3
n <- 1000
n_test <- 500
R <- 200

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && nzchar(args[[1L]])) R <- as.integer(args[[1L]])
if (length(args) >= 2L && nzchar(args[[2L]])) n <- as.integer(args[[2L]])
if (length(args) >= 3L && nzchar(args[[3L]])) n_test <- as.integer(args[[3L]])
if (!is.finite(R) || R < 1L) stop("R must be a positive integer.")
if (!is.finite(n) || n < 10L) stop("n must be at least 10.")
if (!is.finite(n_test) || n_test < 10L) stop("n_test must be at least 10.")

# variable names
start.name   <- "start"
stop.name    <- "stop"
event.name   <- "event"
event.C.name <- "event.C"
id.name      <- "id"

# covariates
covname.pred.baseline    <- NULL
covname.pred.timevarying <- c("L")
covname.C.baseline       <- NULL
covname.C.timevarying    <- c("L")

# model choices
# Paper-focused loc0 rerun for cox_cox_Scox_Xixgb_multiclass, Gtauone.
model.pred <- "cox"             # "cox", "rsf", "hal", "xgb_cox", "xgb_aft"
model.C    <- "cox"             # "km", "cox", "rsf", "hal", "xgb_cox", "xgb_aft"
model.S    <- "cox"             # "cox", "rsf", "hal", "xgb_cox", "xgb_aft"
model.Xi   <- "xgb_multiclass"  # "lm", "rsf", "xgb_reg", "xgb_multiclass", "hal_reg"

# Optional presets for quick switching:
# model.pred <- "cox";     model.C <- "cox";     model.S <- "cox";     model.Xi <- "lm"
# model.pred <- "cox";     model.C <- "rsf";     model.S <- "cox";     model.Xi <- "lm"
# model.pred <- "cox";     model.C <- "hal";     model.S <- "cox";     model.Xi <- "lm"
# model.pred <- "cox";     model.C <- "hal";     model.S <- "hal";     model.Xi <- "hal_reg"



Gtau_mode <- "one"      # "one", "estimated", or "tilted"
Gtau_delta <- 0         # tilt parameter used when Gtau_mode == "tilted"

# Evaluation method for the at-risk estimand {T>tau, C_tilde>tau} when Gtau_mode
# is "estimated" or "tilted" ("one" is unaffected):
#   "weighted" (default): oracle-weighted coverage on uncensored T>tau survivors
#       with weights w_i = G_tilde(tau | path_i) (delta = 0 under "estimated") ->
#       column coverage_test_gtau. Lowest-variance, deterministic given the test
#       set (Rao-Blackwellization of subgroup sampling).
#   "subgroup": physical at-risk subgroup construction.
#       - estimated: LEGACY behavior — censored test data, subgroup X>tau;
#         coverage_test_true is then the at-risk coverage (pre-existing outputs
#         are reproduced exactly under this setting).
#       - tilted: uncensored test data; only the membership indicator
#         I(C_tilde>tau) ~ Bernoulli(w_i) is drawn (its exact law), with the RNG
#         stream isolated so calibration is untouched; plain subgroup mean ->
#         coverage_test_gtau, weighted counterpart -> coverage_test_gtau_wt.
gtau_eval_method <- "weighted"   # "weighted" or "subgroup"
if (!gtau_eval_method %in% c("weighted", "subgroup")) {
    stop("gtau_eval_method must be 'weighted' or 'subgroup'.")
}

# AIPCW controls
localization <- FALSE
censoring_method_AIPCW <- "piecewise"   # "piecewise" currently supported for AIPCW path
pred_use_history <- FALSE
pred_history_mode <- "wide"
G_use_history <- FALSE
G_history_mode <- "wide"
S_use_history <- FALSE
S_history_mode <- "wide"
# Extra RSF args for model.C = "rsf" in dynamicCP_AIPCW_split -> fit_Gk_list.
# Example tuned setting for current DGM:
# rsf_args.C <- list(ntree = 1200L, nodesize = 30L, nsplit = 10L)
rsf_args.C <- list()

# Extra RSF args for model.Xi = "rsf" in dynamicCP_AIPCW_split -> fit_Yk_condmean_at_tk.
# Example:
# rsf_args.Xi <- list(ntree = 1000L, nodesize = 5L, nsplit = 10L)
rsf_args.Xi <- list()

# XGBoost CV/tuning controls (used by model.pred/model.C/model.S/model.Xi when xgb_* is selected)
xgb_nrounds <- 200L
xgb_use_cv <- FALSE
xgb_cv_nfold <- 3L
xgb_cv_early_stopping_rounds <- 20L
xgb_cv_param_grid <- list(
    max_depth = c(2L, 4L),
    min_child_weight = c(1, 5),
    eta = c(0.03, 0.05, 0.1)
)
xgb_cv_seed <- 123L
xgb_cv_verbose <- 0L

# legacy IPCW implementation only supports km/cox censoring models
legacy_ipcw_supported_C <- c("km", "cox")
if (!(model.C %in% legacy_ipcw_supported_C) && ("dCP_IPCW" %in% method.list)) {
    method.list <- setdiff(method.list, "dCP_IPCW")
    message(
        "Skipping dCP_IPCW: model.C='", model.C,
        "' is unsupported by legacy IPCW (supported: ",
        paste(legacy_ipcw_supported_C, collapse = ", "), ")."
    )
}


format_tag_num <- function(x) {
    x_chr <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
    x_chr <- gsub("-", "m", x_chr, fixed = TRUE)
    x_chr <- gsub(".", "p", x_chr, fixed = TRUE)
    x_chr
}
Gtau_delta_tag <- format_tag_num(Gtau_delta)

# Test-data handling for coverage evaluation (see gtau_eval_method above):
#  - "one":                 uncensored test; coverage on the survivor population
#                           {T > tau} via covered_T (coverage_test_true).
#  - "estimated"+weighted:  uncensored test; at-risk coverage on {T>tau, C>tau} by
#                           the oracle weighting identity, weights G_true(tau|path)
#                           -> coverage_test_gtau. coverage_test_true then means
#                           survivor-population coverage (as in "one").
#  - "estimated"+subgroup:  LEGACY — censored test (true DGM law); at-risk coverage
#                           on the X>tau subgroup via covered_T (coverage_test_true).
#  - "tilted" (either):     uncensored test; weights from the exp(delta t)-tilted
#                           TRUE censoring law; weighted -> coverage_test_gtau, or
#                           Bernoulli(w) subgroup mean -> coverage_test_gtau (with
#                           the weighted value in coverage_test_gtau_wt).
# The full covariate path (return_L_full) is retained whenever weights are needed.
needs_gtau_weights <- (Gtau_mode == "tilted") ||
    (Gtau_mode == "estimated" && gtau_eval_method == "weighted")
no_censoring_test <- switch(
    Gtau_mode,
    one = TRUE,
    estimated = (gtau_eval_method == "weighted"),
    tilted = TRUE,
    stop("Unknown Gtau_mode: ", Gtau_mode)
)
if (needs_gtau_weights && dgm_name != "linear_weibull") {
    stop("Gtau evaluation weights are only implemented for dgm_name='linear_weibull' ",
         "(use gtau_eval_method='subgroup' with Gtau_mode='estimated' for other DGMs).")
}

TT.name <- "TT"
train_ratio <- 0.5

# DGM for the current standard simulation setting
change_times <- c(3, 6)
if (dgm_name == "linear_weibull") {
    dgm_parK <- list(
        shape_T = c(4, 4, 5),
        shape_C = c(3, 3, 3),
        beta_T0 = c(-8, -8, -5),
        beta_TL = c(1, 2, 3),
        beta_C0 = c(-6, -6, -5),
        beta_CL = c(2, 2, 2)
    )
} else if (dgm_name == "mixC_uniform_cox") {
    dgm_parK <- list(
        shape_T = c(4, 4, 5),
        shape_C = c(3, 3, 3),
        beta_T0 = c(-8, -8, -5),
        beta_TL = c(1, 2, 3),
        beta_TB = c(0.7, 0.9, 1.1),
        beta_C0 = c(-8.0, -7.5, -7.0),
        beta_CL = c(0.5, 0.6, 0.7),
        beta_CL2 = c(0.4, 0.5, 0.6),
        prob_cox = 0.5,
        uniform_cmax = 12
    )
    
    # Expose the subgroup indicator B to the nuisance and prediction models.
    if (is.null(covname.pred.baseline)) covname.pred.baseline <- c("B")
    if (is.null(covname.C.baseline)) covname.C.baseline <- c("B")
} else {
    stop("Unknown dgm_name: ", dgm_name)
}

tau_grid <- c(0, change_times)

# IPCW baseline implementation in src/dynamicCP.R expects a symmetric full grid.
theta_grid_IPCW <- seq(0, 1, by = 0.01)
# AIPCW implementation in src/dynamicCP_AIPCW.R currently enforces theta in [0, 0.5].
theta_grid_AIPCW <- seq(0, 0.5, by = 0.01)
alpha <- 0.1
trim.C <- 0.05


## -------------------------- Wrappers ---------------------------------
run_one_tau_IPCW <- function(dat_long, dat_test, tau,
                             start.name, stop.name, event.name, event.C.name, id.name,
                             covname.pred.baseline, covname.pred.timevarying,
                             model.pred, model.C, covname.C.baseline, covname.C.timevarying,
                             model.S,
                             theta_grid, alpha, trim.C, TT.name,
                             visit_times,
                             seed, train_ratio) {
    out <- tryCatch({
        dynamicCP_IPCW_split(
            dat = dat_long, dat_test = dat_test,
            start.name = start.name, stop.name = stop.name,
            event.name = event.name, event.C.name = event.C.name, id.name = id.name,
            tau = tau,
            model.pred = model.pred,
            covname.pred.baseline = covname.pred.baseline,
            covname.pred.timevarying = covname.pred.timevarying,
            model.C = model.C,
            covname.C.baseline = covname.C.baseline,
            covname.C.timevarying = covname.C.timevarying,
            TT.name = TT.name,
            seed = seed, train_ratio = train_ratio,
            theta_grid = theta_grid,
            alpha = alpha, trim.C = trim.C
        )
    }, error = function(e) {
        list(error = TRUE, msg = conditionMessage(e))
    })
    out
}

run_one_tau_AIPCW <- function(dat_long, dat_test, tau,
                              start.name, stop.name, event.name, event.C.name, id.name,
                              covname.pred.baseline, covname.pred.timevarying,
                              model.pred, model.C, covname.C.baseline, covname.C.timevarying,
                              model.S, model.Xi,
                              theta_grid, alpha, trim.C, TT.name,
                              visit_times, seed, train_ratio,
                              localization,
                              pred_use_history, pred_history_mode,
                              censoring_method,
                              G_use_history, G_history_mode,
                              S_use_history, S_history_mode,
                              rsf_args.C,
                              rsf_args.Xi,
                              xgb_nrounds,
                              xgb_use_cv,
                              xgb_cv_nfold,
                              xgb_cv_early_stopping_rounds,
                              xgb_cv_param_grid,
                              xgb_cv_seed,
                              xgb_cv_verbose,
                              Gtau_mode, Gtau_delta) {
    out <- tryCatch({
        dynamicCP_AIPCW_split(
            dat = dat_long, dat_test = dat_test,
            start.name = start.name, stop.name = stop.name,
            event.name = event.name, event.C.name = event.C.name, id.name = id.name,
            tau = tau,
            model.pred = model.pred,
            covname.pred.baseline = covname.pred.baseline,
            covname.pred.timevarying = covname.pred.timevarying,
            model.C = model.C,
            covname.C.baseline = covname.C.baseline,
            covname.C.timevarying = covname.C.timevarying,
            rsf_args.C = rsf_args.C,
            rsf_args.Xi = rsf_args.Xi,
            model.S = model.S,
            model.Xi = model.Xi,
            TT.name = TT.name,
            visit_times = visit_times,
            seed = seed, train_ratio = train_ratio,
            theta_grid = theta_grid,
            alpha = alpha, trim.C = trim.C,
            pred_use_history = pred_use_history,
            pred_history_mode = pred_history_mode,
            censoring_method = censoring_method,
            G_use_history = G_use_history,
            G_history_mode = G_history_mode,
            S_use_history = S_use_history,
            S_history_mode = S_history_mode,
            xgb_nrounds = xgb_nrounds,
            xgb_use_cv = xgb_use_cv,
            xgb_cv_nfold = xgb_cv_nfold,
            xgb_cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
            xgb_cv_param_grid = xgb_cv_param_grid,
            xgb_cv_seed = xgb_cv_seed,
            xgb_cv_verbose = xgb_cv_verbose,
            localization = localization,
            Gtau_mode = Gtau_mode,
            Gtau_delta = Gtau_delta
        )
    }, error = function(e) {
        list(error = TRUE, msg = conditionMessage(e))
    })
    out
}


## ------------------------ Replicate seeds -----------------------------
set.seed(123)
seeds <- sample(10^7, R, replace = FALSE)
seeds_test <- sample(10^7, R, replace = FALSE)
seeds_split <- sample(10^7, R, replace = FALSE)
reuse_aipcw_fit <- all(c("dCP_IPCW_new", "dCP_AIPCW") %in% method.list)
aipcw_view_cache <- if (reuse_aipcw_fit) vector("list", R) else NULL


## -------------------------- Simulations -------------------------------
for (method in method.list) {
    cat(paste0(setup, " ", method, " rho=", rho, "\n"))
    method_t0 <- proc.time()[["elapsed"]]
    
    results_summary <- vector("list", R)
    results_bounds <- vector("list", R)
    results_survC <- vector("list", R)
    rep_runtime_sec <- rep(NA_real_, R)
    
    cat("Replications: ")
    for (r in seq_len(R)) {
        cat(r, "..")
        rep_t0 <- proc.time()[["elapsed"]]
        
        dat_long_r <- simulate_dataset_long(
            n = n, seed = seeds[r],
            change_times = change_times,
            tau_max = Inf,
            no_censoring = FALSE,
            par = dgm_parK,
            rho = rho,
            dgm_name = dgm_name
        )
        
        dat_test_r <- simulate_dataset_long(
            n = n_test, seed = seeds_test[r],
            change_times = change_times,
            tau_max = Inf,
            no_censoring = no_censoring_test,
            par = dgm_parK,
            rho = rho,
            dgm_name = dgm_name,
            return_L_full = needs_gtau_weights   # full path needed for the Gtau weights
        )
        L_full_test_r <- if (needs_gtau_weights) attr(dat_test_r, "L_full") else NULL
        
        if (!is.null(TT.name) && !(TT.name %in% colnames(dat_test_r))) {
            stop("TT.name=", TT.name, " not found in dat_test_r.")
        }
        
        tab_r <- data.frame(
            rep = r,
            tau = tau_grid,
            ok = NA,
            msg = NA_character_,
            runtime_sec = NA_real_,
            theta_hat = NA_real_,
            coverage_cal_ipcw = NA_real_,
            coverage_test_ipcw = NA_real_,
            coverage_cal_true = NA_real_,
            coverage_test_true = NA_real_,
            mean_len_cal = NA_real_,
            mean_len_test = NA_real_,
            q25_len_test = NA_real_,
            q50_len_test = NA_real_,
            q75_len_test = NA_real_
        )
        # Extra at-risk-coverage columns only when Gtau weights are in play (keeps
        # the "one" and legacy estimated+subgroup output schemas unchanged).
        if (needs_gtau_weights) tab_r$coverage_test_gtau <- NA_real_
        if (Gtau_mode == "tilted" && gtau_eval_method == "subgroup") {
            tab_r$coverage_test_gtau_wt <- NA_real_
        }

        bounds_r <- vector("list", length(tau_grid))
        names(bounds_r) <- paste0("tau_", tau_grid)
        
        survC_r <- vector("list", length(tau_grid))
        names(survC_r) <- paste0("tau_", tau_grid)
        
        for (k in seq_along(tau_grid)) {
            tau <- tau_grid[k]
            tau_t0 <- proc.time()[["elapsed"]]
            id_tr_use <- NULL
            
            if (method == "dCP_IPCW") {
                fit <- run_one_tau_IPCW(
                    dat_long = dat_long_r, dat_test = dat_test_r, tau = tau,
                    start.name = start.name, stop.name = stop.name,
                    event.name = event.name, event.C.name = event.C.name, id.name = id.name,
                    covname.pred.baseline = covname.pred.baseline,
                    covname.pred.timevarying = covname.pred.timevarying,
                    model.pred = model.pred,
                    model.C = model.C,
                    covname.C.baseline = covname.C.baseline,
                    covname.C.timevarying = covname.C.timevarying,
                    model.S = model.S,
                    theta_grid = theta_grid_IPCW,
                    alpha = alpha,
                    trim.C = trim.C,
                    TT.name = TT.name,
                    visit_times = tau_grid,
                    seed = seeds_split[r],
                    train_ratio = train_ratio
                )
                
                if (isTRUE(fit$error)) {
                    tab_r$ok[k] <- FALSE
                    tab_r$msg[k] <- fit$msg
                    tab_r$runtime_sec[k] <- proc.time()[["elapsed"]] - tau_t0
                    next
                }
                
                lower_cal <- fit$lower_IPCW_cal
                upper_cal <- fit$upper_IPCW_cal
                lower_test <- fit$lower_IPCW_test
                upper_test <- fit$upper_IPCW_test
                cov_cal <- fit$coverage_cal
                cov_test <- fit$coverage_test
                cov_cal_true <- fit$coverage_cal_T
                cov_test_true <- fit$coverage_test_T
                theta_hat_use <- if ("theta_hat" %in% names(fit)) fit$theta_hat else NA_real_
                survC_cal_use <- fit$survC_cal
                survC_test_use <- fit$survC_test
                id_tr_use <- fit$id_tr
                
            } else if (method == "dCP_IPCW_new") {
                fit <- run_one_tau_AIPCW(
                    dat_long = dat_long_r, dat_test = dat_test_r, tau = tau,
                    start.name = start.name, stop.name = stop.name,
                    event.name = event.name, event.C.name = event.C.name, id.name = id.name,
                    covname.pred.baseline = covname.pred.baseline,
                    covname.pred.timevarying = covname.pred.timevarying,
                    model.pred = model.pred,
                    model.C = model.C,
                    covname.C.baseline = covname.C.baseline,
                    covname.C.timevarying = covname.C.timevarying,
                    model.S = model.S,
                    model.Xi = model.Xi,
                    theta_grid = theta_grid_AIPCW,
                    alpha = alpha,
                    trim.C = trim.C,
                    TT.name = TT.name,
                    visit_times = tau_grid,
                    seed = seeds_split[r],
                    train_ratio = train_ratio,
                    localization = localization,
                    pred_use_history = pred_use_history,
                    pred_history_mode = pred_history_mode,
                    censoring_method = censoring_method_AIPCW,
                    G_use_history = G_use_history,
                    G_history_mode = G_history_mode,
                    S_use_history = S_use_history,
                    S_history_mode = S_history_mode,
                    rsf_args.C = rsf_args.C,
                    rsf_args.Xi = rsf_args.Xi,
                    xgb_nrounds = xgb_nrounds,
                    xgb_use_cv = xgb_use_cv,
                    xgb_cv_nfold = xgb_cv_nfold,
                    xgb_cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
                    xgb_cv_param_grid = xgb_cv_param_grid,
                    xgb_cv_seed = xgb_cv_seed,
                    xgb_cv_verbose = xgb_cv_verbose,
                    Gtau_mode = Gtau_mode,
                    Gtau_delta = Gtau_delta
                )
                
                if (isTRUE(fit$error)) {
                    if (isTRUE(reuse_aipcw_fit)) {
                        if (is.null(aipcw_view_cache[[r]])) {
                            aipcw_view_cache[[r]] <- vector("list", length(tau_grid))
                        }
                        aipcw_view_cache[[r]][[k]] <- list(error = TRUE, msg = fit$msg)
                    }
                    tab_r$ok[k] <- FALSE
                    tab_r$msg[k] <- fit$msg
                    tab_r$runtime_sec[k] <- proc.time()[["elapsed"]] - tau_t0
                    next
                }
                
                # Use IPCW outputs from new implementation as sanity-check baseline
                lower_cal <- fit$lower_IPCW_cal
                upper_cal <- fit$upper_IPCW_cal
                lower_test <- fit$lower_IPCW_test
                upper_test <- fit$upper_IPCW_test
                cov_cal <- fit$coverage_IPCW_cal
                cov_test <- fit$coverage_IPCW_test
                cov_cal_true <- fit$coverage_IPCW_cal_T
                cov_test_true <- fit$coverage_IPCW_test_T
                theta_hat_use <- fit$theta_hat_IPCW
                survC_cal_use <- fit$survC_cal
                survC_test_use <- fit$survC_test
                id_tr_use <- fit$id_tr
                
                if (isTRUE(reuse_aipcw_fit)) {
                    if (is.null(aipcw_view_cache[[r]])) {
                        aipcw_view_cache[[r]] <- vector("list", length(tau_grid))
                    }
                    aipcw_view_cache[[r]][[k]] <- list(
                        error = FALSE,
                        id_tr = fit$id_tr,
                        lower_cal = fit$lower_AIPCW_cal,
                        upper_cal = fit$upper_AIPCW_cal,
                        lower_test = fit$lower_AIPCW_test,
                        upper_test = fit$upper_AIPCW_test,
                        cov_cal = fit$coverage_AIPCW_cal,
                        cov_test = fit$coverage_AIPCW_test,
                        cov_cal_true = fit$coverage_AIPCW_cal_T,
                        cov_test_true = fit$coverage_AIPCW_test_T,
                        theta_hat = fit$theta_hat_AIPCW,
                        survC_cal = fit$survC_cal,
                        survC_test = fit$survC_test
                    )
                }
                
            } else if (method == "dCP_AIPCW") {
                cached <- NULL
                if (isTRUE(reuse_aipcw_fit) &&
                    !is.null(aipcw_view_cache[[r]]) &&
                    length(aipcw_view_cache[[r]]) >= k) {
                    cached <- aipcw_view_cache[[r]][[k]]
                }
                
                if (!is.null(cached)) {
                    if (isTRUE(cached$error)) {
                        tab_r$ok[k] <- FALSE
                        tab_r$msg[k] <- cached$msg
                        tab_r$runtime_sec[k] <- proc.time()[["elapsed"]] - tau_t0
                        next
                    }
                    
                    lower_cal <- cached$lower_cal
                    upper_cal <- cached$upper_cal
                    lower_test <- cached$lower_test
                    upper_test <- cached$upper_test
                    cov_cal <- cached$cov_cal
                    cov_test <- cached$cov_test
                    cov_cal_true <- cached$cov_cal_true
                    cov_test_true <- cached$cov_test_true
                    theta_hat_use <- cached$theta_hat
                    survC_cal_use <- cached$survC_cal
                    survC_test_use <- cached$survC_test
                    id_tr_use <- cached$id_tr
                } else {
                    fit <- run_one_tau_AIPCW(
                        dat_long = dat_long_r, dat_test = dat_test_r, tau = tau,
                        start.name = start.name, stop.name = stop.name,
                        event.name = event.name, event.C.name = event.C.name, id.name = id.name,
                        covname.pred.baseline = covname.pred.baseline,
                        covname.pred.timevarying = covname.pred.timevarying,
                        model.pred = model.pred,
                        model.C = model.C,
                        covname.C.baseline = covname.C.baseline,
                        covname.C.timevarying = covname.C.timevarying,
                        model.S = model.S,
                        model.Xi = model.Xi,
                        theta_grid = theta_grid_AIPCW,
                        alpha = alpha,
                        trim.C = trim.C,
                        TT.name = TT.name,
                        visit_times = tau_grid,
                        seed = seeds_split[r],
                        train_ratio = train_ratio,
                        localization = localization,
                        pred_use_history = pred_use_history,
                        pred_history_mode = pred_history_mode,
                        censoring_method = censoring_method_AIPCW,
                        G_use_history = G_use_history,
                        G_history_mode = G_history_mode,
                        S_use_history = S_use_history,
                        S_history_mode = S_history_mode,
                        rsf_args.C = rsf_args.C,
                        rsf_args.Xi = rsf_args.Xi,
                        xgb_nrounds = xgb_nrounds,
                        xgb_use_cv = xgb_use_cv,
                        xgb_cv_nfold = xgb_cv_nfold,
                        xgb_cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
                        xgb_cv_param_grid = xgb_cv_param_grid,
                        xgb_cv_seed = xgb_cv_seed,
                        xgb_cv_verbose = xgb_cv_verbose,
                        Gtau_mode = Gtau_mode,
                        Gtau_delta = Gtau_delta
                    )
                    
                    if (isTRUE(fit$error)) {
                        tab_r$ok[k] <- FALSE
                        tab_r$msg[k] <- fit$msg
                        tab_r$runtime_sec[k] <- proc.time()[["elapsed"]] - tau_t0
                        next
                    }
                    
                    lower_cal <- fit$lower_AIPCW_cal
                    upper_cal <- fit$upper_AIPCW_cal
                    lower_test <- fit$lower_AIPCW_test
                    upper_test <- fit$upper_AIPCW_test
                    cov_cal <- fit$coverage_AIPCW_cal
                    cov_test <- fit$coverage_AIPCW_test
                    cov_cal_true <- fit$coverage_AIPCW_cal_T
                    cov_test_true <- fit$coverage_AIPCW_test_T
                    theta_hat_use <- fit$theta_hat_AIPCW
                    survC_cal_use <- fit$survC_cal
                    survC_test_use <- fit$survC_test
                    id_tr_use <- fit$id_tr
                }
                
            } else {
                stop("Unknown method: ", method)
            }
            
            if (is.null(id_tr_use)) {
                stop("Internal error: id_tr_use was not set for method=", method, ", rep=", r, ", tau=", tau)
            }
            
            tab_r$ok[k] <- TRUE
            tab_r$theta_hat[k] <- theta_hat_use
            tab_r$coverage_cal_ipcw[k] <- cov_cal
            tab_r$coverage_test_ipcw[k] <- cov_test
            tab_r$coverage_cal_true[k] <- cov_cal_true
            tab_r$coverage_test_true[k] <- cov_test_true
            
            len_cal <- upper_cal - lower_cal
            len_test <- upper_test - lower_test
            tab_r$mean_len_cal[k] <- mean(len_cal, na.rm = TRUE)
            tab_r$mean_len_test[k] <- mean(len_test, na.rm = TRUE)
            tab_r$q25_len_test[k] <- as.numeric(quantile(len_test, 0.25, na.rm = TRUE))
            tab_r$q50_len_test[k] <- as.numeric(quantile(len_test, 0.50, na.rm = TRUE))
            tab_r$q75_len_test[k] <- as.numeric(quantile(len_test, 0.75, na.rm = TRUE))
            
            dat_test_tau_tmp <- prepare_data_tau(dat_test_r, start.name, stop.name, event.name, id.name, tau)
            dat_cal_tau_tmp <- prepare_data_tau(
                dat_long_r[!(dat_long_r[[id.name]] %in% id_tr_use), ],
                start.name, stop.name, event.name, id.name, tau
            )
            
            w_test_tmp <- dat_test_tau_tmp[[event.name]] / pmax(survC_test_use, trim.C)
            w_cal_tmp <- dat_cal_tau_tmp[[event.name]] / pmax(survC_cal_use, trim.C)
            
            covered_X_test <- as.integer(lower_test <= dat_test_tau_tmp[[stop.name]] &
                                             dat_test_tau_tmp[[stop.name]] <= upper_test)
            covered_T_test <- if (!is.null(TT.name) && TT.name %in% names(dat_test_tau_tmp)) {
                as.integer(lower_test <= dat_test_tau_tmp[[TT.name]] &
                               dat_test_tau_tmp[[TT.name]] <= upper_test)
            } else NA_integer_
            
            covered_X_cal <- as.integer(lower_cal <= dat_cal_tau_tmp[[stop.name]] &
                                            dat_cal_tau_tmp[[stop.name]] <= upper_cal)
            covered_T_cal <- if (!is.null(TT.name) && TT.name %in% names(dat_cal_tau_tmp)) {
                as.integer(lower_cal <= dat_cal_tau_tmp[[TT.name]] &
                               dat_cal_tau_tmp[[TT.name]] <= upper_cal)
            } else NA_integer_

            # At-risk estimand {T>tau, C_tilde>tau} evaluation on the uncensored test
            # data (survivors T>tau), under the TRUE censoring law (delta = 0 for
            # "estimated") or its exp(Gtau_delta * t)-tilt ("tilted"):
            #   gtau_eval_method = "weighted": oracle weighting identity.
            #   gtau_eval_method = "subgroup" (tilted only reaches here): draw the
            #     membership indicator I(C_tilde>tau) ~ Bernoulli(w_i) — its exact
            #     law — and take the plain subgroup mean; the weighted counterpart
            #     is recorded alongside as a per-run noise gauge.
            w_gtau_test <- NULL
            at_risk_gtau_test <- NULL
            if (needs_gtau_weights) {
                surv_idx <- which(dat_test_tau_tmp[[TT.name]] > tau)
                w_gtau_test <- rep(NA_real_, nrow(dat_test_tau_tmp))
                if (length(surv_idx) > 0L) {
                    ids_s <- as.character(dat_test_tau_tmp[[id.name]][surv_idx])
                    delta_eval_use <- if (Gtau_mode == "tilted") Gtau_delta else 0
                    w_gtau_test[surv_idx] <- compute_Gtau_true_tilted_linWB1(
                        L_full_test_r[ids_s, , drop = FALSE], tau, delta_eval_use,
                        change_times, dgm_parK)
                    cov_gtau_wt <- weighted_coverage_summary(
                        covered_T_test[surv_idx],
                        (upper_test - lower_test)[surv_idx],
                        w_gtau_test[surv_idx])
                    if (gtau_eval_method == "weighted") {
                        tab_r$coverage_test_gtau[k] <- cov_gtau_wt$coverage
                    } else {
                        # Isolate the Bernoulli draws on a dedicated, deterministic
                        # RNG stream so the calibration/generation streams are
                        # untouched (results identical to a "weighted" run except
                        # for the evaluation columns).
                        if (!exists(".Random.seed", envir = .GlobalEnv)) set.seed(1L)
                        rng_state <- get(".Random.seed", envir = .GlobalEnv)
                        set.seed(seeds_split[r] + 7907L * k + 424243L)
                        p_at_risk <- pmin(pmax(w_gtau_test[surv_idx], 0), 1)
                        at_risk_gtau_test <- rep(NA_integer_, nrow(dat_test_tau_tmp))
                        at_risk_gtau_test[surv_idx] <- rbinom(length(surv_idx), 1L, p_at_risk)
                        assign(".Random.seed", rng_state, envir = .GlobalEnv)

                        in_sub <- surv_idx[at_risk_gtau_test[surv_idx] == 1L]
                        tab_r$coverage_test_gtau[k] <- if (length(in_sub) > 0L) {
                            mean(covered_T_test[in_sub], na.rm = TRUE)
                        } else NA_real_
                        tab_r$coverage_test_gtau_wt[k] <- cov_gtau_wt$coverage
                    }
                }
            }

            test_bounds_df <- data.frame(
                rep = r, tau = tau, id = dat_test_tau_tmp[[id.name]],
                lower = lower_test, upper = upper_test,
                X = dat_test_tau_tmp[[stop.name]],
                T_true = if (!is.null(TT.name) && TT.name %in% names(dat_test_tau_tmp)) dat_test_tau_tmp[[TT.name]] else NA_real_,
                w_ipcw = w_test_tmp,
                covered_X = covered_X_test,
                covered_T = covered_T_test
            )
            # Gtau-weighting extra columns (keep other modes' bounds schema unchanged)
            if (needs_gtau_weights) test_bounds_df$w_gtau <- w_gtau_test
            if (!is.null(at_risk_gtau_test)) test_bounds_df$at_risk_gtau <- at_risk_gtau_test
            bounds_r[[k]] <- list(
                test = test_bounds_df,
                cal = data.frame(
                    rep = r, tau = tau, id = dat_cal_tau_tmp[[id.name]],
                    lower = lower_cal, upper = upper_cal,
                    X = dat_cal_tau_tmp[[stop.name]],
                    T_true = if (!is.null(TT.name) && TT.name %in% names(dat_cal_tau_tmp)) dat_cal_tau_tmp[[TT.name]] else NA_real_,
                    w_ipcw = w_cal_tmp,
                    covered_X = covered_X_cal,
                    covered_T = covered_T_cal
                )
            )
            
            survC_r[[k]] <- list(survC_cal = survC_cal_use, survC_test = survC_test_use)
            tab_r$runtime_sec[k] <- proc.time()[["elapsed"]] - tau_t0
        }
        
        rep_runtime_sec[r] <- proc.time()[["elapsed"]] - rep_t0
        results_summary[[r]] <- tab_r
        results_bounds[[r]] <- bounds_r
        results_survC[[r]] <- survC_r
    }
    cat("\n")
    method_runtime_sec <- proc.time()[["elapsed"]] - method_t0
    cat("Elapsed (sec) for ", method, ": ", round(method_runtime_sec, 2), "\n", sep = "")
    
    res <- do.call(rbind, results_summary)
    print(summary(res$ok))
    res_ok <- res[res$ok %in% TRUE, ]
    
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_prefix <- paste0(
        "results/", setup, "/", method,
        "_", model.pred, "_", model.C,
        "_S", model.S, "_Xi", model.Xi,
        "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
        "_alpha", alpha,
        "_loc", as.integer(localization),
        "_Gtau", Gtau_mode, "_d", Gtau_delta_tag
    )
    out_rds <- paste0(out_prefix, ".rds")
    if (file.exists(out_rds)) {
        warning("Output file already exists and will be overwritten: ", out_rds)
    }
    
    # gtau_eval_method is persisted only when it is meaningful (non-"one" modes),
    # so Gtau_mode="one" result files keep their exact legacy config schema
    # (bit-identical outputs). It is appended right after save_obj below.
    save_obj <- list(
        config = list(
            setup = setup,
            dgm_name = dgm_name,
            method = method,
            source_script = "main_simu_dynamic_linWB1_cox_multiclass_loc0.R",
            change_times = change_times,
            dgm_parK = dgm_parK,
            rho = rho,
            n = n, n_test = n_test, R = R,
            start.name = start.name, stop.name = stop.name,
            event.name = event.name, event.C.name = event.C.name, id.name = id.name,
            covname.pred.baseline = covname.pred.baseline,
            covname.pred.timevarying = covname.pred.timevarying,
            covname.C.baseline = covname.C.baseline,
            covname.C.timevarying = covname.C.timevarying,
            model.pred = model.pred, model.C = model.C, model.S = model.S, model.Xi = model.Xi,
            tau_grid = tau_grid,
            alpha = alpha, trim.C = trim.C,
            train_ratio = train_ratio, TT.name = TT.name,
            theta_grid_IPCW = theta_grid_IPCW,
            theta_grid_AIPCW = theta_grid_AIPCW,
            localization = localization,
            censoring_method_AIPCW = censoring_method_AIPCW,
            pred_use_history = pred_use_history,
            pred_history_mode = pred_history_mode,
            G_use_history = G_use_history,
            G_history_mode = G_history_mode,
            S_use_history = S_use_history,
            S_history_mode = S_history_mode,
            rsf_args.C = rsf_args.C,
            rsf_args.Xi = rsf_args.Xi,
            xgb_nrounds = xgb_nrounds,
            xgb_use_cv = xgb_use_cv,
            xgb_cv_nfold = xgb_cv_nfold,
            xgb_cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
            xgb_cv_param_grid = xgb_cv_param_grid,
            xgb_cv_seed = xgb_cv_seed,
            xgb_cv_verbose = xgb_cv_verbose,
            Gtau_mode = Gtau_mode,
            Gtau_delta = Gtau_delta
        ),
        seeds = list(
            seeds = seeds,
            seeds_test = seeds_test,
            seeds_split = seeds_split
        ),
        results = list(
            results_summary = results_summary,
            results_bounds = results_bounds,
            results_survC = results_survC,
            res = res,
            res_ok = res_ok
        ),
        meta = list(
            timestamp = timestamp,
            runtime = list(
                method_runtime_sec = method_runtime_sec,
                rep_runtime_sec = rep_runtime_sec
            )
        )
    )
    if (Gtau_mode != "one") {
        save_obj$config$gtau_eval_method <- gtau_eval_method
    }

    saveRDS(save_obj, file = out_rds)
    cat("Saved simulation bundle to:\n  ", out_rds, "\n", sep = "")
}




### Quick table in R console

library(dplyr)

setup_q <- setup
methods_q <- method.list

files_all <- list.files(file.path("results", setup_q), full.names = TRUE)
files_all <- files_all[endsWith(files_all, ".rds")]
if (length(files_all) == 0) {
    stop("No .rds files found under results/", setup_q)
}

config_matches_current <- function(cfg) {
    isTRUE(all.equal(cfg$setup, setup, check.attributes = FALSE)) &&
        isTRUE(all.equal(if (!is.null(cfg$dgm_name)) cfg$dgm_name else "linear_weibull", dgm_name, check.attributes = FALSE)) &&
        isTRUE(all.equal(cfg$model.pred, model.pred, check.attributes = FALSE)) &&
        isTRUE(all.equal(cfg$model.C, model.C, check.attributes = FALSE)) &&
        isTRUE(all.equal(cfg$model.S, model.S, check.attributes = FALSE)) &&
        isTRUE(all.equal(cfg$model.Xi, model.Xi, check.attributes = FALSE)) &&
        isTRUE(all.equal(if (!is.null(cfg$rsf_args.C)) cfg$rsf_args.C else list(), rsf_args.C, check.attributes = FALSE)) &&
        isTRUE(all.equal(if (!is.null(cfg$rsf_args.Xi)) cfg$rsf_args.Xi else list(), rsf_args.Xi, check.attributes = FALSE)) &&
        isTRUE(all.equal(if (!is.null(cfg$xgb_nrounds)) cfg$xgb_nrounds else 200L, xgb_nrounds, check.attributes = FALSE)) &&
        isTRUE(all.equal(if (!is.null(cfg$xgb_use_cv)) cfg$xgb_use_cv else FALSE, xgb_use_cv, check.attributes = FALSE)) &&
        isTRUE(all.equal(if (!is.null(cfg$xgb_cv_nfold)) cfg$xgb_cv_nfold else 5L, xgb_cv_nfold, check.attributes = FALSE)) &&
        isTRUE(all.equal(if (!is.null(cfg$xgb_cv_early_stopping_rounds)) cfg$xgb_cv_early_stopping_rounds else 20L, xgb_cv_early_stopping_rounds, check.attributes = FALSE)) &&
        isTRUE(all.equal(if (!is.null(cfg$xgb_cv_param_grid)) cfg$xgb_cv_param_grid else NULL, xgb_cv_param_grid, check.attributes = FALSE)) &&
        isTRUE(all.equal(if (!is.null(cfg$xgb_cv_seed)) cfg$xgb_cv_seed else NULL, xgb_cv_seed, check.attributes = FALSE)) &&
        isTRUE(all.equal(if (!is.null(cfg$xgb_cv_verbose)) cfg$xgb_cv_verbose else 0L, xgb_cv_verbose, check.attributes = FALSE)) &&
        isTRUE(all.equal(cfg$rho, rho, check.attributes = FALSE)) &&
        isTRUE(all.equal(cfg$n, n, check.attributes = FALSE)) &&
        isTRUE(all.equal(cfg$n_test, n_test, check.attributes = FALSE)) &&
        isTRUE(all.equal(cfg$alpha, alpha, check.attributes = FALSE)) &&
        isTRUE(all.equal(cfg$localization, localization, check.attributes = FALSE)) &&
        isTRUE(all.equal(cfg$Gtau_mode, Gtau_mode, check.attributes = FALSE)) &&
        isTRUE(all.equal(cfg$Gtau_delta, Gtau_delta, check.attributes = FALSE)) &&
        # Evaluation method only matters when Gtau weighting is in play; legacy
        # files (no gtau_eval_method) were "subgroup" for estimated, "weighted"
        # for tilted.
        (Gtau_mode == "one" ||
             isTRUE(all.equal(
                 if (!is.null(cfg$gtau_eval_method)) cfg$gtau_eval_method
                 else if (identical(cfg$Gtau_mode, "estimated")) "subgroup" else "weighted",
                 gtau_eval_method, check.attributes = FALSE)))
}

meta <- do.call(rbind, lapply(files_all, function(f) {
    obj_tmp <- readRDS(f)
    cfg <- obj_tmp$config
    data.frame(
        file = f,
        mtime = file.info(f)$mtime,
        method = if (!is.null(cfg$method)) as.character(cfg$method) else NA_character_,
        config_match = if (!is.null(cfg)) config_matches_current(cfg) else FALSE,
        stringsAsFactors = FALSE
    )
}))
meta <- meta[meta$method %in% methods_q & meta$config_match, , drop = FALSE]
if (nrow(meta) == 0) {
    stop("No .rds files found matching current config under results/", setup_q)
}

pick_latest_for_method <- function(m) {
    cand <- meta[meta$method == m, , drop = FALSE]
    if (nrow(cand) == 0) return(NA_character_)
    cand$file[which.max(cand$mtime)]
}

quick_objs <- setNames(vector("list", length(methods_q)), methods_q)
for (m in methods_q) {
    f <- pick_latest_for_method(m)
    if (is.na(f)) {
        warning("No result file found for current config with config$method == ", m)
        next
    }
    quick_objs[[m]] <- readRDS(f)
    cat(m, " -> ", basename(f), "\n", sep = "")
}
if (all(vapply(quick_objs, is.null, logical(1)))) {
    stop("No matching method results were loaded.")
}

summary_all <- bind_rows(lapply(names(quick_objs), function(m) {
    obj <- quick_objs[[m]]
    if (is.null(obj)) return(NULL)
    obj$results$res_ok %>%
        group_by(tau) %>%
        summarize(
            cov_true = mean(coverage_test_true, na.rm = TRUE),
            cov_ipcw = mean(coverage_test_ipcw, na.rm = TRUE),
            mean_len = mean(mean_len_test, na.rm = TRUE),
            theta = mean(theta_hat, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(method = m)
}))

print(summary_all, n = Inf)
