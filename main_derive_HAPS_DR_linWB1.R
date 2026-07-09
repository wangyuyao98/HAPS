rm(list = ls())

# Run this script from the repository root.
if (!file.exists("src/gen_ICML_simu.R")) {
    stop("Please run main_derive_HAPS_DR_linWB1.R from the repository root directory.")
}

source("src/gen_ICML_simu.R")
source("src/helpers.R")
source("src/helpers_AIPCW.R")
source("src/dynamicCP_AIPCW.R")
source("src/gtau_eval_helpers.R")   # weighted_coverage_summary for Gtau_mode="tilted"

library(survival)
library(randomForestSRC)

`%||%` <- function(x, y) if (is.null(x)) y else x

setup <- "linWB1"
out_dir <- file.path("results", setup)
old_dynamic_dir <- file.path("results", "ICML_simu1")
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

rho <- 0.3
n <- 1000L
n_test <- 500L
R <- 200L
alpha <- 0.1

valid_targets <- c("cox_cox_multiclass_loc0", "rsf_rsf_reg_loc0", "rsf_rsf_reg_loc1")
args <- commandArgs(trailingOnly = TRUE)
target_args <- intersect(args, valid_targets)
# NOTE: do not use setdiff() here — it deduplicates, so repeated numeric values
# (e.g. n == n_test, as in "... 4 400 400") would silently collapse and leave
# later parameters at their defaults.
numeric_args <- args[!(args %in% valid_targets)]
bad_numeric_args <- numeric_args[is.na(suppressWarnings(as.integer(numeric_args)))]
if (length(bad_numeric_args) > 0L) {
    stop(
        "Unknown argument(s): ", paste(bad_numeric_args, collapse = ", "),
        ". Targets: ", paste(valid_targets, collapse = ", "),
        ". Optional numeric arguments: R n n_test"
    )
}
if (length(numeric_args) >= 1L && nzchar(numeric_args[[1L]])) R <- as.integer(numeric_args[[1L]])
if (length(numeric_args) >= 2L && nzchar(numeric_args[[2L]])) n <- as.integer(numeric_args[[2L]])
if (length(numeric_args) >= 3L && nzchar(numeric_args[[3L]])) n_test <- as.integer(numeric_args[[3L]])
if (!is.finite(R) || R < 1L) stop("R must be a positive integer.")
if (!is.finite(n) || n < 10L) stop("n must be at least 10.")
if (!is.finite(n_test) || n_test < 10L) stop("n_test must be at least 10.")

format_dynamic_file <- function(method, model.pred, model.C, model.S, model.Xi, loc) {
    file.path(
        out_dir,
        paste0(
            method, "_", model.pred, "_", model.C,
            "_S", model.S, "_Xi", model.Xi,
            "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
            "_alpha", alpha, "_loc", loc, "_Gtauone_d0.rds"
        )
    )
}

format_dr_file <- function(tag) {
    file.path(
        out_dir,
        paste0(
            "derived_DynaCoPS_DR_", tag,
            "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
            "_alpha", alpha, ".rds"
        )
    )
}

targets <- list(
    cox_cox_multiclass_loc0 = list(
        tag = "cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone",
        base_file = format_dynamic_file("dCP_IPCW", "cox", "cox", "cox", "xgb_multiclass", 0),
        out_file = format_dr_file("cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone")
    ),
    rsf_rsf_reg_loc0 = list(
        tag = "rsf_rsf_Srsf_Xixgb_reg_loc0_Gtauone",
        base_file = file.path(
            old_dynamic_dir,
            "dCP_IPCW_new_rsf_rsf_Srsf_Xixgb_reg_rho0.3_R200_n1000_ntest500_alpha0.1_loc0_Gtauone_d0.rds"
        ),
        out_file = file.path(
            out_dir,
            "derived_DynaCoPS_DR_rsf_rsf_Srsf_Xixgb_reg_loc0_Gtauone.rds"
        )
    ),
    rsf_rsf_reg_loc1 = list(
        tag = "rsf_rsf_Srsf_Xixgb_reg_loc1_Gtauone",
        base_file = file.path(
            old_dynamic_dir,
            "dCP_IPCW_new_rsf_rsf_Srsf_Xixgb_reg_rho0.3_R200_n1000_ntest500_alpha0.1_loc1_Gtauone_d0.rds"
        ),
        out_file = file.path(
            out_dir,
            "derived_DynaCoPS_DR_rsf_rsf_Srsf_Xixgb_reg_loc1_Gtauone.rds"
        )
    )
)

if (length(target_args) > 0L) {
    targets <- targets[target_args]
} else {
    targets <- targets["cox_cox_multiclass_loc0"]
}

cat("DynaCoPS-DR targets: ", paste(names(targets), collapse = ", "), "\n", sep = "")

make_dr_bounds <- function(base_bounds, fit, which_set = c("test", "cal")) {
    which_set <- match.arg(which_set)
    b <- base_bounds[[which_set]]
    if (!is.data.frame(b) || nrow(b) == 0L) {
        return(data.frame())
    }

    if (identical(which_set, "test")) {
        ids_fit <- fit$id_test_tau
        lower_model <- fit$lower_model_test
        upper_model <- fit$upper_model_test
    } else {
        ids_fit <- fit$id_cal_tau
        lower_model <- fit$lower_model_cal
        upper_model <- fit$upper_model_cal
    }

    idx <- match(as.character(b$id), as.character(ids_fit))
    if (anyNA(idx)) {
        stop("Failed to align ", which_set, " ids to reconstructed model quantiles.")
    }

    b$lower_DynaCoPS <- b$lower
    b$upper_DynaCoPS <- b$upper
    b$lower_model <- lower_model[idx]
    b$upper_model <- upper_model[idx]
    b$lower <- pmin(b$lower_DynaCoPS, b$lower_model)
    b$upper <- pmax(b$upper_DynaCoPS, b$upper_model)
    b$covered_X <- as.integer(b$lower <= b$X & b$X <= b$upper)
    b$covered_T <- if ("T_true" %in% names(b)) {
        as.integer(b$lower <= b$T_true & b$T_true <= b$upper)
    } else {
        NA_integer_
    }
    b
}

summarize_dr_bounds <- function(bt, bc) {
    len_test <- bt$upper - bt$lower
    len_cal <- bc$upper - bc$lower
    out <- data.frame(
        coverage_cal_ipcw = mean(bc$w_ipcw * bc$covered_X, na.rm = TRUE) / mean(bc$w_ipcw, na.rm = TRUE),
        coverage_test_ipcw = mean(bt$w_ipcw * bt$covered_X, na.rm = TRUE) / mean(bt$w_ipcw, na.rm = TRUE),
        coverage_cal_true = if (all(is.na(bc$covered_T))) NA_real_ else mean(bc$covered_T, na.rm = TRUE),
        coverage_test_true = if (all(is.na(bt$covered_T))) NA_real_ else mean(bt$covered_T, na.rm = TRUE),
        mean_len_cal = mean(len_cal, na.rm = TRUE),
        mean_len_test = mean(len_test, na.rm = TRUE),
        q25_len_test = as.numeric(stats::quantile(len_test, 0.25, na.rm = TRUE)),
        q50_len_test = as.numeric(stats::quantile(len_test, 0.50, na.rm = TRUE)),
        q75_len_test = as.numeric(stats::quantile(len_test, 0.75, na.rm = TRUE))
    )
    # At-risk coverage on {T>tau, C_tilde>tau} using the Gtau evaluation carried in
    # the base file's test bounds: subgroup indicator (at_risk_gtau, from a
    # tilted+subgroup base run — SAME drawn population, no re-drawing) takes
    # precedence for coverage_test_gtau, with the weighted value reported alongside;
    # otherwise the weighted estimator on w_gtau.
    if ("w_gtau" %in% names(bt)) {
        surv <- is.finite(bt$w_gtau)
        wc <- if (any(surv)) {
            weighted_coverage_summary(bt$covered_T[surv], len_test[surv], bt$w_gtau[surv])$coverage
        } else NA_real_
        if ("at_risk_gtau" %in% names(bt)) {
            in_sub <- which(bt$at_risk_gtau %in% 1L)
            out$coverage_test_gtau <- if (length(in_sub) > 0L) {
                mean(bt$covered_T[in_sub], na.rm = TRUE)
            } else NA_real_
            out$coverage_test_gtau_wt <- wc
        } else {
            out$coverage_test_gtau <- wc
        }
    }
    out
}

derive_one_target <- function(target) {
    if (!file.exists(target$base_file)) {
        stop("Missing base DynaCoPS file: ", target$base_file)
    }

    base_obj <- readRDS(target$base_file)
    cfg <- base_obj$config
    seeds <- base_obj$seeds

    R <- cfg$R
    tau_grid <- cfg$tau_grid
    dgm_name <- cfg$dgm_name %||% "linear_weibull"
    # Test-data handling must mirror the base driver's Gtau evaluation setup:
    # censored test data ONLY for the legacy estimated+subgroup configuration;
    # everything else ("one", tilted, estimated+weighted) uses uncensored test
    # data, with at-risk coverage taken from the w_gtau / at_risk_gtau columns
    # carried in the base file's bounds. Legacy base files (no gtau_eval_method
    # in config) were "subgroup" for estimated.
    gtau_mode_cfg <- cfg$Gtau_mode %||% "one"
    eval_method_cfg <- cfg$gtau_eval_method %||%
        (if (identical(gtau_mode_cfg, "estimated")) "subgroup" else "weighted")
    no_censoring_test <- !(gtau_mode_cfg == "estimated" && eval_method_cfg == "subgroup")

    results_summary <- vector("list", R)
    results_bounds <- vector("list", R)
    results_survC <- vector("list", R)
    rep_runtime_sec <- rep(NA_real_, R)

    cat("Deriving DynaCoPS-DR for ", target$tag, "\n", sep = "")
    cat("Replications: ")
    method_t0 <- proc.time()[["elapsed"]]

    for (r in seq_len(R)) {
        cat(r, "..")
        rep_t0 <- proc.time()[["elapsed"]]

        dat_long_r <- simulate_dataset_long(
            n = cfg$n,
            seed = seeds$seeds[r],
            change_times = cfg$change_times,
            tau_max = Inf,
            no_censoring = FALSE,
            par = cfg$dgm_parK,
            rho = cfg$rho,
            dgm_name = dgm_name
        )
        dat_test_r <- simulate_dataset_long(
            n = cfg$n_test,
            seed = seeds$seeds_test[r],
            change_times = cfg$change_times,
            tau_max = Inf,
            no_censoring = no_censoring_test,
            par = cfg$dgm_parK,
            rho = cfg$rho,
            dgm_name = dgm_name
        )

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
        needs_gtau_cols <- gtau_mode_cfg == "tilted" ||
            (gtau_mode_cfg == "estimated" && eval_method_cfg == "weighted")
        if (needs_gtau_cols) tab_r$coverage_test_gtau <- NA_real_
        if (gtau_mode_cfg == "tilted" && eval_method_cfg == "subgroup") {
            tab_r$coverage_test_gtau_wt <- NA_real_
        }
        bounds_r <- vector("list", length(tau_grid))
        names(bounds_r) <- paste0("tau_", tau_grid)
        survC_r <- vector("list", length(tau_grid))
        names(survC_r) <- paste0("tau_", tau_grid)

        for (k in seq_along(tau_grid)) {
            tau <- tau_grid[k]
            tau_t0 <- proc.time()[["elapsed"]]

            fit <- tryCatch({
                dynamicCP_AIPCW_split(
                    dat = dat_long_r,
                    dat_test = dat_test_r,
                    start.name = cfg$start.name,
                    stop.name = cfg$stop.name,
                    event.name = cfg$event.name,
                    event.C.name = cfg$event.C.name,
                    id.name = cfg$id.name,
                    tau = tau,
                    model.pred = cfg$model.pred,
                    covname.pred.baseline = cfg$covname.pred.baseline,
                    covname.pred.timevarying = cfg$covname.pred.timevarying,
                    model.C = cfg$model.C,
                    covname.C.baseline = cfg$covname.C.baseline,
                    covname.C.timevarying = cfg$covname.C.timevarying,
                    rsf_args.C = cfg$rsf_args.C %||% list(),
                    model.S = cfg$model.S %||% "cox",
                    model.Xi = cfg$model.Xi %||% "lm",
                    TT.name = cfg$TT.name,
                    visit_times = tau_grid,
                    seed = seeds$seeds_split[r],
                    train_ratio = cfg$train_ratio,
                    theta_grid = cfg$theta_grid_AIPCW %||% seq(0, 0.5, by = 0.01),
                    alpha = cfg$alpha,
                    trim.C = cfg$trim.C,
                    pred_use_history = cfg$pred_use_history %||% FALSE,
                    pred_history_mode = cfg$pred_history_mode %||% "wide",
                    censoring_method = cfg$censoring_method_AIPCW %||% "piecewise",
                    G_use_history = cfg$G_use_history %||% FALSE,
                    G_history_mode = cfg$G_history_mode %||% "wide",
                    # This ipcw_only call is used ONLY to reconstruct the raw model
                    # quantiles (alpha/2, 1-alpha/2), which do not depend on Gtau;
                    # the Gtau-dependent intervals/weights come from the base file's
                    # bounds. ipcw_only also supports only Gtau_mode="one", so pin it.
                    Gtau_mode = "one",
                    Gtau_delta = 0,
                    ipcw_only = TRUE
                )
            }, error = function(e) {
                list(error = TRUE, msg = conditionMessage(e))
            })

            if (isTRUE(fit$error)) {
                tab_r$ok[k] <- FALSE
                tab_r$msg[k] <- fit$msg
                bounds_r[[k]] <- list(test = data.frame(), cal = data.frame())
                tab_r$runtime_sec[k] <- proc.time()[["elapsed"]] - tau_t0
                next
            }

            base_bounds_k <- base_obj$results$results_bounds[[r]][[k]]
            out_k <- tryCatch({
                bt <- make_dr_bounds(base_bounds_k, fit, which_set = "test")
                bc <- make_dr_bounds(base_bounds_k, fit, which_set = "cal")
                ss <- summarize_dr_bounds(bt, bc)
                list(ok = TRUE, test = bt, cal = bc, summary = ss)
            }, error = function(e) {
                list(ok = FALSE, msg = conditionMessage(e))
            })

            if (isTRUE(out_k$ok)) {
                tab_r$ok[k] <- TRUE
                if ("theta_hat" %in% names(base_obj$results$res)) {
                    base_row <- base_obj$results$res[base_obj$results$res$rep == r & base_obj$results$res$tau == tau, , drop = FALSE]
                    if (nrow(base_row) > 0L) tab_r$theta_hat[k] <- base_row$theta_hat[1]
                }
                tab_r[k, names(out_k$summary)] <- out_k$summary
                bounds_r[[k]] <- list(test = out_k$test, cal = out_k$cal)
                survC_r[[k]] <- list(
                    survC_cal = out_k$cal$survC %||% NA_real_,
                    survC_test = out_k$test$survC %||% NA_real_
                )
            } else {
                tab_r$ok[k] <- FALSE
                tab_r$msg[k] <- out_k$msg
                bounds_r[[k]] <- list(test = data.frame(), cal = data.frame())
            }
            tab_r$runtime_sec[k] <- proc.time()[["elapsed"]] - tau_t0
        }

        results_summary[[r]] <- tab_r
        results_bounds[[r]] <- bounds_r
        results_survC[[r]] <- survC_r
        rep_runtime_sec[r] <- proc.time()[["elapsed"]] - rep_t0
    }
    cat("\n")

    res <- do.call(rbind, results_summary)
    res_ok <- res[res$ok %in% TRUE, , drop = FALSE]
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    cfg_out <- cfg
    cfg_out$setup <- setup
    cfg_out$method <- "derived_DynaCoPS_DR"
    cfg_out$source_script <- "main_derive_HAPS_DR_linWB1.R"
    cfg_out$base_file <- target$base_file
    cfg_out$derivation <- "lower_DR=pmin(lower_DynaCoPS,q_alpha/2); upper_DR=pmax(upper_DynaCoPS,q_1-alpha/2)"

    save_obj <- list(
        config = cfg_out,
        seeds = seeds,
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
                method_runtime_sec = proc.time()[["elapsed"]] - method_t0,
                rep_runtime_sec = rep_runtime_sec
            )
        )
    )

    saveRDS(save_obj, target$out_file)
    cat("Saved derived DynaCoPS-DR bundle to:\n  ", target$out_file, "\n", sep = "")

    if (requireNamespace("dplyr", quietly = TRUE)) {
        print(
            dplyr::summarize(
                dplyr::group_by(res_ok, tau),
                cov_true = mean(coverage_test_true, na.rm = TRUE),
                cov_ipcw = mean(coverage_test_ipcw, na.rm = TRUE),
                mean_len = mean(mean_len_test, na.rm = TRUE),
                .groups = "drop"
            ),
            n = Inf
        )
    }
    invisible(target$out_file)
}

for (nm in names(targets)) {
    derive_one_target(targets[[nm]])
}
