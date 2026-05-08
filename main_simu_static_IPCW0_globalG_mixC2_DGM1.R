rm(list = ls())

# Run this script from the repository root.
if (!file.exists("src/gen_ICML_simu.R")) {
    stop("Please run main_simu_static_IPCW0_globalG_mixC2_DGM1.R from the repository root directory.")
}

source("src/gen_ICML_simu.R")
source("src/mixC_dgm_registry.R")
source("src/helpers.R")
source("src/helpers_AIPCW.R")
source("src/dynamicCP_AIPCW.R")

library(survival)
library(randomForestSRC)

`%||%` <- function(x, y) if (is.null(x)) y else x

## ----------------------------- Setup ---------------------------------
setup <- "mixC2_DGM1"
dgm_name <- "mixC_uniform_cox"
dgm_family <- "mixC2"
dgm_variant <- "DGM1"
method <- "static_IPCW0_globalG"

target_map <- list(
    cox_cox_multiclass_loc0 = list(
        tag = "cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone",
        model.pred = "cox",
        model.C = "cox"
    ),
    rsf_rsf_reg_loc0 = list(
        tag = "rsf_rsf_Srsf_Xixgb_reg_loc0_Gtauone",
        model.pred = "rsf",
        model.C = "rsf"
    ),
    rsf_cox_reg_loc0 = list(
        tag = "rsf_cox_Scox_Xixgb_reg_loc0_Gtauone",
        model.pred = "rsf",
        model.C = "cox"
    )
)

rho <- 0.3
n <- 1000L
n_test <- 500L
R <- 200L

args <- commandArgs(trailingOnly = TRUE)
valid_dgms <- c("mixC2_DGM1", "mixC2_DGM3")
if (length(args) >= 1L && args[[1L]] %in% valid_dgms) {
    setup <- args[[1L]]
    dgm_variant <- sub("^mixC2_", "", setup)
    args <- args[-1L]
}

run_target <- "rsf_rsf_reg_loc0"
if (length(args) >= 1L && nzchar(args[[1]])) run_target <- args[[1]]
if (!run_target %in% names(target_map)) {
    stop(
        "Unknown run_target: ", run_target,
        ". Choose one of: ", paste(names(target_map), collapse = ", ")
    )
}
target_cfg <- target_map[[run_target]]
if (length(args) >= 2L && nzchar(args[[2]])) R <- as.integer(args[[2]])
if (length(args) >= 3L && nzchar(args[[3]])) n <- as.integer(args[[3]])
if (length(args) >= 4L && nzchar(args[[4]])) n_test <- as.integer(args[[4]])

folder <- file.path("results", setup)
if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
}

start.name <- "start"
stop.name <- "stop"
event.name <- "event"
event.C.name <- "event.C"
id.name <- "id"
TT.name <- "TT"

covname.pred.baseline <- c("B", "L0")
covname.pred.timevarying <- NULL
covname.C.baseline <- c("B", "L0")
covname.C.timevarying <- NULL

model.pred <- target_cfg$model.pred
model.C <- target_cfg$model.C

change_times <- c(3, 6)
tau_grid <- c(0, change_times)
global_visit_times <- 0

dgm_variant <- resolve_mixC_dgm_variant(dgm_variant, dgm_family = dgm_family)
dgm_parK <- get_mixC_dgm_par(dgm_variant, dgm_family = dgm_family)

alpha <- 0.1
trim.C <- 0.05
theta_grid <- seq(0, 0.5, by = 0.01)
train_ratio <- 0.5
Gtau_mode <- "one"
Gtau_delta <- 0
censoring_method <- "piecewise"
G_fit_type <- "global_baseline_single_interval"

pred_use_history <- FALSE
pred_history_mode <- "wide"
G_use_history <- FALSE
G_history_mode <- "wide"

rsf_args.C <- list()

empty_bounds <- function() {
    list(test = data.frame(), cal = data.frame())
}

make_tau0_global_data <- function(dat) {
    dat0 <- prepare_data_tau(dat, start.name, stop.name, event.name, id.name, tau = 0)
    summ <- subject_summary(dat, id.name, stop.name, event.name, event.C.name)
    idx <- match(as.character(dat0[[id.name]]), as.character(summ[[id.name]]))
    if (anyNA(idx)) {
        stop("make_tau0_global_data(): failed to align tau-0 rows to subject summaries.")
    }

    dat0[[start.name]] <- 0
    dat0[[stop.name]] <- summ$X[idx]
    dat0[[event.name]] <- summ$Delta[idx]
    dat0[[event.C.name]] <- summ$DeltaC[idx]
    dat0
}

make_eval_bounds <- function(dat_tau, b0, survC_by_id, rep_id, tau,
                             stop.name, event.name, id.name, TT.name, trim.C) {
    idx <- match(as.character(dat_tau[[id.name]]), as.character(b0$id))
    if (anyNA(idx)) {
        stop("Failed to match tau=", tau, " ids to tau-0 static bounds.")
    }

    survC <- as.numeric(survC_by_id[as.character(dat_tau[[id.name]])])
    if (anyNA(survC)) {
        stop("Failed to match tau=", tau, " ids to censoring survival estimates.")
    }

    lower <- b0$lower[idx]
    upper <- b0$upper[idx]
    w_ipcw <- dat_tau[[event.name]] / pmax(survC, trim.C)
    covered_X <- as.integer(lower <= dat_tau[[stop.name]] & dat_tau[[stop.name]] <= upper)
    covered_T <- if (!is.null(TT.name) && TT.name %in% names(dat_tau)) {
        as.integer(lower <= dat_tau[[TT.name]] & dat_tau[[TT.name]] <= upper)
    } else {
        NA_integer_
    }

    data.frame(
        rep = rep_id,
        tau = tau,
        id = dat_tau[[id.name]],
        lower = lower,
        upper = upper,
        X = dat_tau[[stop.name]],
        T_true = if (!is.null(TT.name) && TT.name %in% names(dat_tau)) dat_tau[[TT.name]] else NA_real_,
        survC = survC,
        w_ipcw = w_ipcw,
        covered_X = covered_X,
        covered_T = covered_T
    )
}

summarize_bounds <- function(bt, bc) {
    len_test <- bt$upper - bt$lower
    len_cal <- bc$upper - bc$lower
    data.frame(
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
}

truncate_bounds_at_tau <- function(bounds, tau) {
    bounds$test$lower_untruncated <- bounds$test$lower
    bounds$test$upper_untruncated <- bounds$test$upper
    bounds$cal$lower_untruncated <- bounds$cal$lower
    bounds$cal$upper_untruncated <- bounds$cal$upper

    bounds$test$lower <- pmax(bounds$test$lower, tau)
    bounds$test$upper <- pmax(bounds$test$upper, tau)
    bounds$cal$lower <- pmax(bounds$cal$lower, tau)
    bounds$cal$upper <- pmax(bounds$cal$upper, tau)

    bounds$test$covered_X <- as.integer(bounds$test$lower <= bounds$test$X & bounds$test$X <= bounds$test$upper)
    bounds$cal$covered_X <- as.integer(bounds$cal$lower <= bounds$cal$X & bounds$cal$X <= bounds$cal$upper)
    bounds$test$covered_T <- if ("T_true" %in% names(bounds$test)) {
        as.integer(bounds$test$lower <= bounds$test$T_true & bounds$test$T_true <= bounds$test$upper)
    } else {
        NA_integer_
    }
    bounds$cal$covered_T <- if ("T_true" %in% names(bounds$cal)) {
        as.integer(bounds$cal$lower <= bounds$cal$T_true & bounds$cal$T_true <= bounds$cal$upper)
    } else {
        NA_integer_
    }

    bounds
}

run_tau0_ipcw0_globalG <- function(dat0, dat_test0, seed) {
    dynamicCP_AIPCW_split(
        dat = dat0,
        dat_test = dat_test0,
        start.name = start.name,
        stop.name = stop.name,
        event.name = event.name,
        event.C.name = event.C.name,
        id.name = id.name,
        tau = 0,
        model.pred = model.pred,
        covname.pred.baseline = covname.pred.baseline,
        covname.pred.timevarying = covname.pred.timevarying,
        model.C = model.C,
        covname.C.baseline = covname.C.baseline,
        covname.C.timevarying = covname.C.timevarying,
        rsf_args.C = rsf_args.C,
        model.S = "cox",
        model.Xi = "lm",
        TT.name = TT.name,
        visit_times = global_visit_times,
        seed = seed,
        train_ratio = train_ratio,
        theta_grid = theta_grid,
        alpha = alpha,
        trim.C = trim.C,
        pred_use_history = pred_use_history,
        pred_history_mode = pred_history_mode,
        censoring_method = censoring_method,
        G_use_history = G_use_history,
        G_history_mode = G_history_mode,
        Gtau_mode = Gtau_mode,
        Gtau_delta = Gtau_delta,
        ipcw_only = TRUE,
        return_Gfits = TRUE
    )
}

save_bundle <- function(method_name, label, results_summary, results_bounds, results_survC,
                        rep_runtime_sec, method_runtime_sec) {
    res <- do.call(rbind, results_summary)
    res_ok <- res[res$ok %in% TRUE, , drop = FALSE]
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    out_rds <- file.path(
        folder,
        paste0(
            method_name, "_", target_cfg$tag,
            "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
            "_alpha", alpha, ".rds"
        )
    )

    save_obj <- list(
        config = list(
            setup = setup,
            dgm_name = dgm_name,
            dgm_family = dgm_family,
            dgm_variant = dgm_variant,
            run_target = run_target,
            method = method_name,
            source_script = "main_simu_static_IPCW0_globalG_mixC2_DGM1.R",
            tag = target_cfg$tag,
            label = label,
            change_times = change_times,
            dgm_parK = dgm_parK,
            rho = rho,
            n = n,
            n_test = n_test,
            R = R,
            start.name = start.name,
            stop.name = stop.name,
            event.name = event.name,
            event.C.name = event.C.name,
            id.name = id.name,
            covname.pred.baseline = covname.pred.baseline,
            covname.pred.timevarying = covname.pred.timevarying,
            covname.C.baseline = covname.C.baseline,
            covname.C.timevarying = covname.C.timevarying,
            model.pred = model.pred,
            model.C = model.C,
            tau_grid = tau_grid,
            fit_tau = 0,
            alpha = alpha,
            trim.C = trim.C,
            train_ratio = train_ratio,
            TT.name = TT.name,
            theta_grid = theta_grid,
            Gtau_mode = Gtau_mode,
            Gtau_delta = Gtau_delta,
            censoring_method = censoring_method,
            G_fit_type = G_fit_type,
            visit_times = global_visit_times,
            pred_use_history = pred_use_history,
            pred_history_mode = pred_history_mode,
            G_use_history = G_use_history,
            G_history_mode = G_history_mode,
            rsf_args.C = rsf_args.C
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

    saveRDS(save_obj, file = out_rds)
    cat("Saved ", label, " bundle to:\n  ", out_rds, "\n", sep = "")

    if (requireNamespace("dplyr", quietly = TRUE)) {
        print(
            dplyr::summarize(
                dplyr::group_by(res_ok, tau),
                cov_true = mean(coverage_test_true, na.rm = TRUE),
                cov_ipcw = mean(coverage_test_ipcw, na.rm = TRUE),
                mean_len = mean(mean_len_test, na.rm = TRUE),
                theta = mean(theta_hat, na.rm = TRUE),
                .groups = "drop"
            ),
            n = Inf
        )
    }

    invisible(out_rds)
}

## ------------------------ Replicate seeds -----------------------------
set.seed(123)
seeds <- sample(10^7, R, replace = FALSE)
seeds_test <- sample(10^7, R, replace = FALSE)
seeds_split <- sample(10^7, R, replace = FALSE)

summary_ipcw <- vector("list", R)
bounds_ipcw <- vector("list", R)
survC_ipcw <- vector("list", R)
summary_trunc <- vector("list", R)
bounds_trunc <- vector("list", R)
survC_trunc <- vector("list", R)
rep_runtime_sec <- rep(NA_real_, R)

cat("Running ", method, " for ", setup, " with model.pred=", model.pred,
    ", model.C=", model.C, ", target=", run_target,
    ", R=", R, ", n=", n, ", n_test=", n_test, "\n", sep = "")
cat("Global-G trick: collapse to tau-0 subject rows and use visit_times=0.\n")
cat("Replications: ")

method_t0 <- proc.time()[["elapsed"]]
for (r in seq_len(R)) {
    cat(r, "..")
    rep_t0 <- proc.time()[["elapsed"]]

    dat_long_r <- simulate_dataset_long(
        n = n,
        seed = seeds[r],
        change_times = change_times,
        tau_max = Inf,
        no_censoring = FALSE,
        par = dgm_parK,
        rho = rho,
        dgm_name = dgm_name
    )
    dat_test_r <- simulate_dataset_long(
        n = n_test,
        seed = seeds_test[r],
        change_times = change_times,
        tau_max = Inf,
        no_censoring = TRUE,
        par = dgm_parK,
        rho = rho,
        dgm_name = dgm_name
    )

    dat_long_r <- add_baseline_covariates(dat_long_r, id.name = id.name, tv_covs = "L", suffix = "0")
    dat_test_r <- add_baseline_covariates(dat_test_r, id.name = id.name, tv_covs = "L", suffix = "0")

    dat_global_r <- make_tau0_global_data(dat_long_r)
    dat_test_global_r <- make_tau0_global_data(dat_test_r)

    tab_ipcw <- data.frame(
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
    tab_trunc <- tab_ipcw

    bounds_ipcw_r <- vector("list", length(tau_grid))
    names(bounds_ipcw_r) <- paste0("tau_", tau_grid)
    bounds_trunc_r <- bounds_ipcw_r
    survC_ipcw_r <- vector("list", length(tau_grid))
    names(survC_ipcw_r) <- paste0("tau_", tau_grid)
    survC_trunc_r <- survC_ipcw_r

    fit <- tryCatch(
        run_tau0_ipcw0_globalG(dat_global_r, dat_test_global_r, seed = seeds_split[r]),
        error = function(e) list(error = TRUE, msg = conditionMessage(e))
    )

    if (isTRUE(fit$error)) {
        tab_ipcw$ok <- FALSE
        tab_ipcw$msg <- fit$msg
        tab_ipcw$runtime_sec <- proc.time()[["elapsed"]] - rep_t0
        tab_trunc <- tab_ipcw
        for (k in seq_along(tau_grid)) {
            bounds_ipcw_r[[k]] <- empty_bounds()
            bounds_trunc_r[[k]] <- empty_bounds()
        }
        summary_ipcw[[r]] <- tab_ipcw
        bounds_ipcw[[r]] <- bounds_ipcw_r
        survC_ipcw[[r]] <- survC_ipcw_r
        summary_trunc[[r]] <- tab_trunc
        bounds_trunc[[r]] <- bounds_trunc_r
        survC_trunc[[r]] <- survC_trunc_r
        rep_runtime_sec[r] <- proc.time()[["elapsed"]] - rep_t0
        next
    }

    dat_cal_r <- dat_long_r[!(dat_long_r[[id.name]] %in% fit$id_tr), , drop = FALSE]
    dat_cal_global_r <- dat_global_r[!(dat_global_r[[id.name]] %in% fit$id_tr), , drop = FALSE]

    b0_cal <- data.frame(
        id = dat_cal_global_r[[id.name]],
        lower = fit$lower_IPCW_cal,
        upper = fit$upper_IPCW_cal
    )
    b0_test <- data.frame(
        id = dat_test_global_r[[id.name]],
        lower = fit$lower_IPCW_test,
        upper = fit$upper_IPCW_test
    )

    survC_cal_by_id <- stats::setNames(fit$survC_cal, as.character(dat_cal_global_r[[id.name]]))
    survC_test_by_id <- stats::setNames(fit$survC_test, as.character(dat_test_global_r[[id.name]]))

    for (k in seq_along(tau_grid)) {
        tau <- tau_grid[k]
        tau_t0 <- proc.time()[["elapsed"]]
        out_k <- tryCatch({
            dat_cal_tau <- prepare_data_tau(dat_cal_r, start.name, stop.name, event.name, id.name, tau = tau)
            dat_test_tau <- prepare_data_tau(dat_test_r, start.name, stop.name, event.name, id.name, tau = tau)

            bt <- make_eval_bounds(
                dat_test_tau, b0_test, survC_test_by_id, r, tau,
                stop.name, event.name, id.name, TT.name, trim.C
            )
            bc <- make_eval_bounds(
                dat_cal_tau, b0_cal, survC_cal_by_id, r, tau,
                stop.name, event.name, id.name, TT.name, trim.C
            )
            bounds_untrunc <- list(test = bt, cal = bc)
            bounds_trunc_k <- truncate_bounds_at_tau(bounds_untrunc, tau)

            list(
                ok = TRUE,
                untrunc = bounds_untrunc,
                trunc = bounds_trunc_k,
                summary_untrunc = summarize_bounds(bounds_untrunc$test, bounds_untrunc$cal),
                summary_trunc = summarize_bounds(bounds_trunc_k$test, bounds_trunc_k$cal)
            )
        }, error = function(e) {
            list(ok = FALSE, msg = conditionMessage(e))
        })

        if (isTRUE(out_k$ok)) {
            tab_ipcw$ok[k] <- TRUE
            tab_ipcw$theta_hat[k] <- fit$theta_hat_IPCW
            tab_ipcw[k, names(out_k$summary_untrunc)] <- out_k$summary_untrunc
            bounds_ipcw_r[[k]] <- out_k$untrunc
            survC_ipcw_r[[k]] <- list(
                survC_cal = out_k$untrunc$cal$survC,
                survC_test = out_k$untrunc$test$survC
            )

            tab_trunc$ok[k] <- TRUE
            tab_trunc$theta_hat[k] <- fit$theta_hat_IPCW
            tab_trunc[k, names(out_k$summary_trunc)] <- out_k$summary_trunc
            bounds_trunc_r[[k]] <- out_k$trunc
            survC_trunc_r[[k]] <- list(
                survC_cal = out_k$trunc$cal$survC,
                survC_test = out_k$trunc$test$survC
            )
        } else {
            tab_ipcw$ok[k] <- FALSE
            tab_ipcw$msg[k] <- out_k$msg
            bounds_ipcw_r[[k]] <- empty_bounds()

            tab_trunc$ok[k] <- FALSE
            tab_trunc$msg[k] <- out_k$msg
            bounds_trunc_r[[k]] <- empty_bounds()
        }
        tab_ipcw$runtime_sec[k] <- proc.time()[["elapsed"]] - tau_t0
        tab_trunc$runtime_sec[k] <- tab_ipcw$runtime_sec[k]
    }

    summary_ipcw[[r]] <- tab_ipcw
    bounds_ipcw[[r]] <- bounds_ipcw_r
    survC_ipcw[[r]] <- survC_ipcw_r
    summary_trunc[[r]] <- tab_trunc
    bounds_trunc[[r]] <- bounds_trunc_r
    survC_trunc[[r]] <- survC_trunc_r
    rep_runtime_sec[r] <- proc.time()[["elapsed"]] - rep_t0
}
cat("\n")

method_runtime_sec <- proc.time()[["elapsed"]] - method_t0

save_bundle(
    method_name = "static_IPCW0_globalG",
    label = "IPCW",
    results_summary = summary_ipcw,
    results_bounds = bounds_ipcw,
    results_survC = survC_ipcw,
    rep_runtime_sec = rep_runtime_sec,
    method_runtime_sec = method_runtime_sec
)

save_bundle(
    method_name = "static_IPCW_trunc_globalG",
    label = "IPCW-trunc",
    results_summary = summary_trunc,
    results_bounds = bounds_trunc,
    results_survC = survC_trunc,
    rep_runtime_sec = rep_runtime_sec,
    method_runtime_sec = method_runtime_sec
)
