rm(list = ls())

# Run this script from the repository root.
if (!file.exists("src/gen_ICML_simu.R")) {
    stop("Please run main_simu_static_IPCW0_linWB1.R from the repository root directory.")
}

source("src/gen_ICML_simu.R")
source("src/helpers.R")
source("src/helpers_AIPCW.R")
source("src/dynamicCP_AIPCW.R")

library(survival)
library(randomForestSRC)

## ----------------------------- Setup ---------------------------------
setup <- "linWB1"
dgm_name <- "linear_weibull"
method <- "static_IPCW0"

rho <- 0.3
n <- 1000L
n_test <- 500L
R <- 200L

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
    )
)

args <- commandArgs(trailingOnly = TRUE)
run_target <- "cox_cox_multiclass_loc0"
target_args <- intersect(args, names(target_map))
if (length(target_args) > 1L) {
    stop("Please provide at most one target. Choose one of: ", paste(names(target_map), collapse = ", "))
}
if (length(target_args) == 1L) {
    run_target <- target_args[[1L]]
    args <- setdiff(args, target_args)
}
bad_numeric_args <- args[is.na(suppressWarnings(as.integer(args)))]
if (length(bad_numeric_args) > 0L) {
    stop(
        "Unknown argument(s): ", paste(bad_numeric_args, collapse = ", "),
        ". Targets: ", paste(names(target_map), collapse = ", "),
        ". Optional numeric arguments: R n n_test"
    )
}
if (length(args) >= 1L && nzchar(args[[1L]])) R <- as.integer(args[[1L]])
if (length(args) >= 2L && nzchar(args[[2L]])) n <- as.integer(args[[2L]])
if (length(args) >= 3L && nzchar(args[[3L]])) n_test <- as.integer(args[[3L]])
if (!is.finite(R) || R < 1L) stop("R must be a positive integer.")
if (!is.finite(n) || n < 10L) stop("n must be at least 10.")
if (!is.finite(n_test) || n_test < 10L) stop("n_test must be at least 10.")
target_cfg <- target_map[[run_target]]

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

covname.pred.baseline <- "L0"
covname.pred.timevarying <- NULL
covname.C.baseline <- "L0"
covname.C.timevarying <- NULL

model.pred <- target_cfg$model.pred
model.C <- target_cfg$model.C

change_times <- c(3, 6)
tau_grid <- c(0, change_times)
visit_times <- tau_grid

dgm_parK <- list(
    shape_T = c(4, 4, 5),
    shape_C = c(3, 3, 3),
    beta_T0 = c(-8, -8, -5),
    beta_TL = c(1, 2, 3),
    beta_C0 = c(-6, -6, -5),
    beta_CL = c(2, 2, 2)
)

alpha <- 0.1
trim.C <- 0.05
theta_grid <- seq(0, 0.5, by = 0.01)
train_ratio <- 0.5
Gtau_mode <- "one"
Gtau_delta <- 0
censoring_method <- "piecewise"

pred_use_history <- FALSE
pred_history_mode <- "wide"
G_use_history <- FALSE
G_history_mode <- "wide"

rsf_args.C <- list()

format_tag_num <- function(x) {
    x_chr <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
    x_chr <- gsub("-", "m", x_chr, fixed = TRUE)
    x_chr <- gsub(".", "p", x_chr, fixed = TRUE)
    x_chr
}

empty_bounds <- function(rep_id, tau) {
    list(
        test = data.frame(),
        cal = data.frame()
    )
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
        q25_len_test = as.numeric(quantile(len_test, 0.25, na.rm = TRUE)),
        q50_len_test = as.numeric(quantile(len_test, 0.50, na.rm = TRUE)),
        q75_len_test = as.numeric(quantile(len_test, 0.75, na.rm = TRUE))
    )
}

run_tau0_ipcw0 <- function(dat_long, dat_test, seed) {
    dynamicCP_AIPCW_split(
        dat = dat_long,
        dat_test = dat_test,
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
        visit_times = visit_times,
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

## ------------------------ Replicate seeds -----------------------------
set.seed(123)
seeds <- sample(10^7, R, replace = FALSE)
seeds_test <- sample(10^7, R, replace = FALSE)
seeds_split <- sample(10^7, R, replace = FALSE)

results_summary <- vector("list", R)
results_bounds <- vector("list", R)
results_survC <- vector("list", R)
rep_runtime_sec <- rep(NA_real_, R)

cat("Running ", method, " for ", setup, " with model.pred=", model.pred,
    ", model.C=", model.C, ", R=", R, ", n=", n, ", n_test=", n_test, "\n", sep = "")
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
    bounds_r <- vector("list", length(tau_grid))
    names(bounds_r) <- paste0("tau_", tau_grid)
    survC_r <- vector("list", length(tau_grid))
    names(survC_r) <- paste0("tau_", tau_grid)

    fit <- tryCatch(
        run_tau0_ipcw0(dat_long_r, dat_test_r, seed = seeds_split[r]),
        error = function(e) list(error = TRUE, msg = conditionMessage(e))
    )

    if (isTRUE(fit$error)) {
        tab_r$ok <- FALSE
        tab_r$msg <- fit$msg
        tab_r$runtime_sec <- proc.time()[["elapsed"]] - rep_t0
        for (k in seq_along(tau_grid)) bounds_r[[k]] <- empty_bounds(r, tau_grid[k])
        results_summary[[r]] <- tab_r
        results_bounds[[r]] <- bounds_r
        results_survC[[r]] <- survC_r
        rep_runtime_sec[r] <- proc.time()[["elapsed"]] - rep_t0
        next
    }

    dat_cal_r <- dat_long_r[!(dat_long_r[[id.name]] %in% fit$id_tr), , drop = FALSE]

    dat_cal_tau0 <- prepare_data_tau(dat_cal_r, start.name, stop.name, event.name, id.name, tau = 0)
    dat_test_tau0 <- prepare_data_tau(dat_test_r, start.name, stop.name, event.name, id.name, tau = 0)

    b0_cal <- data.frame(
        id = dat_cal_tau0[[id.name]],
        lower = fit$lower_IPCW_cal,
        upper = fit$upper_IPCW_cal
    )
    b0_test <- data.frame(
        id = dat_test_tau0[[id.name]],
        lower = fit$lower_IPCW_test,
        upper = fit$upper_IPCW_test
    )

    G0X_cal <- compute_Hk_at_tk(
        dat_cal_r, visit_times, k = 1,
        start.name, stop.name, event.name, event.C.name, id.name,
        Gfits = fit$Gfits,
        trim.C = trim.C
    )
    G0X_test <- compute_Hk_at_tk(
        dat_test_r, visit_times, k = 1,
        start.name, stop.name, event.name, event.C.name, id.name,
        Gfits = fit$Gfits,
        trim.C = trim.C
    )
    survC_cal_by_id <- stats::setNames(G0X_cal$Hk, as.character(G0X_cal$id))
    survC_test_by_id <- stats::setNames(G0X_test$Hk, as.character(G0X_test$id))

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

            summary_k <- summarize_bounds(bt, bc)
            list(ok = TRUE, test = bt, cal = bc, summary = summary_k)
        }, error = function(e) {
            list(ok = FALSE, msg = conditionMessage(e))
        })

        if (isTRUE(out_k$ok)) {
            tab_r$ok[k] <- TRUE
            tab_r$theta_hat[k] <- fit$theta_hat_IPCW
            tab_r[k, names(out_k$summary)] <- out_k$summary
            bounds_r[[k]] <- list(test = out_k$test, cal = out_k$cal)
            survC_r[[k]] <- list(
                survC_cal = out_k$cal$survC,
                survC_test = out_k$test$survC
            )
        } else {
            tab_r$ok[k] <- FALSE
            tab_r$msg[k] <- out_k$msg
            bounds_r[[k]] <- empty_bounds(r, tau)
        }
        tab_r$runtime_sec[k] <- proc.time()[["elapsed"]] - tau_t0
    }

    results_summary[[r]] <- tab_r
    results_bounds[[r]] <- bounds_r
    results_survC[[r]] <- survC_r
    rep_runtime_sec[r] <- proc.time()[["elapsed"]] - rep_t0
}
cat("\n")

method_runtime_sec <- proc.time()[["elapsed"]] - method_t0
res <- do.call(rbind, results_summary)
res_ok <- res[res$ok %in% TRUE, , drop = FALSE]
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

out_rds <- file.path(
    folder,
    paste0(
        method, "_", target_cfg$tag,
        "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
        "_alpha", alpha, ".rds"
    )
)

save_obj <- list(
    config = list(
        setup = setup,
        dgm_name = dgm_name,
        method = method,
        source_script = "main_simu_static_IPCW0_linWB1.R",
        target = run_target,
        tag = target_cfg$tag,
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
cat("Saved static IPCW0 bundle to:\n  ", out_rds, "\n", sep = "")

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
