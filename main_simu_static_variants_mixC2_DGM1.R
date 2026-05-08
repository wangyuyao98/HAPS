rm(list = ls())

# Run this script from the repository root.
if (!file.exists("src/gen_ICML_simu.R")) {
    stop("Please run main_simu_static_variants_mixC2_DGM1.R from the repository root directory.")
}

source("src/gen_ICML_simu.R")
source("src/helpers.R")
source("src/helpers_AIPCW.R")
source("src/dynamicCP_AIPCW.R")

library(survival)
library(randomForestSRC)

`%||%` <- function(x, y) if (is.null(x)) y else x

setup <- "mixC2_DGM1"

rho <- 0.3
n <- 1000L
n_test <- 500L
R <- 200L
alpha <- 0.1

format_static_ipcw0_file <- function(tag, prefix = "static_IPCW0") {
    file.path(
        out_dir,
        paste0(
            prefix, "_", tag,
            "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
            "_alpha", alpha, ".rds"
        )
    )
}

format_dynamic_file <- function(method, model.pred, model.C, model.S, model.Xi) {
    file.path(
        out_dir,
        paste0(
            method, "_", model.pred, "_", model.C,
            "_S", model.S, "_Xi", model.Xi,
            "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
            "_alpha", alpha, "_loc0_Gtauone_d0.rds"
        )
    )
}

format_variant_file <- function(method, tag) {
    file.path(
        out_dir,
        paste0(
            method, "_", tag,
            "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
            "_alpha", alpha, ".rds"
        )
    )
}

targets <- list(
    cox_cox_multiclass_loc0 = list(
        tag = "cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone",
        base_method = "dCP_IPCW",
        model.S = "cox",
        model.Xi = "xgb_multiclass",
        base_dynamic_file = NULL,
        static_source = "static_ready",
        static_file = NULL,
        model.pred = "cox",
        model.C = "cox"
    ),
    rsf_rsf_reg_loc0 = list(
        tag = "rsf_rsf_Srsf_Xixgb_reg_loc0_Gtauone",
        base_method = "dCP_IPCW_new",
        model.S = "rsf",
        model.Xi = "xgb_reg",
        base_dynamic_file = NULL,
        static_source = "static_ready",
        static_prefix = "static_IPCW0_globalG",
        static_file = NULL,
        model.pred = "rsf",
        model.C = "rsf"
    ),
    rsf_cox_reg_loc0 = list(
        tag = "rsf_cox_Scox_Xixgb_reg_loc0_Gtauone",
        base_method = "dCP_IPCW_new",
        model.S = "xgb_cox",
        model.Xi = "xgb_multiclass",
        base_dynamic_file = NULL,
        static_source = "static_ready",
        static_prefix = "static_IPCW0_globalG",
        static_file = NULL,
        model.pred = "rsf",
        model.C = "cox"
    )
)

method_names <- c("IPCW_trunc", "LM_IPCW")
args <- commandArgs(trailingOnly = TRUE)
valid_dgms <- c("mixC2_DGM1", "mixC2_DGM3")
dgm_args <- intersect(args, valid_dgms)
if (length(dgm_args) > 1L) {
    stop("Please provide at most one DGM argument. Choose one of: ", paste(valid_dgms, collapse = ", "))
}
if (length(dgm_args) == 1L) {
    setup <- dgm_args[[1L]]
    args <- setdiff(args, dgm_args)
}
out_dir <- file.path("results", setup)
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}
target_args <- intersect(args, names(targets))
method_args <- intersect(args, method_names)
numeric_args <- setdiff(args, c(names(targets), method_names))
bad_numeric_args <- numeric_args[is.na(suppressWarnings(as.integer(numeric_args)))]
unknown <- bad_numeric_args
if (length(unknown) > 0L) {
    stop(
        "Unknown argument(s): ", paste(unknown, collapse = ", "),
        ". Targets: ", paste(names(targets), collapse = ", "),
        ". Methods: ", paste(method_names, collapse = ", "),
        ". Optional numeric arguments: R n n_test"
    )
}
if (length(target_args) > 0L) targets <- targets[target_args]
if (length(method_args) == 0L) method_args <- method_names
if (length(numeric_args) >= 1L && nzchar(numeric_args[[1]])) R <- as.integer(numeric_args[[1]])
if (length(numeric_args) >= 2L && nzchar(numeric_args[[2]])) n <- as.integer(numeric_args[[2]])
if (length(numeric_args) >= 3L && nzchar(numeric_args[[3]])) n_test <- as.integer(numeric_args[[3]])
for (nm in names(targets)) {
    targets[[nm]]$static_file <- format_static_ipcw0_file(
        targets[[nm]]$tag,
        targets[[nm]]$static_prefix %||% "static_IPCW0"
    )
    targets[[nm]]$base_dynamic_file <- format_dynamic_file(
        method = targets[[nm]]$base_method,
        model.pred = targets[[nm]]$model.pred,
        model.C = targets[[nm]]$model.C,
        model.S = targets[[nm]]$model.S,
        model.Xi = targets[[nm]]$model.Xi
    )
}

cat("Static variant targets: ", paste(names(targets), collapse = ", "), "\n", sep = "")
cat("Static variant methods: ", paste(method_args, collapse = ", "), "\n", sep = "")

load_static_ipcw <- function(target) {
    if (!file.exists(target$static_file)) {
        stop("Missing source IPCW file: ", target$static_file)
    }
    obj <- readRDS(target$static_file)
    if (identical(target$static_source, "static_ready")) {
        return(obj)
    }
    if (!identical(target$static_source, "dynamic_tau0")) {
        stop("Unknown static_source: ", target$static_source)
    }

    stat <- build_static_from_dynamic(
        results_bounds = obj$results$results_bounds,
        tau_grid = obj$config$tau_grid,
        TT.name = obj$config$TT.name
    )
    res <- do.call(rbind, stat$results_summary_static)
    obj$config$method <- "static_IPCW0"
    obj$config$fit_tau <- 0
    obj$results <- list(
        results_summary = stat$results_summary_static,
        results_bounds = stat$results_bounds_static,
        results_survC = vector("list", length(stat$results_bounds_static)),
        res = res,
        res_ok = res[res$ok %in% TRUE, , drop = FALSE]
    )
    obj
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

make_empty_bounds <- function() {
    list(test = data.frame(), cal = data.frame())
}

derive_ipcw_trunc <- function(target_name, target) {
    src <- load_static_ipcw(target)
    tau_grid <- src$config$tau_grid
    R_src <- length(src$results$results_bounds)

    results_summary <- vector("list", R_src)
    results_bounds <- vector("list", R_src)
    results_survC <- vector("list", R_src)

    for (r in seq_len(R_src)) {
        tab_r <- src$results$res[src$results$res$rep == r, , drop = FALSE]
        if (nrow(tab_r) == 0L) {
            tab_r <- data.frame(rep = r, tau = tau_grid, ok = FALSE, msg = "missing source rows")
        }
        bounds_r <- vector("list", length(tau_grid))
        names(bounds_r) <- paste0("tau_", tau_grid)
        survC_r <- vector("list", length(tau_grid))
        names(survC_r) <- paste0("tau_", tau_grid)

        for (k in seq_along(tau_grid)) {
            tau <- tau_grid[k]
            bk <- src$results$results_bounds[[r]][[k]]
            out_k <- tryCatch({
                bt <- bk$test
                bc <- bk$cal
                bt$lower_untruncated <- bt$lower
                bt$upper_untruncated <- bt$upper
                bc$lower_untruncated <- bc$lower
                bc$upper_untruncated <- bc$upper
                bt$lower <- pmax(bt$lower, tau)
                bt$upper <- pmax(bt$upper, tau)
                bc$lower <- pmax(bc$lower, tau)
                bc$upper <- pmax(bc$upper, tau)
                bt$covered_X <- as.integer(bt$lower <= bt$X & bt$X <= bt$upper)
                bc$covered_X <- as.integer(bc$lower <= bc$X & bc$X <= bc$upper)
                bt$covered_T <- if ("T_true" %in% names(bt)) as.integer(bt$lower <= bt$T_true & bt$T_true <= bt$upper) else NA_integer_
                bc$covered_T <- if ("T_true" %in% names(bc)) as.integer(bc$lower <= bc$T_true & bc$T_true <= bc$upper) else NA_integer_
                ss <- summarize_bounds(bt, bc)
                list(ok = TRUE, test = bt, cal = bc, summary = ss)
            }, error = function(e) {
                list(ok = FALSE, msg = conditionMessage(e))
            })

            if (isTRUE(out_k$ok)) {
                tab_r$ok[tab_r$tau == tau] <- TRUE
                tab_r[tab_r$tau == tau, names(out_k$summary)] <- out_k$summary
                bounds_r[[k]] <- list(test = out_k$test, cal = out_k$cal)
                survC_r[[k]] <- list(
                    survC_cal = out_k$cal$survC %||% NA_real_,
                    survC_test = out_k$test$survC %||% NA_real_
                )
            } else {
                tab_r$ok[tab_r$tau == tau] <- FALSE
                tab_r$msg[tab_r$tau == tau] <- out_k$msg
                bounds_r[[k]] <- make_empty_bounds()
            }
        }
        results_summary[[r]] <- tab_r
        results_bounds[[r]] <- bounds_r
        results_survC[[r]] <- survC_r
    }

    res <- do.call(rbind, results_summary)
    cfg <- src$config
    cfg$setup <- setup
    cfg$method <- "static_IPCW_trunc"
    cfg$source_script <- "main_simu_static_variants_mixC2_DGM1.R"
    cfg$source_file <- target$static_file
    cfg$target <- target_name
    cfg$label <- "IPCW-trunc"

    save_obj <- list(
        config = cfg,
        seeds = src$seeds,
        results = list(
            results_summary = results_summary,
            results_bounds = results_bounds,
            results_survC = results_survC,
            res = res,
            res_ok = res[res$ok %in% TRUE, , drop = FALSE]
        ),
        meta = list(timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"))
    )

    out_file <- format_variant_file("static_IPCW_trunc", target$tag)
    saveRDS(save_obj, out_file)
    cat("Saved IPCW-trunc to:\n  ", out_file, "\n", sep = "")
    invisible(out_file)
}

subject_ids_at_risk <- function(dat, tau, cfg) {
    summ <- subject_summary(dat, cfg$id.name, cfg$stop.name, cfg$event.name, cfg$event.C.name)
    as.character(summ[[cfg$id.name]][summ$X > tau])
}

make_landmark_residual_data <- function(dat, tau, cfg) {
    dat_tau <- prepare_data_tau(
        dat, cfg$start.name, cfg$stop.name,
        cfg$event.name, cfg$id.name, tau
    )
    if (nrow(dat_tau) == 0L) {
        return(dat_tau)
    }

    summ <- subject_summary(dat, cfg$id.name, cfg$stop.name, cfg$event.name, cfg$event.C.name)
    idx <- match(as.character(dat_tau[[cfg$id.name]]), as.character(summ[[cfg$id.name]]))
    if (anyNA(idx)) {
        stop("make_landmark_residual_data(): failed to align landmark rows to subject summary.")
    }

    dat_tau[[cfg$start.name]] <- 0
    dat_tau[[cfg$stop.name]] <- pmax(summ$X[idx] - tau, 0)
    dat_tau[[cfg$event.name]] <- summ$Delta[idx]
    dat_tau[[cfg$event.C.name]] <- summ$DeltaC[idx]
    if (!is.null(cfg$TT.name) && cfg$TT.name %in% names(dat_tau)) {
        dat_tau[[cfg$TT.name]] <- dat_tau[[cfg$TT.name]] - tau
    }

    dat_tau
}

make_bounds_from_fit <- function(dat_tau, lower, upper, survC, rep_id, tau, cfg) {
    w_ipcw <- dat_tau[[cfg$event.name]] / pmax(survC, cfg$trim.C)
    covered_X <- as.integer(lower <= dat_tau[[cfg$stop.name]] & dat_tau[[cfg$stop.name]] <= upper)
    covered_T <- if (!is.null(cfg$TT.name) && cfg$TT.name %in% names(dat_tau)) {
        as.integer(lower <= dat_tau[[cfg$TT.name]] & dat_tau[[cfg$TT.name]] <= upper)
    } else {
        NA_integer_
    }

    data.frame(
        rep = rep_id,
        tau = tau,
        id = dat_tau[[cfg$id.name]],
        lower = lower,
        upper = upper,
        X = dat_tau[[cfg$stop.name]],
        T_true = if (!is.null(cfg$TT.name) && cfg$TT.name %in% names(dat_tau)) dat_tau[[cfg$TT.name]] else NA_real_,
        survC = survC,
        w_ipcw = w_ipcw,
        covered_X = covered_X,
        covered_T = covered_T
    )
}

run_lm_ipcw <- function(target_name, target) {
    base_file <- target$base_dynamic_file
    if (!file.exists(base_file) && file.exists(target$static_file)) {
        base_file <- target$static_file
    }
    if (!file.exists(base_file)) {
        stop(
            "Missing source file for IPCW-LM. Tried dynamic file: ",
            target$base_dynamic_file,
            " and static file: ",
            target$static_file
        )
    }
    base_obj <- readRDS(base_file)
    cfg_base <- base_obj$config
    seeds <- base_obj$seeds
    tau_grid <- cfg_base$tau_grid
    dgm_name <- cfg_base$dgm_name %||% "mixC_uniform_cox"
    no_censoring_test <- identical(cfg_base$Gtau_mode, "one")

    cfg_lm <- cfg_base
    cfg_lm$setup <- setup
    cfg_lm$method <- "static_LM_IPCW"
    cfg_lm$source_script <- "main_simu_static_variants_mixC2_DGM1.R"
    cfg_lm$source_file <- base_file
    cfg_lm$target <- target_name
    cfg_lm$label <- "IPCW-LM"
    cfg_lm$model.pred <- target$model.pred
    cfg_lm$model.C <- target$model.C
    cfg_lm$covname.pred.baseline <- c("B", "L0")
    cfg_lm$covname.pred.timevarying <- NULL
    cfg_lm$covname.C.baseline <- c("B", "L0")
    cfg_lm$covname.C.timevarying <- NULL
    cfg_lm$fit_tau <- "landmark_residual_time"

    results_summary <- vector("list", cfg_base$R)
    results_bounds <- vector("list", cfg_base$R)
    results_survC <- vector("list", cfg_base$R)
    rep_runtime_sec <- rep(NA_real_, cfg_base$R)

    cat("Running IPCW-LM for ", target_name, "\n", sep = "")
    cat("Replications: ")
    method_t0 <- proc.time()[["elapsed"]]

    for (r in seq_len(cfg_base$R)) {
        cat(r, "..")
        rep_t0 <- proc.time()[["elapsed"]]

        dat_long_r <- simulate_dataset_long(
            n = cfg_base$n,
            seed = seeds$seeds[r],
            change_times = cfg_base$change_times,
            tau_max = Inf,
            no_censoring = FALSE,
            par = cfg_base$dgm_parK,
            rho = cfg_base$rho,
            dgm_name = dgm_name
        )
        dat_test_r <- simulate_dataset_long(
            n = cfg_base$n_test,
            seed = seeds$seeds_test[r],
            change_times = cfg_base$change_times,
            tau_max = Inf,
            no_censoring = no_censoring_test,
            par = cfg_base$dgm_parK,
            rho = cfg_base$rho,
            dgm_name = dgm_name
        )

        dat_long_r <- add_baseline_covariates(dat_long_r, id.name = cfg_base$id.name, tv_covs = "L", suffix = "0")
        dat_test_r <- add_baseline_covariates(dat_test_r, id.name = cfg_base$id.name, tv_covs = "L", suffix = "0")

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

        for (k in seq_along(tau_grid)) {
            tau <- tau_grid[k]
            tau_t0 <- proc.time()[["elapsed"]]

            out_k <- tryCatch({
                ids_long <- subject_ids_at_risk(dat_long_r, tau, cfg_base)
                ids_test <- subject_ids_at_risk(dat_test_r, tau, cfg_base)
                dat_land <- dat_long_r[as.character(dat_long_r[[cfg_base$id.name]]) %in% ids_long, , drop = FALSE]
                dat_test_land <- dat_test_r[as.character(dat_test_r[[cfg_base$id.name]]) %in% ids_test, , drop = FALSE]
                dat_land_res <- make_landmark_residual_data(dat_land, tau, cfg_base)
                dat_test_land_res <- make_landmark_residual_data(dat_test_land, tau, cfg_base)
                residual_visit_times <- 0

                fit <- dynamicCP_AIPCW_split(
                    dat = dat_land_res,
                    dat_test = dat_test_land_res,
                    start.name = cfg_base$start.name,
                    stop.name = cfg_base$stop.name,
                    event.name = cfg_base$event.name,
                    event.C.name = cfg_base$event.C.name,
                    id.name = cfg_base$id.name,
                    tau = 0,
                    model.pred = target$model.pred,
                    covname.pred.baseline = c("B", "L0"),
                    covname.pred.timevarying = NULL,
                    model.C = target$model.C,
                    covname.C.baseline = c("B", "L0"),
                    covname.C.timevarying = NULL,
                    rsf_args.C = cfg_base$rsf_args.C %||% list(),
                    model.S = cfg_base$model.S %||% "cox",
                    model.Xi = cfg_base$model.Xi %||% "lm",
                    TT.name = cfg_base$TT.name,
                    visit_times = residual_visit_times,
                    seed = seeds$seeds_split[r],
                    train_ratio = cfg_base$train_ratio,
                    theta_grid = cfg_base$theta_grid_AIPCW %||% seq(0, 0.5, by = 0.01),
                    alpha = cfg_base$alpha,
                    trim.C = cfg_base$trim.C,
                    pred_use_history = FALSE,
                    pred_history_mode = cfg_base$pred_history_mode %||% "wide",
                    censoring_method = cfg_base$censoring_method_AIPCW %||% "piecewise",
                    G_use_history = FALSE,
                    G_history_mode = cfg_base$G_history_mode %||% "wide",
                    Gtau_mode = cfg_base$Gtau_mode,
                    Gtau_delta = cfg_base$Gtau_delta,
                    ipcw_only = TRUE
                )

                dat_cal_land <- dat_land[!(dat_land[[cfg_base$id.name]] %in% fit$id_tr), , drop = FALSE]
                dat_cal_tau <- prepare_data_tau(
                    dat_cal_land, cfg_base$start.name, cfg_base$stop.name,
                    cfg_base$event.name, cfg_base$id.name, tau
                )
                dat_test_tau <- prepare_data_tau(
                    dat_test_land, cfg_base$start.name, cfg_base$stop.name,
                    cfg_base$event.name, cfg_base$id.name, tau
                )

                lower_cal <- fit$lower_IPCW_cal + tau
                upper_cal <- fit$upper_IPCW_cal + tau
                lower_test <- fit$lower_IPCW_test + tau
                upper_test <- fit$upper_IPCW_test + tau

                bc <- make_bounds_from_fit(dat_cal_tau, lower_cal, upper_cal, fit$survC_cal, r, tau, cfg_base)
                bt <- make_bounds_from_fit(dat_test_tau, lower_test, upper_test, fit$survC_test, r, tau, cfg_base)
                ss <- summarize_bounds(bt, bc)
                list(ok = TRUE, test = bt, cal = bc, summary = ss, theta = fit$theta_hat_IPCW)
            }, error = function(e) {
                list(ok = FALSE, msg = conditionMessage(e))
            })

            if (isTRUE(out_k$ok)) {
                tab_r$ok[k] <- TRUE
                tab_r$theta_hat[k] <- out_k$theta
                tab_r[k, names(out_k$summary)] <- out_k$summary
                bounds_r[[k]] <- list(test = out_k$test, cal = out_k$cal)
                survC_r[[k]] <- list(
                    survC_cal = out_k$cal$survC,
                    survC_test = out_k$test$survC
                )
            } else {
                tab_r$ok[k] <- FALSE
                tab_r$msg[k] <- out_k$msg
                bounds_r[[k]] <- make_empty_bounds()
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
    save_obj <- list(
        config = cfg_lm,
        seeds = seeds,
        results = list(
            results_summary = results_summary,
            results_bounds = results_bounds,
            results_survC = results_survC,
            res = res,
            res_ok = res[res$ok %in% TRUE, , drop = FALSE]
        ),
        meta = list(
            timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
            runtime = list(
                method_runtime_sec = proc.time()[["elapsed"]] - method_t0,
                rep_runtime_sec = rep_runtime_sec
            )
        )
    )

    out_file <- format_variant_file("static_LM_IPCW", target$tag)
    saveRDS(save_obj, out_file)
    cat("Saved IPCW-LM to:\n  ", out_file, "\n", sep = "")
    invisible(out_file)
}

for (target_name in names(targets)) {
    target <- targets[[target_name]]
    if ("IPCW_trunc" %in% method_args) {
        derive_ipcw_trunc(target_name, target)
    }
    if ("LM_IPCW" %in% method_args) {
        run_lm_ipcw(target_name, target)
    }
}
