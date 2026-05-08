rm(list = ls())

# Generic application driver for non-AIPCW HAPS analyses.
# Run this script from the repository root.
if (!file.exists("src/app_analysis_helpers.R")) {
    stop("Please run main_app_generic.R from the repository root directory.")
}

source("src/dynamicCP.R")
source("src/helpers.R")
source("src/helpers_AIPCW.R")
source("src/app_analysis_helpers.R")

library(survival)

args <- commandArgs(trailingOnly = TRUE)

dataset_arg <- if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else "all"
alpha <- if (length(args) >= 2L && nzchar(args[[2L]])) as.numeric(args[[2L]]) else 0.2
R_limit <- if (length(args) >= 3L && nzchar(args[[3L]])) as.integer(args[[3L]]) else 200L
model.pred <- if (length(args) >= 4L && nzchar(args[[4L]])) args[[4L]] else "cox"
model.C <- if (length(args) >= 5L && nzchar(args[[5L]])) args[[5L]] else "cox"
steps_arg <- if (length(args) >= 6L && nzchar(args[[6L]])) args[[6L]] else "all"
overwrite <- if (length(args) >= 7L && nzchar(args[[7L]])) {
    tolower(args[[7L]]) %in% c("true", "t", "yes", "y", "1", "overwrite")
} else {
    FALSE
}

valid_datasets <- c("all", "pbcseq", "colon")
if (!dataset_arg %in% valid_datasets) {
    stop("Unknown dataset: ", dataset_arg, ". Choose one of: ", paste(valid_datasets, collapse = ", "))
}
if (!is.finite(alpha)) stop("alpha must be numeric.")
if (!is.finite(R_limit) || R_limit < 1L) stop("R must be a positive integer.")

model.pred <- match.arg(model.pred, c("cox", "rsf"))
model.C <- match.arg(model.C, c("cox", "rsf"))
steps <- strsplit(steps_arg, ",", fixed = TRUE)[[1]]
steps <- trimws(steps)
valid_steps <- c("all", "dynamic", "static")
bad_steps <- setdiff(steps, valid_steps)
if (length(bad_steps) > 0L) {
    stop("Unknown step(s): ", paste(bad_steps, collapse = ", "),
         ". Choose from: ", paste(valid_steps, collapse = ", "))
}
if ("all" %in% steps) steps <- c("dynamic", "static")

datasets <- if (identical(dataset_arg, "all")) c("pbcseq", "colon") else dataset_arg

model.C.args <- if (identical(model.C, "rsf")) {
    list(ntree = 500L, nodesize = 15L, nsplit = 10L, forest = TRUE, rsf_min_censor_events = 10L)
} else {
    list()
}

model_tag <- function(model.pred, model.C) {
    pred_tag <- switch(model.pred, cox = "coxpred", rsf = "rsfpred")
    C_tag <- switch(model.C, cox = "coxG", rsf = "rsfG")
    paste(pred_tag, C_tag, sep = "_")
}

analysis_tag <- model_tag(model.pred, model.C)

generic_app_file <- function(cfg, method, alpha, R) {
    file.path(
        app_out_dir(cfg),
        paste0(
            method, "_", analysis_tag,
            "_R", as.integer(R), "_n", cfg$n_all,
            "_alpha", alpha, ".rds"
        )
    )
}

assert_can_write <- function(path) {
    if (file.exists(path) && !isTRUE(overwrite)) {
        stop(
            "Refusing to overwrite existing file: ", path, "\n",
            "Pass TRUE as the seventh argument to overwrite."
        )
    }
}

make_app_summary_tab <- function(rep_id, tau_grid) {
    data.frame(
        rep = rep_id,
        tau = tau_grid,
        ok = FALSE,
        msg = NA_character_,
        runtime_sec = NA_real_,
        theta_hat = NA_real_,
        coverage_cal_ipcw = NA_real_,
        coverage_test_ipcw = NA_real_,
        mean_len_cal = NA_real_,
        mean_len_test = NA_real_,
        q25_len_test = NA_real_,
        q50_len_test = NA_real_,
        q75_len_test = NA_real_
    )
}

save_app_generic_bundle <- function(cfg0, base_cfg, seeds, results_summary, results_bounds,
                                    results_survC, method, label, R, fit_tau,
                                    extra_config = list()) {
    res <- do.call(rbind, results_summary)
    cfg_save <- base_cfg
    cfg_save$setup <- cfg0$out_setup
    cfg_save$dataset_name <- cfg0$dataset_name
    cfg_save$n_all <- cfg0$n_all
    cfg_save$method <- method
    cfg_save$label <- label
    cfg_save$source_script <- "main_app_generic.R"
    cfg_save$analysis_tag <- analysis_tag
    cfg_save$model.pred <- model.pred
    cfg_save$model.C <- model.C
    cfg_save$model.C.args <- model.C.args
    cfg_save$R <- R
    cfg_save$fit_tau <- fit_tau
    cfg_save$alpha <- alpha
    for (nm in names(extra_config)) {
        cfg_save[[nm]] <- extra_config[[nm]]
    }

    save_obj <- list(
        config = cfg_save,
        seeds = seeds,
        results = list(
            results_summary = results_summary,
            results_bounds = results_bounds,
            results_survC = results_survC,
            res = res,
            res_ok = res[res$ok %in% TRUE, , drop = FALSE]
        ),
        meta = list(timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"))
    )

    out_file <- generic_app_file(cfg0, method, alpha, R)
    assert_can_write(out_file)
    saveRDS(save_obj, out_file)
    cat("Saved ", label, " to:\n  ", out_file, "\n", sep = "")
    invisible(out_file)
}

load_base_application_run <- function(cfg0) {
    app_base_run(cfg0, alpha = alpha, R = R_limit)
}

split_app_data <- function(dat0, cfg, seeds, r) {
    id_all <- unique(dat0[[cfg$id.name]])
    set.seed(seeds$seeds_test[r])
    id_test <- sample(id_all, floor(length(id_all) * cfg$test_ratio), replace = FALSE)
    list(
        dat_test = dat0[dat0[[cfg$id.name]] %in% id_test, , drop = FALSE],
        dat_long = dat0[!(dat0[[cfg$id.name]] %in% id_test), , drop = FALSE]
    )
}

run_app_dynamic_generic <- function(dataset_name, alpha, R_limit) {
    cfg0 <- app_dataset_config(dataset_name, alpha)
    dat0 <- prepare_app_data(dataset_name)
    base <- load_base_application_run(cfg0)
    base_cfg <- base$config
    seeds <- base$seeds
    R <- min(as.integer(R_limit), base$R_full)
    tau_grid <- cfg0$tau_list_plot

    results <- list(
        dCP_IPCW = list(
            label = "HAPS",
            results_summary = vector("list", R),
            results_bounds = vector("list", R),
            results_survC = vector("list", R)
        ),
        derived_DynaCoPS_DR = list(
            label = "HAPS-DR",
            results_summary = vector("list", R),
            results_bounds = vector("list", R),
            results_survC = vector("list", R)
        )
    )

    cat("\nGeneric application HAPS run for ", dataset_name, "\n", sep = "")
    cat("Analysis tag: ", analysis_tag, "\n", sep = "")
    cat("Replications: ", R, "/", base$R_full, "\n", sep = "")
    cat("Tau grid: ", paste(tau_grid, collapse = ", "), "\n", sep = "")
    cat("Dynamic prediction covariates: baseline + current time-varying covariates.\n")
    cat("Dynamic G covariates: baseline + current time-varying covariates.\n")
    cat("Replications: ")

    for (r in seq_len(R)) {
        cat(r, "..")
        split <- split_app_data(dat0, base_cfg, seeds, r)

        rep_tabs <- lapply(results, function(...) make_app_summary_tab(r, tau_grid))
        rep_bounds <- lapply(results, function(...) {
            x <- vector("list", length(tau_grid))
            names(x) <- paste0("tau_", tau_grid)
            x
        })
        rep_survC <- rep_bounds

        for (k in seq_along(tau_grid)) {
            tau <- tau_grid[k]
            tau_t0 <- proc.time()[["elapsed"]]

            fit <- tryCatch({
                dynamicCP_split(
                    dat = split$dat_long,
                    dat_test = split$dat_test,
                    start.name = base_cfg$start.name,
                    stop.name = base_cfg$stop.name,
                    event.name = base_cfg$event.name,
                    event.C.name = base_cfg$event.C.name,
                    id.name = base_cfg$id.name,
                    tau = tau,
                    model.pred = model.pred,
                    covname.pred.baseline = base_cfg$covname.pred.baseline,
                    covname.pred.timevarying = base_cfg$covname.pred.timevarying,
                    model.C = model.C,
                    model.C.args = model.C.args,
                    covname.C.baseline = base_cfg$covname.C.baseline,
                    covname.C.timevarying = base_cfg$covname.C.timevarying,
                    TT.name = NULL,
                    seed = seeds$seeds_split[r],
                    train_ratio = base_cfg$train_ratio,
                    theta_grid = base_cfg$theta_grid,
                    alpha = alpha,
                    trim.C = base_cfg$trim.C
                )
            }, error = function(e) {
                list(error = TRUE, msg = conditionMessage(e))
            })
            runtime_k <- proc.time()[["elapsed"]] - tau_t0

            if (isTRUE(fit$error)) {
                for (method in names(results)) {
                    rep_tabs[[method]]$msg[k] <- fit$msg
                    rep_tabs[[method]]$runtime_sec[k] <- runtime_k
                    rep_bounds[[method]][[k]] <- empty_app_bounds()
                }
                next
            }

            dat_cal_tau <- prepare_data_tau(
                split$dat_long[!(split$dat_long[[base_cfg$id.name]] %in% fit$id_tr), , drop = FALSE],
                base_cfg$start.name,
                base_cfg$stop.name,
                base_cfg$event.name,
                base_cfg$id.name,
                tau
            )
            dat_test_tau <- prepare_data_tau(
                split$dat_test,
                base_cfg$start.name,
                base_cfg$stop.name,
                base_cfg$event.name,
                base_cfg$id.name,
                tau
            )

            method_bounds <- list(
                dCP_IPCW = list(
                    lower_cal = fit$lower_IPCW_cal,
                    upper_cal = fit$upper_IPCW_cal,
                    lower_test = fit$lower_IPCW_test,
                    upper_test = fit$upper_IPCW_test
                ),
                derived_DynaCoPS_DR = list(
                    lower_cal = fit$lower_IPCW_DR_cal,
                    upper_cal = fit$upper_IPCW_DR_cal,
                    lower_test = fit$lower_IPCW_DR_test,
                    upper_test = fit$upper_IPCW_DR_test
                )
            )

            for (method in names(method_bounds)) {
                mb <- method_bounds[[method]]
                bc <- make_app_bounds_from_fit(
                    dat_tau = dat_cal_tau,
                    lower = mb$lower_cal,
                    upper = mb$upper_cal,
                    survC = fit$survC_cal,
                    rep_id = r,
                    tau = tau,
                    cfg = base_cfg
                )
                bt <- make_app_bounds_from_fit(
                    dat_tau = dat_test_tau,
                    lower = mb$lower_test,
                    upper = mb$upper_test,
                    survC = fit$survC_test,
                    rep_id = r,
                    tau = tau,
                    cfg = base_cfg
                )
                ss <- summarize_app_bounds(bt, bc)

                rep_tabs[[method]]$ok[k] <- TRUE
                rep_tabs[[method]]$theta_hat[k] <- fit$theta_hat_IPCW
                rep_tabs[[method]][k, names(ss)] <- ss
                rep_tabs[[method]]$runtime_sec[k] <- runtime_k
                rep_bounds[[method]][[k]] <- list(test = bt, cal = bc)
                rep_survC[[method]][[k]] <- list(
                    survC_cal = fit$survC_cal,
                    survC_test = fit$survC_test
                )
            }
        }

        for (method in names(results)) {
            results[[method]]$results_summary[[r]] <- rep_tabs[[method]]
            results[[method]]$results_bounds[[r]] <- rep_bounds[[method]]
            results[[method]]$results_survC[[r]] <- rep_survC[[method]]
        }
    }
    cat("\n")

    out_files <- character(0)
    for (method in names(results)) {
        out_files <- c(
            out_files,
            save_app_generic_bundle(
                cfg0 = cfg0,
                base_cfg = base_cfg,
                seeds = seeds,
                results_summary = results[[method]]$results_summary,
                results_bounds = results[[method]]$results_bounds,
                results_survC = results[[method]]$results_survC,
                method = method,
                label = results[[method]]$label,
                R = R,
                fit_tau = "landmark_dynamic",
                extra_config = list(
                    source_file = base$path,
                    covname.pred.baseline = base_cfg$covname.pred.baseline,
                    covname.pred.timevarying = base_cfg$covname.pred.timevarying,
                    covname.C.baseline = base_cfg$covname.C.baseline,
                    covname.C.timevarying = base_cfg$covname.C.timevarying
                )
            )
        )
    }

    invisible(out_files)
}

run_app_static_source_generic <- function(dataset_name, alpha, R_limit) {
    cfg0 <- app_dataset_config(dataset_name, alpha)
    dat0 <- prepare_app_data(dataset_name)
    base <- load_base_application_run(cfg0)
    base_cfg <- base$config
    seeds <- base$seeds
    R <- min(as.integer(R_limit), base$R_full)
    tau_grid <- cfg0$tau_list_plot

    tv_covs_all <- union(base_cfg$covname.pred.timevarying, base_cfg$covname.C.timevarying)
    cov_pred_baseline_static <- c(
        base_cfg$covname.pred.baseline,
        paste0(base_cfg$covname.pred.timevarying, "0")
    )
    cov_C_baseline_static <- c(
        base_cfg$covname.C.baseline,
        paste0(base_cfg$covname.C.timevarying, "0")
    )

    results_summary_source <- vector("list", R)
    results_bounds_source <- vector("list", R)
    results_survC_source <- vector("list", R)

    cat("\nGeneric application static source run for ", dataset_name, "\n", sep = "")
    cat("Analysis tag: ", analysis_tag, "\n", sep = "")
    cat("Static prediction/G covariates: baseline covariates plus time-zero covariates only.\n")
    cat("Replications: ")

    for (r in seq_len(R)) {
        cat(r, "..")
        split <- split_app_data(dat0, base_cfg, seeds, r)
        split$dat_long <- add_baseline_covariates(
            split$dat_long, id.name = base_cfg$id.name, tv_covs = tv_covs_all, suffix = "0"
        )
        split$dat_test <- add_baseline_covariates(
            split$dat_test, id.name = base_cfg$id.name, tv_covs = tv_covs_all, suffix = "0"
        )

        tab_r <- make_app_summary_tab(r, tau_grid)
        bounds_r <- vector("list", length(tau_grid))
        names(bounds_r) <- paste0("tau_", tau_grid)
        survC_r <- vector("list", length(tau_grid))
        names(survC_r) <- paste0("tau_", tau_grid)

        for (k in seq_along(tau_grid)) {
            tau <- tau_grid[k]
            tau_t0 <- proc.time()[["elapsed"]]

            fit <- tryCatch({
                dynamicCP_split(
                    dat = split$dat_long,
                    dat_test = split$dat_test,
                    start.name = base_cfg$start.name,
                    stop.name = base_cfg$stop.name,
                    event.name = base_cfg$event.name,
                    event.C.name = base_cfg$event.C.name,
                    id.name = base_cfg$id.name,
                    tau = tau,
                    model.pred = model.pred,
                    covname.pred.baseline = cov_pred_baseline_static,
                    covname.pred.timevarying = NULL,
                    model.C = model.C,
                    model.C.args = model.C.args,
                    covname.C.baseline = cov_C_baseline_static,
                    covname.C.timevarying = NULL,
                    TT.name = NULL,
                    seed = seeds$seeds_split[r],
                    train_ratio = base_cfg$train_ratio,
                    theta_grid = base_cfg$theta_grid,
                    alpha = alpha,
                    trim.C = base_cfg$trim.C
                )
            }, error = function(e) {
                list(error = TRUE, msg = conditionMessage(e))
            })

            runtime_k <- proc.time()[["elapsed"]] - tau_t0
            tab_r$runtime_sec[k] <- runtime_k

            if (isTRUE(fit$error)) {
                tab_r$msg[k] <- fit$msg
                bounds_r[[k]] <- empty_app_bounds()
                next
            }

            dat_cal_tau <- prepare_data_tau(
                split$dat_long[!(split$dat_long[[base_cfg$id.name]] %in% fit$id_tr), , drop = FALSE],
                base_cfg$start.name,
                base_cfg$stop.name,
                base_cfg$event.name,
                base_cfg$id.name,
                tau
            )
            dat_test_tau <- prepare_data_tau(
                split$dat_test,
                base_cfg$start.name,
                base_cfg$stop.name,
                base_cfg$event.name,
                base_cfg$id.name,
                tau
            )

            bc <- make_app_bounds_from_fit(
                dat_cal_tau, fit$lower_IPCW_cal, fit$upper_IPCW_cal, fit$survC_cal, r, tau, base_cfg
            )
            bt <- make_app_bounds_from_fit(
                dat_test_tau, fit$lower_IPCW_test, fit$upper_IPCW_test, fit$survC_test, r, tau, base_cfg
            )
            ss <- summarize_app_bounds(bt, bc)

            tab_r$ok[k] <- TRUE
            tab_r$theta_hat[k] <- fit$theta_hat_IPCW
            tab_r[k, names(ss)] <- ss
            bounds_r[[k]] <- list(test = bt, cal = bc)
            survC_r[[k]] <- list(survC_cal = fit$survC_cal, survC_test = fit$survC_test)
        }

        results_summary_source[[r]] <- tab_r
        results_bounds_source[[r]] <- bounds_r
        results_survC_source[[r]] <- survC_r
    }
    cat("\n")

    save_app_generic_bundle(
        cfg0 = cfg0,
        base_cfg = base_cfg,
        seeds = seeds,
        results_summary = results_summary_source,
        results_bounds = results_bounds_source,
        results_survC = results_survC_source,
        method = "dCP0_IPCW0_source",
        label = "dCP0-IPCW0 source",
        R = R,
        fit_tau = "all_tau_source_for_static_tau0_bounds",
        extra_config = list(
            source_file = base$path,
            covname.pred.baseline = cov_pred_baseline_static,
            covname.pred.timevarying = NULL,
            covname.C.baseline = cov_C_baseline_static,
            covname.C.timevarying = NULL
        )
    )

    stat <- build_static_from_dynamic(
        results_bounds = results_bounds_source,
        tau_grid = tau_grid,
        TT.name = NULL
    )
    res_static <- do.call(rbind, stat$results_summary_static)
    source_static <- list(
        results_summary = stat$results_summary_static,
        results_bounds = stat$results_bounds_static,
        results_survC = vector("list", length(stat$results_bounds_static)),
        res = res_static,
        res_ok = res_static[res_static$ok %in% TRUE, , drop = FALSE]
    )

    save_app_generic_bundle(
        cfg0 = cfg0,
        base_cfg = base_cfg,
        seeds = seeds,
        results_summary = source_static$results_summary,
        results_bounds = source_static$results_bounds,
        results_survC = source_static$results_survC,
        method = "static_IPCW",
        label = "IPCW",
        R = R,
        fit_tau = 0,
        extra_config = list(
            source_file = base$path,
            source_method = "dCP0_IPCW0_source",
            covname.pred.baseline = cov_pred_baseline_static,
            covname.pred.timevarying = NULL,
            covname.C.baseline = cov_C_baseline_static,
            covname.C.timevarying = NULL
        )
    )

    list(
        cfg0 = cfg0,
        base_cfg = base_cfg,
        seeds = seeds,
        R = R,
        source_static = source_static,
        cov_pred_baseline_static = cov_pred_baseline_static,
        cov_C_baseline_static = cov_C_baseline_static
    )
}

derive_app_static_trunc_generic <- function(source) {
    cfg0 <- source$cfg0
    base_cfg <- source$base_cfg
    seeds <- source$seeds
    R <- source$R
    src <- source$source_static
    tau_grid <- cfg0$tau_list_plot

    results_summary <- vector("list", R)
    results_bounds <- vector("list", R)
    results_survC <- vector("list", R)

    for (r in seq_len(R)) {
        tab_r <- src$res[src$res$rep == r, , drop = FALSE]
        bounds_r <- vector("list", length(tau_grid))
        names(bounds_r) <- paste0("tau_", tau_grid)
        survC_r <- vector("list", length(tau_grid))
        names(survC_r) <- paste0("tau_", tau_grid)

        for (k in seq_along(tau_grid)) {
            tau <- tau_grid[k]
            out_k <- tryCatch({
                bk <- src$results_bounds[[r]][[k]]
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
                ss <- summarize_app_bounds(bt, bc)
                list(ok = TRUE, test = bt, cal = bc, summary = ss)
            }, error = function(e) {
                list(ok = FALSE, msg = conditionMessage(e))
            })

            if (isTRUE(out_k$ok)) {
                tab_r$ok[tab_r$tau == tau] <- TRUE
                tab_r[tab_r$tau == tau, names(out_k$summary)] <- out_k$summary
                bounds_r[[k]] <- list(test = out_k$test, cal = out_k$cal)
                survC_r[[k]] <- list(survC_cal = out_k$cal$w_ipcw, survC_test = out_k$test$w_ipcw)
            } else {
                tab_r$ok[tab_r$tau == tau] <- FALSE
                tab_r$msg[tab_r$tau == tau] <- out_k$msg
                bounds_r[[k]] <- empty_app_bounds()
            }
        }

        results_summary[[r]] <- tab_r
        results_bounds[[r]] <- bounds_r
        results_survC[[r]] <- survC_r
    }

    save_app_generic_bundle(
        cfg0 = cfg0,
        base_cfg = base_cfg,
        seeds = seeds,
        results_summary = results_summary,
        results_bounds = results_bounds,
        results_survC = results_survC,
        method = "static_IPCW_trunc",
        label = "IPCW-trunc",
        R = R,
        fit_tau = 0,
        extra_config = list(
            source_method = "static_IPCW",
            covname.pred.baseline = source$cov_pred_baseline_static,
            covname.pred.timevarying = NULL,
            covname.C.baseline = source$cov_C_baseline_static,
            covname.C.timevarying = NULL
        )
    )
}

run_app_lm_ipcw_generic <- function(dataset_name, alpha, source) {
    cfg0 <- source$cfg0
    base_cfg <- source$base_cfg
    seeds <- source$seeds
    R <- source$R
    tau_grid <- cfg0$tau_list_plot
    dat0 <- prepare_app_data(dataset_name)

    tv_covs_all <- union(base_cfg$covname.pred.timevarying, base_cfg$covname.C.timevarying)
    cov_pred_baseline_lm <- c(
        base_cfg$covname.pred.baseline,
        paste0(base_cfg$covname.pred.timevarying, "0")
    )
    cov_C_baseline_lm <- c(
        base_cfg$covname.C.baseline,
        paste0(base_cfg$covname.C.timevarying, "0")
    )

    results_summary <- vector("list", R)
    results_bounds <- vector("list", R)
    results_survC <- vector("list", R)

    cat("Running IPCW-LM for ", dataset_name, " / ", analysis_tag, "\n", sep = "")
    cat("IPCW-LM uses residual time after each landmark and baseline/time-zero covariates.\n")
    cat("Replications: ")

    for (r in seq_len(R)) {
        cat(r, "..")
        split <- split_app_data(dat0, base_cfg, seeds, r)
        split$dat_long <- add_baseline_covariates(
            split$dat_long, id.name = base_cfg$id.name, tv_covs = tv_covs_all, suffix = "0"
        )
        split$dat_test <- add_baseline_covariates(
            split$dat_test, id.name = base_cfg$id.name, tv_covs = tv_covs_all, suffix = "0"
        )

        tab_r <- make_app_summary_tab(r, tau_grid)
        bounds_r <- vector("list", length(tau_grid))
        names(bounds_r) <- paste0("tau_", tau_grid)
        survC_r <- vector("list", length(tau_grid))
        names(survC_r) <- paste0("tau_", tau_grid)

        for (k in seq_along(tau_grid)) {
            tau <- tau_grid[k]
            tau_t0 <- proc.time()[["elapsed"]]

            out_k <- tryCatch({
                ids_long <- subject_ids_at_risk_app(split$dat_long, tau, base_cfg)
                ids_test <- subject_ids_at_risk_app(split$dat_test, tau, base_cfg)
                dat_land <- split$dat_long[
                    as.character(split$dat_long[[base_cfg$id.name]]) %in% ids_long, , drop = FALSE
                ]
                dat_test_land <- split$dat_test[
                    as.character(split$dat_test[[base_cfg$id.name]]) %in% ids_test, , drop = FALSE
                ]

                dat_land_res <- make_landmark_residual_data_app(dat_land, tau, base_cfg)
                dat_test_land_res <- make_landmark_residual_data_app(dat_test_land, tau, base_cfg)

                fit <- dynamicCP_split(
                    dat = dat_land_res,
                    dat_test = dat_test_land_res,
                    start.name = base_cfg$start.name,
                    stop.name = base_cfg$stop.name,
                    event.name = base_cfg$event.name,
                    event.C.name = base_cfg$event.C.name,
                    id.name = base_cfg$id.name,
                    tau = 0,
                    model.pred = model.pred,
                    covname.pred.baseline = cov_pred_baseline_lm,
                    covname.pred.timevarying = NULL,
                    model.C = model.C,
                    model.C.args = model.C.args,
                    covname.C.baseline = cov_C_baseline_lm,
                    covname.C.timevarying = NULL,
                    TT.name = NULL,
                    seed = seeds$seeds_split[r],
                    train_ratio = base_cfg$train_ratio,
                    theta_grid = base_cfg$theta_grid,
                    alpha = alpha,
                    trim.C = base_cfg$trim.C
                )

                dat_cal_land <- dat_land[!(dat_land[[base_cfg$id.name]] %in% fit$id_tr), , drop = FALSE]
                dat_cal_tau <- prepare_data_tau(
                    dat_cal_land,
                    base_cfg$start.name,
                    base_cfg$stop.name,
                    base_cfg$event.name,
                    base_cfg$id.name,
                    tau
                )
                dat_test_tau <- prepare_data_tau(
                    dat_test_land,
                    base_cfg$start.name,
                    base_cfg$stop.name,
                    base_cfg$event.name,
                    base_cfg$id.name,
                    tau
                )

                lower_cal <- fit$lower_IPCW_cal + tau
                upper_cal <- fit$upper_IPCW_cal + tau
                lower_test <- fit$lower_IPCW_test + tau
                upper_test <- fit$upper_IPCW_test + tau

                bc <- make_app_bounds_from_fit(dat_cal_tau, lower_cal, upper_cal, fit$survC_cal, r, tau, base_cfg)
                bt <- make_app_bounds_from_fit(dat_test_tau, lower_test, upper_test, fit$survC_test, r, tau, base_cfg)
                ss <- summarize_app_bounds(bt, bc)
                list(ok = TRUE, test = bt, cal = bc, summary = ss, theta = fit$theta_hat_IPCW)
            }, error = function(e) {
                list(ok = FALSE, msg = conditionMessage(e))
            })

            tab_r$runtime_sec[k] <- proc.time()[["elapsed"]] - tau_t0
            if (isTRUE(out_k$ok)) {
                tab_r$ok[k] <- TRUE
                tab_r$theta_hat[k] <- out_k$theta
                tab_r[k, names(out_k$summary)] <- out_k$summary
                bounds_r[[k]] <- list(test = out_k$test, cal = out_k$cal)
                survC_r[[k]] <- list(survC_cal = out_k$cal$w_ipcw, survC_test = out_k$test$w_ipcw)
            } else {
                tab_r$ok[k] <- FALSE
                tab_r$msg[k] <- out_k$msg
                bounds_r[[k]] <- empty_app_bounds()
            }
        }

        results_summary[[r]] <- tab_r
        results_bounds[[r]] <- bounds_r
        results_survC[[r]] <- survC_r
    }
    cat("\n")

    save_app_generic_bundle(
        cfg0 = cfg0,
        base_cfg = base_cfg,
        seeds = seeds,
        results_summary = results_summary,
        results_bounds = results_bounds,
        results_survC = results_survC,
        method = "static_LM_IPCW",
        label = "IPCW-LM",
        R = R,
        fit_tau = "landmark_residual_time",
        extra_config = list(
            covname.pred.baseline = cov_pred_baseline_lm,
            covname.pred.timevarying = NULL,
            covname.C.baseline = cov_C_baseline_lm,
            covname.C.timevarying = NULL
        )
    )
}

cat("Generic application workflow\n")
cat("Datasets: ", paste(datasets, collapse = ", "), "\n", sep = "")
cat("Alpha: ", alpha, "\n", sep = "")
cat("R limit: ", R_limit, "\n", sep = "")
cat("Prediction model: ", model.pred, "\n", sep = "")
cat("G model: ", model.C, "\n", sep = "")
cat("Analysis tag: ", analysis_tag, "\n", sep = "")
cat("Steps: ", paste(steps, collapse = ", "), "\n", sep = "")
cat("Overwrite existing outputs: ", overwrite, "\n", sep = "")

for (dataset_name in datasets) {
    if ("dynamic" %in% steps) {
        run_app_dynamic_generic(dataset_name, alpha = alpha, R_limit = R_limit)
    }
    if ("static" %in% steps) {
        src <- run_app_static_source_generic(dataset_name, alpha = alpha, R_limit = R_limit)
        derive_app_static_trunc_generic(src)
        run_app_lm_ipcw_generic(dataset_name, alpha = alpha, source = src)
    }
}
