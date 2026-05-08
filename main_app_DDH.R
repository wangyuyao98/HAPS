rm(list = ls())

# Run this script from the repository root.
if (!file.exists("src/app_analysis_helpers.R")) {
    stop("Please run main_app_DDH.R from the repository root directory.")
}

# Keep reticulate/PyTorch conservative on laptops, matching the DDH simulation runner.
thread_env_defaults <- c(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    KMP_DUPLICATE_LIB_OK = "TRUE"
)
for (env_name in names(thread_env_defaults)) {
    if (!nzchar(Sys.getenv(env_name, unset = ""))) {
        do.call(Sys.setenv, as.list(setNames(thread_env_defaults[[env_name]], env_name)))
    }
}

source("src/dynamicCP.R")
source("src/helpers.R")
source("src/ddh_tuning_presets.R")
source("src/ddh_bridge.R")
source("src/app_analysis_helpers.R")

library(survival)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && args[[1L]] %in% c("list", "presets")) {
    print(ddh_tuning_preset_table(), row.names = FALSE)
    quit(save = "no", status = 0)
}

dataset_arg <- if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else "all"
R_limit_arg <- if (length(args) >= 2L && nzchar(args[[2L]])) as.integer(args[[2L]]) else 200L
ddh_arg <- if (length(args) >= 3L && nzchar(args[[3L]])) args[[3L]] else "20"
time_feature_arg <- if (length(args) >= 4L && nzchar(args[[4L]])) args[[4L]] else "time"
model.C.arg <- if (length(args) >= 5L && nzchar(args[[5L]])) args[[5L]] else "cox"
checkpoint_every <- if (length(args) >= 6L && nzchar(args[[6L]])) as.integer(args[[6L]]) else 0L

ddh_tuning_preset <- NULL
ddh_tuning_info <- NULL
ddh_epochs <- suppressWarnings(as.integer(ddh_arg))
if (is.na(ddh_epochs)) {
    ddh_tuning_preset <- ddh_arg
    ddh_tuning_info <- ddh_get_tuning_preset(ddh_tuning_preset)
    ddh_epochs <- as.integer(ddh_tuning_info$args$epochs)
}

valid_datasets <- c("all", "pbcseq", "colon")
if (!dataset_arg %in% valid_datasets) {
    stop("Unknown dataset: ", dataset_arg, ". Choose one of: ", paste(valid_datasets, collapse = ", "))
}
if (!is.finite(R_limit_arg) || R_limit_arg < 1L) stop("R must be a positive integer.")
if (!is.finite(ddh_epochs) || ddh_epochs < 1L) stop("ddh_epochs must be a positive integer.")
if (!is.finite(checkpoint_every) || checkpoint_every < 0L) {
    stop("checkpoint_every must be a non-negative integer.")
}
checkpoint_every <- as.integer(checkpoint_every)
time_feature_arg <- match.arg(
    time_feature_arg,
    c("time", "notime", "TRUE", "FALSE", "true", "false")
)
include_time_features <- time_feature_arg %in% c("time", "TRUE", "true")
time_feature_tag <- if (include_time_features) "time" else "notime"
model.C <- match.arg(model.C.arg, c("cox", "rsf"))

if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for model.pred='ddh'.")
}
if (!reticulate::py_module_available("torch")) {
    stop(
        "Python module 'torch' is not available in the reticulate Python environment.\n",
        "Run with an environment that has torch, e.g.\n",
        "  RETICULATE_PYTHON_ENV=r-ddh Rscript main_app_DDH.R pbcseq 5 wide64_epochs50 time cox"
    )
}
torch_mod <- reticulate::import("torch", delay_load = FALSE)
torch_threads <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1")))
if (is.na(torch_threads) || torch_threads < 1L) torch_threads <- 1L
torch_mod$set_num_threads(torch_threads)

datasets <- if (identical(dataset_arg, "all")) c("pbcseq", "colon") else dataset_arg

alpha <- 0.2
model.pred <- "ddh"
TT.name <- NULL

model.pred.args.base <- list(
    ddh_bridge_path = "src/ddh_bridge.py",
    backend = "external",
    external_module = "ddh_backend",
    external_module_path = "src/ddh_backend.py",
    external_function = "fit_predict",
    seed = 123L,
    epochs = ddh_epochs,
    hidden_size = 32L,
    num_layers = 1L,
    dropout = 0.1,
    lr = 1e-3,
    weight_decay = 1e-5,
    batch_size = 64L,
    num_time_bins = 40L,
    include_time_features = include_time_features,
    use_gpu = FALSE
)
if (is.null(ddh_tuning_preset)) {
    model.pred.args.base$epochs <- ddh_epochs
} else {
    model.pred.args.base <- ddh_apply_tuning_preset(model.pred.args.base, ddh_tuning_preset)
}
ddh_epochs <- as.integer(model.pred.args.base$epochs)

model.C.args <- if (identical(model.C, "rsf")) {
    list(ntree = 500L, nodesize = 15L, nsplit = 10L, forest = TRUE, rsf_min_censor_events = 10L)
} else {
    list()
}

app_ddh_tune_part <- function(ddh_tuning_preset) {
    if (is.null(ddh_tuning_preset)) return("")
    paste0("_tune", ddh_tuning_tag(ddh_tuning_preset))
}

app_ddh_file <- function(cfg, method, alpha, R, ddh_epochs, time_feature_tag,
                         model.C, ddh_tuning_preset = NULL) {
    file.path(
        app_out_dir(cfg),
        paste0(
            method,
            "_ddh_", model.C,
            "_", time_feature_tag,
            app_ddh_tune_part(ddh_tuning_preset),
            "_R", as.integer(R),
            "_n", cfg$n_all,
            "_alpha", alpha,
            "_epochs", as.integer(ddh_epochs),
            ".rds"
        )
    )
}

app_ddh_checkpoint_file <- function(cfg, R, ddh_epochs, time_feature_tag,
                                    model.C, ddh_tuning_preset = NULL) {
    file.path(
        app_out_dir(cfg),
        paste0(
            "progress_app_DDH_ddh_", model.C,
            "_", time_feature_tag,
            app_ddh_tune_part(ddh_tuning_preset),
            "_R", as.integer(R),
            "_n", cfg$n_all,
            "_alpha", alpha,
            "_epochs", as.integer(ddh_epochs),
            ".csv"
        )
    )
}

new_app_store <- function(R) {
    list(
        results_summary = vector("list", R),
        results_bounds = vector("list", R),
        results_survC = vector("list", R),
        rep_runtime_sec = rep(NA_real_, R)
    )
}

summarise_app_ddh_progress <- function(stores, method_labels, rep_done, elapsed_min) {
    rows <- lapply(names(stores), function(method) {
        tabs <- stores[[method]]$results_summary[seq_len(rep_done)]
        tabs <- tabs[!vapply(tabs, is.null, logical(1))]
        if (length(tabs) == 0L) return(NULL)

        df <- do.call(rbind, tabs)
        out <- dplyr::summarise(
            dplyr::group_by(df, tau),
            ok_n = sum(ok %in% TRUE),
            fail_n = sum(!(ok %in% TRUE)),
            cov_ipcw = mean(coverage_test_ipcw, na.rm = TRUE),
            mean_len = mean(mean_len_test, na.rm = TRUE),
            median_len = stats::median(mean_len_test, na.rm = TRUE),
            theta = mean(theta_hat, na.rm = TRUE),
            .groups = "drop"
        )
        out$method <- method
        out$label <- method_labels[[method]]
        out$rep_done <- rep_done
        out$elapsed_min <- elapsed_min
        out[, c(
            "rep_done", "elapsed_min", "method", "label", "tau",
            "ok_n", "fail_n", "cov_ipcw", "mean_len", "median_len", "theta"
        )]
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0L) return(NULL)
    do.call(rbind, rows)
}

write_app_ddh_checkpoint <- function(cfg0, stores, method_labels, R, ddh_epochs,
                                     time_feature_tag, model.C, ddh_tuning_preset,
                                     rep_done, elapsed_min) {
    progress <- summarise_app_ddh_progress(stores, method_labels, rep_done, elapsed_min)
    if (is.null(progress)) return(invisible(NULL))

    out_file <- app_ddh_checkpoint_file(
        cfg = cfg0,
        R = R,
        ddh_epochs = ddh_epochs,
        time_feature_tag = time_feature_tag,
        model.C = model.C,
        ddh_tuning_preset = ddh_tuning_preset
    )
    write.csv(progress, out_file, row.names = FALSE)
    cat("\nCheckpoint saved to:\n  ", out_file, "\n", sep = "")
    invisible(out_file)
}

save_app_ddh_bundle <- function(cfg0, base_cfg, seeds, store, method, label,
                                tau_grid, R, ddh_epochs, model.pred.args) {
    res <- do.call(rbind, store$results_summary)
    out_file <- app_ddh_file(
        cfg = cfg0,
        method = method,
        alpha = alpha,
        R = R,
        ddh_epochs = ddh_epochs,
        time_feature_tag = time_feature_tag,
        model.C = model.C,
        ddh_tuning_preset = ddh_tuning_preset
    )
    save_obj <- list(
        config = list(
            setup = cfg0$out_setup,
            method = method,
            label = label,
            source_script = "main_app_DDH.R",
            source_file = "generated_from_survival_package_data",
            dataset_name = cfg0$dataset_name,
            n_all = cfg0$n_all,
            tau_grid = tau_grid,
            model.pred = model.pred,
            model.pred.args = model.pred.args,
            ddh_tuning_preset = ddh_tuning_preset,
            ddh_tuning_description = if (!is.null(ddh_tuning_info)) ddh_tuning_info$description else NULL,
            time_feature_tag = time_feature_tag,
            model.C = model.C,
            model.C.args = model.C.args,
            covname.pred.baseline = base_cfg$covname.pred.baseline,
            covname.pred.timevarying = base_cfg$covname.pred.timevarying,
            covname.C.baseline = base_cfg$covname.C.baseline,
            covname.C.timevarying = base_cfg$covname.C.timevarying,
            theta_grid = base_cfg$theta_grid,
            alpha = alpha,
            trim.C = base_cfg$trim.C,
            train_ratio = base_cfg$train_ratio,
            test_ratio = base_cfg$test_ratio,
            TT.name = TT.name,
            start.name = base_cfg$start.name,
            stop.name = base_cfg$stop.name,
            event.name = base_cfg$event.name,
            event.C.name = base_cfg$event.C.name,
            id.name = base_cfg$id.name,
            checkpoint_every = checkpoint_every
        ),
        seeds = seeds,
        results = list(
            results_summary = store$results_summary,
            results_bounds = store$results_bounds,
            results_survC = store$results_survC,
            res = res,
            res_ok = res[res$ok %in% TRUE, , drop = FALSE]
        ),
        meta = list(
            timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
            rep_runtime_sec = store$rep_runtime_sec
        )
    )
    saveRDS(save_obj, out_file)
    cat("Saved ", label, " to:\n  ", out_file, "\n", sep = "")
    invisible(out_file)
}

run_one_app_dataset_ddh <- function(dataset_name) {
    cfg0 <- app_dataset_config(dataset_name, alpha)

    base <- app_base_run(cfg0, alpha = alpha, R = R_limit_arg)
    base_cfg <- base$config
    seeds <- base$seeds
    R_full <- base$R_full
    R <- as.integer(R_limit_arg)
    tau_grid <- cfg0$tau_list_plot

    dat0 <- prepare_app_data(dataset_name)
    id_all <- unique(dat0[[base_cfg$id.name]])
    n_all <- length(id_all)

    stores <- list(
        dCP_IPCW = new_app_store(R),
        derived_DynaCoPS_DR = new_app_store(R)
    )
    method_labels <- c(
        dCP_IPCW = "HAPS-DDH",
        derived_DynaCoPS_DR = "HAPS-DDH-DR"
    )

    cat("\nApplication DDH workflow for ", dataset_name, "\n", sep = "")
    cat("Replications: ", R, "/", R_full, "\n", sep = "")
    cat("Tau grid: ", paste(tau_grid, collapse = ", "), "\n", sep = "")
    cat("DDH epochs: ", ddh_epochs, "\n", sep = "")
    cat("DDH tuning preset: ", ddh_tuning_preset %||% "manual", "\n", sep = "")
    cat("DDH time features: ", time_feature_tag, "\n", sep = "")
    cat("G model: ", model.C, "\n", sep = "")
    if (checkpoint_every > 0L) cat("Checkpoint every: ", checkpoint_every, " reps\n", sep = "")
    cat("Python: ", reticulate::py_config()$python, "\n", sep = "")
    cat("Replications: ")

    run_t0 <- proc.time()[["elapsed"]]
    for (r in seq_len(R)) {
        cat(r, "..")
        rep_t0 <- proc.time()[["elapsed"]]

        set.seed(seeds$seeds_test[r])
        id_test <- sample(id_all, floor(n_all * base_cfg$test_ratio), replace = FALSE)
        dat_test_r <- dat0[dat0[[base_cfg$id.name]] %in% id_test, , drop = FALSE]
        dat_long_r <- dat0[!(dat0[[base_cfg$id.name]] %in% id_test), , drop = FALSE]

        rep_tabs <- lapply(seq_along(stores), function(...) {
            data.frame(
                rep = r,
                tau = tau_grid,
                ok = FALSE,
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
        })
        names(rep_tabs) <- names(stores)
        rep_bounds <- lapply(seq_along(stores), function(...) {
            x <- vector("list", length(tau_grid))
            names(x) <- paste0("tau_", tau_grid)
            x
        })
        names(rep_bounds) <- names(stores)
        rep_survC <- rep_bounds

        for (k in seq_along(tau_grid)) {
            tau <- tau_grid[k]
            tau_t0 <- proc.time()[["elapsed"]]
            ddh_args_k <- model.pred.args.base
            ddh_args_k$seed <- as.integer(seeds$seeds_split[r] + 1000L * k)

            fit <- tryCatch({
                dynamicCP_split(
                    dat = dat_long_r,
                    dat_test = dat_test_r,
                    start.name = base_cfg$start.name,
                    stop.name = base_cfg$stop.name,
                    event.name = base_cfg$event.name,
                    event.C.name = base_cfg$event.C.name,
                    id.name = base_cfg$id.name,
                    tau = tau,
                    model.pred = model.pred,
                    model.pred.args = ddh_args_k,
                    covname.pred.baseline = base_cfg$covname.pred.baseline,
                    covname.pred.timevarying = base_cfg$covname.pred.timevarying,
                    model.C = model.C,
                    model.C.args = model.C.args,
                    covname.C.baseline = base_cfg$covname.C.baseline,
                    covname.C.timevarying = base_cfg$covname.C.timevarying,
                    TT.name = TT.name,
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
                for (method in names(stores)) {
                    rep_tabs[[method]]$msg[k] <- fit$msg
                    rep_tabs[[method]]$runtime_sec[k] <- runtime_k
                    rep_bounds[[method]][[k]] <- empty_app_bounds()
                }
                next
            }

            dat_cal_tau <- prepare_data_tau(
                dat_long_r[!(dat_long_r[[base_cfg$id.name]] %in% fit$id_tr), , drop = FALSE],
                base_cfg$start.name,
                base_cfg$stop.name,
                base_cfg$event.name,
                base_cfg$id.name,
                tau
            )
            dat_test_tau <- prepare_data_tau(
                dat_test_r,
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
                rep_survC[[method]][[k]] <- list(survC_cal = fit$survC_cal, survC_test = fit$survC_test)
            }
        }

        for (method in names(stores)) {
            stores[[method]]$results_summary[[r]] <- rep_tabs[[method]]
            stores[[method]]$results_bounds[[r]] <- rep_bounds[[method]]
            stores[[method]]$results_survC[[r]] <- rep_survC[[method]]
            stores[[method]]$rep_runtime_sec[r] <- proc.time()[["elapsed"]] - rep_t0
        }
        if (checkpoint_every > 0L && (r %% checkpoint_every == 0L || r == R)) {
            write_app_ddh_checkpoint(
                cfg0 = cfg0,
                stores = stores,
                method_labels = method_labels,
                R = R,
                ddh_epochs = ddh_epochs,
                time_feature_tag = time_feature_tag,
                model.C = model.C,
                ddh_tuning_preset = ddh_tuning_preset,
                rep_done = r,
                elapsed_min = (proc.time()[["elapsed"]] - run_t0) / 60
            )
        }
    }
    cat("\n")

    out_files <- character(0)
    for (method in names(stores)) {
        out_files <- c(
            out_files,
            save_app_ddh_bundle(
                cfg0 = cfg0,
                base_cfg = base_cfg,
                seeds = seeds,
                store = stores[[method]],
                method = method,
                label = method_labels[[method]],
                tau_grid = tau_grid,
                R = R,
                ddh_epochs = ddh_epochs,
                model.pred.args = model.pred.args.base
            )
        )
    }

    for (method in names(stores)) {
        res_ok <- do.call(rbind, stores[[method]]$results_summary)
        res_ok <- res_ok[res_ok$ok %in% TRUE, , drop = FALSE]
        cat("\nSummary for ", method_labels[[method]], ":\n", sep = "")
        print(
            dplyr::summarise(
                dplyr::group_by(res_ok, tau),
                cov_ipcw = mean(coverage_test_ipcw, na.rm = TRUE),
                mean_len = mean(mean_len_test, na.rm = TRUE),
                theta = mean(theta_hat, na.rm = TRUE),
                .groups = "drop"
            ),
            n = Inf
        )
    }

    invisible(out_files)
}

cat("Application DDH workflow\n")
cat("Datasets: ", paste(datasets, collapse = ", "), "\n", sep = "")
cat("Alpha: ", alpha, "\n", sep = "")
cat("DDH epochs: ", ddh_epochs, "\n", sep = "")
cat("DDH tuning preset: ", ddh_tuning_preset %||% "manual", "\n", sep = "")
cat("DDH time features: ", time_feature_tag, "\n", sep = "")
cat("G model: ", model.C, "\n", sep = "")
if (checkpoint_every > 0L) cat("Checkpoint every: ", checkpoint_every, " reps\n", sep = "")

for (dataset_name in datasets) {
    run_one_app_dataset_ddh(dataset_name)
}
