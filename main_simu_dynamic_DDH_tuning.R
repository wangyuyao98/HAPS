rm(list = ls())

# Separate DDH tuning runner. This intentionally writes to results/<setup>/ddh_tuning/
# so exploratory fits do not overwrite or get mixed with the main DDH simulations.
if (!file.exists("src/gen_ICML_simu.R")) {
    stop("Please run main_simu_dynamic_DDH_tuning.R from the repository root directory.")
}

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
if (!nzchar(Sys.getenv("RETICULATE_PYTHON_ENV", unset = ""))) {
    Sys.setenv(RETICULATE_PYTHON_ENV = "r-ddh")
}

source("src/gen_ICML_simu.R")
source("src/helpers.R")
source("src/helpers_AIPCW.R")
source("src/ddh_tuning_presets.R")
source("src/ddh_bridge.R")
source("src/dynamicCP_AIPCW.R")

library(survival)
library(dplyr)

if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required for model.Xi='xgb_multiclass'.")
}
if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for model.pred='ddh'.")
}
if (!reticulate::py_module_available("torch")) {
    stop(
        "Python module 'torch' is not available in the reticulate Python environment.\n",
        "Create/select an environment with torch before running this script."
    )
}
torch_mod <- reticulate::import("torch", delay_load = FALSE)
torch_threads <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "1")))
if (is.na(torch_threads) || torch_threads < 1L) torch_threads <- 1L
torch_mod$set_num_threads(torch_threads)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && identical(args[[1L]], "list")) {
    print(ddh_tuning_preset_table(), row.names = FALSE)
    quit(save = "no")
}

dgm_arg <- if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else "mixC2_DGM1"
preset_name <- if (length(args) >= 2L && nzchar(args[[2L]])) args[[2L]] else "wide64_epochs50"
R <- if (length(args) >= 3L && nzchar(args[[3L]])) as.integer(args[[3L]]) else 20L
n <- if (length(args) >= 4L && nzchar(args[[4L]])) as.integer(args[[4L]]) else 1000L
n_test <- if (length(args) >= 5L && nzchar(args[[5L]])) as.integer(args[[5L]]) else 500L
model.C.arg <- if (length(args) >= 6L && nzchar(args[[6L]])) args[[6L]] else "rsf"
model.S.arg <- if (length(args) >= 7L && nzchar(args[[7L]])) args[[7L]] else "xgb_cox"
model.Xi.arg <- if (length(args) >= 8L && nzchar(args[[8L]])) args[[8L]] else "xgb_multiclass"
checkpoint_every <- if (length(args) >= 9L && nzchar(args[[9L]])) as.integer(args[[9L]]) else 1L

valid_dgms <- c("linWB1", "mixC2_DGM1", "mixC2_DGM3")
if (!dgm_arg %in% valid_dgms) {
    stop("Unknown DGM: ", dgm_arg, ". Choose one of: ", paste(valid_dgms, collapse = ", "))
}
if (!is.finite(R) || R < 1L) stop("R must be a positive integer.")
if (!is.finite(n) || n < 10L) stop("n must be at least 10.")
if (!is.finite(n_test) || n_test < 10L) stop("n_test must be at least 10.")
if (!is.finite(checkpoint_every) || checkpoint_every < 1L) {
    stop("checkpoint_every must be a positive integer.")
}

rho <- 0.3
alpha <- 0.1
trim.C <- 0.05
train_ratio <- 0.5
change_times <- c(3, 6)
tau_grid <- c(0, change_times)
theta_grid_AIPCW <- seq(0, 0.5, by = 0.01)

start.name <- "start"
stop.name <- "stop"
event.name <- "event"
event.C.name <- "event.C"
id.name <- "id"
TT.name <- "TT"

model.pred <- "ddh"
model.C <- match.arg(model.C.arg, c("km", "cox", "rsf", "hal", "grf", "xgb_cox", "xgb_aft"))
model.S <- match.arg(model.S.arg, c("km", "cox", "rsf", "hal", "grf", "xgb_cox", "xgb_aft"))
model.Xi <- match.arg(model.Xi.arg, c("lm", "rsf", "xgb_reg", "xgb_multiclass", "hal_reg"))
localization <- FALSE
Gtau_mode <- "one"
Gtau_delta <- 0
Gtau_delta_tag <- "0"
censoring_method_AIPCW <- "piecewise"

if (identical(dgm_arg, "linWB1")) {
    setup <- "linWB1"
    dgm_name <- "linear_weibull"
    dgm_family <- NULL
    dgm_variant <- NULL
    covname.pred.baseline <- NULL
    covname.pred.timevarying <- c("L")
    covname.C.baseline <- NULL
    covname.C.timevarying <- c("L")
    dgm_parK <- list(
        shape_T = c(4, 4, 5),
        shape_C = c(3, 3, 3),
        beta_T0 = c(-8, -8, -5),
        beta_TL = c(1, 2, 3),
        beta_C0 = c(-6, -6, -5),
        beta_CL = c(2, 2, 2)
    )
} else if (startsWith(dgm_arg, "mixC2_DGM")) {
    source("src/mixC_dgm_registry.R")
    setup <- dgm_arg
    dgm_name <- "mixC_uniform_cox"
    dgm_family <- "mixC2"
    dgm_variant <- sub("^mixC2_", "", dgm_arg)
    dgm_variant <- resolve_mixC_dgm_variant(dgm_variant, dgm_family = dgm_family)
    covname.pred.baseline <- c("B")
    covname.pred.timevarying <- c("L")
    covname.C.baseline <- c("B")
    covname.C.timevarying <- c("L")
    dgm_parK <- get_mixC_dgm_par(dgm_variant, dgm_family = dgm_family)
} else {
    stop("Unsupported DGM dispatch for: ", dgm_arg)
}

folder <- file.path("results", setup, "ddh_tuning")
if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)

base_model.pred.args <- list(
    ddh_bridge_path = "src/ddh_bridge.py",
    backend = "external",
    external_module = "ddh_backend",
    external_module_path = "src/ddh_backend.py",
    external_function = "fit_predict",
    seed = 123L,
    epochs = 20L,
    hidden_size = 32L,
    num_layers = 1L,
    dropout = 0.1,
    lr = 1e-3,
    weight_decay = 1e-5,
    batch_size = 64L,
    num_time_bins = 40L,
    include_time_features = FALSE,
    use_gpu = FALSE
)
model.pred.args <- ddh_apply_tuning_preset(base_model.pred.args, preset_name)
tuning_preset <- ddh_get_tuning_preset(preset_name)
tuning_tag <- ddh_tuning_tag(preset_name)
run_tag <- paste0(
    model.pred, "_", model.C,
    "_S", model.S, "_Xi", model.Xi,
    "_tune", tuning_tag,
    "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
    "_alpha", alpha,
    "_loc", as.integer(localization),
    "_Gtau", Gtau_mode, "_d", Gtau_delta_tag
)
progress_csv <- file.path(folder, paste0("progress_", run_tag, ".csv"))

empty_bounds <- function() {
    list(test = data.frame(), cal = data.frame())
}

make_bounds_from_fit <- function(dat_tau, lower, upper, survC, rep_id, tau) {
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

make_empty_tab <- function(rep_id) {
    data.frame(
        rep = rep_id,
        tau = tau_grid,
        ok = FALSE,
        msg = NA_character_,
        theta_hat = NA_real_,
        coverage_cal_ipcw = NA_real_,
        coverage_test_ipcw = NA_real_,
        coverage_cal_true = NA_real_,
        coverage_test_true = NA_real_,
        mean_len_cal = NA_real_,
        mean_len_test = NA_real_,
        q25_len_test = NA_real_,
        q50_len_test = NA_real_,
        q75_len_test = NA_real_,
        runtime_sec = NA_real_
    )
}

new_method_store <- function() {
    list(
        results_summary = vector("list", R),
        results_bounds = vector("list", R),
        results_survC = vector("list", R),
        rep_runtime_sec = rep(NA_real_, R)
    )
}

summarize_progress <- function(stores, rep_done, elapsed_sec) {
    rows <- lapply(names(stores), function(method) {
        completed <- stores[[method]]$results_summary[seq_len(rep_done)]
        has_result <- !vapply(completed, is.null, logical(1))
        completed <- completed[has_result]
        if (length(completed) == 0L) return(NULL)

        res <- do.call(rbind, completed)
        total_by_tau <- res |>
            dplyr::group_by(tau) |>
            dplyr::summarise(total_n = dplyr::n(), .groups = "drop")

        res_ok <- res[res$ok %in% TRUE, , drop = FALSE]
        if (nrow(res_ok) > 0L) {
            ok_sum <- res_ok |>
                dplyr::group_by(tau) |>
                dplyr::summarise(
                    ok_n = dplyr::n(),
                    coverage_test_true = mean(coverage_test_true, na.rm = TRUE),
                    coverage_test_ipcw = mean(coverage_test_ipcw, na.rm = TRUE),
                    mean_len_test = mean(mean_len_test, na.rm = TRUE),
                    median_len_test = stats::median(mean_len_test, na.rm = TRUE),
                    theta_hat = mean(theta_hat, na.rm = TRUE),
                    runtime_sec = mean(runtime_sec, na.rm = TRUE),
                    .groups = "drop"
                )
            out <- dplyr::left_join(total_by_tau, ok_sum, by = "tau")
            out$ok_n[is.na(out$ok_n)] <- 0L
        } else {
            out <- total_by_tau
            out$ok_n <- 0L
            out$coverage_test_true <- NA_real_
            out$coverage_test_ipcw <- NA_real_
            out$mean_len_test <- NA_real_
            out$median_len_test <- NA_real_
            out$theta_hat <- NA_real_
            out$runtime_sec <- NA_real_
        }

        out$fail_n <- out$total_n - out$ok_n
        data.frame(
            setup = setup,
            preset = preset_name,
            method = method,
            label = unname(method_labels[[method]]),
            model.C = model.C,
            model.S = model.S,
            model.Xi = model.Xi,
            R_target = R,
            rep_done = rep_done,
            elapsed_min = elapsed_sec / 60,
            out,
            stringsAsFactors = FALSE
        )
    })
    dplyr::bind_rows(rows)
}

write_progress_checkpoint <- function(stores, rep_done, elapsed_sec, quiet = FALSE) {
    progress <- summarize_progress(stores, rep_done = rep_done, elapsed_sec = elapsed_sec)
    utils::write.csv(progress, progress_csv, row.names = FALSE)
    if (!quiet) {
        cat("\nProgress checkpoint after rep ", rep_done, "/", R, ":\n", sep = "")
        show_methods <- c("DDH_model_quantile", "dCP_IPCW_new")
        show_cols <- c(
            "method", "tau", "ok_n", "fail_n",
            "coverage_test_true", "mean_len_test", "theta_hat"
        )
        print(
            progress[progress$method %in% show_methods, show_cols, drop = FALSE],
            row.names = FALSE
        )
        cat("Progress CSV:\n  ", progress_csv, "\n", sep = "")
    }
    invisible(progress)
}

save_method_bundle <- function(method, store, method_label) {
    res <- do.call(rbind, store$results_summary)
    res_ok <- res[res$ok %in% TRUE, , drop = FALSE]
    out_rds <- file.path(
        folder,
        paste0(
            method, "_", model.pred, "_", model.C,
            "_S", model.S, "_Xi", model.Xi,
            "_tune", tuning_tag,
            "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
            "_alpha", alpha,
            "_loc", as.integer(localization),
            "_Gtau", Gtau_mode, "_d", Gtau_delta_tag,
            ".rds"
        )
    )

    save_obj <- list(
        config = list(
            setup = setup,
            dgm_name = dgm_name,
            dgm_family = dgm_family,
            dgm_variant = dgm_variant,
            method = method,
            label = method_label,
            source_script = "main_simu_dynamic_DDH_tuning.R",
            tuning_preset = preset_name,
            tuning_preset_details = tuning_preset,
            tuning_tag = tuning_tag,
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
            model.pred.args = model.pred.args,
            model.C = model.C,
            model.S = model.S,
            model.Xi = model.Xi,
            tau_grid = tau_grid,
            alpha = alpha,
            trim.C = trim.C,
            train_ratio = train_ratio,
            TT.name = TT.name,
            theta_grid_AIPCW = theta_grid_AIPCW,
            localization = localization,
            checkpoint_every = checkpoint_every,
            progress_csv = progress_csv,
            censoring_method_AIPCW = censoring_method_AIPCW,
            pred_use_history = FALSE,
            pred_history_mode = "wide",
            G_use_history = FALSE,
            G_history_mode = "wide",
            S_use_history = FALSE,
            S_history_mode = "wide",
            Gtau_mode = Gtau_mode,
            Gtau_delta = Gtau_delta
        ),
        seeds = list(
            seeds = seeds,
            seeds_test = seeds_test,
            seeds_split = seeds_split
        ),
        results = list(
            results_summary = store$results_summary,
            results_bounds = store$results_bounds,
            results_survC = store$results_survC,
            res = res,
            res_ok = res_ok
        ),
        meta = list(
            timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
            runtime = list(rep_runtime_sec = store$rep_runtime_sec)
        )
    )

    saveRDS(save_obj, out_rds)
    cat("Saved ", method_label, " to:\n  ", out_rds, "\n", sep = "")
    invisible(out_rds)
}

set.seed(123)
seeds <- sample(10^7, R, replace = FALSE)
seeds_test <- sample(10^7, R, replace = FALSE)
seeds_split <- sample(10^7, R, replace = FALSE)

stores <- list(
    DDH_model_quantile = new_method_store(),
    dCP_IPCW_new = new_method_store(),
    derived_DynaCoPS_DR = new_method_store(),
    dCP_AIPCW = new_method_store()
)
method_labels <- c(
    DDH_model_quantile = paste0("DDH model ", 100 * (1 - alpha), "% quantile"),
    dCP_IPCW_new = "HAPS-DDH",
    derived_DynaCoPS_DR = "HAPS-DDH-DR",
    dCP_AIPCW = "HAPS-DDH-A"
)

cat("Dynamic DDH tuning simulation\n")
cat("DGM: ", dgm_arg, "\n", sep = "")
cat("R/n/n_test: ", R, "/", n, "/", n_test, "\n", sep = "")
cat("Preset: ", preset_name, "\n", sep = "")
print(ddh_tuning_preset_table()[ddh_tuning_preset_table()$preset == preset_name, ], row.names = FALSE)
cat("Nuisance models: C=", model.C, ", S=", model.S, ", Xi=", model.Xi, "\n", sep = "")
cat("Output folder: ", folder, "\n", sep = "")
cat("Checkpoint every: ", checkpoint_every, " rep(s)\n", sep = "")
cat("Progress CSV: ", progress_csv, "\n", sep = "")
cat("Python: ", reticulate::py_config()$python, "\n", sep = "")

run_t0 <- proc.time()[["elapsed"]]
for (r in seq_len(R)) {
    cat("rep ", r, "..", sep = "")
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

    rep_tabs <- lapply(seq_along(stores), function(...) make_empty_tab(r))
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
        ddh_args_k <- model.pred.args
        ddh_args_k$seed <- as.integer(seeds_split[r] + 1000L * k)

        fit <- tryCatch({
            dynamicCP_AIPCW_split(
                dat = dat_long_r,
                dat_test = dat_test_r,
                start.name = start.name,
                stop.name = stop.name,
                event.name = event.name,
                event.C.name = event.C.name,
                id.name = id.name,
                tau = tau,
                model.pred = model.pred,
                model.pred.args = ddh_args_k,
                covname.pred.baseline = covname.pred.baseline,
                covname.pred.timevarying = covname.pred.timevarying,
                model.C = model.C,
                covname.C.baseline = covname.C.baseline,
                covname.C.timevarying = covname.C.timevarying,
                model.S = model.S,
                model.Xi = model.Xi,
                TT.name = TT.name,
                visit_times = tau_grid,
                seed = seeds_split[r],
                train_ratio = train_ratio,
                theta_grid = theta_grid_AIPCW,
                alpha = alpha,
                trim.C = trim.C,
                pred_use_history = FALSE,
                pred_history_mode = "wide",
                censoring_method = censoring_method_AIPCW,
                G_use_history = FALSE,
                G_history_mode = "wide",
                S_use_history = FALSE,
                S_history_mode = "wide",
                localization = localization,
                Gtau_mode = Gtau_mode,
                Gtau_delta = Gtau_delta,
                return_on_aipcw_fail = TRUE
            )
        }, error = function(e) {
            list(error = TRUE, msg = conditionMessage(e))
        })

        runtime_k <- proc.time()[["elapsed"]] - tau_t0
        method_msg <- if (is.null(fit$msg)) NA_character_ else fit$msg

        if (is.null(fit$lower_IPCW_test)) {
            for (method in names(stores)) {
                rep_tabs[[method]]$msg[k] <- method_msg
                rep_tabs[[method]]$runtime_sec[k] <- runtime_k
                rep_bounds[[method]][[k]] <- empty_bounds()
            }
            next
        }

        dat_cal_tau <- prepare_data_tau(
            dat_long_r[!(dat_long_r[[id.name]] %in% fit$id_tr), , drop = FALSE],
            start.name, stop.name, event.name, id.name, tau
        )
        dat_test_tau <- prepare_data_tau(
            dat_test_r, start.name, stop.name, event.name, id.name, tau
        )

        method_bounds <- list(
            DDH_model_quantile = list(
                lower_cal = fit$lower_model_cal,
                upper_cal = fit$upper_model_cal,
                lower_test = fit$lower_model_test,
                upper_test = fit$upper_model_test,
                theta = alpha / 2,
                ok = TRUE,
                msg = NA_character_
            ),
            dCP_IPCW_new = list(
                lower_cal = fit$lower_IPCW_cal,
                upper_cal = fit$upper_IPCW_cal,
                lower_test = fit$lower_IPCW_test,
                upper_test = fit$upper_IPCW_test,
                theta = fit$theta_hat_IPCW,
                ok = TRUE,
                msg = NA_character_
            ),
            derived_DynaCoPS_DR = list(
                lower_cal = fit$lower_IPCW_DR_cal,
                upper_cal = fit$upper_IPCW_DR_cal,
                lower_test = fit$lower_IPCW_DR_test,
                upper_test = fit$upper_IPCW_DR_test,
                theta = fit$theta_hat_IPCW,
                ok = TRUE,
                msg = NA_character_
            ),
            dCP_AIPCW = list(
                lower_cal = fit$lower_AIPCW_cal,
                upper_cal = fit$upper_AIPCW_cal,
                lower_test = fit$lower_AIPCW_test,
                upper_test = fit$upper_AIPCW_test,
                theta = fit$theta_hat_AIPCW,
                ok = !isTRUE(fit$error),
                msg = if (isTRUE(fit$error)) method_msg else NA_character_
            )
        )

        for (method in names(method_bounds)) {
            mb <- method_bounds[[method]]
            if (!isTRUE(mb$ok) || is.null(mb$lower_test)) {
                rep_tabs[[method]]$msg[k] <- mb$msg
                rep_tabs[[method]]$runtime_sec[k] <- runtime_k
                rep_bounds[[method]][[k]] <- empty_bounds()
                next
            }

            bt <- make_bounds_from_fit(dat_test_tau, mb$lower_test, mb$upper_test, fit$survC_test, r, tau)
            bc <- make_bounds_from_fit(dat_cal_tau, mb$lower_cal, mb$upper_cal, fit$survC_cal, r, tau)
            ss <- summarize_bounds(bt, bc)

            rep_tabs[[method]]$ok[k] <- TRUE
            rep_tabs[[method]]$theta_hat[k] <- mb$theta
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

    if (r %% checkpoint_every == 0L || r == R) {
        write_progress_checkpoint(
            stores = stores,
            rep_done = r,
            elapsed_sec = proc.time()[["elapsed"]] - run_t0
        )
    }
}
cat("\n")

for (method in names(stores)) {
    save_method_bundle(method, stores[[method]], method_labels[[method]])
}

for (method in names(stores)) {
    res_ok <- do.call(rbind, stores[[method]]$results_summary)
    res_ok <- res_ok[res_ok$ok %in% TRUE, , drop = FALSE]
    cat("\nSummary for ", method_labels[[method]], ":\n", sep = "")
    if (nrow(res_ok) == 0L) {
        cat("  no successful tau/rep rows\n")
    } else {
        print(
            dplyr::summarise(
                dplyr::group_by(res_ok, tau),
                cov_true = mean(coverage_test_true, na.rm = TRUE),
                cov_ipcw = mean(coverage_test_ipcw, na.rm = TRUE),
                mean_len = mean(mean_len_test, na.rm = TRUE),
                median_len = stats::median(mean_len_test, na.rm = TRUE),
                theta = mean(theta_hat, na.rm = TRUE),
                ok_n = dplyr::n(),
                .groups = "drop"
            ),
            n = Inf
        )
    }
}
