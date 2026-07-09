## =====================================================================
## Shared per-replication core of the Gtau tilt sensitivity study.
##
## Extracted verbatim from main_simu_gtau_tilt_sensitivity.R so that the
## local serial driver and the OSG worker (osg/gtau_job.R) run EXACTLY the
## same code — one source of truth. Given the same three per-rep seeds,
## run_one_gtau_rep() reproduces the corresponding replication of the
## serial driver bit-for-bit (same machine + package versions).
##
## Requires (sourced by the caller, in this order as in the driver):
##   src/gen_ICML_simu.R, src/helpers.R, src/helpers_AIPCW.R,
##   src/dynamicCP.R, src/dynamicCP_AIPCW.R, src/gtau_eval_helpers.R,
##   src/linWB_dgm_registry.R
## =====================================================================

## Single source of truth for every study constant. Both the local driver
## and the OSG worker build their configuration through this function; any
## change here changes both.
gtau_study_cfg <- function(setup, n, n_test, alpha = 0.1) {
    dgm <- get_linWB_dgm_par(setup)    # validates setup
    list(
        setup = setup, dgm_name = "linear_weibull",
        n = n, n_test = n_test, alpha = alpha,
        rho = dgm$rho, T_max = dgm$T_max, change_times = dgm$change_times,
        dgm_parK = dgm$par,
        tau_grid = c(0, dgm$change_times),           # 0, 3, 6
        delta_grid = dgm$delta_grid,                 # per-setup tilt grid
        # models: correctly-specified Cox everywhere; regression Xi (multiclass
        # is invalid once g_tau != 1, so xgb_reg uniformly across ALL arms)
        model.pred = "cox", model.C = "cox", model.S = "cox", model.Xi = "xgb_reg",
        theta_grid_AIPCW = seq(0, 0.5, by = 0.01),
        train_ratio = 0.5, trim.C = 0.05,
        # variable names
        start.name = "start", stop.name = "stop", event.name = "event",
        event.C.name = "event.C", id.name = "id", TT.name = "TT",
        # covariates
        covname.pred.timevarying = c("L"),
        covname.C.timevarying    = c("L")
    )
}

.gtau_cfg_required <- c(
    "setup", "dgm_name", "n", "n_test", "alpha", "rho", "T_max", "change_times",
    "dgm_parK", "tau_grid", "delta_grid", "model.pred", "model.C", "model.S",
    "model.Xi", "theta_grid_AIPCW", "train_ratio", "trim.C",
    "start.name", "stop.name", "event.name", "event.C.name", "id.name", "TT.name",
    "covname.pred.timevarying", "covname.C.timevarying")

gtau_cfg_check <- function(cfg) {
    missing <- setdiff(.gtau_cfg_required, names(cfg))
    if (length(missing)) {
        stop("gtau study cfg is missing fields: ", paste(missing, collapse = ", "))
    }
    invisible(cfg)
}

## Canonical output path of the study (relative to `results_root`), exactly the
## sprintf the driver has always used. The collector reuses it under
## results/osg/<experiment>/collected/.
gtau_sensitivity_filename <- function(setup, R, n, n_test, alpha,
                                      results_root = "results") {
    file.path(results_root, setup, "gtau_tilt",
              sprintf("gtau_tilt_sensitivity_R%d_n%d_ntest%d_alpha%s.rds",
                      R, n, n_test, format(alpha)))
}

## ------------------------- Calibration wrapper -----------------------
.gtau_calibrate <- function(dat_long, dat_test, tau, gtau_mode, gtau_delta,
                            seed_split, cfg) {
    tryCatch(
        dynamicCP_AIPCW_split(
            dat = dat_long, dat_test = dat_test, tau = tau,
            start.name = cfg$start.name, stop.name = cfg$stop.name,
            event.name = cfg$event.name, event.C.name = cfg$event.C.name,
            id.name = cfg$id.name,
            covname.pred.timevarying = cfg$covname.pred.timevarying,
            model.pred = cfg$model.pred,
            model.C = cfg$model.C, covname.C.timevarying = cfg$covname.C.timevarying,
            model.S = cfg$model.S, model.Xi = cfg$model.Xi,
            theta_grid = cfg$theta_grid_AIPCW, alpha = cfg$alpha, trim.C = cfg$trim.C,
            TT.name = cfg$TT.name, visit_times = cfg$tau_grid,
            seed = seed_split, train_ratio = cfg$train_ratio,
            localization = FALSE, censoring_method = "piecewise",
            Gtau_mode = gtau_mode, Gtau_delta = gtau_delta,
            return_on_aipcw_fail = TRUE,
            support_upper = cfg$T_max   # Inf for linWB1 (legacy candidate sets);
                                        # (tau, T_max] widest member for linWB2
        ),
        error = function(e) list(error = TRUE, msg = conditionMessage(e))
    )
}

# Turn one calibrated fit into per-delta_eval evaluation rows (weighted coverage
# on {T>tau, C_tilde>tau} for each delta_eval). Bounds/covered_T are computed once;
# only the weight vector changes with delta_eval.
# NOTE: this study uses the WEIGHTED evaluator for every arm (the consistent,
# lowest-variance choice; see gtau_eval_method in the shared dynamic driver and
# validation/compare_gtau_evaluators.R for the empirical comparison).
.gtau_eval_rows <- function(fit, method, delta_cal, tau, r, L_full, TT_by_id, cfg) {
    delta_grid <- cfg$delta_grid
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
        w <- compute_Gtau_true_tilted_linWB1(Lsub, tau, de, cfg$change_times, cfg$dgm_parK)
        s <- weighted_coverage_summary(covered, len, w)
        data.frame(rep = r, tau = tau, method = method, delta_cal = delta_cal,
                   delta_eval = de, coverage = s$coverage, med_len = s$med_len,
                   ess = s$ess, theta_hat = fit$theta_hat_AIPCW, ok = TRUE,
                   stringsAsFactors = FALSE)
    })
}

## One full replication: simulate train + (uncensored) test data from the rep's
## seeds, then calibrate "one" + the tilt family (delta_cal = 0 is the
## "estimated" anchor) at every tau and evaluate each fit at every delta_eval.
## Returns a list of one data.frame per (tau, arm) — 18 elements — in the same
## order the serial driver appended them, so do.call(rbind, ...) over the
## concatenated per-rep lists reproduces the serial `results` table exactly.
run_one_gtau_rep <- function(r, seed_dat, seed_test, seed_split, cfg) {
    gtau_cfg_check(cfg)

    dat_tr <- simulate_dataset_long(
        n = cfg$n, seed = seed_dat, change_times = cfg$change_times, tau_max = Inf,
        no_censoring = FALSE, par = cfg$dgm_parK, rho = cfg$rho, dgm_name = cfg$dgm_name,
        T_max = cfg$T_max)
    dat_te <- simulate_dataset_long(
        n = cfg$n_test, seed = seed_test, change_times = cfg$change_times, tau_max = Inf,
        no_censoring = TRUE, par = cfg$dgm_parK, rho = cfg$rho, dgm_name = cfg$dgm_name,
        return_L_full = TRUE, T_max = cfg$T_max)
    L_full   <- attr(dat_te, "L_full")
    rownames(L_full) <- as.character(seq_len(cfg$n_test))
    TT_by_id <- setNames(
        dat_te[[cfg$TT.name]][match(seq_len(cfg$n_test), dat_te[[cfg$id.name]])],
        as.character(seq_len(cfg$n_test)))

    rep_rows <- list()
    for (tau in cfg$tau_grid) {
        # calibrations: "one" + tilt family (delta_cal = 0 is the "estimated" anchor)
        fit_one <- .gtau_calibrate(dat_tr, dat_te, tau, "one", 0, seed_split, cfg)
        rep_rows[[length(rep_rows) + 1L]] <- do.call(rbind,
            .gtau_eval_rows(fit_one, "one", NA_real_, tau, r, L_full, TT_by_id, cfg))

        for (dcal in cfg$delta_grid) {
            method <- if (dcal == 0) "estimated" else "tilted"
            fit <- .gtau_calibrate(dat_tr, dat_te, tau, "tilted", dcal, seed_split, cfg)
            rep_rows[[length(rep_rows) + 1L]] <- do.call(rbind,
                .gtau_eval_rows(fit, method, dcal, tau, r, L_full, TT_by_id, cfg))
        }
    }
    rep_rows
}

## The exact result object the driver has always saved (results/config/seeds/
## meta, with T_max/support_upper added ONLY for finite-support setups so linWB1
## files keep the schema of the existing runs). `meta_extra` lets the OSG
## collector append provenance fields (after the legacy fields, so driver-made
## and collector-made objects agree on everything the plot script reads).
build_gtau_result_object <- function(rows_list, cfg, R,
                                     seeds, seeds_test, seeds_split,
                                     elapsed_sec, completed_reps,
                                     script = "main_simu_gtau_tilt_sensitivity.R",
                                     meta_extra = NULL) {
    out <- list(
        results = do.call(rbind, rows_list),
        config = list(setup = cfg$setup, dgm_name = cfg$dgm_name, rho = cfg$rho,
                      n = cfg$n, n_test = cfg$n_test,
                      R = R, alpha = cfg$alpha, tau_grid = cfg$tau_grid,
                      delta_grid = cfg$delta_grid,
                      model.pred = cfg$model.pred, model.C = cfg$model.C,
                      model.S = cfg$model.S, model.Xi = cfg$model.Xi,
                      change_times = cfg$change_times, dgm_parK = cfg$dgm_parK),
        seeds = list(seeds = seeds, seeds_test = seeds_test, seeds_split = seeds_split),
        meta = list(script = script,
                    elapsed_sec = elapsed_sec,
                    completed_reps = completed_reps,
                    complete = (completed_reps == R))
    )
    # setup-specific fields added conditionally so linWB1 result files keep the
    # exact schema of the existing runs (bit-comparable across driver versions)
    if (is.finite(cfg$T_max)) {
        out$config$T_max <- cfg$T_max
        out$config$support_upper <- cfg$T_max
    }
    if (!is.null(meta_extra)) out$meta <- c(out$meta, meta_extra)
    out
}
