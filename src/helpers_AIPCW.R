# helpers for AIPCW approach 


# subject-level summary (needed for interval endpoints)
subject_summary <- function(dat, id.name, stop.name, event.name, event.C.name, start.name = NULL) {
    # assumes long format, last stop per id is X = min(T,C,admin)
    dat <- dat[order(dat[[id.name]], dat[[stop.name]]), ]
    X     <- ave(dat[[stop.name]], dat[[id.name]], FUN = function(x) x[length(x)])
    Delta <- ave(dat[[event.name]], dat[[id.name]], FUN = function(x) as.integer(any(x == 1)))
    DeltaC <- ave(dat[[event.C.name]], dat[[id.name]], FUN = function(x) as.integer(any(x == 1)))
    
    out <- dat[!duplicated(dat[[id.name]]), id.name, drop = FALSE]
    out$X <- X[!duplicated(dat[[id.name]])]
    out$Delta <- Delta[!duplicated(dat[[id.name]])]
    out$DeltaC <- DeltaC[!duplicated(dat[[id.name]])]
    
    out
}

# fit: coxph object
# newdata: data.frame with N rows
# t_rel: numeric vector length N (relative times since t_k)
# matrix: whether to output a matrix of the survival probabilities; if FALSE output a vector of S(t_rel[i] | newdata[i,])
# returns: numeric vector length N of S(t_rel[i] | newdata[i,]) or a matrix size N*length(t_rel) of S(t_rel[j] | newdata[i,])
predict_surv_coxph_fast <- function(fit, S0_step, newdata, t_rel, matrix = FALSE) {
    t_rel <- pmax(0, t_rel)

    X <- stats::model.matrix(fit, newdata)      # N x p
    # model.matrix() silently DROPS rows with NA covariates; a shorter eta would
    # then be recycled against S0 below, silently assigning survival probabilities
    # to the wrong subjects. Fail loudly instead (the xgb/grf paths already do).
    if (nrow(X) != nrow(newdata)) {
        stop("predict_surv_coxph_fast(): model.matrix dropped ",
             nrow(newdata) - nrow(X), " of ", nrow(newdata),
             " rows (missing covariate values in newdata are not supported).")
    }
    eta <- drop(X %*% stats::coef(fit))         # N
    r <- exp(eta)                               
    S0 <- S0_step(t_rel)                         # N
    
    if(matrix){
        # out[i, j] = S0[j] ^ r[i]
        out <- outer(r, S0, FUN = function(ri, s0j) s0j^ri)
        rownames(out) <- rownames(newdata)
        colnames(out) <- t_rel
    }else{
        out <- S0^r
        names(out) <- rownames(newdata)
    }
    
    out
}

# helper: choose best test metric column from xgb.cv evaluation_log
.xgb_pick_test_metric_col <- function(eval_log) {
    cols <- grep("^test_.*_mean$", names(eval_log), value = TRUE)
    if (length(cols) == 0L) {
        stop(".xgb_pick_test_metric_col(): no test metric column found in xgb.cv evaluation_log.")
    }
    cols[[1]]
}

# helper: count the fitted covariate dimension after factor expansion
.covariate_design_ncoef <- function(dat, cov_use) {
    cov_use <- cov_use[!is.na(cov_use) & nzchar(cov_use)]
    cov_use <- intersect(cov_use, names(dat))
    if (length(cov_use) == 0L || nrow(dat) == 0L) return(0L)
    
    dat_cov <- dat[, cov_use, drop = FALSE]
    keep <- stats::complete.cases(dat_cov)
    dat_cc <- dat[keep, , drop = FALSE]
    if (nrow(dat_cc) == 0L) return(as.integer(length(cov_use)))
    
    form <- stats::as.formula(paste0("~ ", paste(cov_use, collapse = " + ")))
    X <- tryCatch(
        stats::model.matrix(form, data = dat_cc),
        error = function(e) NULL
    )
    if (is.null(X)) return(as.integer(length(cov_use)))
    
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    as.integer(ncol(X))
}

# helper: require both an absolute minimum and enough events per coefficient
.event_threshold_by_dim <- function(base_min_events, min_events_per_coef, n_coef) {
    base_min_events <- as.numeric(base_min_events[[1L]])
    if (!is.finite(base_min_events) || base_min_events < 0) {
        stop(".event_threshold_by_dim(): base_min_events must be nonnegative.")
    }
    min_events_per_coef <- as.numeric(min_events_per_coef[[1L]])
    if (!is.finite(min_events_per_coef) || min_events_per_coef < 0) {
        stop(".event_threshold_by_dim(): min_events_per_coef must be nonnegative.")
    }
    if (length(n_coef) == 0L) n_coef <- 0L
    n_coef <- as.integer(n_coef[[1L]])
    if (!is.finite(n_coef) || n_coef < 0L) n_coef <- 0L
    
    as.integer(ceiling(max(base_min_events, min_events_per_coef * n_coef)))
}

# helper: effective sample size for positive finite weights
.weighted_ess <- function(w) {
    w <- as.numeric(w)
    w <- w[is.finite(w) & w > 0]
    if (length(w) == 0L) return(0)
    (sum(w)^2) / sum(w^2)
}

# helper: build candidate parameter sets for optional CV search
.xgb_param_candidates <- function(base_params, cv_param_grid = NULL) {
    if (is.null(cv_param_grid)) return(list(base_params))
    
    # data.frame grid
    if (is.data.frame(cv_param_grid)) {
        grid_df <- cv_param_grid
        out <- vector("list", nrow(grid_df))
        for (i in seq_len(nrow(grid_df))) {
            out[[i]] <- utils::modifyList(base_params, as.list(grid_df[i, , drop = FALSE]))
        }
        return(out)
    }
    
    # named list -> expand.grid
    if (is.list(cv_param_grid) && !is.null(names(cv_param_grid)) && length(cv_param_grid) > 0L) {
        grid_df <- expand.grid(cv_param_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        out <- vector("list", nrow(grid_df))
        for (i in seq_len(nrow(grid_df))) {
            out[[i]] <- utils::modifyList(base_params, as.list(grid_df[i, , drop = FALSE]))
        }
        return(out)
    }
    
    # list-of-lists
    if (is.list(cv_param_grid) && all(vapply(cv_param_grid, is.list, logical(1)))) {
        return(lapply(cv_param_grid, function(g) utils::modifyList(base_params, g)))
    }
    
    stop(".xgb_param_candidates(): cv_param_grid must be NULL, data.frame, named list, or list-of-lists.")
}

# helper: train xgboost with optional cross-validation tuning
.xgb_train_optional_cv <- function(
        dtrain, params, nrounds = 200L,
        use_cv = FALSE,
        cv_nfold = 5L,
        cv_early_stopping_rounds = 20L,
        cv_param_grid = NULL,
        cv_seed = NULL,
        cv_verbose = 0,
        maximize = FALSE,
        train_verbose = 0
) {
    nrounds <- as.integer(nrounds)
    if (!is.finite(nrounds) || nrounds < 1L) {
        stop(".xgb_train_optional_cv(): nrounds must be >= 1.")
    }
    
    obj <- if (!is.null(params$objective)) as.character(params$objective) else ""
    if (isTRUE(use_cv) && identical(obj, "survival:aft")) {
        warning(
            ".xgb_train_optional_cv(): xgb.cv is not supported in this pipeline for objective='survival:aft'. ",
            "Falling back to xgb.train without CV for this fit."
        )
        use_cv <- FALSE
    }
    
    if (!isTRUE(use_cv)) {
        fit <- xgboost::xgb.train(
            params = params,
            data = dtrain,
            nrounds = nrounds,
            verbose = train_verbose
        )
        return(list(
            fit = fit,
            tuning = list(
                use_cv = FALSE,
                params = params,
                best_nrounds = nrounds
            )
        ))
    }
    
    cand_params <- .xgb_param_candidates(params, cv_param_grid = cv_param_grid)
    best <- NULL
    
    for (ii in seq_along(cand_params)) {
        if (!is.null(cv_seed)) set.seed(as.integer(cv_seed) + ii - 1L)
        
        cv_obj <- xgboost::xgb.cv(
            params = cand_params[[ii]],
            data = dtrain,
            nrounds = nrounds,
            nfold = as.integer(cv_nfold),
            early_stopping_rounds = as.integer(cv_early_stopping_rounds),
            maximize = isTRUE(maximize),
            verbose = cv_verbose
        )
        
        metric_col <- .xgb_pick_test_metric_col(cv_obj$evaluation_log)
        it <- cv_obj$best_iteration
        if (is.null(it) || !is.finite(it) || it < 1L) {
            it <- nrow(cv_obj$evaluation_log)
        }
        score <- cv_obj$evaluation_log[[metric_col]][it]
        
        better <- is.null(best) ||
            ((!maximize) && score < best$score) ||
            (isTRUE(maximize) && score > best$score)
        
        if (better) {
            best <- list(
                params = cand_params[[ii]],
                score = score,
                metric = metric_col,
                best_nrounds = as.integer(it)
            )
        }
    }
    
    fit <- xgboost::xgb.train(
        params = best$params,
        data = dtrain,
        nrounds = best$best_nrounds,
        verbose = train_verbose
    )
    
    list(
        fit = fit,
        tuning = c(
            list(use_cv = TRUE),
            best
        )
    )
}

.fit_hal_cox_bundle <- function(dat, cov_use, Y_surv, weights = NULL) {
    if (!requireNamespace("hal9001", quietly = TRUE)) {
        stop(".fit_hal_cox_bundle(): package 'hal9001' is required for model='hal'.")
    }
    
    dm_fit <- .fit_design_matrix(dat, cov_use = cov_use)
    if (is.null(dm_fit$X) || nrow(dm_fit$X) == 0) {
        stop(".fit_hal_cox_bundle(): no complete-case rows left after design matrix construction.")
    }
    
    Y_fit <- Y_surv[dm_fit$keep]
    fit_args <- list(
        X = dm_fit$X,
        Y = Y_fit,
        family = "cox",
        yolo = FALSE
    )
    if (!is.null(weights)) {
        fit_args$weights <- as.numeric(weights[dm_fit$keep])
    }
    fit_hal <- do.call(hal9001::fit_hal, fit_args)
    
    lp_fit <- as.numeric(stats::predict(fit_hal, new_data = dm_fit$X, type = "link"))
    lp_fit[!is.finite(lp_fit)] <- 0
    lp_fit <- pmax(pmin(lp_fit, 30), -30)
    lp_center <- mean(lp_fit)
    dat_anchor <- data.frame(lp_off = lp_fit - lp_center)
    dat_anchor$Y_surv <- Y_fit
    fit_anchor <- survival::coxph(Y_surv ~ offset(lp_off), data = dat_anchor, ties = "breslow", x = TRUE)
    bh <- survival::basehaz(fit_anchor, centered = FALSE)
    
    list(
        fit_hal = fit_hal,
        design_info = dm_fit$design_info,
        cov_use = cov_use,
        keep_train = dm_fit$keep,
        lp_center = lp_center,
        bh_time = bh$time,
        bh_cumhaz = bh$hazard
    )
}

.predict_hal_cox_lp <- function(hal_bundle, newdata) {
    if (!is.null(hal_bundle$cov_use)) {
        miss <- setdiff(hal_bundle$cov_use, names(newdata))
        if (length(miss) > 0) {
            stop(".predict_hal_cox_lp(): newdata is missing covariates: ", paste(miss, collapse = ", "))
        }
    }
    
    dm_new <- .predict_design_matrix(
        dat = newdata,
        cov_use = hal_bundle$cov_use,
        design_info = hal_bundle$design_info
    )
    if (any(!dm_new$keep)) {
        stop(".predict_hal_cox_lp(): missing values in covariates are not supported for HAL prediction.")
    }
    
    lp <- as.numeric(stats::predict(hal_bundle$fit_hal, new_data = dm_new$X, type = "link"))
    if (length(lp) != nrow(newdata)) {
        stop(".predict_hal_cox_lp(): prediction length mismatch.")
    }
    lp
}

.predict_surv_hal_cox <- function(hal_bundle, lp, t_stop, t_start = NULL) {
    if (length(lp) != length(t_stop)) {
        stop(".predict_surv_hal_cox(): lp and t_stop must have the same length.")
    }
    if (is.null(t_start)) t_start <- rep(0, length(t_stop))
    if (length(t_start) != length(t_stop)) {
        stop(".predict_surv_hal_cox(): t_start and t_stop must have the same length.")
    }
    
    if (length(hal_bundle$bh_time) == 0 || length(hal_bundle$bh_cumhaz) == 0) {
        return(rep(1.0, length(t_stop)))
    }
    
    H0_step <- stats::stepfun(hal_bundle$bh_time, c(0, hal_bundle$bh_cumhaz))
    H0_stop <- H0_step(t_stop)
    H0_start <- H0_step(t_start)
    dH0 <- pmax(0, H0_stop - H0_start)
    lp_center <- if (!is.null(hal_bundle$lp_center) && is.finite(hal_bundle$lp_center)) hal_bundle$lp_center else 0
    lp_adj <- lp - lp_center
    exp(-dH0 * exp(lp_adj))
}


# grid: length m increasing
# S: n x m survival values at grid times
# t_rel: length n query times per subject
# returns: length n survival evaluated as right-continuous step, with S(t<grid[1])=1
eval_step_surv <- function(S, grid, t_rel) {
    stopifnot(is.matrix(S), length(grid) == ncol(S), length(t_rel) == nrow(S))
    
    m <- length(grid)
    n <- nrow(S)
    
    # j = max { j : grid[j] <= tt }  (0 if tt < grid[1])
    j0 <- findInterval(t_rel, grid, left.open = FALSE, rightmost.closed = FALSE)
    # j0 in {0,1,...,m}
    
    out <- rep(1.0, n)
    ok  <- j0 > 0
    out[ok] <- S[cbind(which(ok), j0[ok])]
    out
}


# Build the landmark dataset at tau, optionally with covariate history up to tau
# - If use_history=FALSE: uses prepare_data_tau()
# - If use_history=TRUE:  uses prepare_data_tau_with_history() with k_tau = match(tau, visit_times)
make_tau_df <- function(dat_long, tau, visit_times,
                        start.name, stop.name, event.name, id.name,
                        cov_tv = NULL,
                        use_history = FALSE,
                        history_mode = c("wide", "summary")) {
    history_mode <- match.arg(history_mode)
    
    k_tau <- match(tau, visit_times)
    if (is.na(k_tau)) stop("make_tau_df(): tau must be one of visit_times.")
    
    if (use_history && !is.null(cov_tv) && length(cov_tv) > 0) {
        prepare_data_tau_with_history(
            dat_long, visit_times, k_tau,
            start.name, stop.name, event.name, id.name,
            covars = cov_tv,
            mode = history_mode
        )
    } else {
        prepare_data_tau(dat_long, start.name, stop.name, event.name, id.name, tau)
    }
}

# helper to expand covariate names
expand_hist_covars <- function(cov_tv, k, mode = c("wide","summary")) {
    mode <- match.arg(mode)
    if (is.null(cov_tv) || length(cov_tv) == 0) return(character(0))
    if (mode == "wide") {
        unlist(lapply(seq_len(k), function(j) paste0(cov_tv, "_t", j)))
    } else {
        unlist(lapply(cov_tv, function(v) paste0(v, c("_hist_mean","_hist_last","_hist_slope"))))
    }
}


# Build features at tau = visit_times[k] that include covariate history at visit_times[1:k]
# mode = "wide": append each covariate at each past visit time as var_t1, var_t2, ...
# mode = "summary": append summaries over the history (mean, last, slope) for each covariate
prepare_data_tau_with_history <- function(dat, visit_times, k,
                                          start.name, stop.name, event.name, id.name,
                                          covars,
                                          mode = c("wide", "summary")
                                          ) {
    mode <- match.arg(mode)
    stopifnot(k >= 1, k <= length(visit_times))
    
    t_k <- visit_times[k]
    
    # Base snapshot at t_k (this also filters to the correct interval for tau)
    dat_tk <- prepare_data_tau(dat, start.name, stop.name, event.name, id.name, tau = t_k)
    
    # Restrict to subjects at risk at t_k (using subject-level stop that prepare_data_tau writes)
    dat_tk <- dat_tk[dat_tk[[stop.name]] > t_k, , drop = FALSE]
    if (nrow(dat_tk) == 0) return(dat_tk)
    
    ids <- dat_tk[[id.name]]
    
    # Collect covariates at each visit time t_1,...,t_k for these ids
    hist_list <- vector("list", k)
    for (j in seq_len(k)) {
        tj <- visit_times[j]
        dj <- prepare_data_tau(dat[dat[[id.name]] %in% ids, , drop = FALSE],
                               start.name, stop.name, event.name, id.name, tau = tj)
        dj <- dj[, c(id.name, covars), drop = FALSE]
        # keep only 1 row per id (prepare_data_tau should already give 1 per id)
        dj <- dj[!duplicated(dj[[id.name]]), , drop = FALSE]
        
        if (mode == "wide") {
            suf <- paste0("_t", j) 
            names(dj)[names(dj) %in% covars] <- paste0(names(dj)[names(dj) %in% covars], suf)
            hist_list[[j]] <- dj
        } else if (mode == "summary"){
            # in "summary" mode, keep raw per-time values for later summarization
            names(dj)[names(dj) %in% covars] <- paste0(names(dj)[names(dj) %in% covars], "_raw")
            dj$._time_index <- j
            hist_list[[j]] <- dj
        } else {
            stop(paste0("Unknown mode: ", mode))
        }
    }
    
    if (mode == "wide") {
        hist_wide <- Reduce(function(a, b) merge(a, b, by = id.name, all = TRUE, sort = FALSE), hist_list)
        out <- merge(dat_tk, hist_wide, by = id.name, all.x = TRUE, sort = FALSE)
        return(out)
    }
    
    # summary mode
    hist_long <- do.call(rbind, hist_list)
    # ensure aligned order
    hist_long <- hist_long[order(hist_long[[id.name]], hist_long$._time_index), , drop = FALSE]
    
    # build summaries per covariate per id
    out_sum <- hist_long[!duplicated(hist_long[[id.name]]), id.name, drop = FALSE]
    
    for (v in covars) {
        vv <- paste0(v, "_raw")
        x <- hist_long[[vv]]
        idv <- hist_long[[id.name]]
        tv <- hist_long$._time_index
        
        # mean and last
        meanv <- ave(x, idv, FUN = function(z) mean(z, na.rm = TRUE))
        lastv <- ave(x, idv, FUN = function(z) z[length(z)])
        
        # simple slope vs time index (needs >=2 points)
        slopev <- ave(seq_along(x), idv, FUN = function(idx) {
            # idx indexes into x/idv/tv for that id
            # We'll compute in-place using the global vectors
            NA_real_
        })
        # compute slope properly per id
        slope_by_id <- tapply(seq_len(nrow(hist_long)), idv, function(ii) {
            if (length(ii) < 2) return(NA_real_)
            zz <- x[ii]
            tt <- tv[ii]
            if (all(is.na(zz))) return(NA_real_)
            coef(lm(zz ~ tt))[2]
        })
        slopev2 <- as.numeric(slope_by_id[as.character(idv)])
        
        tmp <- hist_long[!duplicated(idv), id.name, drop = FALSE]
        tmp[[paste0(v, "_hist_mean")]]  <- meanv[!duplicated(idv)]
        tmp[[paste0(v, "_hist_last")]]  <- lastv[!duplicated(idv)]
        tmp[[paste0(v, "_hist_slope")]] <- slopev2[!duplicated(idv)]
        
        out_sum <- merge(out_sum, tmp, by = id.name, all.x = TRUE, sort = FALSE)
    }
    
    out <- merge(dat_tk, out_sum, by = id.name, all.x = TRUE, sort = FALSE)
    out
}




fit_Gk_list <- function(dat_tr, visit_times,
                        start.name, stop.name, event.name, event.C.name, id.name,
                        model.C = c("km", "cox", "rsf", "hal", "xgb_cox", "xgb_aft"),
                        covname.C.baseline = NULL, covname.C.timevarying = NULL,
                        use_history = FALSE,
                        history_mode = c("wide", "summary"),
                        cox_min_censor_events = 10L,
                        cox_min_censor_events_per_coef = 5L,
                        cox_warn_fallback = TRUE,
                        rsf_min_censor_events = 10L,
                        rsf_min_censor_events_per_coef = 5L,
                        rsf_args = list(),
                        xgb_nrounds = 200L,
                        xgb_use_cv = FALSE,
                        xgb_cv_nfold = 5L,
                        xgb_cv_early_stopping_rounds = 20L,
                        xgb_cv_param_grid = NULL,
                        xgb_cv_seed = NULL,
                        xgb_cv_verbose = 0
                        ) {
    
    cox_min_censor_events_per_coef_missing <- missing(cox_min_censor_events_per_coef)
    rsf_min_censor_events_per_coef_missing <- missing(rsf_min_censor_events_per_coef)
    model.C <- match.arg(model.C)
    history_mode <- match.arg(history_mode)
    if (!is.list(rsf_args)) {
        stop("fit_Gk_list(): rsf_args must be a list.")
    }
    cox_base_tmp <- suppressWarnings(as.integer(cox_min_censor_events[[1L]]))
    rsf_base_tmp <- suppressWarnings(as.integer(rsf_min_censor_events[[1L]]))
    if (is.finite(cox_base_tmp) && cox_base_tmp == 0L && cox_min_censor_events_per_coef_missing) {
        cox_min_censor_events_per_coef <- 0L
    }
    if (is.finite(rsf_base_tmp) && rsf_base_tmp == 0L && rsf_min_censor_events_per_coef_missing) {
        rsf_min_censor_events_per_coef <- 0L
    }
    
    cov_base <- covname.C.baseline
    cov_tv   <- covname.C.timevarying
    
    K <- length(visit_times)
    stopifnot(K >= 1)
    
    summ <- subject_summary(dat_tr, id.name, stop.name, event.name, event.C.name)
    
    Gfits <- vector("list", K)
    names(Gfits) <- paste0("k=", seq_len(K))
    
    for (k in 1:K) {
        t_k  <- visit_times[k]
        t_k1 <- if (k < K) visit_times[k + 1] else Inf
        
        if (use_history && length(cov_tv) > 0) {
            dat_tk <- prepare_data_tau_with_history(
                dat_tr, visit_times, k,
                start.name, stop.name, event.name, id.name,
                covars = cov_tv,
                mode = history_mode
            )
        } else {
            dat_tk <- prepare_data_tau(dat_tr, start.name, stop.name, event.name, id.name, tau = t_k)
        }
        
        # avoid X/Delta duplication if present
        dat_tk <- dat_tk[, setdiff(names(dat_tk), c("X", "Delta", "DeltaC")), drop = FALSE]
        dat_tk <- merge(dat_tk, summ, by = id.name, all.x = TRUE, sort = FALSE)
        
        # at risk
        dat_tk <- dat_tk[dat_tk$X > t_k, , drop = FALSE]
        if (nrow(dat_tk) == 0) {
            Gfits[k] <- list(NULL)
            next
        }
        
        # interval censoring outcome on relative time
        dat_tk$timeC_rel <- pmin(dat_tk$X, t_k1) - t_k
        dat_tk$eventC_k  <- as.integer(dat_tk$DeltaC == 1 & dat_tk$X < t_k1)
        n_censor_events_k <- sum(dat_tk$eventC_k == 1L, na.rm = TRUE)
        
        if (!use_history) {
            cov_C_use <- c(cov_base, cov_tv)
        } else {
            # baseline enters once; only time-varying gets history features
            if (history_mode == "wide") {
                cov_hist <- unlist(lapply(seq_len(k), function(j) paste0(cov_tv, "_t", j)))
            } else if (history_mode == "summary") { # "summary"
                cov_hist <- unlist(lapply(cov_tv, function(v) paste0(v, c("_hist_mean","_hist_last","_hist_slope"))))
            } else {
                stop(paste0("Unknown history_mode: ", history_mode))
            }
            cov_C_use <- c(cov_base, cov_hist)
        }
        
        cov_C_use <- keep_nonconstant_covars(dat_tk, cov_C_use)
        n_coef_C <- .covariate_design_ncoef(dat_tk, cov_C_use)
        cox_min_censor_events_k <- .event_threshold_by_dim(
            cox_min_censor_events,
            cox_min_censor_events_per_coef,
            n_coef_C
        )
        rsf_min_censor_events_k <- .event_threshold_by_dim(
            rsf_min_censor_events,
            rsf_min_censor_events_per_coef,
            n_coef_C
        )
        
        use_km_fallback <- FALSE
        km_reason <- NULL
        if (model.C == "km" || length(cov_C_use) == 0) {
            use_km_fallback <- TRUE
            if (model.C != "km" && length(cov_C_use) == 0) {
                km_reason <- "no_usable_covariates"
            }
        } else if (model.C == "cox" && n_censor_events_k < cox_min_censor_events_k) {
            use_km_fallback <- TRUE
            km_reason <- paste0("too_few_censor_events(<", cox_min_censor_events_k, ")")
            if (isTRUE(cox_warn_fallback)) {
                warning(
                    "fit_Gk_list(): model.C='cox' fallback to KM at interval [",
                    t_k, ", ", ifelse(is.finite(t_k1), as.character(t_k1), "Inf"),
                    ") because censoring events=", n_censor_events_k,
                    " < required=", cox_min_censor_events_k,
                    " based on ", n_coef_C, " coefficient(s)."
                )
            }
        } else if (model.C == "rsf" && n_censor_events_k < rsf_min_censor_events_k) {
            use_km_fallback <- TRUE
            km_reason <- paste0("too_few_censor_events(<", rsf_min_censor_events_k, ")")
        }
        
        if (use_km_fallback) {
            fit_km <- survival::survfit(stats::as.formula("Surv(timeC_rel, eventC_k) ~ 1"),
                                        data = dat_tk)
            ss <- summary(fit_km)
            if (length(ss$time) == 0L) {
                # No censoring jumps in this interval => estimated G_k(t)=1 over interval.
                S_step <- function(tt) rep(1, length(tt))
                jump_times_rel <- numeric(0)
            } else {
                S_step <- stats::stepfun(ss$time, c(1, ss$surv))
                jump_times_rel <- ss$time
            }
            
            Gfits[[k]] <- list(
                model = "km",
                t_k = t_k, t_k1 = t_k1,
                S_step = S_step,
                jump_times_rel = jump_times_rel,
                cov_use = character(0),
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base,
                n_censor_events = n_censor_events_k,
                n_coef = n_coef_C,
                min_censor_events_required = if (model.C == "rsf") rsf_min_censor_events_k else cox_min_censor_events_k,
                requested_model = model.C,
                fallback_reason = km_reason
            )
            
        } else if (model.C == "cox"){
            rhs <- paste(cov_C_use, collapse = " + ")
            form <- stats::as.formula(paste0("Surv(timeC_rel, eventC_k) ~ ", rhs))
            fit_cox <- survival::coxph(form, data = dat_tk, x = TRUE)
            
            # baseline survival step function
            bh <- survival::basehaz(fit_cox, centered = FALSE)
            H0 <- bh$hazard
            tt <- bh$time
            S0_step <- stats::stepfun(tt, c(1, exp(-H0)))
            
            Gfits[[k]] <- list(
                model = "cox",
                t_k = t_k, t_k1 = t_k1,
                fit = fit_cox,
                S0_step = S0_step,
                jump_times_rel = tt, 
                cov_use = cov_C_use,
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base,
                n_censor_events = n_censor_events_k,
                n_coef = n_coef_C,
                min_censor_events_required = cox_min_censor_events_k
            )
            
        } else if (model.C == "rsf") {
            
            if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
                stop("fit_Gk_list(): package 'randomForestSRC' required for model.C='rsf'.")
            }
            
            # If no usable covariates, RSF reduces to ~1 anyway; you can either:
            # (a) fit ~1 via RSF, or (b) fall back to KM. I'll fit ~1 via RSF.
            rhs <- if (length(cov_C_use) > 0) paste(cov_C_use, collapse = " + ") else "1"
            form <- stats::as.formula(paste0("Surv(timeC_rel, eventC_k) ~ ", rhs))
            
            rsf_fit_args <- utils::modifyList(
                list(
                    formula  = form,
                    data     = dat_tk,
                    ntree    = 1000L,
                    nodesize = 30L,
                    nsplit   = 10L,
                    forest   = TRUE
                ),
                rsf_args
            )
            # keep formula/data tied to the current interval data
            rsf_fit_args$formula <- form
            rsf_fit_args$data <- dat_tk
            fit_rsf <- do.call(randomForestSRC::rfsrc, rsf_fit_args)
            
            Gfits[[k]] <- list(
                model = "rsf",
                t_k = t_k, t_k1 = t_k1,
                fit = fit_rsf,
                jump_times_rel = fit_rsf$time.interest,
                cov_use = cov_C_use,
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base,
                n_censor_events = n_censor_events_k,
                n_coef = n_coef_C,
                min_censor_events_required = rsf_min_censor_events_k
            )
            
        } else if (model.C == "grf") {
            if (!requireNamespace("grf", quietly = TRUE)) {
                stop("fit_Gk_list(): package 'grf' required for model.C='grf'.")
            }
            
            dm_fit <- .fit_design_matrix(dat_tk, cov_use = cov_C_use)
            if (is.null(dm_fit$X) || nrow(dm_fit$X) == 0) {
                stop("fit_Gk_list(): no complete-case rows left after design matrix construction for model.C='grf'.")
            }
            dat_fit_ml <- dat_tk[dm_fit$keep, , drop = FALSE]
            
            fit_grf <- grf::survival_forest(
                X = dm_fit$X,
                Y = as.numeric(dat_fit_ml$timeC_rel),
                D = as.numeric(dat_fit_ml$eventC_k)
            )
            
            Gfits[[k]] <- list(
                model = "grf",
                t_k = t_k, t_k1 = t_k1,
                fit = fit_grf,
                jump_times_rel = fit_grf$failure.times,
                cov_use = cov_C_use,
                design_info = dm_fit$design_info,
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base
            )
            
        } else if (model.C == "hal") {
            
            Y_surv <- survival::Surv(dat_tk$timeC_rel, dat_tk$eventC_k)
            hal_bundle <- .fit_hal_cox_bundle(
                dat = dat_tk,
                cov_use = cov_C_use,
                Y_surv = Y_surv
            )
            
            Gfits[[k]] <- list(
                model = "hal",
                t_k = t_k, t_k1 = t_k1,
                fit = hal_bundle$fit_hal,
                lp_center = hal_bundle$lp_center,
                bh_time = hal_bundle$bh_time,
                bh_cumhaz = hal_bundle$bh_cumhaz,
                jump_times_rel = hal_bundle$bh_time,
                cov_use = cov_C_use,
                design_info = hal_bundle$design_info,
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base
            )
            
        } else if (model.C %in% c("xgb_cox", "xgb_aft")) {
            if (!requireNamespace("xgboost", quietly = TRUE)) {
                stop("fit_Gk_list(): package 'xgboost' required for model.C='", model.C, "'.")
            }
            
            dm_fit <- .fit_design_matrix(dat_tk, cov_use = cov_C_use)
            if (is.null(dm_fit$X) || nrow(dm_fit$X) == 0) {
                stop("fit_Gk_list(): no complete-case rows left after design matrix construction for model.C='", model.C, "'.")
            }
            dat_fit_ml <- dat_tk[dm_fit$keep, , drop = FALSE]
            X_ml <- dm_fit$X
            
            xgb_params_common <- list(
                eta = 0.05,
                max_depth = 4L,
                min_child_weight = 1,
                subsample = 0.8,
                colsample_bytree = 0.8,
                lambda = 1,
                gamma = 0
            )
            
            if (model.C == "xgb_cox") {
                # xgboost Cox uses sign(label): censored observations have negative labels.
                y_cox <- ifelse(dat_fit_ml$eventC_k == 1L,
                                as.numeric(dat_fit_ml$timeC_rel),
                                -as.numeric(dat_fit_ml$timeC_rel))
                dtrain <- xgboost::xgb.DMatrix(data = X_ml, label = y_cox)
                params <- c(list(
                    objective = "survival:cox",
                    eval_metric = "cox-nloglik"
                ), xgb_params_common)
                xgb_fit <- .xgb_train_optional_cv(
                    dtrain = dtrain,
                    params = params,
                    nrounds = as.integer(xgb_nrounds),
                    use_cv = xgb_use_cv,
                    cv_nfold = xgb_cv_nfold,
                    cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
                    cv_param_grid = xgb_cv_param_grid,
                    cv_seed = xgb_cv_seed,
                    cv_verbose = xgb_cv_verbose,
                    maximize = FALSE,
                    train_verbose = 0
                )
                fit_xgb <- xgb_fit$fit
                
                lp_fit <- as.numeric(stats::predict(fit_xgb, newdata = X_ml, outputmargin = TRUE))
                lp_fit[!is.finite(lp_fit)] <- 0
                lp_fit <- pmax(pmin(lp_fit, 30), -30)
                lp_center <- mean(lp_fit)
                dat_anchor <- data.frame(
                    timeC_rel = dat_fit_ml$timeC_rel,
                    eventC_k = dat_fit_ml$eventC_k,
                    lp_off = lp_fit - lp_center
                )
                fit_anchor <- survival::coxph(
                    survival::Surv(timeC_rel, eventC_k) ~ offset(lp_off),
                    data = dat_anchor,
                    ties = "breslow",
                    x = TRUE
                )
                bh <- survival::basehaz(fit_anchor, centered = FALSE)
                
                Gfits[[k]] <- list(
                    model = "xgb_cox",
                    t_k = t_k, t_k1 = t_k1,
                    fit = fit_xgb,
                    tuning = xgb_fit$tuning,
                    lp_center = lp_center,
                    bh_time = bh$time,
                    bh_cumhaz = bh$hazard,
                    jump_times_rel = bh$time,
                    cov_use = cov_C_use,
                    design_info = dm_fit$design_info,
                    use_history = use_history,
                    history_mode = if (use_history) history_mode else NULL,
                    cov_tv = cov_tv,
                    cov_base = cov_base
                )
                
            } else {
                dtrain <- xgboost::xgb.DMatrix(data = X_ml)
                xgboost::setinfo(dtrain, "label_lower_bound", as.numeric(dat_fit_ml$timeC_rel))
                xgboost::setinfo(
                    dtrain, "label_upper_bound",
                    ifelse(dat_fit_ml$eventC_k == 1L, as.numeric(dat_fit_ml$timeC_rel), Inf)
                )
                
                aft_dist <- "normal"
                aft_scale <- 1
                params <- c(list(
                    objective = "survival:aft",
                    eval_metric = "aft-nloglik",
                    aft_loss_distribution = aft_dist,
                    aft_loss_distribution_scale = aft_scale
                ), xgb_params_common)
                xgb_fit <- .xgb_train_optional_cv(
                    dtrain = dtrain,
                    params = params,
                    nrounds = as.integer(xgb_nrounds),
                    use_cv = xgb_use_cv,
                    cv_nfold = xgb_cv_nfold,
                    cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
                    cv_param_grid = xgb_cv_param_grid,
                    cv_seed = xgb_cv_seed,
                    cv_verbose = xgb_cv_verbose,
                    maximize = FALSE,
                    train_verbose = 0
                )
                fit_xgb <- xgb_fit$fit
                
                Gfits[[k]] <- list(
                    model = "xgb_aft",
                    t_k = t_k, t_k1 = t_k1,
                    fit = fit_xgb,
                    tuning = xgb_fit$tuning,
                    aft_dist = aft_dist,
                    aft_scale = aft_scale,
                    cov_use = cov_C_use,
                    design_info = dm_fit$design_info,
                    use_history = use_history,
                    history_mode = if (use_history) history_mode else NULL,
                    cov_tv = cov_tv,
                    cov_base = cov_base
                )
            }
            
        }else{
            stop("Unknown model.C: ", model.C)
        }
    }
    
    Gfits
}


# Return a vector for G(t_i|Z_i)
predict_Gk <- function(Gfit_k, newdata_tk, t_abs) {
    if (is.null(Gfit_k)) return(rep(NA_real_, nrow(newdata_tk)))
    
    t_k  <- Gfit_k$t_k
    t_k1 <- Gfit_k$t_k1
    
    # relative time within interval (clamped)
    t_rel <- pmax(0, pmin(t_abs, t_k1) - t_k)
    
    if (Gfit_k$model == "km") {
        return(Gfit_k$S_step(t_rel))
        
    } else if (Gfit_k$model == "cox") {
        # Ensure the covariates needed by the Cox model exist in newdata_tk
        if (!is.null(Gfit_k$cov_use)) {
            miss <- setdiff(Gfit_k$cov_use, names(newdata_tk))
            if (length(miss) > 0) {
                stop("predict_Gk(): newdata_tk is missing covariates: ", paste(miss, collapse = ", "))
            }
        }
        
        return(
            predict_surv_coxph_fast(
                fit   = Gfit_k$fit,
                S0_step = Gfit_k$S0_step,
                newdata = newdata_tk,
                t_rel = t_rel
            )
        )
        
    }else if (Gfit_k$model == "rsf") {
        
        if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
            stop("predict_Gk(): package 'randomForestSRC' required for model='rsf'.")
        }
        
        # Ensure covariates exist
        if (!is.null(Gfit_k$cov_use)) {
            miss <- setdiff(Gfit_k$cov_use, names(newdata_tk))
            if (length(miss) > 0) {
                stop("predict_Gk(): newdata_tk is missing covariates: ", paste(miss, collapse = ", "))
            }
        }
        
        # Predict survival curve for censoring time at RSF internal grid
        pred <- predict(Gfit_k$fit, newdata = newdata_tk)
        
        grid <- pred$time.interest
        S    <- pred$survival  # n x length(grid); this is Ghat(t | Lk)
        
        t_k  <- Gfit_k$t_k
        t_k1 <- Gfit_k$t_k1
        # relative time within interval (clamped)
        t_rel <- pmax(0, pmin(t_abs, t_k1) - t_k)
        out <- eval_step_surv(S, grid, t_rel = t_rel)
        
        return(out)
        
    } else if (Gfit_k$model == "grf") {
        
        if (!requireNamespace("grf", quietly = TRUE)) {
            stop("predict_Gk(): package 'grf' required for model='grf'.")
        }
        
        # Ensure covariates exist
        if (!is.null(Gfit_k$cov_use)) {
            miss <- setdiff(Gfit_k$cov_use, names(newdata_tk))
            if (length(miss) > 0) {
                stop("predict_Gk(): newdata_tk is missing covariates: ", paste(miss, collapse = ", "))
            }
        }
        
        dm_new <- .predict_design_matrix(
            dat = newdata_tk,
            cov_use = Gfit_k$cov_use,
            design_info = Gfit_k$design_info
        )
        if (any(!dm_new$keep)) {
            stop("predict_Gk(): missing values in covariates are not supported for GRF prediction.")
        }
        
        pred <- predict(Gfit_k$fit, newdata = dm_new$X, failure.times = Gfit_k$jump_times_rel)
        
        grid <- pred$failure.times
        S <- pred$predictions
        
        t_rel <- pmax(0, pmin(t_abs, t_k1) - t_k)
        out <- eval_step_surv(S, grid, t_rel = t_rel)
        
        return(out)
        
    } else if (Gfit_k$model == "hal") {
        
        hal_bundle <- list(
            fit_hal = Gfit_k$fit,
            design_info = Gfit_k$design_info,
            cov_use = Gfit_k$cov_use,
            lp_center = Gfit_k$lp_center,
            bh_time = Gfit_k$bh_time,
            bh_cumhaz = Gfit_k$bh_cumhaz
        )
        lp_new <- .predict_hal_cox_lp(hal_bundle, newdata_tk)
        
        out <- .predict_surv_hal_cox(
            hal_bundle = hal_bundle,
            lp = lp_new,
            t_start = rep(0, length(t_rel)),
            t_stop = t_rel
        )
        
        return(out)
        
    } else if (Gfit_k$model == "xgb_cox") {
        
        if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("predict_Gk(): package 'xgboost' required for model='xgb_cox'.")
        }
        
        if (!is.null(Gfit_k$cov_use)) {
            miss <- setdiff(Gfit_k$cov_use, names(newdata_tk))
            if (length(miss) > 0) {
                stop("predict_Gk(): newdata_tk is missing covariates: ", paste(miss, collapse = ", "))
            }
        }
        
        dm_new <- .predict_design_matrix(
            dat = newdata_tk,
            cov_use = Gfit_k$cov_use,
            design_info = Gfit_k$design_info
        )
        if (any(!dm_new$keep)) {
            stop("predict_Gk(): missing values in covariates are not supported for xgb_cox prediction.")
        }
        
        lp <- as.numeric(stats::predict(Gfit_k$fit, newdata = dm_new$X, outputmargin = TRUE))
        lp[!is.finite(lp)] <- 0
        lp <- pmax(pmin(lp, 30), -30)
        lp_center <- if (!is.null(Gfit_k$lp_center) && is.finite(Gfit_k$lp_center)) Gfit_k$lp_center else 0
        lp <- lp - lp_center
        
        t_rel <- pmax(0, pmin(t_abs, t_k1) - t_k)
        H0_t <- approx(
            x = Gfit_k$bh_time, y = Gfit_k$bh_cumhaz, xout = t_rel,
            method = "constant", f = 0, yleft = 0, yright = max(Gfit_k$bh_cumhaz, 0),
            ties = "ordered"
        )$y
        out <- exp(-H0_t * exp(lp))
        out[out < 0] <- 0
        out[out > 1] <- 1
        
        return(out)
        
    } else if (Gfit_k$model == "xgb_aft") {
        
        if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("predict_Gk(): package 'xgboost' required for model='xgb_aft'.")
        }
        
        if (!is.null(Gfit_k$cov_use)) {
            miss <- setdiff(Gfit_k$cov_use, names(newdata_tk))
            if (length(miss) > 0) {
                stop("predict_Gk(): newdata_tk is missing covariates: ", paste(miss, collapse = ", "))
            }
        }
        
        dm_new <- .predict_design_matrix(
            dat = newdata_tk,
            cov_use = Gfit_k$cov_use,
            design_info = Gfit_k$design_info
        )
        if (any(!dm_new$keep)) {
            stop("predict_Gk(): missing values in covariates are not supported for xgb_aft prediction.")
        }
        
        mu <- as.numeric(stats::predict(Gfit_k$fit, newdata = dm_new$X))
        tiny <- .Machine$double.eps
        t_rel <- pmax(0, pmin(t_abs, t_k1) - t_k)
        t_pos <- pmax(t_rel, tiny)
        z <- (log(t_pos) - mu) / as.numeric(Gfit_k$aft_scale)
        
        Fz <- switch(
            Gfit_k$aft_dist,
            normal = stats::pnorm(z),
            logistic = stats::plogis(z),
            extreme = exp(-exp(-z)),
            stop("predict_Gk(): unsupported aft_dist = ", Gfit_k$aft_dist)
        )
        out <- 1 - Fz
        out[t_rel <= 0] <- 1
        out[out < 0] <- 0
        out[out > 1] <- 1
        
        return(out)
        
    } else {
        stop("Unknown Gfit_k$model")
    }
}


# Build global absolute-time grid from all fitted G_k's
get_G_jump_times_from_Gfits <- function(Gfits, include_zero = FALSE) {
    times_abs <- numeric(0)
    
    for (k in seq_along(Gfits)) {
        gk <- Gfits[[k]]
        if (is.null(gk)) next
        
        if (is.null(gk$t_k)) {
            stop("get_G_jump_times_from_Gfits(): Gfits[[", k, "]] is missing t_k.")
        }
        
        # Prefer explicitly stored jump times
        jr <- gk$jump_times_rel
        
        # Fallbacks
        if (is.null(jr)) {
            if (!is.null(gk$S_step)) {
                jr <- knots(gk$S_step)
            } else if (!is.null(gk$fit) && gk$model == "rsf" && !is.null(gk$fit$time.interest)) {
                jr <- gk$fit$time.interest
            } else {
                next
            }
        }
        
        jr <- jr[is.finite(jr) & jr > 0]
        if (length(jr) == 0) next
        
        ta <- gk$t_k + jr
        
        # Optional clip to interval end
        if (!is.null(gk$t_k1) && is.finite(gk$t_k1)) {
            ta <- ta[ta <= gk$t_k1 + 1e-12]
        }
        
        times_abs <- c(times_abs, ta)
    }
    
    times_abs <- sort(unique(times_abs))
    if (include_zero) times_abs <- sort(unique(c(0, times_abs)))
    
    times_abs
}



predict_G_from_Gk_matrix <- function(
        dat_new,
        time_vec = NULL,              # if NULL, use all jump times from fitted Gk's
        start.name, stop.name, event.name, event.C.name, id.name,
        Gfits,
        trim.G = 1e-7,
        return_subject_info = TRUE,   # return ids and X along with matrix
        fast_rsf_predict = FALSE      # opt-in: one rsf prediction per interval instead of
                                      # per evaluation time (identical values, ~big speedup,
                                      # but changes the R RNG draw count => downstream
                                      # randomized fits differ bitwise; default off)
) {
    if (is.null(Gfits) || length(Gfits) == 0) {
        stop("Gfits must be a non-empty list.")
    }
    
    K <- length(Gfits)
    visit_times <- vapply(Gfits, function(g) if (is.null(g)) NA_real_ else g$t_k, numeric(1))
    if (anyNA(visit_times)) {
        stop("predict_G_from_Gk_matrix(): every non-NULL Gfits[[k]] should store t_k.")
    }
    
    # Default evaluation times = union of jump times from fitted Gk's (absolute scale)
    if (is.null(time_vec)) {
        time_vec <- get_G_jump_times_from_Gfits(Gfits)
    }
    if (length(time_vec) == 0) {
        stop("time_vec is empty (and no jump times found in Gfits).")
    }
    
    # Sort times for interval-wise processing; restore original order at the end
    ord <- order(time_vec)
    time_sorted <- time_vec[ord]
    inv_ord <- order(ord)
    
    # Subject-level summary (X = observed/censored/admin time)
    summ <- subject_summary(dat_new, id.name, stop.name, event.name, event.C.name)
    ids <- summ[[id.name]]
    X   <- summ$X
    n   <- length(ids)
    m   <- length(time_sorted)
    
    # Output on sorted time grid first; initialize NA (requested for t > X_i)
    G_sorted <- matrix(NA_real_, nrow = n, ncol = m)
    rownames(G_sorted) <- as.character(ids)
    colnames(G_sorted) <- as.character(time_sorted)
    
    # ------------------------------------------------------------------
    # Cache covariate rows at each visit time t_k for subjects with X > t_k
    # ------------------------------------------------------------------
    cache <- vector("list", K)
    
    for (k in seq_len(K)) {
        t_k <- Gfits[[k]]$t_k
        idxk <- which(X > t_k)  # subjects at risk at t_k in new data
        
        if (length(idxk) == 0) {
            cache[[k]] <- list(idx = integer(0), dat = NULL)
            next
        }
        
        dat_slice <- dat_new[dat_new[[id.name]] %in% ids[idxk], , drop = FALSE]
        
        use_hist_k  <- isTRUE(Gfits[[k]]$use_history)
        cov_tv_k    <- Gfits[[k]]$cov_tv
        hist_mode_k <- Gfits[[k]]$history_mode
        
        if (use_hist_k && !is.null(cov_tv_k) && length(cov_tv_k) > 0) {
            dat_tk <- prepare_data_tau_with_history(
                dat = dat_slice, visit_times = visit_times, k = k,
                start.name = start.name, stop.name = stop.name, event.name = event.name, id.name = id.name,
                covars = cov_tv_k,
                mode = hist_mode_k
            )
        } else {
            dat_tk <- prepare_data_tau(
                dat = dat_slice,
                start.name = start.name, stop.name = stop.name, event.name = event.name, id.name = id.name,
                tau = t_k
            )
        }
        
        # Expect one row per id; check rather than silently dropping
        if (anyDuplicated(dat_tk[[id.name]])) {
            stop("predict_G_from_Gk_matrix(): duplicated ids found in cached dat_tk at k=", k)
        }
        
        # Align rows to ids[idxk]
        dat_tk <- dat_tk[match(ids[idxk], dat_tk[[id.name]]), , drop = FALSE]
        if (anyNA(dat_tk[[id.name]])) {
            stop("predict_G_from_Gk_matrix(): failed to align cached rows at k=", k)
        }
        
        cache[[k]] <- list(idx = idxk, dat = dat_tk)
    }
    
    # ------------------------------------------------------------------
    # Map each evaluation time to its interval k:
    #   k = max{j : visit_times[j] <= t}
    # For t in [t_k, t_{k+1}), k is the interval index.
    # rightmost.closed = FALSE is important so t = visit_times[K] maps to K.
    # ------------------------------------------------------------------
    k_of_t <- findInterval(
        time_sorted, visit_times,
        left.open = FALSE,
        rightmost.closed = FALSE
    )
    # k_of_t in {0,1,...,K}; 0 means t < visit_times[1]
    
    # If there are times before the first visit, define empty-product = 1 for t <= X_i
    cols_pre <- which(k_of_t == 0L)
    if (length(cols_pre) > 0) {
        for (cc in cols_pre) {
            tt <- time_sorted[cc]
            valid <- (tt <= X)
            G_sorted[valid, cc] <- 1
        }
    }
    
    # ------------------------------------------------------------------
    # Hoisted per-interval curve evaluators: for the km/cox/rsf backends the
    # subject-level censoring-survival curve does not depend on the evaluation
    # time, so run the model prediction ONCE per interval and evaluate the
    # resulting step functions at each requested time. Values are identical to
    # per-column predict_Gk() calls (same per-row computations); backends not
    # handled here fall back to predict_Gk() per column as before.
    # eval(rows, t_abs): survival at scalar time t_abs for datk[rows, ].
    make_Gk_evaluator <- function(gk, datk) {
        if (is.null(gk) || is.null(datk) || nrow(datk) == 0) return(NULL)
        t_k <- gk$t_k
        t_k1 <- gk$t_k1
        if (gk$model == "km") {
            function(rows, t_abs) {
                t_rel <- pmax(0, pmin(t_abs, t_k1) - t_k)
                gk$S_step(rep(t_rel, length(rows)))
            }
        } else if (gk$model == "cox") {
            if (!is.null(gk$cov_use)) {
                miss <- setdiff(gk$cov_use, names(datk))
                if (length(miss) > 0) {
                    stop("predict_Gk(): newdata_tk is missing covariates: ", paste(miss, collapse = ", "))
                }
            }
            X <- stats::model.matrix(gk$fit, datk)
            if (nrow(X) != nrow(datk)) {
                stop("predict_G_from_Gk_matrix(): model.matrix dropped ",
                     nrow(datk) - nrow(X), " of ", nrow(datk),
                     " cached rows at k with t_k=", t_k,
                     " (missing covariate values are not supported).")
            }
            r <- exp(drop(X %*% stats::coef(gk$fit)))
            function(rows, t_abs) {
                t_rel <- pmax(0, pmin(t_abs, t_k1) - t_k)
                gk$S0_step(rep(t_rel, length(rows)))^r[rows]
            }
        } else if (gk$model == "rsf" && isTRUE(fast_rsf_predict)) {
            # OPT-IN fast path (see fast_rsf_predict below / HAPS_FAST_RSF_PREDICT):
            # one forest prediction per interval instead of one per evaluation time.
            # The predicted values are identical to the per-column path, but
            # predict.rfsrc() consumes the R RNG stream on every call, so changing
            # the NUMBER of predict calls shifts the seeds of downstream randomized
            # fits => bitwise-different (statistically equivalent) results. Hence
            # DEFAULT OFF to keep outputs bit-identical to the original code.
            if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
                stop("predict_Gk(): package 'randomForestSRC' required for model='rsf'.")
            }
            if (!is.null(gk$cov_use)) {
                miss <- setdiff(gk$cov_use, names(datk))
                if (length(miss) > 0) {
                    stop("predict_Gk(): newdata_tk is missing covariates: ", paste(miss, collapse = ", "))
                }
            }
            pred <- predict(gk$fit, newdata = datk)
            grid <- pred$time.interest
            S <- pred$survival
            if (NROW(S) != nrow(datk)) {
                stop("predict_G_from_Gk_matrix(): rsf prediction returned ", NROW(S),
                     " rows for ", nrow(datk),
                     " cached rows (missing covariate values are not supported).")
            }
            function(rows, t_abs) {
                t_rel <- pmax(0, pmin(t_abs, t_k1) - t_k)
                eval_step_surv(S[rows, , drop = FALSE], grid, t_rel = rep(t_rel, length(rows)))
            }
        } else {
            # No hoisted evaluator for 'rsf' (unless fast_rsf_predict=TRUE, above)
            # or other backends: rsf's predict.rfsrc() consumes the R RNG stream on
            # every call, so hoisting would change the draw count and shift all
            # downstream randomized fits (bitwise-different results, statistically
            # equivalent). The per-column predict_Gk() fallback keeps outputs
            # bit-identical. The cox/km evaluators above consume no RNG (verified).
            NULL
        }
    }

    # ------------------------------------------------------------------
    # Prefix product entering interval k:
    # prefix[i] = product_{j<k} G_j(t_{j+1}-t_j | L_j) for subject i
    # We update this recursively as k increases.
    # ------------------------------------------------------------------
    prefix <- rep(1, n)  # entering interval 1

    # Process intervals in order
    for (k in seq_len(K)) {
        cols_k <- which(k_of_t == k)
        gk <- Gfits[[k]]
        idxk <- cache[[k]]$idx
        datk <- cache[[k]]$dat

        # One model prediction per interval (NULL for backends without a
        # hoisted evaluator; those use the per-column predict_Gk fallback)
        eval_k <- if (length(idxk) > 0 && !is.null(gk) && (length(cols_k) > 0 || k < K)) {
            make_Gk_evaluator(gk, datk)
        } else {
            NULL
        }

        # Fill columns in interval k using:
        # G(t) = prefix * G_k(t - t_k)
        if (length(cols_k) > 0) {
            if (length(idxk) > 0 && !is.null(gk)) {
                for (cc in cols_k) {
                    tt <- time_sorted[cc]

                    # valid among subjects cached for this interval
                    keep_local <- (tt <= X[idxk])   # still require t <= X_i
                    if (!any(keep_local)) next

                    subj_idx <- idxk[keep_local]

                    if (!is.null(eval_k)) {
                        Gj <- eval_k(which(keep_local), tt)
                    } else {
                        dat_use  <- datk[keep_local, , drop = FALSE]
                        Gj <- predict_Gk(gk, dat_use, t_abs = rep(tt, length(subj_idx)))
                    }
                    if (!is.null(trim.G)) Gj <- pmax(Gj, trim.G)

                    G_sorted[subj_idx, cc] <- prefix[subj_idx] * Gj
                }
            }
            # else: remain NA (no available model / no subjects at risk)
        }

        # Update prefix for entering interval k+1 using full interval factor of G_k
        if (k < K) {
            t_next <- visit_times[k + 1]

            # Subjects who enter the next interval satisfy X > t_{k+1}
            if (length(idxk) > 0) {
                enter_next_local <- (X[idxk] > t_next)

                if (any(enter_next_local)) {
                    subj_next <- idxk[enter_next_local]

                    if (is.null(gk)) {
                        prefix[subj_next] <- NA_real_
                    } else if (!is.null(eval_k)) {
                        Gj_end <- eval_k(which(enter_next_local), t_next)
                        if (!is.null(trim.G)) Gj_end <- pmax(Gj_end, trim.G)

                        prefix[subj_next] <- prefix[subj_next] * Gj_end
                    } else {
                        Gj_end <- predict_Gk(
                            gk,
                            datk[enter_next_local, , drop = FALSE],
                            t_abs = rep(t_next, length(subj_next))
                        )
                        if (!is.null(trim.G)) Gj_end <- pmax(Gj_end, trim.G)

                        prefix[subj_next] <- prefix[subj_next] * Gj_end
                    }
                }
            }
        }
    }
    
    # Restore original column order
    Gmat <- G_sorted[, inv_ord, drop = FALSE]
    colnames(Gmat) <- as.character(time_vec)
    
    if (return_subject_info) {
        return(list(
            G = Gmat,
            ids = ids,
            X = X,
            time_vec = time_vec
        ))
    } else {
        return(Gmat)
    }
}





# Compute the conditional censoring survival probabilities H_k(X) for subjects at risk at t_k,
# where H_k(t) = \PP(C>t|C>t_k, \bar L_k)
compute_Hk_at_tk <- function(dat_long, visit_times, k,
                             start.name, stop.name, event.name, event.C.name, id.name,
                             Gfits, trim.C = 1e-7) {
    
    K <- length(visit_times)
    t_k <- visit_times[k]
    
    # subjects at risk at t_k with covariates at t_k and subject-level X, Delta, DeltaC
    summ <- subject_summary(dat_long, id.name, stop.name, event.name, event.C.name)
    
    dat_tk <- prepare_data_tau(dat_long, start.name, stop.name, event.name, id.name, tau = t_k)
    dat_tk <- dat_tk[, setdiff(names(dat_tk), c("X","Delta","DeltaC")), drop=FALSE]  # if they exist
    dat_tk <- merge(dat_tk, summ, by = id.name, all.x = TRUE, sort = FALSE)
    
    dat_tk <- dat_tk[dat_tk$X > t_k, , drop = FALSE]
    if (nrow(dat_tk) == 0) return(data.frame())
    
    ids <- dat_tk[[id.name]]
    X   <- dat_tk$X
    
    H <- rep(1, length(ids))
    
    for (j in k:K) {
        t_j <- visit_times[j]
        
        # only those with X > t_j contribute
        idx <- which(X > t_j)
        
        if (length(idx) == 0) next
        
        if (is.null(Gfits[[j]])) {
            stop("compute_Hk_at_tk(): Gfits[[", j, "]] is NULL but there are subjects with X > t_j. ",
                 "Cannot evaluate G_j on this dataset.")
        }
        
        dat_slice <- dat_long[dat_long[[id.name]] %in% ids[idx], , drop = FALSE]
        
        if (!is.null(Gfits[[j]]) && isTRUE(Gfits[[j]]$use_history)) {
            cov_tv <- Gfits[[j]]$cov_tv
            if (is.null(cov_tv) || length(cov_tv) == 0) {
                stop("compute_Hk_at_tk(): Gfits[[", j, "]] uses history but does not store base cov_tv. ",
                     "Store cov_tv in fit_Gk_list() or pass it into compute_Hk_*().")
            }
            dat_tj <- prepare_data_tau_with_history(
                dat_slice, visit_times, j,
                start.name, stop.name, event.name, id.name,
                covars = cov_tv,
                mode = Gfits[[j]]$history_mode
            )
        } else {
            dat_tj <- prepare_data_tau(dat_slice, start.name, stop.name, event.name, id.name, tau = t_j)
        }
        
        
        dat_tj <- dat_tj[, setdiff(names(dat_tj), c("X","Delta","DeltaC")), drop=FALSE]  # if they exist
        dat_tj <- merge(dat_tj, summ, by = id.name, all.x = TRUE, sort = FALSE)

        # align to ids[idx]
        idx_align <- match(ids[idx], dat_tj[[id.name]])
        if (anyNA(idx_align)) {
            miss_ids <- ids[idx][is.na(idx_align)]
            stop("compute_Hk_at_tk(): no interval-", j, " covariate row for some at-risk ids ",
                 "(an unmatched id would inject an all-NA row and silently corrupt H_k). Example: ",
                 paste(utils::head(miss_ids, 10), collapse = ", "))
        }
        dat_tj <- dat_tj[idx_align, , drop = FALSE]

        t_eval_abs <- pmin(X[idx], if (j < K) visit_times[j + 1] else X[idx])
        Gj <- predict_Gk(Gfits[[j]], dat_tj, t_eval_abs)
        
        H[idx] <- H[idx] * pmax(Gj, trim.C)
    }
    
    out <- data.frame(id = ids, X = X, Hk = H)
    out
}



# return a list of H_k for all k
compute_Hk_list <- function(dat_long, visit_times, start.name, stop.name, event.name, event.C.name, id.name,
                            Gfits, trim.C = 1e-6) {
    K <- length(visit_times)
    out <- vector("list", K)
    for (k in 1:K) {
        out[[k]] <- compute_Hk_at_tk(dat_long, visit_times, k,
                                     start.name, stop.name, event.name, event.C.name, id.name,
                                     Gfits, trim.C = trim.C)
    }
    out
}


## Fit S_k(t|\bar L_k) on training data.
## Default uses the current event-only IPCW fit; optional all-at-risk mode
## keeps censored subjects in the risk set and relies on right-censoring.
## The hybrid mode uses all-at-risk fitting on (t_k, t_{k+1}] and an IPCW
## tail fit beyond t_{k+1} while keeping predictors fixed at L_k.
fit_Sk_list <- function(dat_tr, visit_times,
                        start.name, stop.name, event.name, event.C.name, id.name,
                        covname.S.baseline = NULL, covname.S.timevarying = NULL,
                        Hk_list_tr,  # output from compute_Hk_list(dat_tr,...)
                        model.S = c("cox", "km", "rsf", "hal", "xgb_cox", "xgb_aft"),
                        trim.H = 1e-7,
                        fit_sample_mode = c("event_only_ipcw", "all_at_risk", "hybrid_interval_ipcw"),
                        small_event_threshold = 10L,
                        small_event_min_events_per_coef = 5L,
                        small_event_fallback_model = "km",
                        small_event_metric = c("min_count_ess", "count", "ess"),
                        small_event_guard_models = c("cox", "rsf", "xgb_cox", "xgb_aft"),
                        coxph_args = list(),   # extra args passed to survival::coxph
                        rsf_args   = list(ntree = 1000L, nodesize = 30L, nsplit = 10L, forest = TRUE),
                        use_history = FALSE,
                        history_mode = c("wide", "summary"),
                        xgb_nrounds = 200L,
                        xgb_use_cv = FALSE,
                        xgb_cv_nfold = 5L,
                        xgb_cv_early_stopping_rounds = 20L,
                        xgb_cv_param_grid = NULL,
                        xgb_cv_seed = NULL,
                        xgb_cv_verbose = 0
) {

    small_event_min_events_per_coef_missing <- missing(small_event_min_events_per_coef)
    model.S <- match.arg(model.S)
    fit_sample_mode <- match.arg(fit_sample_mode)
    history_mode <- match.arg(history_mode)
    if (!is.list(rsf_args)) stop("fit_Sk_list(): rsf_args must be a list.")
    small_event_threshold <- as.integer(small_event_threshold[[1L]])
    if (!is.finite(small_event_threshold) || small_event_threshold < 0L) {
        stop("fit_Sk_list(): small_event_threshold must be a nonnegative integer.")
    }
    if (small_event_threshold == 0L && small_event_min_events_per_coef_missing) {
        small_event_min_events_per_coef <- 0L
    }
    small_event_min_events_per_coef <- as.numeric(small_event_min_events_per_coef[[1L]])
    if (!is.finite(small_event_min_events_per_coef) || small_event_min_events_per_coef < 0) {
        stop("fit_Sk_list(): small_event_min_events_per_coef must be nonnegative.")
    }
    small_event_fallback_model <- as.character(small_event_fallback_model[[1L]])
    if (!small_event_fallback_model %in% c("cox", "km")) {
        stop("fit_Sk_list(): unsupported small_event_fallback_model. Supported: 'cox', 'km'.")
    }
    small_event_metric <- match.arg(small_event_metric)
    small_event_guard_models <- unique(as.character(small_event_guard_models))
    bad_guard_models <- setdiff(small_event_guard_models, c("cox", "rsf", "hal", "grf", "xgb_cox", "xgb_aft"))
    if (length(bad_guard_models) > 0L) {
        stop(
            "fit_Sk_list(): unsupported small_event_guard_models: ",
            paste(bad_guard_models, collapse = ", ")
        )
    }

    cov_base <- covname.S.baseline
    cov_tv   <- covname.S.timevarying

    K <- length(visit_times)
    summ <- subject_summary(dat_tr, id.name, stop.name, event.name, event.C.name)

    Sfits <- vector("list", K)
    names(Sfits) <- paste0("k=", seq_len(K))

    fit_single_S_model <- function(dat_landmark,
                                   cov_use,
                                   t_origin,
                                   fit_mode = c("event_only_ipcw", "all_at_risk"),
                                   fit_weight_scheme,
                                   n_all_reference = nrow(dat_landmark)) {
        fit_mode <- match.arg(fit_mode)
        rhs <- if (length(cov_use) == 0L) "1" else paste(cov_use, collapse = "+")
        form <- stats::as.formula(paste0("Surv(timeT_rel, Delta) ~ ", rhs))

        dat_fit_base <- dat_landmark
        dat_fit_base$timeT_rel <- dat_fit_base$X - t_origin

        if (fit_mode == "event_only_ipcw") {
            dat_fit <- dat_fit_base[dat_fit_base$Delta == 1, , drop = FALSE]
            fit_weights <- dat_fit$w_Sk
        } else {
            dat_fit <- dat_fit_base
            fit_weights <- NULL
        }
        if (nrow(dat_fit) == 0L) return(NULL)
        n_events_fit <- sum(dat_fit$Delta == 1)
        n_coef_fit <- .covariate_design_ncoef(dat_fit, cov_use)
        if (is.null(fit_weights)) {
            w_events <- rep(1, n_events_fit)
        } else {
            w_events <- as.numeric(fit_weights[dat_fit$Delta == 1])
            w_events <- w_events[is.finite(w_events) & w_events > 0]
        }
        event_ess_fit <- if (length(w_events) == 0L) {
            0
        } else {
            (sum(w_events)^2) / sum(w_events^2)
        }
        small_event_stat <- switch(
            small_event_metric,
            count = as.numeric(n_events_fit),
            ess = as.numeric(event_ess_fit),
            min_count_ess = min(as.numeric(n_events_fit), as.numeric(event_ess_fit))
        )
        small_event_threshold_used <- .event_threshold_by_dim(
            small_event_threshold,
            small_event_min_events_per_coef,
            n_coef_fit
        )
        effective_model <- model.S
        fallback_reason <- NULL
        if (small_event_threshold_used > 0L &&
            model.S %in% small_event_guard_models &&
            small_event_stat < small_event_threshold_used &&
            model.S != small_event_fallback_model) {
            effective_model <- small_event_fallback_model
            fallback_reason <- paste0(
                "few_events_",
                small_event_metric,
                "_lt_",
                small_event_threshold_used
            )
        }

        if (effective_model == "km") {
            km_form <- stats::as.formula("Surv(timeT_rel, Delta) ~ 1")
            km_args <- list(formula = km_form, data = dat_fit)
            if (!is.null(fit_weights)) km_args$weights <- fit_weights
            fit <- do.call(survival::survfit, km_args)
            ss <- summary(fit)
            out <- list(
                model = "km",
                t_k = t_origin,
                fit = fit,
                grid = ss$time,
                S_grid = ss$surv,
                cov_use = character(0),
                requested_cov_use = cov_use,
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base,
                n_all = n_all_reference,
                n_used = nrow(dat_fit),
                n_events_used = n_events_fit,
                event_ess_used = as.numeric(event_ess_fit),
                n_coef = n_coef_fit,
                small_event_metric = small_event_metric,
                small_event_stat = as.numeric(small_event_stat),
                small_event_threshold = small_event_threshold_used,
                small_event_min_events_per_coef = small_event_min_events_per_coef,
                weight_scheme = fit_weight_scheme
            )
            if (!identical(effective_model, model.S)) out$requested_model <- model.S
            if (!is.null(fallback_reason)) out$fallback_reason <- fallback_reason
            return(out)

        } else if (effective_model == "cox") {
            cox_args <- list(formula = form, data = dat_fit, x = TRUE)
            if (!is.null(fit_weights)) cox_args$weights <- fit_weights
            fit <- do.call(survival::coxph, c(cox_args, coxph_args))

            out <- list(
                model = "cox",
                t_k = t_origin,
                fit = fit,
                cov_use = cov_use,
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base,
                n_all = n_all_reference,
                n_used = nrow(dat_fit),
                n_events_used = n_events_fit,
                event_ess_used = as.numeric(event_ess_fit),
                n_coef = n_coef_fit,
                small_event_metric = small_event_metric,
                small_event_stat = as.numeric(small_event_stat),
                small_event_threshold = small_event_threshold_used,
                small_event_min_events_per_coef = small_event_min_events_per_coef,
                weight_scheme = fit_weight_scheme
            )
            if (!identical(effective_model, model.S)) out$requested_model <- model.S
            if (!is.null(fallback_reason)) out$fallback_reason <- fallback_reason
            return(out)

        } else if (effective_model == "rsf") {
            if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
                stop("Package 'randomForestSRC' is required for model.S='rsf'. Install via install.packages('randomForestSRC').")
            }

            rsf_fit_args <- utils::modifyList(
                list(
                    formula = form,
                    data = dat_fit,
                    ntree = 1000L,
                    nodesize = 30L,
                    nsplit = 10L,
                    forest = TRUE
                ),
                rsf_args
            )
            rsf_fit_args$formula <- form
            rsf_fit_args$data <- dat_fit
            if (!is.null(fit_weights)) {
                rsf_fit_args$case.wt <- fit_weights
            } else if (!is.null(rsf_fit_args$case.wt)) {
                rsf_fit_args$case.wt <- NULL
            }
            rsf_fit_args$forest <- TRUE
            fit <- do.call(randomForestSRC::rfsrc, rsf_fit_args)

            out <- list(
                model = "rsf",
                t_k = t_origin,
                fit = fit,
                cov_use = cov_use,
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base,
                n_all = n_all_reference,
                n_used = nrow(dat_fit),
                n_events_used = n_events_fit,
                event_ess_used = as.numeric(event_ess_fit),
                n_coef = n_coef_fit,
                small_event_metric = small_event_metric,
                small_event_stat = as.numeric(small_event_stat),
                small_event_threshold = small_event_threshold_used,
                small_event_min_events_per_coef = small_event_min_events_per_coef,
                weight_scheme = fit_weight_scheme
            )
            if (!identical(effective_model, model.S)) out$requested_model <- model.S
            if (!is.null(fallback_reason)) out$fallback_reason <- fallback_reason
            return(out)

        } else if (effective_model == "hal") {
            if (length(cov_use) == 0L) {
                fit <- do.call(
                    survival::coxph,
                    c({
                        cox_args <- list(formula = form, data = dat_fit, x = TRUE)
                        if (!is.null(fit_weights)) cox_args$weights <- fit_weights
                        cox_args
                    }, coxph_args)
                )

                return(list(
                    model = "cox",
                    t_k = t_origin,
                    fit = fit,
                    cov_use = cov_use,
                    use_history = use_history,
                    history_mode = if (use_history) history_mode else NULL,
                    cov_tv = cov_tv,
                    cov_base = cov_base,
                    n_all = n_all_reference,
                    n_used = nrow(dat_fit),
                    n_events_used = n_events_fit,
                    event_ess_used = as.numeric(event_ess_fit),
                    n_coef = n_coef_fit,
                    small_event_metric = small_event_metric,
                    small_event_stat = as.numeric(small_event_stat),
                    small_event_threshold = small_event_threshold_used,
                    small_event_min_events_per_coef = small_event_min_events_per_coef,
                    weight_scheme = fit_weight_scheme,
                    requested_model = model.S,
                    fallback_reason = if (is.null(fallback_reason)) "no_usable_covariates" else fallback_reason
                ))
            }

            Y_surv <- survival::Surv(dat_fit$timeT_rel, dat_fit$Delta)
            hal_bundle <- .fit_hal_cox_bundle(
                dat = dat_fit,
                cov_use = cov_use,
                Y_surv = Y_surv,
                weights = fit_weights
            )

            out <- list(
                model = "hal",
                t_k = t_origin,
                fit = hal_bundle$fit_hal,
                lp_center = hal_bundle$lp_center,
                bh_time = hal_bundle$bh_time,
                bh_cumhaz = hal_bundle$bh_cumhaz,
                design_info = hal_bundle$design_info,
                cov_use = cov_use,
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base,
                n_all = n_all_reference,
                n_used = nrow(dat_fit),
                n_events_used = n_events_fit,
                event_ess_used = as.numeric(event_ess_fit),
                n_coef = n_coef_fit,
                small_event_metric = small_event_metric,
                small_event_stat = as.numeric(small_event_stat),
                small_event_threshold = small_event_threshold_used,
                small_event_min_events_per_coef = small_event_min_events_per_coef,
                weight_scheme = fit_weight_scheme
            )
            if (!identical(effective_model, model.S)) out$requested_model <- model.S
            if (!is.null(fallback_reason)) out$fallback_reason <- fallback_reason
            return(out)

        } else if (effective_model == "grf") {
            if (!requireNamespace("grf", quietly = TRUE)) {
                stop("Package 'grf' is required for model.S='grf'. Install via install.packages('grf').")
            }

            if (length(cov_use) == 0L) {
                fit <- do.call(
                    survival::coxph,
                    c({
                        cox_args <- list(formula = form, data = dat_fit, x = TRUE)
                        if (!is.null(fit_weights)) cox_args$weights <- fit_weights
                        cox_args
                    }, coxph_args)
                )

                return(list(
                    model = "cox",
                    t_k = t_origin,
                    fit = fit,
                    cov_use = cov_use,
                    use_history = use_history,
                    history_mode = if (use_history) history_mode else NULL,
                    cov_tv = cov_tv,
                    cov_base = cov_base,
                    n_all = n_all_reference,
                    n_used = nrow(dat_fit),
                    n_events_used = n_events_fit,
                    event_ess_used = as.numeric(event_ess_fit),
                    n_coef = n_coef_fit,
                    small_event_metric = small_event_metric,
                    small_event_stat = as.numeric(small_event_stat),
                    small_event_threshold = small_event_threshold_used,
                    small_event_min_events_per_coef = small_event_min_events_per_coef,
                    weight_scheme = fit_weight_scheme,
                    requested_model = model.S,
                    fallback_reason = if (is.null(fallback_reason)) "no_usable_covariates" else fallback_reason
                ))
            }

            dm_fit <- .fit_design_matrix(dat_fit, cov_use = cov_use)
            if (is.null(dm_fit$X) || nrow(dm_fit$X) == 0) {
                stop("fit_Sk_list(): no complete-case rows left after design matrix construction for model.S='grf'.")
            }

            dat_fit_ml <- dat_fit[dm_fit$keep, , drop = FALSE]
            X_ml <- dm_fit$X
            w_ml <- fit_weights
            if (!is.null(w_ml)) w_ml <- as.numeric(w_ml[dm_fit$keep])
            if (!is.null(w_ml) && all(is.finite(w_ml)) && sum(w_ml) > 0) {
                w_ml <- w_ml / mean(w_ml)
            }

            grf_args <- list(
                X = X_ml,
                Y = as.numeric(dat_fit_ml$timeT_rel),
                D = as.numeric(dat_fit_ml$Delta)
            )
            if (!is.null(w_ml)) grf_args$sample.weights <- w_ml
            fit_grf <- do.call(grf::survival_forest, grf_args)

            out <- list(
                model = "grf",
                t_k = t_origin,
                fit = fit_grf,
                cov_use = cov_use,
                design_info = dm_fit$design_info,
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base,
                n_all = n_all_reference,
                n_used = nrow(dat_fit_ml),
                n_events_used = sum(dat_fit_ml$Delta == 1),
                event_ess_used = as.numeric(event_ess_fit),
                n_coef = n_coef_fit,
                small_event_metric = small_event_metric,
                small_event_stat = as.numeric(small_event_stat),
                small_event_threshold = small_event_threshold_used,
                small_event_min_events_per_coef = small_event_min_events_per_coef,
                weight_scheme = fit_weight_scheme
            )
            if (!identical(effective_model, model.S)) out$requested_model <- model.S
            if (!is.null(fallback_reason)) out$fallback_reason <- fallback_reason
            return(out)

        } else if (effective_model %in% c("xgb_cox", "xgb_aft")) {
            if (!requireNamespace("xgboost", quietly = TRUE)) {
                stop("Package 'xgboost' is required for model.S='", effective_model, "'.")
            }

            if (length(cov_use) == 0L) {
                fit <- do.call(
                    survival::coxph,
                    c({
                        cox_args <- list(formula = form, data = dat_fit, x = TRUE)
                        if (!is.null(fit_weights)) cox_args$weights <- fit_weights
                        cox_args
                    }, coxph_args)
                )

                return(list(
                    model = "cox",
                    t_k = t_origin,
                    fit = fit,
                    cov_use = cov_use,
                    use_history = use_history,
                    history_mode = if (use_history) history_mode else NULL,
                    cov_tv = cov_tv,
                    cov_base = cov_base,
                    n_all = n_all_reference,
                    n_used = nrow(dat_fit),
                    n_events_used = n_events_fit,
                    event_ess_used = as.numeric(event_ess_fit),
                    n_coef = n_coef_fit,
                    small_event_metric = small_event_metric,
                    small_event_stat = as.numeric(small_event_stat),
                    small_event_threshold = small_event_threshold_used,
                    small_event_min_events_per_coef = small_event_min_events_per_coef,
                    weight_scheme = fit_weight_scheme,
                    requested_model = model.S,
                    fallback_reason = if (is.null(fallback_reason)) "no_usable_covariates" else fallback_reason
                ))
            }

            dm_fit <- .fit_design_matrix(dat_fit, cov_use = cov_use)
            if (is.null(dm_fit$X) || nrow(dm_fit$X) == 0) {
                stop("fit_Sk_list(): no complete-case rows left after design matrix construction for model.S='", effective_model, "'.")
            }

            dat_fit_ml <- dat_fit[dm_fit$keep, , drop = FALSE]
            X_ml <- dm_fit$X
            w_ml <- fit_weights
            if (!is.null(w_ml)) w_ml <- as.numeric(w_ml[dm_fit$keep])
            if (!is.null(w_ml) && all(is.finite(w_ml)) && sum(w_ml) > 0) {
                w_ml <- w_ml / mean(w_ml)
            }

            xgb_params_common <- list(
                eta = 0.05,
                max_depth = 4L,
                min_child_weight = 1,
                subsample = 0.8,
                colsample_bytree = 0.8,
                lambda = 1,
                gamma = 0
            )

            if (effective_model == "xgb_cox") {
                y_cox <- ifelse(dat_fit_ml$Delta == 1L,
                                as.numeric(dat_fit_ml$timeT_rel),
                                -as.numeric(dat_fit_ml$timeT_rel))
                dtrain <- xgboost::xgb.DMatrix(
                    data = X_ml,
                    label = y_cox
                )
                if (!is.null(w_ml)) xgboost::setinfo(dtrain, "weight", w_ml)
                params <- c(list(
                    objective = "survival:cox",
                    eval_metric = "cox-nloglik"
                ), xgb_params_common)
                xgb_fit <- .xgb_train_optional_cv(
                    dtrain = dtrain,
                    params = params,
                    nrounds = as.integer(xgb_nrounds),
                    use_cv = xgb_use_cv,
                    cv_nfold = xgb_cv_nfold,
                    cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
                    cv_param_grid = xgb_cv_param_grid,
                    cv_seed = xgb_cv_seed,
                    cv_verbose = xgb_cv_verbose,
                    maximize = FALSE,
                    train_verbose = 0
                )
                fit_xgb <- xgb_fit$fit

                lp_fit <- as.numeric(stats::predict(fit_xgb, newdata = X_ml, outputmargin = TRUE))
                lp_fit[!is.finite(lp_fit)] <- 0
                lp_fit <- pmax(pmin(lp_fit, 30), -30)
                lp_center <- mean(lp_fit)
                dat_anchor <- data.frame(
                    timeT_rel = dat_fit_ml$timeT_rel,
                    Delta = dat_fit_ml$Delta,
                    lp_off = lp_fit - lp_center
                )
                anchor_args <- list(
                    formula = survival::Surv(timeT_rel, Delta) ~ offset(lp_off),
                    data = dat_anchor,
                    ties = "breslow",
                    x = TRUE
                )
                if (!is.null(w_ml)) anchor_args$weights <- w_ml
                fit_anchor <- do.call(survival::coxph, anchor_args)
                bh <- survival::basehaz(fit_anchor, centered = FALSE)

                out <- list(
                    model = "xgb_cox",
                    t_k = t_origin,
                    fit = fit_xgb,
                    tuning = xgb_fit$tuning,
                    lp_center = lp_center,
                    bh_time = bh$time,
                    bh_cumhaz = bh$hazard,
                    cov_use = cov_use,
                    design_info = dm_fit$design_info,
                    use_history = use_history,
                    history_mode = if (use_history) history_mode else NULL,
                    cov_tv = cov_tv,
                    cov_base = cov_base,
                    n_all = n_all_reference,
                    n_used = nrow(dat_fit_ml),
                    n_events_used = sum(dat_fit_ml$Delta == 1),
                    event_ess_used = as.numeric(event_ess_fit),
                    n_coef = n_coef_fit,
                    small_event_metric = small_event_metric,
                    small_event_stat = as.numeric(small_event_stat),
                    small_event_threshold = small_event_threshold_used,
                    small_event_min_events_per_coef = small_event_min_events_per_coef,
                    weight_scheme = fit_weight_scheme
                )
                if (!identical(effective_model, model.S)) out$requested_model <- model.S
                if (!is.null(fallback_reason)) out$fallback_reason <- fallback_reason
                return(out)
            }

            dtrain <- xgboost::xgb.DMatrix(data = X_ml)
            if (!is.null(w_ml)) xgboost::setinfo(dtrain, "weight", w_ml)
            xgboost::setinfo(dtrain, "label_lower_bound", as.numeric(dat_fit_ml$timeT_rel))
            xgboost::setinfo(
                dtrain,
                "label_upper_bound",
                ifelse(dat_fit_ml$Delta == 1L, as.numeric(dat_fit_ml$timeT_rel), Inf)
            )

            aft_dist <- "normal"
            aft_scale <- 1
            params <- c(list(
                objective = "survival:aft",
                eval_metric = "aft-nloglik",
                aft_loss_distribution = aft_dist,
                aft_loss_distribution_scale = aft_scale
            ), xgb_params_common)

            xgb_fit <- .xgb_train_optional_cv(
                dtrain = dtrain,
                params = params,
                nrounds = as.integer(xgb_nrounds),
                use_cv = xgb_use_cv,
                cv_nfold = xgb_cv_nfold,
                cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
                cv_param_grid = xgb_cv_param_grid,
                cv_seed = xgb_cv_seed,
                cv_verbose = xgb_cv_verbose,
                maximize = FALSE,
                train_verbose = 0
            )
            fit_xgb <- xgb_fit$fit

            out <- list(
                model = "xgb_aft",
                t_k = t_origin,
                fit = fit_xgb,
                tuning = xgb_fit$tuning,
                aft_dist = aft_dist,
                aft_scale = aft_scale,
                cov_use = cov_use,
                design_info = dm_fit$design_info,
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base,
                n_all = n_all_reference,
                n_used = nrow(dat_fit_ml),
                n_events_used = sum(dat_fit_ml$Delta == 1),
                event_ess_used = as.numeric(event_ess_fit),
                n_coef = n_coef_fit,
                small_event_metric = small_event_metric,
                small_event_stat = as.numeric(small_event_stat),
                small_event_threshold = small_event_threshold_used,
                small_event_min_events_per_coef = small_event_min_events_per_coef,
                weight_scheme = fit_weight_scheme
            )
            if (!identical(effective_model, model.S)) out$requested_model <- model.S
            if (!is.null(fallback_reason)) out$fallback_reason <- fallback_reason
            return(out)
        }

        stop("Unsupported model.S: ", effective_model)
    }

    for (k in seq_len(K)) {
        t_k <- visit_times[k]

        if (use_history && length(cov_tv) > 0) {
            dat_tk <- prepare_data_tau_with_history(
                dat_tr, visit_times, k,
                start.name, stop.name, event.name, id.name,
                covars = cov_tv,
                mode = history_mode
            )
        } else {
            dat_tk <- prepare_data_tau(dat_tr, start.name, stop.name, event.name, id.name, tau = t_k)
        }

        dat_tk <- dat_tk[, setdiff(names(dat_tk), c("X", "Delta", "DeltaC")), drop = FALSE]
        dat_tk <- merge(dat_tk, summ, by = id.name, all.x = TRUE, sort = FALSE)
        dat_tk <- dat_tk[dat_tk$X > t_k, , drop = FALSE]
        if (nrow(dat_tk) == 0L) {
            Sfits[k] <- list(NULL)
            next
        }

        if (!use_history) {
            cov_use <- c(cov_base, cov_tv)
        } else {
            if (history_mode == "wide") {
                cov_hist <- unlist(lapply(seq_len(k), function(j) paste0(cov_tv, "_t", j)))
            } else if (history_mode == "summary") {
                cov_hist <- unlist(lapply(cov_tv, function(v) paste0(v, c("_hist_mean", "_hist_last", "_hist_slope"))))
            } else {
                stop(paste0("Unknown history_mode: ", history_mode))
            }
            cov_use <- c(cov_base, cov_hist)
        }
        cov_use <- keep_nonconstant_covars(dat_tk, cov_use)

        if (fit_sample_mode == "hybrid_interval_ipcw") {
            bridge_time <- if (k < K) visit_times[k + 1L] else Inf

            dat_short <- dat_tk
            if (is.finite(bridge_time)) {
                X_orig <- dat_short$X
                Delta_orig <- dat_short$Delta
                dat_short$X <- pmin(X_orig, bridge_time)
                dat_short$Delta <- as.integer(Delta_orig == 1L & X_orig <= bridge_time)
            }

            fit_short <- fit_single_S_model(
                dat_landmark = dat_short,
                cov_use = cov_use,
                t_origin = t_k,
                fit_mode = "all_at_risk",
                fit_weight_scheme = if (is.finite(bridge_time)) "all_at_risk_short" else "all_at_risk",
                n_all_reference = nrow(dat_tk)
            )
            if (is.null(fit_short)) {
                Sfits[k] <- list(NULL)
                next
            }

            fit_tail <- NULL
            n_all_tail <- 0L
            if (is.finite(bridge_time)) {
                dat_tail <- dat_tk[dat_tk$X > bridge_time, , drop = FALSE]
                n_all_tail <- nrow(dat_tail)
                if (n_all_tail > 0L) {
                    Hdf_tail <- Hk_list_tr[[k + 1L]]
                    dat_tail <- merge(
                        dat_tail,
                        Hdf_tail[, c("id", "Hk")],
                        by.x = id.name, by.y = "id", all.x = TRUE, sort = FALSE
                    )
                    if (anyNA(dat_tail$Hk)) {
                        stop("fit_Sk_list(): missing Hk after merge for hybrid tail fit at k=", k)
                    }
                    dat_tail$w_Sk <- dat_tail$Delta / pmax(dat_tail$Hk, trim.H)
                    fit_tail <- fit_single_S_model(
                        dat_landmark = dat_tail,
                        cov_use = cov_use,
                        t_origin = bridge_time,
                        fit_mode = "event_only_ipcw",
                        fit_weight_scheme = "event_only_ipcw_tail",
                        n_all_reference = n_all_tail
                    )
                }
            }

            Sfits[[k]] <- list(
                model = "hybrid_interval_ipcw",
                component_model = model.S,
                t_k = t_k,
                bridge_time = bridge_time,
                short_fit = fit_short,
                tail_fit = fit_tail,
                cov_use = cov_use,
                use_history = use_history,
                history_mode = if (use_history) history_mode else NULL,
                cov_tv = cov_tv,
                cov_base = cov_base,
                n_all = nrow(dat_tk),
                n_used = fit_short$n_used,
                n_events_used = fit_short$n_events_used,
                n_all_tail = n_all_tail,
                n_used_tail = if (!is.null(fit_tail)) fit_tail$n_used else 0L,
                n_events_used_tail = if (!is.null(fit_tail)) fit_tail$n_events_used else 0L,
                weight_scheme = "hybrid_interval_ipcw"
            )
            next
        }

        Hdf <- Hk_list_tr[[k]]
        dat_tk <- merge(dat_tk, Hdf[, c("id", "Hk")],
                        by.x = id.name, by.y = "id", all.x = TRUE, sort = FALSE)
        if (anyNA(dat_tk$Hk)) stop("fit_Sk_list(): missing Hk after merge at k=", k)
        dat_tk$w_Sk <- dat_tk$Delta / pmax(dat_tk$Hk, trim.H)

        fit_k <- fit_single_S_model(
            dat_landmark = dat_tk,
            cov_use = cov_use,
            t_origin = t_k,
            fit_mode = fit_sample_mode,
            fit_weight_scheme = fit_sample_mode,
            n_all_reference = nrow(dat_tk)
        )
        Sfits[k] <- list(fit_k)
    }

    Sfits
}







## ------------------------------------------------------------------
## Helpers for ML design matrix construction (xgboost / grf)
## ------------------------------------------------------------------

.make_cov_formula_no_intercept <- function(cov_use) {
    if (length(cov_use) == 0) return(NULL)
    stats::as.formula(paste0("~ 0 + ", paste(cov_use, collapse = " + ")))
}

.fit_design_matrix <- function(dat, cov_use) {
    # Returns numeric design matrix + metadata to rebuild on newdata.
    # Drops rows with missing raw covariates first (conservative / robust).
    if (length(cov_use) == 0) {
        return(list(
            X = NULL,
            keep = rep(TRUE, nrow(dat)),
            design_info = NULL
        ))
    }
    
    keep <- stats::complete.cases(dat[, cov_use, drop = FALSE])
    dat_cc <- dat[keep, , drop = FALSE]
    if (nrow(dat_cc) == 0) {
        return(list(
            X = matrix(numeric(0), nrow = 0, ncol = 0),
            keep = keep,
            design_info = NULL
        ))
    }
    
    f_rhs <- .make_cov_formula_no_intercept(cov_use)
    
    # model.frame + xlev to handle factor levels consistently at prediction time
    mf <- stats::model.frame(f_rhs, data = dat_cc, na.action = stats::na.pass,
                             drop.unused.levels = FALSE)
    tt <- attr(mf, "terms")
    xlev <- stats::.getXlevels(tt, mf)
    X <- stats::model.matrix(tt, data = mf)
    
    # Ensure numeric matrix
    storage.mode(X) <- "double"
    
    list(
        X = X,
        keep = keep,
        design_info = list(
            terms = tt,
            xlev = xlev,
            colnames = colnames(X)
        )
    )
}

.predict_design_matrix <- function(dat, cov_use, design_info) {
    # Rebuild design matrix using stored terms/xlev, align columns to training matrix.
    if (length(cov_use) == 0) {
        return(list(
            X = NULL,
            keep = rep(TRUE, nrow(dat))
        ))
    }
    
    keep <- stats::complete.cases(dat[, cov_use, drop = FALSE])
    dat_cc <- dat[keep, , drop = FALSE]
    if (nrow(dat_cc) == 0) {
        return(list(
            X = matrix(numeric(0), nrow = 0, ncol = length(design_info$colnames),
                       dimnames = list(NULL, design_info$colnames)),
            keep = keep
        ))
    }
    
    # If new factor levels appear, model.frame may error.
    # Catch and return informative message.
    mf <- tryCatch(
        stats::model.frame(design_info$terms, data = dat_cc, xlev = design_info$xlev,
                           na.action = stats::na.pass, drop.unused.levels = FALSE),
        error = function(e) {
            stop("predict_Yk_condmean_at_tk(): failed to build model frame on new data. ",
                 "Possible unseen factor levels in covariates. Original error: ", e$message)
        }
    )
    
    X_new_raw <- stats::model.matrix(design_info$terms, data = mf)
    storage.mode(X_new_raw) <- "double"
    
    # Align columns exactly to training design columns
    train_cols <- design_info$colnames
    X_new <- matrix(0, nrow = nrow(X_new_raw), ncol = length(train_cols),
                    dimnames = list(NULL, train_cols))
    common <- intersect(colnames(X_new_raw), train_cols)
    if (length(common) > 0) {
        X_new[, common] <- X_new_raw[, common, drop = FALSE]
    }
    
    list(X = X_new, keep = keep)
}

.coerce_xgb_multiclass_pred <- function(
        pred, n, num_class,
        check_prob = TRUE,
        prob_tol = 1e-6
) {
    # Robust to xgboost versions returning vector or matrix/array.
    # Returns an n x num_class numeric matrix with rows = observations.
    
    if (is.null(pred)) stop("xgboost prediction returned NULL.")
    
    # Case 1: already a matrix
    if (is.matrix(pred)) {
        if (nrow(pred) == n && ncol(pred) == num_class) {
            pmat <- pred
        } else if (ncol(pred) == n && nrow(pred) == num_class) {
            pmat <- t(pred)
        } else {
            # fall through to vector reshape
            pred <- as.vector(pred)
            pmat <- NULL
        }
    } else {
        pmat <- NULL
    }
    
    # Case 2: array-like or vector -> flatten
    if (is.null(pmat)) {
        if (length(dim(pred)) >= 2) {
            pred <- as.vector(pred)
        }
        if (!is.numeric(pred)) pred <- as.numeric(pred)
        
        if (length(pred) != n * num_class) {
            stop("Unexpected xgboost multiclass prediction length: got ", length(pred),
                 ", expected ", n * num_class, ".")
        }
        
        # xgboost softprob is typically interpreted row-wise:
        # (p_11, p_12, ..., p_1K, p_21, ..., p_nK)
        pmat <- matrix(pred, nrow = n, ncol = num_class, byrow = TRUE)
    }
    
    storage.mode(pmat) <- "double"
    
    # Optional sanity checks for probabilities
    if (check_prob) {
        if (any(!is.finite(pmat))) {
            stop(".coerce_xgb_multiclass_pred(): non-finite values found in probability matrix.")
        }
        # allow tiny numerical noise but flag clearly impossible values
        if (any(pmat < -prob_tol) || any(pmat > 1 + prob_tol)) {
            rng <- range(pmat, na.rm = TRUE)
            stop(".coerce_xgb_multiclass_pred(): probabilities outside [0,1] beyond tolerance. ",
                 "Range = [", signif(rng[1], 6), ", ", signif(rng[2], 6), "].")
        }
        
        rs <- rowSums(pmat)
        bad <- which(abs(rs - 1) > prob_tol)
        
        if (length(bad) > 0) {
            show_idx <- utils::head(bad, 5)
            stop(
                ".coerce_xgb_multiclass_pred(): row sums are not ~1 for ",
                length(bad), " row(s). Example row sums: ",
                paste(signif(rs[show_idx], 6), collapse = ", "),
                ". Consider checking xgboost objective/prediction format."
            )
        }
    }
    
    pmat
}


## ------------------------------------------------------------------
## Fit m_k(L_k) = E[Y_k(theta, A, alpha) | L_k, X >= t_k]
## using IPCW pseudo-outcome regression/classification
## ------------------------------------------------------------------

fit_Yk_condmean_at_tk <- function(
        dat_long, visit_times, k, tau,
        start.name, stop.name, event.name, event.C.name, id.name,
        Hk_df,
        Ctau_df, q_lower.name = "q_lower", q_upper.name = "q_upper",
        g_tau.name = "g_tau",
        alpha = 0.1,
        covname.baseline = NULL,
        covname.timevarying = NULL,
        use_history = FALSE,
        history_mode = c("wide", "summary"),
        model = c("lm", "rsf", "xgb_reg", "xgb_multiclass", "hal_reg"),
        trim.H = 1e-7,
        interval_open = TRUE,

        ## randomForestSRC tuning (used by rsf)
        rsf_args = list(),
        
        ## XGBoost tuning (used by xgb_reg / xgb_multiclass)
        xgb_nrounds = 200L,
        xgb_use_cv = FALSE,
        xgb_cv_nfold = 5L,
        xgb_cv_early_stopping_rounds = 20L,
        xgb_cv_param_grid = NULL,
        xgb_cv_seed = NULL,
        xgb_cv_verbose = 0,
        xgb_eta = 0.05,
        xgb_max_depth = 4L,
        xgb_min_child_weight = 1,
        xgb_subsample = 0.8,
        xgb_colsample_bytree = 0.8,
        xgb_lambda = 1,
        xgb_gamma = 0,
        xgb_nthread = NULL,
        xgb_seed = NULL,
        xgb_verbose = 0,
        xgb_params_extra = list(),
        
        ## GRF tuning (used by grf)
        grf_num.trees = 2000L,
        grf_mtry = NULL,
        grf_min.node.size = 5L,
        grf_sample.fraction = 0.5,
        grf_honesty = TRUE,
        grf_seed = NULL,
        grf_tune.parameters = "none",

        ## Direct fallback to weighted mean for sparse Xi_k fits
        Xi_fallback_to_mean = TRUE,
        Xi_min_n = 20L,
        Xi_min_ess = 20,
        Xi_min_obs_per_coef = 5,
        Xi_min_class_n = 5L,
        Xi_min_class_ess = 5
) {
    history_mode <- match.arg(history_mode)
    model <- match.arg(model)
    if (!is.list(rsf_args)) {
        stop("fit_Yk_condmean_at_tk(): rsf_args must be a list.")
    }
    Xi_min_n <- as.integer(Xi_min_n[[1L]])
    Xi_min_ess <- as.numeric(Xi_min_ess[[1L]])
    Xi_min_obs_per_coef <- as.numeric(Xi_min_obs_per_coef[[1L]])
    Xi_min_class_n <- as.integer(Xi_min_class_n[[1L]])
    Xi_min_class_ess <- as.numeric(Xi_min_class_ess[[1L]])
    if (!is.finite(Xi_min_n) || Xi_min_n < 0L) {
        stop("fit_Yk_condmean_at_tk(): Xi_min_n must be a nonnegative integer.")
    }
    if (!is.finite(Xi_min_ess) || Xi_min_ess < 0) {
        stop("fit_Yk_condmean_at_tk(): Xi_min_ess must be nonnegative.")
    }
    if (!is.finite(Xi_min_obs_per_coef) || Xi_min_obs_per_coef < 0) {
        stop("fit_Yk_condmean_at_tk(): Xi_min_obs_per_coef must be nonnegative.")
    }
    if (!is.finite(Xi_min_class_n) || Xi_min_class_n < 0L) {
        stop("fit_Yk_condmean_at_tk(): Xi_min_class_n must be a nonnegative integer.")
    }
    if (!is.finite(Xi_min_class_ess) || Xi_min_class_ess < 0) {
        stop("fit_Yk_condmean_at_tk(): Xi_min_class_ess must be nonnegative.")
    }
    
    K <- length(visit_times)
    if (K == 0) stop("fit_Yk_condmean_at_tk(): visit_times must be nonempty.")
    if (length(k) != 1 || is.na(k) || k < 1 || k > K) {
        stop("fit_Yk_condmean_at_tk(): k must be in 1:length(visit_times).")
    }
    k <- as.integer(k)
    t_k <- visit_times[k]
    
    ## -------- Subject summary (X, Delta, DeltaC) --------
    summ <- subject_summary(dat_long, id.name, stop.name, event.name, event.C.name)
    if (nrow(summ) == 0) stop("fit_Yk_condmean_at_tk(): no subjects in dat_long.")
    
    ## -------- Build L_k data --------
    if (use_history) {
        if (is.null(covname.timevarying) || length(covname.timevarying) == 0) {
            stop("fit_Yk_condmean_at_tk(): use_history=TRUE requires covname.timevarying.")
        }
        dat_k <- prepare_data_tau_with_history(
            dat = dat_long, visit_times = visit_times, k = k,
            start.name = start.name, stop.name = stop.name, event.name = event.name, id.name = id.name,
            covars = covname.timevarying, mode = history_mode
        )
    } else {
        dat_k <- prepare_data_tau(
            dat = dat_long,
            start.name = start.name, stop.name = stop.name, event.name = event.name, id.name = id.name,
            tau = t_k
        )
    }
    
    if (nrow(dat_k) == 0) {
        stop("fit_Yk_condmean_at_tk(): no rows returned by prepare_data_tau* at t_k.")
    }
    
    dat_k <- dat_k[, setdiff(names(dat_k), c("X", "Delta", "DeltaC")), drop = FALSE]
    dat_k <- merge(dat_k, summ, by = id.name, all.x = TRUE, sort = FALSE)
    
    ## Condition on X >= t_k; since time is continuous/interval data and your code uses >, keep >.
    dat_k <- dat_k[dat_k$X > t_k, , drop = FALSE]
    if (nrow(dat_k) == 0) stop("fit_Yk_condmean_at_tk(): no subjects at risk at t_k.")
    
    ## -------- Merge H_k(X) and interval endpoints --------
    if (!all(c("id", "Hk") %in% names(Hk_df))) {
        stop("fit_Yk_condmean_at_tk(): Hk_df must contain columns 'id' and 'Hk'.")
    }
    if (!all(c("id", q_lower.name, q_upper.name) %in% names(Ctau_df))) {
        stop("fit_Yk_condmean_at_tk(): Ctau_df must contain columns 'id', '", q_lower.name,
             "', and '", q_upper.name, "'.")
    }
    
    Hk_use <- Hk_df[, c("id", "Hk"), drop = FALSE]
    use_g_tau_col <- (!is.null(g_tau.name) && nzchar(g_tau.name) && (g_tau.name %in% names(Ctau_df)))
    C_cols <- c("id", q_lower.name, q_upper.name)
    if (use_g_tau_col) C_cols <- c(C_cols, g_tau.name)
    C_use <- Ctau_df[, C_cols, drop = FALSE]
    
    dat_k <- merge(dat_k, Hk_use, by.x = id.name, by.y = "id", all.x = TRUE, sort = FALSE)
    dat_k <- merge(dat_k, C_use,  by.x = id.name, by.y = "id", all.x = TRUE, sort = FALSE)
    
    if (anyNA(dat_k$Hk)) {
        miss_ids <- unique(dat_k[[id.name]][is.na(dat_k$Hk)])
        stop("fit_Yk_condmean_at_tk(): missing Hk for some ids. Example: ",
             paste(utils::head(miss_ids, 10), collapse = ", "))
    }
    
    # q-limits are only needed when X > tau (because I(X > tau) gates Y_tilde)
    need_q <- dat_k$X > tau
    bad_q <- need_q & (is.na(dat_k[[q_lower.name]]) | is.na(dat_k[[q_upper.name]]))
    if (any(bad_q)) {
        miss_ids <- unique(dat_k[[id.name]][bad_q])
        stop("fit_Yk_condmean_at_tk(): missing interval endpoints for some ids. Example: ",
             paste(utils::head(miss_ids, 10), collapse = ", "))
    }
    
    ## -------- Construct pseudo-outcome Y_tilde --------
    X  <- dat_k$X
    ql <- dat_k[[q_lower.name]]
    qu <- dat_k[[q_upper.name]]
    
    if (use_g_tau_col) {
        bad_g <- need_q & (!is.finite(dat_k[[g_tau.name]]) | dat_k[[g_tau.name]] < 0)
        if (any(bad_g)) {
            miss_ids <- unique(dat_k[[id.name]][bad_g])
            stop("fit_Yk_condmean_at_tk(): missing/invalid ", g_tau.name, " for some ids. Example: ",
                 paste(utils::head(miss_ids, 10), collapse = ", "))
        }
        g_tau_vec <- as.numeric(dat_k[[g_tau.name]])
    } else {
        g_tau_vec <- rep(1.0, nrow(dat_k))
    }
    
    in_interval <- integer(length(X))
    if (any(need_q)) {
        if (interval_open) {
            in_interval[need_q] <- as.integer(X[need_q] > ql[need_q] & X[need_q] < qu[need_q])
        } else {
            in_interval[need_q] <- as.integer(X[need_q] >= ql[need_q] & X[need_q] <= qu[need_q])
        }
    }
    
    # Y_tilde = \tilde G(tau|L_nu) * I(X > tau) * ( I{X in C_tau} - (1-alpha) )
    dat_k$Y_tilde <- numeric(length(X)) # Initiate to be all zeros
    dat_k$Y_tilde[need_q] <- g_tau_vec[need_q] * (in_interval[need_q] - (1 - alpha))
    
    # 3-class label (sign pattern of Y_tilde states):
    #   zero : X <= tau
    #   out  : X > tau and X not in interval
    #   in   : X > tau and X in interval
    dat_k$Y_class <- ifelse(!need_q, "zero", ifelse(in_interval == 1L, "in", "out"))
    dat_k$Y_class <- factor(dat_k$Y_class, levels = c("zero", "out", "in"))
    
    ## -------- IPCW weights and fitting sample --------
    dat_k$w_ipcw <- dat_k$Delta / pmax(dat_k$Hk, trim.H)
    
    fit_idx <- which(is.finite(dat_k$w_ipcw) & dat_k$w_ipcw > 0)
    if (length(fit_idx) == 0) {
        stop("fit_Yk_condmean_at_tk(): no uncensored weighted observations to fit regression.")
    }
    dat_fit <- dat_k[fit_idx, , drop = FALSE]

    # The 3-class construction of 'xgb_multiclass' identifies E[Y_tilde | L] only
    # because Y_tilde takes the three values {0, -(1-alpha), alpha}. When a
    # nonunit g_tau multiplies the pseudo-outcome (Gtau_mode = "estimated" or
    # "tilted"), Y_tilde is no longer 3-valued and a classification model is
    # conceptually inapplicable: computing pmat %*% value_map would silently
    # drop the g_tau factor and target a different estimand than the regression
    # branches. Fail loudly instead. Inert when g_tau is identically 1
    # (Gtau_mode = "one", the paper's setting).
    if (model == "xgb_multiclass" && use_g_tau_col) {
        g_fit <- g_tau_vec[fit_idx]
        g_active <- g_fit[dat_fit$Y_class != "zero" & is.finite(g_fit)]
        if (length(g_active) > 0 && any(abs(g_active - 1) > 1e-12)) {
            stop("fit_Yk_condmean_at_tk(): model='xgb_multiclass' requires the ",
                 "pseudo-outcome to take exactly the three values {0, -(1-alpha), alpha}, ",
                 "but a nonunit '", g_tau.name, "' factor is present (Gtau_mode != 'one'), ",
                 "so Y_tilde = g_tau * (...) is no longer 3-valued. ",
                 "Use a regression Xi model instead (e.g. model = 'xgb_reg', 'lm', 'rsf', or 'hal_reg').")
        }
    }

    ## -------- Covariates for L_k --------
    if (!use_history) {
        cov_use <- c(covname.baseline, covname.timevarying)
    } else {
        cov_use <- c(covname.baseline, expand_hist_covars(covname.timevarying, k, history_mode))
    }
    cov_use <- keep_nonconstant_covars(dat_fit, cov_use)
    n_coef_fit <- .covariate_design_ncoef(dat_fit, cov_use)
    fit_n_used <- nrow(dat_fit)
    fit_ess <- .weighted_ess(dat_fit$w_ipcw)
    Xi_min_n_required <- .event_threshold_by_dim(Xi_min_n, Xi_min_obs_per_coef, n_coef_fit)
    Xi_min_ess_required <- .event_threshold_by_dim(Xi_min_ess, Xi_min_obs_per_coef, n_coef_fit)
    class_counts <- table(dat_fit$Y_class)
    class_ess <- stats::setNames(
        vapply(levels(dat_k$Y_class), function(cl) {
            idx_cl <- dat_fit$Y_class == cl
            .weighted_ess(dat_fit$w_ipcw[idx_cl])
        }, numeric(1)),
        levels(dat_k$Y_class)
    )
    class_counts_full <- stats::setNames(rep(0L, length(levels(dat_k$Y_class))), levels(dat_k$Y_class))
    class_counts_full[names(class_counts)] <- as.integer(class_counts)
    
    make_mean_output <- function(fallback_reason, train_df = dat_fit, cov_use_out = cov_use) {
        mu0 <- stats::weighted.mean(train_df$Y_tilde, w = train_df$w_ipcw)
        fit_obj <- list(
            type = "mean",
            mu = unname(mu0),
            requested_type = model,
            fallback_reason = fallback_reason
        )
        out <- list(
            k = k, t_k = t_k, tau = tau, alpha = alpha,
            id.name = id.name, start.name = start.name, stop.name = stop.name,
            event.name = event.name, event.C.name = event.C.name,
            visit_times = visit_times,
            use_history = use_history, history_mode = history_mode,
            covname.baseline = covname.baseline,
            covname.timevarying = covname.timevarying,
            cov_use = cov_use_out,
            q_lower.name = q_lower.name, q_upper.name = q_upper.name,
            interval_open = interval_open,
            model = fit_obj,
            train_df = train_df,
            requested_model = model,
            fallback_reason = fallback_reason,
            n_used = nrow(train_df),
            fit_ess = .weighted_ess(train_df$w_ipcw),
            n_coef = n_coef_fit,
            Xi_min_n_required = Xi_min_n_required,
            Xi_min_ess_required = Xi_min_ess_required,
            Xi_min_class_n = Xi_min_class_n,
            Xi_min_class_ess = Xi_min_class_ess,
            class_counts = class_counts_full,
            class_ess = class_ess
        )
        class(out) <- "YkCondMeanFit"
        out
    }
    
    ## no-covariate fallback (constant weighted mean)
    if (length(cov_use) == 0) {
        return(make_mean_output("no_usable_covariates"))
    }

    fallback_reason <- NULL
    if (isTRUE(Xi_fallback_to_mean)) {
        if (fit_n_used < Xi_min_n_required) {
            fallback_reason <- paste0("few_rows_lt_", Xi_min_n_required)
        } else if (fit_ess < Xi_min_ess_required) {
            fallback_reason <- paste0("low_ess_lt_", Xi_min_ess_required)
        }
    }
    if (!is.null(fallback_reason)) {
        return(make_mean_output(fallback_reason))
    }
    
    ## -------- LM branch --------
    if (model == "lm") {
        # complete-case for raw covariates
        cc <- stats::complete.cases(dat_fit[, cov_use, drop = FALSE])
        dat_fit_cc <- dat_fit[cc, , drop = FALSE]
        if (nrow(dat_fit_cc) == 0) {
            stop("fit_Yk_condmean_at_tk(): no complete-case rows left for lm.")
        }
        if (isTRUE(Xi_fallback_to_mean)) {
            cc_ess <- .weighted_ess(dat_fit_cc$w_ipcw)
            if (nrow(dat_fit_cc) < Xi_min_n_required) {
                return(make_mean_output(paste0("complete_case_few_rows_lt_", Xi_min_n_required), dat_fit_cc))
            }
            if (cc_ess < Xi_min_ess_required) {
                return(make_mean_output(paste0("complete_case_low_ess_lt_", Xi_min_ess_required), dat_fit_cc))
            }
        }
        
        rhs <- paste(cov_use, collapse = " + ")
        fml <- stats::as.formula(paste0("Y_tilde ~ ", rhs))
        fit_lm <- stats::lm(fml, data = dat_fit_cc, weights = dat_fit_cc$w_ipcw)
        
        fit_obj <- list(type = "lm", fit = fit_lm)
        
        out <- list(
            k = k, t_k = t_k, tau = tau, alpha = alpha,
            id.name = id.name, start.name = start.name, stop.name = stop.name,
            event.name = event.name, event.C.name = event.C.name,
            visit_times = visit_times,
            use_history = use_history, history_mode = history_mode,
            covname.baseline = covname.baseline,
            covname.timevarying = covname.timevarying,
            cov_use = cov_use,
            q_lower.name = q_lower.name, q_upper.name = q_upper.name,
            interval_open = interval_open,
            model = fit_obj,
            train_df = dat_fit_cc,
            requested_model = model,
            fallback_reason = NA_character_,
            n_used = nrow(dat_fit_cc),
            fit_ess = .weighted_ess(dat_fit_cc$w_ipcw),
            n_coef = n_coef_fit,
            Xi_min_n_required = Xi_min_n_required,
            Xi_min_ess_required = Xi_min_ess_required,
            Xi_min_class_n = Xi_min_class_n,
            Xi_min_class_ess = Xi_min_class_ess,
            class_counts = class_counts_full,
            class_ess = class_ess
        )
        class(out) <- "YkCondMeanFit"
        return(out)
    }
    
    ## -------- Build numeric design matrix for tree/boosting methods --------
    dm_fit <- .fit_design_matrix(dat_fit, cov_use = cov_use)
    if (is.null(dm_fit$X) || nrow(dm_fit$X) == 0) {
        stop("fit_Yk_condmean_at_tk(): no complete-case rows left after design matrix construction.")
    }
    dat_fit_ml <- dat_fit[dm_fit$keep, , drop = FALSE]
    X_ml <- dm_fit$X
    w_ml <- as.numeric(dat_fit_ml$w_ipcw)
    if (isTRUE(Xi_fallback_to_mean)) {
        ml_ess <- .weighted_ess(dat_fit_ml$w_ipcw)
        if (nrow(dat_fit_ml) < Xi_min_n_required) {
            return(make_mean_output(paste0("complete_case_few_rows_lt_", Xi_min_n_required), dat_fit_ml))
        }
        if (ml_ess < Xi_min_ess_required) {
            return(make_mean_output(paste0("complete_case_low_ess_lt_", Xi_min_ess_required), dat_fit_ml))
        }
        if (identical(model, "xgb_multiclass")) {
            ml_class_counts <- table(dat_fit_ml$Y_class)
            ml_class_counts_full <- stats::setNames(rep(0L, length(levels(dat_k$Y_class))), levels(dat_k$Y_class))
            ml_class_counts_full[names(ml_class_counts)] <- as.integer(ml_class_counts)
            ml_class_ess <- stats::setNames(
                vapply(levels(dat_k$Y_class), function(cl) {
                    idx_cl <- dat_fit_ml$Y_class == cl
                    .weighted_ess(dat_fit_ml$w_ipcw[idx_cl])
                }, numeric(1)),
                levels(dat_k$Y_class)
            )
            nonzero_classes <- ml_class_counts_full[ml_class_counts_full > 0L]
            positive_class_ess <- ml_class_ess[ml_class_counts_full > 0L]
            if (length(nonzero_classes) < 2L) {
                return(make_mean_output("fewer_than_two_classes", dat_fit_ml))
            }
            if (min(nonzero_classes) < Xi_min_class_n) {
                return(make_mean_output(paste0("small_class_count_lt_", Xi_min_class_n), dat_fit_ml))
            }
            if (min(positive_class_ess) < Xi_min_class_ess) {
                return(make_mean_output(paste0("small_class_ess_lt_", Xi_min_class_ess), dat_fit_ml))
            }
        }
    }
    
    # Optional weight normalization for numerical stability (does not change target scale)
    if (all(is.finite(w_ml)) && sum(w_ml) > 0) {
        w_ml <- w_ml / mean(w_ml)
    }
    
    ## -------- Random forest regression (randomForestSRC) --------
    if (model == "rsf") {
        if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
            stop("Package 'randomForestSRC' is required for model='rsf'.")
        }

        y_ml <- as.numeric(dat_fit_ml$Y_tilde)
        rsf_df <- as.data.frame(X_ml)
        rsf_df$Y_tilde <- y_ml

        rsf_fit_args <- utils::modifyList(
            list(
                formula = stats::as.formula("Y_tilde ~ ."),
                data = rsf_df,
                case.wt = w_ml,
                ntree = 1000L,
                nodesize = 5L,
                forest = TRUE
            ),
            rsf_args
        )
        rsf_fit_args$formula <- stats::as.formula("Y_tilde ~ .")
        rsf_fit_args$data <- rsf_df
        rsf_fit_args$case.wt <- w_ml
        rsf_fit_args$forest <- TRUE

        fit_rsf <- do.call(randomForestSRC::rfsrc, rsf_fit_args)

        fit_obj <- list(
            type = "rsf",
            fit = fit_rsf,
            design_info = dm_fit$design_info,
            tuning = rsf_fit_args
        )

        ## -------- XGBoost regression --------
    } else if (model == "xgb_reg") {
        y_ml <- as.numeric(dat_fit_ml$Y_tilde)
        
        dtrain <- xgboost::xgb.DMatrix(data = X_ml, label = y_ml, weight = w_ml)
        
        params <- c(list(
            objective = "reg:squarederror",
            eta = xgb_eta,
            max_depth = as.integer(xgb_max_depth),
            min_child_weight = xgb_min_child_weight,
            subsample = xgb_subsample,
            colsample_bytree = xgb_colsample_bytree,
            lambda = xgb_lambda,
            gamma = xgb_gamma
        ), xgb_params_extra)
        
        if (!is.null(xgb_nthread)) params$nthread <- as.integer(xgb_nthread)
        
        xgb_fit <- .xgb_train_optional_cv(
            dtrain = dtrain,
            params = params,
            nrounds = as.integer(xgb_nrounds),
            use_cv = xgb_use_cv,
            cv_nfold = xgb_cv_nfold,
            cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
            cv_param_grid = xgb_cv_param_grid,
            cv_seed = xgb_cv_seed,
            cv_verbose = xgb_cv_verbose,
            maximize = FALSE,
            train_verbose = xgb_verbose
        )
        fit_xgb <- xgb_fit$fit
        
        fit_obj <- list(
            type = "xgb_reg",
            fit = fit_xgb,
            design_info = dm_fit$design_info,
            tuning = xgb_fit$tuning
            # xgb_version = as.character(utils::packageVersion("xgboost"))
        )
        
        ## -------- XGBoost multiclass probability -> mean --------
    } else if (model == "xgb_multiclass") {
        
        # If only one class among fitting rows, fallback to constant weighted mean
        if (length(unique(dat_fit_ml$Y_class)) < 2) {
            mu0 <- stats::weighted.mean(dat_fit_ml$Y_tilde, w = w_ml)
            fit_obj <- list(
                type = "mean",
                mu = unname(mu0),
                requested_type = model,
                fallback_reason = "fewer_than_two_classes"
            )
        } else {
            # Map factor levels c("zero","out","in") -> 0,1,2
            class_levels <- c("zero", "out", "in")
            y_fac <- factor(dat_fit_ml$Y_class, levels = class_levels)
            y_int0 <- as.integer(y_fac) - 1L
            
            dtrain <- xgboost::xgb.DMatrix(data = X_ml, label = y_int0, weight = w_ml)
            
            params <- c(list(
                objective = "multi:softprob",
                num_class = 3L,
                eta = xgb_eta,
                max_depth = as.integer(xgb_max_depth),
                min_child_weight = xgb_min_child_weight,
                subsample = xgb_subsample,
                colsample_bytree = xgb_colsample_bytree,
                lambda = xgb_lambda,
                gamma = xgb_gamma,
                eval_metric = "mlogloss"
            ), xgb_params_extra)
            
            if (!is.null(xgb_nthread)) params$nthread <- as.integer(xgb_nthread)
            
            xgb_fit <- .xgb_train_optional_cv(
                dtrain = dtrain,
                params = params,
                nrounds = as.integer(xgb_nrounds),
                use_cv = xgb_use_cv,
                cv_nfold = xgb_cv_nfold,
                cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
                cv_param_grid = xgb_cv_param_grid,
                cv_seed = xgb_cv_seed,
                cv_verbose = xgb_cv_verbose,
                maximize = FALSE,
                train_verbose = xgb_verbose
            )
            fit_xgb <- xgb_fit$fit
            
            # E[Y_tilde | L] = 0 * p_zero + (-(1-alpha))*p_out + alpha*p_in
            value_map <- c(zero = 0, out = -(1 - alpha), `in` = alpha)
            
            fit_obj <- list(
                type = "xgb_multiclass",
                fit = fit_xgb,
                design_info = dm_fit$design_info,
                class_levels = class_levels,
                value_map = value_map,
                tuning = xgb_fit$tuning
                # xgb_version = as.character(utils::packageVersion("xgboost"))
            )
        }
        
        ## -------- GRF regression forest --------
    } else if (model == "grf") {
        y_ml <- as.numeric(dat_fit_ml$Y_tilde)
        
        grf_args <- list(
            X = X_ml,
            Y = y_ml,
            sample.weights = w_ml,
            num.trees = as.integer(grf_num.trees),
            min.node.size = as.integer(grf_min.node.size),
            sample.fraction = grf_sample.fraction,
            honesty = isTRUE(grf_honesty),
            tune.parameters = grf_tune.parameters
        )
        if (!is.null(grf_mtry)) grf_args$mtry <- as.integer(grf_mtry)
        if (!is.null(grf_seed)) grf_args$seed <- as.integer(grf_seed)
        
        fit_grf <- do.call(grf::regression_forest, grf_args)
        
        fit_obj <- list(
            type = "grf",
            fit = fit_grf,
            design_info = dm_fit$design_info
        )
        
        ## -------- HAL gaussian regression --------
    } else if (model == "hal_reg") {
        if (!requireNamespace("hal9001", quietly = TRUE)) {
            stop("Package 'hal9001' is required for model='hal_reg'.")
        }
        
        y_ml <- as.numeric(dat_fit_ml$Y_tilde)
        fit_hal <- hal9001::fit_hal(
            X = X_ml,
            Y = y_ml,
            family = "gaussian",
            weights = w_ml,
            yolo = FALSE
        )
        
        fit_obj <- list(
            type = "hal_reg",
            fit = fit_hal,
            design_info = dm_fit$design_info
        )
        
    } else {
        stop("Unsupported model in fit_Yk_condmean_at_tk().")
    }
    
    ## -------- Pack result --------
    out <- list(
        k = k,
        t_k = t_k,
        tau = tau,
        alpha = alpha,
        id.name = id.name,
        start.name = start.name,
        stop.name = stop.name,
        event.name = event.name,
        event.C.name = event.C.name,
        visit_times = visit_times,
        use_history = use_history,
        history_mode = history_mode,
        covname.baseline = covname.baseline,
        covname.timevarying = covname.timevarying,
        cov_use = cov_use,
        q_lower.name = q_lower.name,
        q_upper.name = q_upper.name,
        interval_open = interval_open,
        model = fit_obj,
        train_df = dat_fit_ml,
        requested_model = model,
        fallback_reason = if (!is.null(fit_obj$fallback_reason)) fit_obj$fallback_reason else NA_character_,
        n_used = nrow(dat_fit_ml),
        fit_ess = .weighted_ess(dat_fit_ml$w_ipcw),
        n_coef = n_coef_fit,
        Xi_min_n_required = Xi_min_n_required,
        Xi_min_ess_required = Xi_min_ess_required,
        Xi_min_class_n = Xi_min_class_n,
        Xi_min_class_ess = Xi_min_class_ess,
        class_counts = class_counts_full,
        class_ess = class_ess
    )
    class(out) <- "YkCondMeanFit"
    out
}


## ------------------------------------------------------------------
## Predict m_k(L_k) on new long-format data (returns at-risk subjects only)
## ------------------------------------------------------------------

predict_Yk_condmean_at_tk <- function(
        object,
        dat_long_new,
        clip_to_support = TRUE
) {
    if (!inherits(object, "YkCondMeanFit")) {
        stop("predict_Yk_condmean_at_tk(): object must be of class 'YkCondMeanFit'.")
    }
    
    k <- object$k
    t_k <- object$t_k
    
    ## -------- Build L_k rows --------
    if (object$use_history) {
        dat_k_new <- prepare_data_tau_with_history(
            dat = dat_long_new, visit_times = object$visit_times, k = k,
            start.name = object$start.name, stop.name = object$stop.name,
            event.name = object$event.name, id.name = object$id.name,
            covars = object$covname.timevarying, mode = object$history_mode
        )
    } else {
        dat_k_new <- prepare_data_tau(
            dat = dat_long_new,
            start.name = object$start.name, stop.name = object$stop.name,
            event.name = object$event.name, id.name = object$id.name,
            tau = t_k
        )
    }
    
    if (nrow(dat_k_new) == 0) return(data.frame())
    
    ## -------- Restrict to at-risk at t_k --------
    summ_new <- subject_summary(
        dat_long_new, object$id.name, object$stop.name, object$event.name, object$event.C.name
    )
    dat_k_new <- dat_k_new[, setdiff(names(dat_k_new), c("X", "Delta", "DeltaC")), drop = FALSE]
    dat_k_new <- merge(dat_k_new, summ_new, by = object$id.name, all.x = TRUE, sort = FALSE)
    dat_k_new <- dat_k_new[dat_k_new$X > t_k, , drop = FALSE]
    if (nrow(dat_k_new) == 0) return(data.frame())
    
    pred <- rep(NA_real_, nrow(dat_k_new))
    
    ## -------- Predict by model type --------
    if (object$model$type == "mean") {
        pred[] <- object$model$mu
        
    } else if (object$model$type == "lm") {
        pred <- as.numeric(stats::predict(object$model$fit, newdata = dat_k_new))
        
    } else if (object$model$type %in% c("rsf", "xgb_reg", "xgb_multiclass", "grf", "hal_reg")) {
        # Build aligned numeric design matrix
        dm_new <- .predict_design_matrix(
            dat = dat_k_new,
            cov_use = object$cov_use,
            design_info = object$model$design_info
        )
        
        if (any(dm_new$keep)) {
            X_new <- dm_new$X
            idx_keep <- which(dm_new$keep)
            
            if (object$model$type == "rsf") {
                pr <- predict(object$model$fit, newdata = as.data.frame(X_new))
                pred[idx_keep] <- as.numeric(pr$predicted)
            } else if (object$model$type == "xgb_reg") {
                dnew <- xgboost::xgb.DMatrix(data = X_new)
                pr <- predict(object$model$fit, newdata = dnew)
                pred[idx_keep] <- as.numeric(pr)
                
            } else if (object$model$type == "xgb_multiclass") {
                dnew <- xgboost::xgb.DMatrix(data = X_new)
                
                # Prefer strict_shape if supported (more explicit output shape)
                pr <- tryCatch(
                    predict(object$model$fit, newdata = dnew, strict_shape = TRUE),
                    error = function(e) predict(object$model$fit, newdata = dnew)
                )
                
                pmat <- .coerce_xgb_multiclass_pred(
                    pred = pr,
                    n = nrow(X_new),
                    num_class = length(object$model$class_levels),
                    check_prob = TRUE
                )
                
                # xgboost class order should match labels 0,1,2 used in training
                vals <- object$model$value_map[object$model$class_levels]
                pred[idx_keep] <- as.numeric(pmat %*% vals)
                
            } else if (object$model$type == "grf") {
                pr <- predict(object$model$fit, newdata = X_new)
                # grf predict returns a list with $predictions
                pred[idx_keep] <- as.numeric(pr$predictions)
            } else if (object$model$type == "hal_reg") {
                pr <- predict(object$model$fit, new_data = X_new, type = "response")
                pred[idx_keep] <- as.numeric(pr)
            }
        }
        
    } else {
        stop("predict_Yk_condmean_at_tk(): unsupported model type in object.")
    }
    
    ## -------- Optional clipping to support of Y_tilde --------
    if (clip_to_support) {
        lo <- -(1 - object$alpha)
        hi <- object$alpha
        pred <- pmin(hi, pmax(lo, pred))
    }
    
    data.frame(
        id = dat_k_new[[object$id.name]],
        mYk_hat = pred
    )
}



### Helpers and the function for h_{\tau} -----------------------------------------------------------------------------------

## ============================================================
## External helpers for h_tau matrix construction
## ============================================================

# Row-wise survival quantile from predicted survival curves at pred_time
# beta in [0,1], pred_surv rows aligned with subjects, columns aligned with pred_time
# Returns q_beta = inf{t : 1 - S(t) >= beta}, approximated on pred_time grid
surv_quantile_vec_from_pred <- function(pred_time, pred_surv, beta) {
    if (!is.matrix(pred_surv)) stop("surv_quantile_vec_from_pred(): pred_surv must be a matrix.")
    if (length(pred_time) != ncol(pred_surv)) {
        stop("surv_quantile_vec_from_pred(): length(pred_time) must equal ncol(pred_surv).")
    }
    if (length(beta) != 1L || !is.finite(beta) || beta < 0 || beta > 1) {
        stop("surv_quantile_vec_from_pred(): beta must be a scalar in [0,1].")
    }
    
    n <- nrow(pred_surv)
    out <- rep(NA_real_, n)
    tmax <- max(pred_time)
    
    # final attainable CDF on grid
    cdf_end <- 1 - pred_surv[, ncol(pred_surv)]
    
    for (i in seq_len(n)) {
        if (!is.finite(cdf_end[i])) next
        
        # if threshold not reached on grid, use max(pred_time)
        if (beta > cdf_end[i]) {
            out[i] <- tmax
        } else {
            idx <- which((1 - pred_surv[i, ]) >= beta)
            out[i] <- if (length(idx) == 0L) tmax else pred_time[min(idx)]
        }
    }
    out
}


# Evaluate a single right-continuous step survival curve S(t) at vector t
# grid increasing; S_grid[j] = S(grid[j])
# Returns 1 for t < grid[1]
eval_step_surv_vec_right <- function(S_grid, grid, t) {
    stopifnot(length(S_grid) == length(grid))
    j <- findInterval(t, grid, left.open = FALSE, rightmost.closed = FALSE)
    out <- rep(1.0, length(t))
    ok <- (j > 0L)
    out[ok] <- S_grid[j[ok]]
    out
}


# Evaluate a single step survival curve left-limit S(t-) at vector t
# If t equals a jump time exactly, returns the value just before the jump
eval_step_surv_vec_left <- function(S_grid, grid, t, tol = 1e-12) {
    stopifnot(length(S_grid) == length(grid))
    j_le <- findInterval(t, grid, left.open = FALSE, rightmost.closed = FALSE)
    
    out <- rep(1.0, length(t))
    ok <- (j_le > 0L)
    if (!any(ok)) return(out)
    
    j_use <- j_le
    idx_ok <- which(ok)
    eq_jump <- rep(FALSE, length(t))
    eq_jump[idx_ok] <- abs(t[idx_ok] - grid[j_le[idx_ok]]) <= tol
    j_use[eq_jump] <- j_use[eq_jump] - 1L
    
    ok2 <- (j_use > 0L)
    out[ok2] <- S_grid[j_use[ok2]]
    out
}


# Row-wise right-continuous evaluation for matrix survival curves
# S: n x m, grid length m, t_rel length n (one query time per row)
# Returns vector length n with S_i(t_rel[i])
eval_step_surv_row_right <- function(S, grid, t_rel) {
    stopifnot(is.matrix(S), length(grid) == ncol(S), length(t_rel) == nrow(S))
    
    j <- findInterval(t_rel, grid, left.open = FALSE, rightmost.closed = FALSE)
    out <- rep(1.0, nrow(S))
    ok <- (j > 0L)
    out[ok] <- S[cbind(which(ok), j[ok])]
    out
}


# Row-wise left-limit evaluation for matrix survival curves
# S: n x m, grid length m, t_rel length n (one query time per row)
# Returns vector length n with S_i(t_rel[i]-)
eval_step_surv_row_left <- function(S, grid, t_rel, tol = 1e-12) {
    stopifnot(is.matrix(S), length(grid) == ncol(S), length(t_rel) == nrow(S))
    
    j_le <- findInterval(t_rel, grid, left.open = FALSE, rightmost.closed = FALSE)
    out <- rep(1.0, nrow(S))
    
    ok <- (j_le > 0L)
    if (!any(ok)) return(out)
    
    j_use <- j_le
    idx_ok <- which(ok)
    eq_jump <- rep(FALSE, nrow(S))
    eq_jump[idx_ok] <- abs(t_rel[idx_ok] - grid[j_le[idx_ok]]) <= tol
    j_use[eq_jump] <- j_use[eq_jump] - 1L
    
    ok2 <- (j_use > 0L)
    out[ok2] <- S[cbind(which(ok2), j_use[ok2])]
    out
}


# Build landmark rows at visit_times[k] aligned to a target id order (subset with X > t_k)
# Uses Sfit_k metadata to decide whether to build history features
build_Sk_landmark_rows_aligned <- function(
        dat_long, ids_target, X_target, visit_times, k,
        start.name, stop.name, event.name, id.name,
        Sfit_k
) {
    t_k <- visit_times[k]
    idxk <- which(X_target > t_k)
    if (length(idxk) == 0L) return(list(idx = integer(0), dat = NULL))
    
    ids_k <- ids_target[idxk]
    dat_slice <- dat_long[dat_long[[id.name]] %in% ids_k, , drop = FALSE]
    
    use_hist_k  <- isTRUE(Sfit_k$use_history)
    cov_tv_k    <- Sfit_k$cov_tv
    hist_mode_k <- Sfit_k$history_mode
    
    if (use_hist_k && !is.null(cov_tv_k) && length(cov_tv_k) > 0) {
        dat_tk <- prepare_data_tau_with_history(
            dat = dat_slice, visit_times = visit_times, k = k,
            start.name = start.name, stop.name = stop.name, event.name = event.name, id.name = id.name,
            covars = cov_tv_k, mode = hist_mode_k
        )
    } else {
        dat_tk <- prepare_data_tau(
            dat = dat_slice,
            start.name = start.name, stop.name = stop.name, event.name = event.name, id.name = id.name,
            tau = t_k
        )
    }
    
    if (anyDuplicated(dat_tk[[id.name]])) {
        stop("build_Sk_landmark_rows_aligned(): duplicated ids at k=", k)
    }
    
    dat_tk <- dat_tk[match(ids_k, dat_tk[[id.name]]), , drop = FALSE]
    if (anyNA(dat_tk[[id.name]])) {
        stop("build_Sk_landmark_rows_aligned(): failed to align landmark rows at k=", k)
    }
    
    list(idx = idxk, dat = dat_tk)
}


# Precompute/caches needed to evaluate S_k(t | L_k) repeatedly for the same newdata_tk
# Returns an object used by predict_Sk_cache_* helpers below
make_Sk_predictor_cache <- function(Sfit_k, newdata_tk) {
    if (is.null(Sfit_k)) stop("make_Sk_predictor_cache(): Sfit_k is NULL.")
    n <- nrow(newdata_tk)
    
    # Covariate presence check
    if (!is.null(Sfit_k$cov_use)) {
        miss <- setdiff(Sfit_k$cov_use, names(newdata_tk))
        if (length(miss) > 0) {
            stop("make_Sk_predictor_cache(): newdata_tk missing covariates: ",
                 paste(miss, collapse = ", "))
        }
    }
    
    out <- list(
        model = Sfit_k$model,
        t_k   = Sfit_k$t_k,
        n     = n
    )

    if (Sfit_k$model == "hybrid_interval_ipcw") {
        short_cache <- make_Sk_predictor_cache(Sfit_k$short_fit, newdata_tk)
        tail_cache <- if (is.null(Sfit_k$tail_fit)) NULL else make_Sk_predictor_cache(Sfit_k$tail_fit, newdata_tk)
        bridge_time <- Sfit_k$bridge_time
        S_bridge_right <- if (is.finite(bridge_time)) {
            as.numeric(
                predict_Sk_cache_common_times(
                    short_cache,
                    t_abs_vec = bridge_time,
                    left_limit = FALSE
                )[, 1L]
            )
        } else {
            rep(NA_real_, n)
        }

        out$bridge_time <- bridge_time
        out$short_cache <- short_cache
        out$tail_cache <- tail_cache
        out$S_bridge_right <- S_bridge_right
        return(out)
    }
    
    if (Sfit_k$model == "km") {
        if (is.null(Sfit_k$grid) || is.null(Sfit_k$S_grid)) {
            ss <- summary(Sfit_k$fit)
            out$grid <- ss$time
            out$S_grid <- ss$surv
        } else {
            out$grid <- Sfit_k$grid
            out$S_grid <- Sfit_k$S_grid
        }
        
    } else if (Sfit_k$model == "cox") {
        fit <- Sfit_k$fit
        
        # baseline survival grid
        bh <- survival::basehaz(fit, centered = FALSE)
        out$grid <- bh$time
        out$S0_grid <- exp(-bh$hazard)
        
        beta <- stats::coef(fit)
        if (length(beta) == 0L) {
            out$risk_mult <- rep(1, n)
        } else {
            Xmat <- stats::model.matrix(fit, newdata_tk)
            # model.matrix() silently drops NA-covariate rows; a short risk_mult
            # would be recycled by sweep()/^ downstream, attributing survival
            # probabilities to the wrong subjects. Fail loudly (as xgb/grf do).
            if (nrow(Xmat) != n) {
                stop("make_Sk_predictor_cache(): model.matrix dropped ",
                     n - nrow(Xmat), " of ", n,
                     " rows (missing covariate values are not supported for cox prediction).")
            }
            eta  <- drop(Xmat %*% beta)
            out$risk_mult <- exp(eta)
        }

    } else if (Sfit_k$model == "rsf") {
        if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
            stop("make_Sk_predictor_cache(): package 'randomForestSRC' required for RSF.")
        }

        pr <- predict(Sfit_k$fit, newdata = newdata_tk)
        out$grid <- pr$time.interest
        out$Smat <- pr$survival  # n x length(grid)
        # predict.rfsrc() na.action default omits NA-covariate rows; a short Smat
        # would be silently recycled by sweep() downstream. Fail loudly instead.
        if (is.null(out$Smat) || NROW(out$Smat) != n) {
            stop("make_Sk_predictor_cache(): rsf prediction returned ",
                 NROW(out$Smat), " rows for ", n,
                 " subjects (missing covariate values are not supported for rsf prediction).")
        }
        
    } else if (Sfit_k$model == "hal") {
        hal_bundle <- list(
            fit_hal = Sfit_k$fit,
            design_info = Sfit_k$design_info,
            cov_use = Sfit_k$cov_use,
            lp_center = Sfit_k$lp_center,
            bh_time = Sfit_k$bh_time,
            bh_cumhaz = Sfit_k$bh_cumhaz
        )
        
        out$grid <- Sfit_k$bh_time
        out$S0_grid <- exp(-Sfit_k$bh_cumhaz)
        lp <- .predict_hal_cox_lp(hal_bundle, newdata_tk)
        lp_center <- if (!is.null(Sfit_k$lp_center) && is.finite(Sfit_k$lp_center)) Sfit_k$lp_center else 0
        out$risk_mult <- exp(lp - lp_center)
        
    } else if (Sfit_k$model == "xgb_cox") {
        if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("make_Sk_predictor_cache(): package 'xgboost' required for xgb_cox.")
        }
        
        dm_new <- .predict_design_matrix(
            dat = newdata_tk,
            cov_use = Sfit_k$cov_use,
            design_info = Sfit_k$design_info
        )
        if (any(!dm_new$keep)) {
            stop("make_Sk_predictor_cache(): missing values in covariates are not supported for xgb_cox prediction.")
        }
        
        lp <- as.numeric(stats::predict(Sfit_k$fit, newdata = dm_new$X, outputmargin = TRUE))
        lp[!is.finite(lp)] <- 0
        lp <- pmax(pmin(lp, 30), -30)
        lp_center <- if (!is.null(Sfit_k$lp_center) && is.finite(Sfit_k$lp_center)) Sfit_k$lp_center else 0
        out$grid <- Sfit_k$bh_time
        out$S0_grid <- exp(-Sfit_k$bh_cumhaz)
        out$risk_mult <- exp(lp - lp_center)
        
    } else if (Sfit_k$model == "xgb_aft") {
        if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("make_Sk_predictor_cache(): package 'xgboost' required for xgb_aft.")
        }
        
        dm_new <- .predict_design_matrix(
            dat = newdata_tk,
            cov_use = Sfit_k$cov_use,
            design_info = Sfit_k$design_info
        )
        if (any(!dm_new$keep)) {
            stop("make_Sk_predictor_cache(): missing values in covariates are not supported for xgb_aft prediction.")
        }
        
        out$aft_mu <- as.numeric(stats::predict(Sfit_k$fit, newdata = dm_new$X))
        out$aft_dist <- Sfit_k$aft_dist
        out$aft_scale <- as.numeric(Sfit_k$aft_scale)
        
    } else if (Sfit_k$model == "grf") {
        if (!requireNamespace("grf", quietly = TRUE)) {
            stop("make_Sk_predictor_cache(): package 'grf' required for GRF survival forest.")
        }
        
        dm_new <- .predict_design_matrix(
            dat = newdata_tk,
            cov_use = Sfit_k$cov_use,
            design_info = Sfit_k$design_info
        )
        if (any(!dm_new$keep)) {
            stop("make_Sk_predictor_cache(): missing values in covariates are not supported for GRF prediction.")
        }
        
        pr <- predict(Sfit_k$fit, newdata = dm_new$X)
        out$grid <- pr$failure.times
        out$Smat <- pr$predictions
        
    } else {
        stop("make_Sk_predictor_cache(): unsupported Sfit_k$model = ", Sfit_k$model)
    }
    
    out
}


# Evaluate cached S_k at common absolute times for all rows
# Returns n x length(t_abs_vec) matrix
predict_Sk_cache_common_times <- function(cache, t_abs_vec, left_limit = FALSE, tol = 1e-12) {
    n <- cache$n
    m <- length(t_abs_vec)
    if (n == 0L || m == 0L) return(matrix(numeric(0), nrow = n, ncol = m))

    if (cache$model == "hybrid_interval_ipcw") {
        bridge_time <- cache$bridge_time
        if (!is.finite(bridge_time)) {
            return(predict_Sk_cache_common_times(cache$short_cache, t_abs_vec, left_limit = left_limit, tol = tol))
        }

        out <- matrix(NA_real_, nrow = n, ncol = m)
        use_short <- t_abs_vec <= (bridge_time + tol)
        if (any(use_short)) {
            out[, use_short] <- predict_Sk_cache_common_times(
                cache$short_cache,
                t_abs_vec[use_short],
                left_limit = left_limit,
                tol = tol
            )
        }
        if (any(!use_short)) {
            if (is.null(cache$tail_cache)) {
                out[, !use_short] <- matrix(cache$S_bridge_right, nrow = n, ncol = sum(!use_short))
            } else {
                tail_part <- predict_Sk_cache_common_times(
                    cache$tail_cache,
                    t_abs_vec[!use_short],
                    left_limit = left_limit,
                    tol = tol
                )
                out[, !use_short] <- sweep(tail_part, 1, cache$S_bridge_right, FUN = "*")
            }
        }
        colnames(out) <- as.character(t_abs_vec)
        return(out)
    }
    
    t_rel <- pmax(0, t_abs_vec - cache$t_k)
    
    if (cache$model == "km") {
        if (left_limit) {
            S_eval <- eval_step_surv_vec_left(cache$S_grid, cache$grid, t_rel, tol = tol)
        } else {
            S_eval <- eval_step_surv_vec_right(cache$S_grid, cache$grid, t_rel)
        }
        out <- matrix(rep(S_eval, each = n), nrow = n, ncol = m)
        colnames(out) <- as.character(t_abs_vec)
        return(out)

    } else if (cache$model %in% c("cox", "hal", "xgb_cox")) {
        if (left_limit) {
            S0_eval <- eval_step_surv_vec_left(cache$S0_grid, cache$grid, t_rel, tol = tol)
        } else {
            S0_eval <- eval_step_surv_vec_right(cache$S0_grid, cache$grid, t_rel)
        }
        
        out <- outer(cache$risk_mult, S0_eval, FUN = function(r, s) s^r)
        colnames(out) <- as.character(t_abs_vec)
        return(out)
        
    } else if (cache$model %in% c("rsf", "grf")) {
        out <- matrix(NA_real_, nrow = n, ncol = m)
        for (j in seq_len(m)) {
            tj <- rep(t_rel[j], n)
            if (left_limit) {
                out[, j] <- eval_step_surv_row_left(cache$Smat, cache$grid, tj, tol = tol)
            } else {
                out[, j] <- eval_step_surv_row_right(cache$Smat, cache$grid, tj)
            }
        }
        colnames(out) <- as.character(t_abs_vec)
        return(out)
        
    } else if (cache$model == "xgb_aft") {
        tiny <- .Machine$double.eps
        t_pos <- pmax(t_rel, tiny)
        log_t <- log(t_pos)
        z <- outer(cache$aft_mu, log_t, FUN = function(mu_i, lt_j) (lt_j - mu_i) / cache$aft_scale)
        
        Fz <- switch(
            cache$aft_dist,
            normal = stats::pnorm(z),
            logistic = stats::plogis(z),
            extreme = exp(-exp(-z)),
            stop("predict_Sk_cache_common_times(): unsupported aft_dist = ", cache$aft_dist)
        )
        out <- 1 - Fz
        if (any(t_rel <= 0)) out[, t_rel <= 0] <- 1
        out[out < 0] <- 0
        out[out > 1] <- 1
        colnames(out) <- as.character(t_abs_vec)
        return(out)
        
    } else {
        stop("predict_Sk_cache_common_times(): unsupported cache$model")
    }
}


# Evaluate cached S_k at row-specific absolute times arranged in a matrix
# t_abs_mat: n x c matrix (entry [i,j] is time for subject i, column j)
# Returns n x c matrix with S_k(t_abs_mat[i,j] | L_k_i) (or left-limit)
predict_Sk_cache_rowtime_matrix <- function(cache, t_abs_mat, left_limit = FALSE, tol = 1e-12) {
    if (!is.matrix(t_abs_mat)) stop("predict_Sk_cache_rowtime_matrix(): t_abs_mat must be a matrix.")
    n <- cache$n
    if (nrow(t_abs_mat) != n) {
        stop("predict_Sk_cache_rowtime_matrix(): nrow(t_abs_mat) must equal cache$n.")
    }
    
    c <- ncol(t_abs_mat)
    if (n == 0L || c == 0L) return(matrix(numeric(0), nrow = n, ncol = c))

    if (cache$model == "hybrid_interval_ipcw") {
        bridge_time <- cache$bridge_time
        if (!is.finite(bridge_time)) {
            return(predict_Sk_cache_rowtime_matrix(cache$short_cache, t_abs_mat, left_limit = left_limit, tol = tol))
        }

        out <- matrix(NA_real_, nrow = n, ncol = c)
        use_short <- t_abs_mat <= (bridge_time + tol)
        if (any(use_short)) {
            short_all <- predict_Sk_cache_rowtime_matrix(
                cache$short_cache,
                t_abs_mat,
                left_limit = left_limit,
                tol = tol
            )
            out[use_short] <- short_all[use_short]
        }
        if (any(!use_short)) {
            if (is.null(cache$tail_cache)) {
                bridge_mat <- matrix(cache$S_bridge_right, nrow = n, ncol = c)
                out[!use_short] <- bridge_mat[!use_short]
            } else {
                tail_all <- predict_Sk_cache_rowtime_matrix(
                    cache$tail_cache,
                    t_abs_mat,
                    left_limit = left_limit,
                    tol = tol
                )
                tail_all <- sweep(tail_all, 1, cache$S_bridge_right, FUN = "*")
                out[!use_short] <- tail_all[!use_short]
            }
        }
        return(out)
    }
    
    # Keep matrix shape even when ncol(t_abs_mat) == 1.
    t_rel_mat <- t_abs_mat - cache$t_k
    t_rel_mat[t_rel_mat < 0] <- 0
    
    if (cache$model == "km") {
        t_vec <- as.vector(t_rel_mat)
        if (left_limit) {
            S_vec <- eval_step_surv_vec_left(cache$S_grid, cache$grid, t_vec, tol = tol)
        } else {
            S_vec <- eval_step_surv_vec_right(cache$S_grid, cache$grid, t_vec)
        }
        out <- matrix(S_vec, nrow = n, ncol = c)
        return(out)

    } else if (cache$model %in% c("cox", "hal", "xgb_cox")) {
        t_vec <- as.vector(t_rel_mat)
        if (left_limit) {
            S0_vec <- eval_step_surv_vec_left(cache$S0_grid, cache$grid, t_vec, tol = tol)
        } else {
            S0_vec <- eval_step_surv_vec_right(cache$S0_grid, cache$grid, t_vec)
        }
        S0_mat <- matrix(S0_vec, nrow = n, ncol = c)
        
        # row-wise exponentiation by risk multiplier
        out <- sweep(S0_mat, 1, cache$risk_mult, FUN = "^")
        return(out)
        
    } else if (cache$model %in% c("rsf", "grf")) {
        out <- matrix(NA_real_, nrow = n, ncol = c)
        for (j in seq_len(c)) {
            tj <- t_rel_mat[, j]
            if (left_limit) {
                out[, j] <- eval_step_surv_row_left(cache$Smat, cache$grid, tj, tol = tol)
            } else {
                out[, j] <- eval_step_surv_row_right(cache$Smat, cache$grid, tj)
            }
        }
        return(out)
        
    } else if (cache$model == "xgb_aft") {
        tiny <- .Machine$double.eps
        t_pos <- pmax(t_rel_mat, tiny)
        z <- (log(t_pos) - matrix(cache$aft_mu, nrow = n, ncol = c)) / cache$aft_scale
        
        Fz <- switch(
            cache$aft_dist,
            normal = stats::pnorm(z),
            logistic = stats::plogis(z),
            extreme = exp(-exp(-z)),
            stop("predict_Sk_cache_rowtime_matrix(): unsupported aft_dist = ", cache$aft_dist)
        )
        out <- 1 - Fz
        out[t_rel_mat <= 0] <- 1
        out[out < 0] <- 0
        out[out > 1] <- 1
        return(out)
        
    } else {
        stop("predict_Sk_cache_rowtime_matrix(): unsupported cache$model")
    }
}


## ============================================================
## Main function: h_tau matrix on calibration data
## ============================================================
# Rows correspond to dat_cal_tau rows (must align with pred_surv_cal rows)
# Columns correspond to time_vec (default = all jump times of G, absolute time)
#
# Xi_fits:
#   list of length >= length(visit_times), where Xi_fits[[k]] is a YkCondMeanFit object
#   used only for k < nu (u < tau branch)
#
# Gtau_vec:
#   optional vector length nrow(dat_cal_tau) or nrow(output rows):
#   \tilde G(tau | L_nu). If NULL, computed internally from Gfits + dat_cal.
build_h_tau_matrix_cal <- function(
        dat_cal, dat_cal_tau,
        pred_time, pred_surv_cal,
        visit_times, tau, theta = NULL, alpha,
        start.name, stop.name, event.name, event.C.name, id.name,
        Gfits, Sfits,
        Xi_fits = NULL,
        time_vec = NULL,
        Gtau_vec = NULL,
        trim.S = 0.05,
        trim.Gtau = 0,
        row_data = c("all", "tau"),
        mask_after_X = TRUE,
        after_X_value = NA_real_,
        tol = 1e-12,
        precomp = NULL,
        precomp_only = FALSE,
        fast_rsf_predict = FALSE,
        support_upper = Inf   # known event-support bound: at theta == 0 the
                              # candidate set is (tau, support_upper] (see
                              # .candidate_bounds_at_theta in dynamicCP_AIPCW.R)
) {
    # precomp: an optional list of theta-INDEPENDENT objects returned by a prior
    #   call with precomp_only = TRUE. Passing it back avoids recomputing the
    #   subject summary, G(tau|.), the interval map, the S_k predictor caches,
    #   and the S_k(u-|L_k) denominators on every theta of the AIPCW grid search
    #   (these do not depend on theta). Results are bit-identical either way.
    # precomp_only: build and return that cache without needing theta.
    row_data <- match.arg(row_data)
    ## ---------------- checks ----------------
    if (!is.data.frame(dat_cal) || !is.data.frame(dat_cal_tau)) {
        stop("build_h_tau_matrix_cal(): dat_cal and dat_cal_tau must be data.frames.")
    }
    if (!is.matrix(pred_surv_cal)) {
        stop("build_h_tau_matrix_cal(): pred_surv_cal must be a matrix.")
    }
    if (nrow(pred_surv_cal) != nrow(dat_cal_tau)) {
        stop("build_h_tau_matrix_cal(): nrow(pred_surv_cal) must equal nrow(dat_cal_tau).")
    }
    if (length(pred_time) != ncol(pred_surv_cal)) {
        stop("build_h_tau_matrix_cal(): length(pred_time) must equal ncol(pred_surv_cal).")
    }
    if (!isTRUE(precomp_only)) {
        if (length(theta) != 1L || !is.finite(theta) || theta < 0 || theta > 0.5) {
            stop("build_h_tau_matrix_cal(): theta must be a scalar in [0, 0.5].")
        }
    }
    if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
        stop("build_h_tau_matrix_cal(): alpha must be a scalar in (0, 1).")
    }

    K <- length(visit_times)
    if (K == 0L) stop("build_h_tau_matrix_cal(): visit_times must be nonempty.")
    nu <- match(tau, visit_times)
    if (is.na(nu)) stop("build_h_tau_matrix_cal(): tau must be one of visit_times.")
    
    if (length(Sfits) < K) stop("build_h_tau_matrix_cal(): Sfits length < length(visit_times).")
    if (length(Gfits) < K) stop("build_h_tau_matrix_cal(): Gfits length < length(visit_times).")
    
    ## ============ theta-INDEPENDENT precompute (reused across the theta grid) ============
    ## Built once when precomp is NULL; otherwise the caller passes it back to skip
    ## the subject summary, G(tau), interval map, S_k predictor caches and S_k(u-)
    ## denominators (none depend on theta). local({}) keeps the scratch variables
    ## contained and returns only the cache list.
    pc <- if (is.null(precomp)) local({

    # default evaluation grid = all absolute jump times from Gfits
    if (is.null(time_vec)) {
        time_vec <- get_G_jump_times_from_Gfits(Gfits)
    }
    if (length(time_vec) == 0L) stop("build_h_tau_matrix_cal(): time_vec is empty.")
    time_vec <- sort(unique(as.numeric(time_vec)))

    # dat_cal_tau ids define the prediction-model row space for q_lo/q_hi
    ids_tau <- dat_cal_tau[[id.name]]
    if (anyDuplicated(ids_tau)) {
        stop("build_h_tau_matrix_cal(): dat_cal_tau must have one row per id (duplicated ids found).")
    }
    
    # output row ids
    summ <- subject_summary(dat_cal, id.name, stop.name, event.name, event.C.name)
    if (row_data == "all") {
        ids_row <- summ[[id.name]]
    } else {
        ids_row <- ids_tau
    }
    if (anyDuplicated(ids_row)) {
        stop("build_h_tau_matrix_cal(): output row ids have duplicates.")
    }
    n <- length(ids_row)
    m <- length(time_vec)
    
    ## ---------------- subject summary and X ----------------
    idx_summ <- match(as.character(ids_row), as.character(summ[[id.name]]))
    if (anyNA(idx_summ)) {
        miss_ids <- ids_row[is.na(idx_summ)]
        stop("build_h_tau_matrix_cal(): some output-row ids not found in dat_cal summary. Example: ",
             paste(utils::head(miss_ids, 10), collapse = ", "))
    }
    X_row <- summ$X[idx_summ]

    # id alignment from output rows to the tau (prediction-model) row space
    # (theta-independent; the quantile VALUES are filled in after the cache)
    idx_tau_in_row <- match(as.character(ids_row), as.character(ids_tau))
    has_tau <- !is.na(idx_tau_in_row)

    ## ---------------- \tilde G(tau | L_nu) ----------------
    if (is.null(Gtau_vec)) {
        G_tau_obj <- predict_G_from_Gk_matrix(
            dat_new = dat_cal,
            time_vec = tau,
            start.name = start.name, stop.name = stop.name, event.name = event.name,
            event.C.name = event.C.name, id.name = id.name,
            Gfits = Gfits,
            trim.G = NULL,
            return_subject_info = TRUE,
            fast_rsf_predict = fast_rsf_predict
        )
        idx_g_tau <- match(as.character(ids_tau), as.character(G_tau_obj$ids))
        if (anyNA(idx_g_tau)) {
            stop("build_h_tau_matrix_cal(): failed to align internally computed G(tau|.) to dat_cal_tau ids.")
        }
        Gtau_tau <- as.numeric(G_tau_obj$G[idx_g_tau, 1])
    } else {
        if (length(Gtau_vec) == nrow(dat_cal_tau)) {
            Gtau_tau <- as.numeric(Gtau_vec)
        } else if (length(Gtau_vec) == n) {
            Gtau_full <- as.numeric(Gtau_vec)
            Gtau_tau <- NULL
        } else {
            stop("build_h_tau_matrix_cal(): Gtau_vec must have length nrow(dat_cal_tau) or nrow(output rows).")
        }
    }
    if (!exists("Gtau_full", inherits = FALSE)) {
        gt_map <- stats::setNames(Gtau_tau, as.character(ids_tau))
        Gtau_full <- as.numeric(gt_map[as.character(ids_row)])
    }
    Gtau_vec <- pmax(Gtau_full, trim.Gtau)
    
    ## ---------------- map u to interval index k under (t_k, t_{k+1}] ----------------
    # Returns k in {0,1,...,K}; k=0 means u <= visit_times[1]
    k_of_u <- findInterval(time_vec, visit_times, left.open = TRUE, rightmost.closed = FALSE)
    
    # Formula uses the u >= tau branch at u = tau, i.e., force k = nu for that column
    is_tau_col <- abs(time_vec - tau) <= tol
    k_of_u[is_tau_col] <- nu
    
    k_of_u[k_of_u < 0L] <- 0L
    k_of_u[k_of_u > K]  <- K
    
    ## ---------------- cache L_k rows + S_k predictor cache ----------------
    Sk_row_cache  <- vector("list", K)   # stores idx + dat
    Sk_pred_cache <- vector("list", K)   # stores precomputed predictor cache
    
    for (k in seq_len(K)) {
        if (is.null(Sfits[[k]])) {
            Sk_row_cache[[k]]  <- list(idx = integer(0), dat = NULL)
            Sk_pred_cache[k] <- list(NULL)
            next
        }
        
        rc <- build_Sk_landmark_rows_aligned(
            dat_long = dat_cal,
            ids_target = ids_row,
            X_target = X_row,
            visit_times = visit_times,
            k = k,
            start.name = start.name, stop.name = stop.name, event.name = event.name, id.name = id.name,
            Sfit_k = Sfits[[k]]
        )
        Sk_row_cache[[k]] <- rc
        
        if (length(rc$idx) > 0L && !is.null(rc$dat)) {
            Sk_pred_cache[[k]] <- make_Sk_predictor_cache(Sfits[[k]], rc$dat)
        } else {
            Sk_pred_cache[k] <- list(NULL)
        }
    }

    ## ---------------- denominator S_k(u- | L_k) per interval (theta-independent) ----------------
    S_u_left_list <- vector("list", K)
    for (k in seq_len(K)) {
        cols_k <- which(k_of_u == k)
        if (length(cols_k) == 0L) next
        cache_k <- Sk_pred_cache[[k]]
        idxk    <- Sk_row_cache[[k]]$idx
        if (is.null(cache_k) || length(idxk) == 0L) next
        S_u_left_list[[k]] <- pmax(
            predict_Sk_cache_common_times(cache_k, time_vec[cols_k], left_limit = TRUE, tol = tol),
            trim.S
        )
    }

    list(
        nu = nu, ids_tau = ids_tau, ids_row = ids_row, n = n, m = m,
        time_vec = time_vec, X_row = X_row,
        idx_tau_in_row = idx_tau_in_row, has_tau = has_tau,
        Gtau_vec = Gtau_vec, k_of_u = k_of_u,
        Sk_row_cache = Sk_row_cache, Sk_pred_cache = Sk_pred_cache,
        S_u_left_list = S_u_left_list
    )
    }) else precomp   # end theta-independent precompute

    if (isTRUE(precomp_only)) return(pc)

    ## ---------------- unpack the theta-independent cache ----------------
    nu             <- pc$nu
    ids_tau        <- pc$ids_tau
    ids_row        <- pc$ids_row
    n              <- pc$n
    m              <- pc$m
    time_vec       <- pc$time_vec
    X_row          <- pc$X_row
    idx_tau_in_row <- pc$idx_tau_in_row
    has_tau        <- pc$has_tau
    Gtau_vec       <- pc$Gtau_vec
    k_of_u         <- pc$k_of_u
    Sk_row_cache   <- pc$Sk_row_cache
    Sk_pred_cache  <- pc$Sk_pred_cache
    S_u_left_list  <- pc$S_u_left_list

    ## ---------------- prediction quantiles q_theta, q_{1-theta} (theta-DEPENDENT) ----------------
    if (is.finite(support_upper) && theta == 0) {
        # widest candidate member under a known bounded support: (tau, support_upper]
        q_lo_tau <- rep(tau, nrow(pred_surv_cal))
        q_hi_tau <- rep(support_upper, nrow(pred_surv_cal))
    } else {
        q_lo_tau <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = theta)
        q_hi_tau <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = 1 - theta)
    }

    # expand to output rows by id; non-tau ids remain NA (expected when row_data="all")
    q_lo <- rep(NA_real_, n)
    q_hi <- rep(NA_real_, n)
    if (any(has_tau)) {
        q_lo[has_tau] <- q_lo_tau[idx_tau_in_row[has_tau]]
        q_hi[has_tau] <- q_hi_tau[idx_tau_in_row[has_tau]]
    }

    ## ---------------- cache xi_k predictions for k < nu ----------------
    xi_cache <- vector("list", K)
    if (nu > 1L) {
        if (is.null(Xi_fits)) {
            stop("build_h_tau_matrix_cal(): Xi_fits must be provided when tau is not the first visit time.")
        }
        if (length(Xi_fits) < (nu - 1L)) {
            stop("build_h_tau_matrix_cal(): Xi_fits length must be at least nu-1 for u < tau branch.")
        }
        
        for (k in seq_len(nu - 1L)) {
            if (is.null(Xi_fits[[k]])) {
                xi_cache[k] <- list(NULL)
                next
            }
            xi_df <- predict_Yk_condmean_at_tk(Xi_fits[[k]], dat_long_new = dat_cal, clip_to_support = TRUE)
            if (!is.data.frame(xi_df) || !all(c("id", "mYk_hat") %in% names(xi_df))) {
                stop("build_h_tau_matrix_cal(): Xi_fits[[", k, "]] prediction must return columns id,mYk_hat.")
            }
            xi_cache[[k]] <- xi_df
        }
    }
    
    ## ---------------- build h_tau matrix ----------------
    h_mat <- matrix(NA_real_, nrow = n, ncol = m,
                    dimnames = list(as.character(ids_row), as.character(time_vec)))
    
    # Process each interval k
    for (k in seq_len(K)) {
        cols_k <- which(k_of_u == k)
        if (length(cols_k) == 0L) next
        
        cache_k <- Sk_pred_cache[[k]]
        idxk    <- Sk_row_cache[[k]]$idx
        
        if (is.null(cache_k) || length(idxk) == 0L) next
        
        # absolute times in this interval-block
        u_k <- time_vec[cols_k]

        # Denominator S_k(u- | L_k), matrix n_k x c (theta-independent; precomputed)
        S_u_left <- S_u_left_list[[k]]
        
        if (k < nu) {
            ## ---------- u < tau branch ----------
            if (is.null(xi_cache[[k]])) {
                stop("build_h_tau_matrix_cal(): Xi_fits[[", k, "]] (or xi_cache[[k]]) is missing for u < tau branch.")
            }
            
            xi_df <- xi_cache[[k]]
            idx_xi <- match(as.character(ids_row[idxk]), as.character(xi_df$id))
            if (anyNA(idx_xi)) {
                miss_ids <- ids_row[idxk][is.na(idx_xi)]
                stop("build_h_tau_matrix_cal(): failed to align xi_k to calibration rows at k=", k,
                     ". Example ids: ", paste(utils::head(miss_ids, 10), collapse = ", "))
            }
            xi_vec <- as.numeric(xi_df$mYk_hat[idx_xi])   # length n_k
            
            # h = xi_k / S_k(u- | L_k)
            h_block <- sweep(S_u_left, 1, xi_vec, FUN = function(s, xi) xi / s)
            
        } else {
            ## ---------- u >= tau branch ----------
            qlo_k  <- q_lo[idxk]
            qhi_k  <- q_hi[idxk]
            Gtau_k <- Gtau_vec[idxk]
            if (anyNA(qlo_k) || anyNA(qhi_k) || anyNA(Gtau_k)) {
                stop("build_h_tau_matrix_cal(): missing tau-level quantities on u>=tau branch at k=", k,
                     ". Check dat_cal_tau alignment and row_data setting.")
            }
            
            # v_lo = u \/ q_theta(L_nu), v_hi = u \/ q_{1-theta}(L_nu)
            # each is n_k x c
            v_lo <- outer(qlo_k, u_k, FUN = pmax)
            v_hi <- outer(qhi_k, u_k, FUN = pmax)
            
            # S_k(v_lo | L_k), S_k((v_hi)- | L_k)
            S_lo <- predict_Sk_cache_rowtime_matrix(cache_k, v_lo, left_limit = FALSE, tol = tol)
            S_hi_left <- predict_Sk_cache_rowtime_matrix(cache_k, v_hi, left_limit = TRUE, tol = tol)
            
            numer <- S_lo - S_hi_left
            
            # h = Gtau * ( numer / S_k(u- | L_k) - (1 - alpha) )
            h_block <- sweep(numer / S_u_left - (1 - alpha), 1, Gtau_k, FUN = "*")
        }
        
        h_mat[idxk, cols_k] <- h_block
    }
    
    ## ---------------- optional mask u > X_i ----------------
    if (isTRUE(mask_after_X)) {
        mask_mat <- outer(X_row, time_vec, FUN = "<")
        h_mat[mask_mat] <- after_X_value
    }
    
    list(
        h = h_mat,
        row_data = row_data,
        ids = ids_row,
        X = X_row,
        time_vec = time_vec,
        q_lo = q_lo,
        q_hi = q_hi,
        Gtau = Gtau_vec,
        k_of_u = k_of_u
    )
}
