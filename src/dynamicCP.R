## Implementation for dynamic conformal prediction (CP) algorithm with sample splitting - general time-varying covariates setting
##  -- 2 splits: one for training the prediction model and estimate the censoring model (nuisance parameter); 
##               one for computing the miscoverage rate

#' Construct the prediction intervals with split conformal prediction algorithm for survivors at time tau
#' @param dat a data frame in long format with time-varying covariates
#' @param dat_test a data frame in long format with time-varying covariates representing the test data; if is NULL, the calibration data is used
#' @param start.name the variable name for the start of the counting process interval
#' @param stop.name the variable name for the stop of the counting process interval
#' @param event.name the variable name for the event indicator
#' @param event.C.name the variable name for the censoring indicator
#' @param id.name the variable name for the id of subjects
#' @param tau the time point at which prediction intervals are constructed
#' @param model.pred the prediction model
#' @param model.C the censoring model
#' @param covname.pred.baseline names for the baseline covariates to be used in the prediction
#' @param covname.pred.timevarying names for the time-varying covariates to be used in the prediction
#' @param TT.name the variable name for the true event time
#' @param seed seed for randomly split the data
#' @param train_ratio the proportion of data used for the training set
#' @param theta_grid a vector of theta values to be selected from
#' @param alpha the miscoverage tolerance level
dynamicCP_split <- function(dat, dat_test = NULL, start.name, stop.name, event.name, event.C.name, id.name, tau,
                            model.pred = "cox", covname.pred.baseline = NULL, covname.pred.timevarying = NULL, # pred_tvtime = FALSE,
                            model.pred.args = list(),
                            model.C = "km", covname.C.baseline = NULL, covname.C.timevarying = NULL,
                            TT.name = NULL, 
                            seed = NULL, train_ratio = 0.5,
                            theta_grid, 
                            alpha = 0.1, trim.C = 0.05,
                            model.C.args = list()){
    
    
    # split into training the calibration data
    if(!is.null(seed)){
        set.seed(seed)
    }
    
    id = dat[[id.name]]
    id_unique <- unique(id)
    id_tr = sample(id_unique, floor(length(id_unique)*train_ratio), replace = FALSE)
    dat_tr = dat[id %in% id_tr, ] 
    dat_cal = dat[!(id %in% id_tr), ]
    
    if(is.null(dat_test) || nrow(dat_test) == 0){
        dat_test <- dat_cal
    }
    
    
    ## On the training data --------------------------------------------------------
    # Train the prediction model 
    if(model.pred == "cox"){
        dat_tr_tau <- prepare_data_tau(dat_tr, start.name, stop.name, event.name, id.name, tau)
        
        cov_pred_use <- c(covname.pred.baseline, covname.pred.timevarying)
        cov_pred_use <- keep_nonconstant_covars(dat_tr_tau, cov_pred_use)
        
        rhs <- if (length(cov_pred_use) == 0) "1" else paste(cov_pred_use, collapse = "+")
        formula.pred <- as.formula(paste0("Surv(", start.name, ",", stop.name, ",", event.name, ") ~ ", rhs))
        
        fit.pred <- coxph(formula.pred, data = dat_tr_tau)
        
        # prediction on the calibration data and test data
        dat_cal_tau <- prepare_data_tau(dat_cal, start.name, stop.name, event.name, id.name, tau)
        dat_test_tau <- prepare_data_tau(dat_test, start.name, stop.name, event.name, id.name, tau)
        
        # se.fit/conf.type disabled: only $time/$surv are used below, and skipping
        # the standard-error/CI computation leaves $surv identical while making
        # survfit() several times faster on large newdata.
        sf_cal <- survfit(fit.pred, newdata = dat_cal_tau, se.fit = FALSE, conf.type = "none")
        pred_time <- sf_cal$time
        pred_surv_cal <- t(sf_cal$surv)  # row - individual, col - times points

        sf_test <- survfit(fit.pred, newdata = dat_test_tau, se.fit = FALSE, conf.type = "none")
        pred_surv_test <- t(sf_test$surv)  # row - individual, col - times points
        
    } else if (model.pred == "rsf") {
        
        dat_tr_tau <- prepare_data_tau(dat_tr, start.name, stop.name, event.name, id.name, tau)
        
        cov_pred_use <- c(covname.pred.baseline, covname.pred.timevarying)
        cov_pred_use <- keep_nonconstant_covars(dat_tr_tau, cov_pred_use)
        
        make_rsf_df <- function(dat_tau) {
            # remaining time after tau (relative time scale for landmarking)
            dat_tau$.time_rsf <- dat_tau[[stop.name]] - tau
            dat_tau
        }
        
        dat_tr_rsf <- make_rsf_df(dat_tr_tau)
        
        dat_cal_tau  <- prepare_data_tau(dat_cal,  start.name, stop.name, event.name, id.name, tau)
        dat_test_tau <- prepare_data_tau(dat_test, start.name, stop.name, event.name, id.name, tau)
        
        dat_cal_rsf  <- make_rsf_df(dat_cal_tau)
        dat_test_rsf <- make_rsf_df(dat_test_tau)
        
        if (length(cov_pred_use) == 0) {
            # KM on remaining time after tau, then shift time grid back to original scale
            km_fit <- survfit(Surv(.time_rsf, dat_tr_rsf[[event.name]]) ~ 1, data = dat_tr_rsf)
            
            pred_time_rel <- km_fit$time                 # time since tau
            pred_time     <- tau + pred_time_rel         # ORIGINAL time scale
            S_km <- km_fit$surv
            
            pred_surv_cal  <- matrix(rep(S_km, nrow(dat_cal_rsf)),  nrow = nrow(dat_cal_rsf),  byrow = TRUE)
            pred_surv_test <- matrix(rep(S_km, nrow(dat_test_rsf)), nrow = nrow(dat_test_rsf), byrow = TRUE)
            
        } else {
            
            rhs <- paste(cov_pred_use, collapse = "+")
            rsf_formula <- as.formula(paste0("Surv(.time_rsf, ", event.name, ") ~ ", rhs))
            
            fit.pred <- randomForestSRC::rfsrc(
                formula  = rsf_formula,
                data     = dat_tr_rsf,
                ntree    = 1000,
                nodesize = 15, # default for survival forest
                forest   = TRUE  # Save forest object for future prediction
            )
            
            # Predict on calibration to define a RELATIVE time grid
            pr_cal <- predict(fit.pred, newdata = dat_cal_rsf)
            pred_time_rel <- pr_cal$time.interest        # time since tau
            
            pred_time     <- tau + pred_time_rel         # ORIGINAL time scale (what you want)
            pred_surv_cal <- pr_cal$survival             # rows = subjects, cols aligned with pred_time_rel
            
            # Predict on test using the SAME RELATIVE grid
            pr_test <- predict(
                fit.pred,
                newdata = dat_test_rsf,
                time.interest = pred_time_rel
            )
            pred_surv_test <- pr_test$survival
        }
    } else if (model.pred %in% c("ddh", "ddh_py")) {
        
        if (!exists("predict_surv_ddh", mode = "function") &&
            !exists("predict_surv_ddh_py", mode = "function")) {
            stop("model.pred='ddh' requires source('src/ddh_bridge.R') before running.")
        }
        predict_fn <- if (exists("predict_surv_ddh", mode = "function")) {
            predict_surv_ddh
        } else {
            predict_surv_ddh_py
        }
        
        dat_tr_tau <- prepare_data_tau(dat_tr, start.name, stop.name, event.name, id.name, tau)
        dat_cal_tau_pred <- prepare_data_tau(dat_cal, start.name, stop.name, event.name, id.name, tau)
        dat_test_tau_pred <- prepare_data_tau(dat_test, start.name, stop.name, event.name, id.name, tau)
        
        cov_pred_use <- c(covname.pred.baseline, covname.pred.timevarying)
        cov_pred_use <- keep_nonconstant_covars(dat_tr_tau, cov_pred_use)
        
        ddh_pred <- predict_fn(
            dat_tr = dat_tr,
            dat_cal = dat_cal,
            dat_test = dat_test,
            dat_tr_tau = dat_tr_tau,
            dat_cal_tau = dat_cal_tau_pred,
            dat_test_tau = dat_test_tau_pred,
            start.name = start.name,
            stop.name = stop.name,
            event.name = event.name,
            id.name = id.name,
            tau = tau,
            cov_pred_use = cov_pred_use,
            model.pred.args = model.pred.args
        )
        
        pred_time <- ddh_pred$pred_time
        pred_surv_cal <- ddh_pred$pred_surv_cal
        pred_surv_test <- ddh_pred$pred_surv_test
        
    } else {
        stop("This prediction model is not implemented in the function")
    }
    
    
    # Fit the censoring model on training data, compute the IPCW weights for the calibration data and test data
    if(model.C == "km"){
        dat_tr_0 <- prepare_data_tau(dat_tr, start.name, stop.name, event.C.name, id.name, tau = 0)
        formula_C <- as.formula(paste0("Surv(", stop.name, ", ", event.C.name, ") ~ 1"))
        kmfit_C <- survfit(formula_C, data = dat_tr_0)
        km_C <- summary(kmfit_C)
        Sc_stepfun <- stepfun(km_C$time, c(1, km_C$surv))
        
        dat_cal$surv.C.subj <- ave(dat_cal[, stop.name], dat_cal[, id.name],
                                   FUN = function(x) Sc_stepfun(max(x)))
        dat_test$surv.C.subj <- ave(dat_test[, stop.name], dat_test[, id.name],
                                    FUN = function(x) Sc_stepfun(max(x)))
        
        dat_cal_tau <- prepare_data_tau(dat_cal, start.name, stop.name, event.name, id.name, tau)
        dat_test_tau <- prepare_data_tau(dat_test, start.name, stop.name, event.name, id.name, tau)
        
    }else if(model.C == "cox"){
        
        cov_C_use <- c(covname.C.baseline, covname.C.timevarying)
        cov_C_use <- keep_nonconstant_covars(dat_tr, cov_C_use)
        
        rhsC <- if (length(cov_C_use) == 0) "1" else paste(cov_C_use, collapse = "+")
        formula.C <- as.formula(paste0("Surv(", start.name, ",", stop.name, ",", event.C.name, ") ~ ", rhsC))
        
        fit.C <- coxph(formula.C, data = dat_tr)
        
        dat_cal$surv.C <- predict(fit.C, newdata = dat_cal, type = "survival")
        dat_test$surv.C <- predict(fit.C, newdata = dat_test, type = "survival")
        
        dat_cal$surv.C.subj <- ave(dat_cal[, "surv.C"], dat_cal[, id.name],
                                   FUN = function(x) prod(x))
        dat_test$surv.C.subj <- ave(dat_test[, "surv.C"], dat_test[, id.name],
                                    FUN = function(x) prod(x))
        
        dat_cal_tau <- prepare_data_tau(dat_cal, start.name, stop.name, event.name, id.name, tau)
        dat_test_tau <- prepare_data_tau(dat_test, start.name, stop.name, event.name, id.name, tau)
        
    } else if (model.C == "rsf") {
        
        if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
            stop("Package 'randomForestSRC' is required for model.C='rsf'.")
        }
        
        dat_tr_C0 <- prepare_data_tau(dat_tr, start.name, stop.name, event.C.name, id.name, tau = 0)
        dat_cal_C0 <- prepare_data_tau(dat_cal, start.name, stop.name, event.C.name, id.name, tau = 0)
        dat_test_C0 <- prepare_data_tau(dat_test, start.name, stop.name, event.C.name, id.name, tau = 0)
        
        cov_C_use <- c(covname.C.baseline, covname.C.timevarying)
        cov_C_use <- keep_nonconstant_covars(dat_tr_C0, cov_C_use)
        min_censor_events <- as.integer(model.C.args[["rsf_min_censor_events"]] %||% 10L)
        use_km_C <- length(cov_C_use) == 0 || sum(dat_tr_C0[[event.C.name]] == 1, na.rm = TRUE) < min_censor_events
        
        if (use_km_C) {
            formula_C <- as.formula(paste0("Surv(", stop.name, ", ", event.C.name, ") ~ 1"))
            kmfit_C <- survfit(formula_C, data = dat_tr_C0)
            km_C <- summary(kmfit_C)
            Sc_stepfun <- stepfun(km_C$time, c(1, km_C$surv))
            dat_cal_C0$surv.C.subj <- Sc_stepfun(dat_cal_C0[[stop.name]])
            dat_test_C0$surv.C.subj <- Sc_stepfun(dat_test_C0[[stop.name]])
        } else {
            dat_tr_C0$.time_C_rsf <- dat_tr_C0[[stop.name]]
            dat_cal_C0$.time_C_rsf <- dat_cal_C0[[stop.name]]
            dat_test_C0$.time_C_rsf <- dat_test_C0[[stop.name]]
            
            rhsC <- paste(cov_C_use, collapse = "+")
            formula.C <- as.formula(paste0("Surv(.time_C_rsf, ", event.C.name, ") ~ ", rhsC))
            rsf_defaults <- list(ntree = 500L, nodesize = 15L, nsplit = 10L, forest = TRUE)
            rsf_args <- modifyList(rsf_defaults, model.C.args)
            rsf_args$rsf_min_censor_events <- NULL
            
            fit.C <- do.call(
                randomForestSRC::rfsrc,
                c(list(formula = formula.C, data = dat_tr_C0), rsf_args)
            )
            
            predict_rsf_surv_at <- function(fit, newdata, time_vec) {
                pr <- predict(fit, newdata = newdata)
                pred_time <- pr$time.interest
                pred_surv <- pr$survival
                idx <- findInterval(time_vec, pred_time)
                out <- rep(1, length(time_vec))
                use <- idx > 0
                out[use] <- pred_surv[cbind(which(use), idx[use])]
                out
            }
            
            dat_cal_C0$surv.C.subj <- predict_rsf_surv_at(fit.C, dat_cal_C0, dat_cal_C0[[stop.name]])
            dat_test_C0$surv.C.subj <- predict_rsf_surv_at(fit.C, dat_test_C0, dat_test_C0[[stop.name]])
        }
        
        dat_cal_tau <- prepare_data_tau(dat_cal, start.name, stop.name, event.name, id.name, tau)
        dat_test_tau <- prepare_data_tau(dat_test, start.name, stop.name, event.name, id.name, tau)
        
        idx_cal <- match(as.character(dat_cal_tau[[id.name]]), as.character(dat_cal_C0[[id.name]]))
        idx_test <- match(as.character(dat_test_tau[[id.name]]), as.character(dat_test_C0[[id.name]]))
        if (anyNA(idx_cal)) {
            miss_ids <- unique(dat_cal_tau[[id.name]][is.na(idx_cal)])
            stop("rsf censoring: missing baseline G prediction for calibration ids. Example: ",
                 paste(utils::head(miss_ids, 10), collapse = ", "))
        }
        if (anyNA(idx_test)) {
            miss_ids <- unique(dat_test_tau[[id.name]][is.na(idx_test)])
            stop("rsf censoring: missing baseline G prediction for test ids. Example: ",
                 paste(utils::head(miss_ids, 10), collapse = ", "))
        }
        
        dat_cal_tau$surv.C.subj <- pmax(dat_cal_C0$surv.C.subj[idx_cal], trim.C)
        dat_test_tau$surv.C.subj <- pmax(dat_test_C0$surv.C.subj[idx_test], trim.C)
        
    } else if (model.C == "xgb_global") {
        
        if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("Package 'xgboost' is required for model.C='xgb_global'.")
        }
        
        cov_C_use <- c(covname.C.baseline, covname.C.timevarying)
        cov_C_use <- keep_nonconstant_covars(dat_tr, cov_C_use)
        
        if (length(cov_C_use) == 0) {
            dat_tr_0 <- prepare_data_tau(dat_tr, start.name, stop.name, event.C.name, id.name, tau = 0)
            formula_C <- as.formula(paste0("Surv(", stop.name, ", ", event.C.name, ") ~ 1"))
            kmfit_C <- survfit(formula_C, data = dat_tr_0)
            km_C <- summary(kmfit_C)
            Sc_stepfun <- stepfun(km_C$time, c(1, km_C$surv))
            
            dat_cal$surv.C.subj <- ave(dat_cal[, stop.name], dat_cal[, id.name],
                                       FUN = function(x) Sc_stepfun(max(x)))
            dat_test$surv.C.subj <- ave(dat_test[, stop.name], dat_test[, id.name],
                                        FUN = function(x) Sc_stepfun(max(x)))
        } else {
            X_tr <- xgb_build_design_matrix(dat_tr, cov_C_use)
            time_cox <- pmax(as.numeric(dat_tr[[stop.name]]), 1e-8)
            y_cox <- ifelse(dat_tr[[event.C.name]] == 1L, time_cox, -time_cox)
            
            dtrain <- xgboost::xgb.DMatrix(data = X_tr$X, label = y_cox)
            
            xgb_defaults <- list(
                eta = 0.05,
                max_depth = 4L,
                min_child_weight = 1,
                subsample = 0.8,
                colsample_bytree = 0.8,
                lambda = 1,
                gamma = 0
            )
            
            xgb_nrounds <- as.integer(model.C.args[["xgb_nrounds"]] %||% 200L)
            xgb_verbose <- as.integer(model.C.args[["xgb_verbose"]] %||% 0L)
            xgb_seed <- model.C.args[["xgb_seed"]]
            xgb_nthread <- model.C.args[["xgb_nthread"]]
            
            xgb_params <- c(
                list(objective = "survival:cox", eval_metric = "cox-nloglik"),
                xgb_defaults,
                model.C.args[["xgb_params_extra"]] %||% list()
            )
            for (nm in names(xgb_defaults)) {
                arg_nm <- paste0("xgb_", nm)
                if (!is.null(model.C.args[[arg_nm]])) {
                    xgb_params[[nm]] <- model.C.args[[arg_nm]]
                }
            }
            # xgboost R ignores params$seed; use set.seed() before training instead.
            if (!is.null(xgb_params$seed)) xgb_params$seed <- NULL
            if (!is.null(xgb_nthread)) xgb_params$nthread <- as.integer(xgb_nthread)
            if (!is.null(xgb_seed)) set.seed(as.integer(xgb_seed))
            
            fit_xgb <- xgboost::xgb.train(
                params = xgb_params,
                data = dtrain,
                nrounds = xgb_nrounds,
                verbose = xgb_verbose
            )
            
            lp_tr <- as.numeric(stats::predict(fit_xgb, newdata = X_tr$X, outputmargin = TRUE))
            lp_tr[!is.finite(lp_tr)] <- 0
            lp_tr <- pmax(pmin(lp_tr, 30), -30)

            # Center the offset before the anchor fit: basehaz(fit, centered = FALSE)
            # on an offset-only coxph returns the cumulative hazard at the MEAN offset,
            # so an uncentered lp would scale G by exp(mean(lp)) once exp(lp) is
            # multiplied back in at prediction time. Centering here (and subtracting
            # lp_center again at prediction) anchors the baseline at offset 0,
            # mirroring the xgb_cox nuisance paths in helpers_AIPCW.R.
            lp_center <- mean(lp_tr)
            dat_anchor <- dat_tr
            dat_anchor$lp_off <- lp_tr - lp_center
            form_anchor <- as.formula(
                paste0("Surv(", start.name, ",", stop.name, ",", event.C.name, ") ~ offset(lp_off)")
            )
            fit_anchor <- coxph(form_anchor, data = dat_anchor, ties = "breslow", x = TRUE)
            bh <- survival::basehaz(fit_anchor, centered = FALSE)
            
            dat_cal$surv.C <- predict_surv_xgb_cox_interval(
                fit_xgb = fit_xgb,
                bh_time = bh$time,
                bh_cumhaz = bh$hazard,
                newdata = dat_cal,
                cov_C_use = cov_C_use,
                design_info = X_tr$design_info,
                start.name = start.name,
                stop.name = stop.name,
                lp_center = lp_center
            )
            dat_test$surv.C <- predict_surv_xgb_cox_interval(
                fit_xgb = fit_xgb,
                bh_time = bh$time,
                bh_cumhaz = bh$hazard,
                newdata = dat_test,
                cov_C_use = cov_C_use,
                design_info = X_tr$design_info,
                start.name = start.name,
                stop.name = stop.name,
                lp_center = lp_center
            )
            
            dat_cal$surv.C.subj <- ave(dat_cal[, "surv.C"], dat_cal[, id.name],
                                       FUN = function(x) prod(x))
            dat_test$surv.C.subj <- ave(dat_test[, "surv.C"], dat_test[, id.name],
                                        FUN = function(x) prod(x))
        }
        
        dat_cal_tau <- prepare_data_tau(dat_cal, start.name, stop.name, event.name, id.name, tau)
        dat_test_tau <- prepare_data_tau(dat_test, start.name, stop.name, event.name, id.name, tau)
        
    } else if (model.C == "xgb_piecewise") {
        
        if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("Package 'xgboost' is required for model.C='xgb_piecewise'.")
        }
        
        if (!exists("fit_Gk_list_xgb_piecewise", mode = "function") ||
            !exists("compute_Hk_at_tk_xgb_piecewise", mode = "function")) {
            stop("model.C='xgb_piecewise' requires source('src/helpers.R') ",
                 "with fit_Gk_list_xgb_piecewise() and compute_Hk_at_tk_xgb_piecewise().")
        }
        
        visit_times <- model.C.args[["visit_times"]]
        if (is.null(visit_times)) {
            visit_times <- sort(unique(as.numeric(dat_tr[[start.name]])))
        }
        visit_times <- sort(unique(as.numeric(visit_times[is.finite(visit_times)])))
        if (length(visit_times) == 0) {
            stop("model.C='xgb_piecewise': visit_times is empty.")
        }
        if (abs(visit_times[1] - 0) > 1e-12) {
            visit_times <- sort(unique(c(0, visit_times)))
        }
        
        use_history <- isTRUE(model.C.args[["use_history"]])
        history_mode <- model.C.args[["history_mode"]] %||% "wide"
        trim_piecewise <- as.numeric(model.C.args[["trim_piecewise"]] %||% 1e-7)
        
        xgb_piece_args <- list(
            nrounds = as.integer(model.C.args[["xgb_nrounds"]] %||% 200L),
            verbose = as.integer(model.C.args[["xgb_verbose"]] %||% 0L),
            eta = as.numeric(model.C.args[["xgb_eta"]] %||% 0.05),
            max_depth = as.integer(model.C.args[["xgb_max_depth"]] %||% 4L),
            min_child_weight = as.numeric(model.C.args[["xgb_min_child_weight"]] %||% 1),
            subsample = as.numeric(model.C.args[["xgb_subsample"]] %||% 0.8),
            colsample_bytree = as.numeric(model.C.args[["xgb_colsample_bytree"]] %||% 0.8),
            lambda = as.numeric(model.C.args[["xgb_lambda"]] %||% 1),
            gamma = as.numeric(model.C.args[["xgb_gamma"]] %||% 0),
            seed = model.C.args[["xgb_seed"]],
            nthread = model.C.args[["xgb_nthread"]],
            params_extra = model.C.args[["xgb_params_extra"]] %||% list(),
            min_censor_events = as.integer(model.C.args[["xgb_min_censor_events"]] %||% 1L)
        )
        if (!is.list(xgb_piece_args$params_extra)) {
            stop("model.C.args[['xgb_params_extra']] must be a list for model.C='xgb_piecewise'.")
        }
        
        Gfits <- fit_Gk_list_xgb_piecewise(
            dat_tr = dat_tr,
            visit_times = visit_times,
            start.name = start.name,
            stop.name = stop.name,
            event.name = event.name,
            event.C.name = event.C.name,
            id.name = id.name,
            covname.C.baseline = covname.C.baseline,
            covname.C.timevarying = covname.C.timevarying,
            use_history = use_history,
            history_mode = history_mode,
            xgb_args = xgb_piece_args
        )
        
        G0X_cal_df <- compute_Hk_at_tk_xgb_piecewise(
            dat_long = dat_cal,
            visit_times = visit_times,
            k = 1,
            start.name = start.name,
            stop.name = stop.name,
            event.name = event.name,
            event.C.name = event.C.name,
            id.name = id.name,
            Gfits = Gfits,
            trim.C = trim_piecewise
        )
        G0X_test_df <- compute_Hk_at_tk_xgb_piecewise(
            dat_long = dat_test,
            visit_times = visit_times,
            k = 1,
            start.name = start.name,
            stop.name = stop.name,
            event.name = event.name,
            event.C.name = event.C.name,
            id.name = id.name,
            Gfits = Gfits,
            trim.C = trim_piecewise
        )
        
        dat_cal_tau <- prepare_data_tau(dat_cal, start.name, stop.name, event.name, id.name, tau)
        dat_test_tau <- prepare_data_tau(dat_test, start.name, stop.name, event.name, id.name, tau)
        
        idx_cal <- match(as.character(dat_cal_tau[[id.name]]), as.character(G0X_cal_df$id))
        idx_test <- match(as.character(dat_test_tau[[id.name]]), as.character(G0X_test_df$id))
        
        dat_cal_tau$surv.C.subj <- G0X_cal_df$Hk[idx_cal]
        dat_test_tau$surv.C.subj <- G0X_test_df$Hk[idx_test]
        
        if (anyNA(dat_cal_tau$surv.C.subj)) {
            miss_ids <- unique(dat_cal_tau[[id.name]][is.na(dat_cal_tau$surv.C.subj)])
            stop("xgb_piecewise censoring: missing G(0->X) for some calibration ids. Example: ",
                 paste(utils::head(miss_ids, 10), collapse = ", "))
        }
        if (anyNA(dat_test_tau$surv.C.subj)) {
            miss_ids <- unique(dat_test_tau[[id.name]][is.na(dat_test_tau$surv.C.subj)])
            stop("xgb_piecewise censoring: missing G(0->X) for some test ids. Example: ",
                 paste(utils::head(miss_ids, 10), collapse = ", "))
        }
        
        dat_cal_tau$surv.C.subj <- pmax(dat_cal_tau$surv.C.subj, trim.C)
        dat_test_tau$surv.C.subj <- pmax(dat_test_tau$surv.C.subj, trim.C)
        
    } else {
        stop("This censoring model is not implemented in the function")
    }


    
    ## On the calibration data -----------------------------------------------------
    
    # search for the smallest theta
    n_cal_tau = nrow(dat_cal_tau)
    n_grid <- length(theta_grid)
    # Precompute all candidate interval endpoints in one vectorized pass
    # (theta-independent work hoisted out of the search loop; the values are
    # identical to the former per-(i, j) surv_quantile() calls).
    Q_cal <- surv_quantile_grid(pred_time, pred_surv_cal, theta_grid)
    X_cal <- dat_cal_tau[[stop.name]]
    event_cal <- dat_cal_tau[[event.name]]
    G_floor_cal <- pmax(dat_cal_tau[["surv.C.subj"]], trim.C)

    hat_W_IPCW = 1
    j = 1
    while(hat_W_IPCW[j] >= 0 & j <= n_grid/2){

        # hat_W_IPCW
        lo_j <- Q_cal[, j]
        hi_j <- Q_cal[, n_grid + 1 - j]
        inside_j <- as.numeric(X_cal >= lo_j & X_cal <= hi_j)
        U_vec <- ifelse(event_cal == 0, 0, (inside_j - (1 - alpha)) / G_floor_cal)

        hat_W_IPCW <- c(hat_W_IPCW, mean(U_vec))

        j <- j+1
    } # End while(hat_W_IPCW[j] >= 0 & j <= length(theta_grid))
    
    hat_W_IPCW = hat_W_IPCW[-1] # remove the place holder.
    # Select the best theta
    # Bounded scan: if the estimating function never goes negative on the searched
    # grid (e.g., all uncensored survivors fall inside even the narrowest candidate
    # interval), stop at the last evaluated index instead of reading past the end
    # of hat_W_IPCW (which errors with "missing value where TRUE/FALSE needed").
    j <- 1
    while (j <= length(hat_W_IPCW) && hat_W_IPCW[j] >= 0) {
        j <- j + 1
    }
    theta_index_lower_IPCW <- max(j - 1, 1) # Ensure at least index 1
    theta_index_upper_IPCW <- length(theta_grid) + 1 - theta_index_lower_IPCW
    theta_hat_IPCW <- theta_grid[theta_index_lower_IPCW]
    
    # Raw model quantile intervals at alpha/2 and 1 - alpha/2
    Q_model_cal <- surv_quantile_grid(pred_time, pred_surv_cal, c(alpha / 2, 1 - alpha / 2))
    lower_model_cal <- Q_model_cal[, 1]
    upper_model_cal <- Q_model_cal[, 2]
    Q_model_test <- surv_quantile_grid(pred_time, pred_surv_test, c(alpha / 2, 1 - alpha / 2))
    lower_model_test <- Q_model_test[, 1]
    upper_model_test <- Q_model_test[, 2]


    # Compute the IPCW lower/upper bounds at the selected grid indices
    # (calibration side reuses the precomputed Q_cal matrix)
    n_test_tau = nrow(dat_test_tau)
    lower_IPCW_cal <- Q_cal[, theta_index_lower_IPCW]
    upper_IPCW_cal <- Q_cal[, theta_index_upper_IPCW]
    Q_test_sel <- surv_quantile_grid(pred_time, pred_surv_test,
                                     theta_grid[c(theta_index_lower_IPCW, theta_index_upper_IPCW)])
    lower_IPCW_test <- Q_test_sel[, 1]
    upper_IPCW_test <- Q_test_sel[, 2]
    
    
    # Estimate the empirical coverage in the calibration data and test data using IPCW
    w_cal  <- dat_cal_tau[[event.name]] / pmax(dat_cal_tau[["surv.C.subj"]], trim.C)
    w_test <- dat_test_tau[[event.name]] / pmax(dat_test_tau[["surv.C.subj"]], trim.C)
    
    coverage_cal  <- mean(w_cal  * (lower_IPCW_cal  <= dat_cal_tau[[stop.name]]  & dat_cal_tau[[stop.name]]  <= upper_IPCW_cal)) / mean(w_cal)
    coverage_test <- mean(w_test * (lower_IPCW_test <= dat_test_tau[[stop.name]] & dat_test_tau[[stop.name]] <= upper_IPCW_test)) / mean(w_test)
    
    if(!is.null(TT.name)){
        coverage_cal_T  <- mean((lower_IPCW_cal  <= dat_cal_tau[[TT.name]]  & dat_cal_tau[[TT.name]]  <= upper_IPCW_cal)) 
        coverage_test_T <- mean((lower_IPCW_test <= dat_test_tau[[TT.name]] & dat_test_tau[[TT.name]] <= upper_IPCW_test)) 
    }else{
        coverage_cal_T  <- NA
        coverage_test_T <- NA
    }
    
    
    lower_IPCW_DR_cal <- pmin(lower_IPCW_cal, lower_model_cal)
    upper_IPCW_DR_cal <- pmax(upper_IPCW_cal, upper_model_cal)
    lower_IPCW_DR_test <- pmin(lower_IPCW_test, lower_model_test)
    upper_IPCW_DR_test <- pmax(upper_IPCW_test, upper_model_test)
    
    return(list(lower_IPCW_cal = lower_IPCW_cal, upper_IPCW_cal = upper_IPCW_cal, coverage_cal = coverage_cal, 
                lower_IPCW_test = lower_IPCW_test, upper_IPCW_test = upper_IPCW_test, coverage_test = coverage_test, 
                lower_model_cal = lower_model_cal, upper_model_cal = upper_model_cal,
                lower_model_test = lower_model_test, upper_model_test = upper_model_test,
                lower_IPCW_DR_cal = lower_IPCW_DR_cal, upper_IPCW_DR_cal = upper_IPCW_DR_cal,
                lower_IPCW_DR_test = lower_IPCW_DR_test, upper_IPCW_DR_test = upper_IPCW_DR_test,
                coverage_cal_T = coverage_cal_T, coverage_test_T = coverage_test_T, 
                theta_hat_IPCW = theta_hat_IPCW,
                theta_index_IPCW = theta_index_lower_IPCW,
                hat_W_IPCW = hat_W_IPCW,
                id_tr = id_tr,
                survC_cal = dat_cal_tau[["surv.C.subj"]], survC_test = dat_test_tau[["surv.C.subj"]]))

}



# NULL-coalescing helper used for defaults
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Build a design matrix for xgboost and keep train-time columns for prediction-time alignment
xgb_build_design_matrix <- function(dat, cov_use, design_info = NULL) {
    if (length(cov_use) == 0) {
        X0 <- matrix(0, nrow = nrow(dat), ncol = 0)
        return(list(X = X0, design_info = list(cov_use = cov_use, ref_cols = character(0))))
    }
    
    fmla <- as.formula(paste0("~ ", paste(cov_use, collapse = "+")))
    mm <- stats::model.matrix(fmla, data = dat, na.action = na.pass)
    if ("(Intercept)" %in% colnames(mm)) {
        mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    }
    storage.mode(mm) <- "double"
    
    if (is.null(design_info)) {
        design_info <- list(cov_use = cov_use, ref_cols = colnames(mm))
        return(list(X = mm, design_info = design_info))
    }
    
    ref_cols <- design_info$ref_cols %||% character(0)
    if (length(ref_cols) == 0) {
        X0 <- matrix(0, nrow = nrow(dat), ncol = 0)
        return(list(X = X0, design_info = design_info))
    }
    
    if (ncol(mm) == 0) {
        mm <- matrix(0, nrow = nrow(dat), ncol = 0)
    }
    
    miss_cols <- setdiff(ref_cols, colnames(mm))
    if (length(miss_cols) > 0) {
        add <- matrix(0, nrow = nrow(mm), ncol = length(miss_cols))
        colnames(add) <- miss_cols
        mm <- cbind(mm, add)
    }
    extra_cols <- setdiff(colnames(mm), ref_cols)
    if (length(extra_cols) > 0) {
        mm <- mm[, setdiff(colnames(mm), extra_cols), drop = FALSE]
    }
    mm <- mm[, ref_cols, drop = FALSE]
    storage.mode(mm) <- "double"
    
    list(X = mm, design_info = design_info)
}

# Predict interval survival P(C > stop | C > start, Z) using xgb Cox + anchored baseline hazard
# lp_center must match the centering used when fitting the anchored baseline hazard
# (see the xgb_global branch above), so that the baseline corresponds to offset 0.
predict_surv_xgb_cox_interval <- function(fit_xgb, bh_time, bh_cumhaz, newdata,
                                          cov_C_use, design_info, start.name, stop.name,
                                          lp_center = 0) {
    if (nrow(newdata) == 0) return(numeric(0))
    if (length(bh_time) == 0 || length(bh_cumhaz) == 0) return(rep(1, nrow(newdata)))

    X_new <- xgb_build_design_matrix(newdata, cov_C_use, design_info = design_info)
    lp <- as.numeric(stats::predict(fit_xgb, newdata = X_new$X, outputmargin = TRUE))
    lp[!is.finite(lp)] <- 0
    lp <- pmax(pmin(lp, 30), -30)
    lp <- lp - lp_center
    
    t_start <- as.numeric(newdata[[start.name]])
    t_stop <- as.numeric(newdata[[stop.name]])
    
    max_h <- max(bh_cumhaz, na.rm = TRUE)
    H0_start <- stats::approx(
        x = bh_time, y = bh_cumhaz, xout = pmax(t_start, 0),
        method = "constant", f = 0, yleft = 0, yright = max_h, ties = "ordered"
    )$y
    H0_stop <- stats::approx(
        x = bh_time, y = bh_cumhaz, xout = pmax(t_stop, 0),
        method = "constant", f = 0, yleft = 0, yright = max_h, ties = "ordered"
    )$y
    
    dH <- pmax(H0_stop - H0_start, 0)
    S <- exp(-dH * exp(lp))
    S[S < 0] <- 0
    S[S > 1] <- 1
    S
}


# -------- helper: only keep bariables that are not constant --------
keep_nonconstant_covars <- function(df, vars) {
    vars <- vars[!is.na(vars) & nzchar(vars)]
    if (length(vars) == 0) return(character(0))
    
    keep <- vapply(vars, function(v) {
        if (!v %in% names(df)) return(FALSE)
        x <- df[[v]]
        if (is.factor(x)) {
            nlevels(droplevels(x)) >= 2
        } else {
            ux <- unique(x[is.finite(x)])
            length(ux) >= 2
        }
    }, logical(1))
    
    vars[keep]
}




#' Prepare the data set for prediction at time tau
#' @param dat a data frame in long format with time varying covariates
#' @param start.name the variable name for the start of the counting process interval
#' @param stop.name the variable name for the stop of the counting process interval
#' @param event.name the variable name for the event indicator
#' @param id.name the variable name for the id of subjects
#' @param tau the time point at which prediction intervals are constructed
#' @return A data set that contains subjects at risk with the time-varying covariate value at time tau
prepare_data_tau <- function(dat, start.name, stop.name, event.name, id.name, tau){
    
    dat <- dat[order(dat[[id.name]], dat[[start.name]]), ]
    select_tau = (dat[, start.name] <= tau & tau < dat[, stop.name])
    
    # The censored event time and event indicator for each subject
    stop_subj <- ave(dat[, stop.name], dat[, id.name], FUN = function(x) x[length(x)])
    event_subj <- ave(dat[, event.name], dat[, id.name], FUN = function(x) as.integer(any(x == 1)))
    
    dat[, stop.name] <- stop_subj
    dat[, event.name] <- event_subj
    
    dat_tau = dat[select_tau, , drop = FALSE]
    
    dat_tau[, start.name] <- tau  # set the start time to tau
    
    return(dat_tau)
}











## ----------------------------------------------------------------------------------

#' Compute the estimated survival quantile for a given probability level
#'
#' Finds the smallest time `t` such that the survival probability is at most `1 - beta_grid[j]`
#' for an individual `i`. This corresponds to the estimated quantile of the survival distribution.
#'
#' @param i Index of the individual in the calibration dataset
#' @param j Index in the beta_grid, corresponding to a probability threshold
#' @param pred_time Vector of time points for which survival probabilities were estimated
#' @param pred_surv Matrix of survival probabilities (rows correspond to individuals, columns to time points)
#' @param beta_grid Vector of probability thresholds for quantile computation
#' @return Estimated survival quantile for individual `i` at probability level `beta_grid[j]`
#'
surv_quantile <- function(
        i, j, pred_time, pred_surv, beta_grid
) {
    # If the probability threshold is greater than the final survival probability, return the max time
    if (beta_grid[j] > 1 - pred_surv[i, length(pred_time)]) {
        return(max(pred_time))
    } else {
        # Find the smallest time index where survival probability is at most 1 - beta_grid[j]
        return(pred_time[min(which(1 - pred_surv[i, ] >= beta_grid[j]))])
    }
}


#' Vectorized version of surv_quantile() over a grid of probability levels
#'
#' Returns the n x length(beta_grid) matrix Q with
#'   Q[i, j] = surv_quantile(i, j, pred_time, pred_surv, beta_grid),
#' i.e., exactly reproducing the per-(i, j) lookups:
#'   if (beta > 1 - pred_surv[i, m]) max(pred_time)
#'   else pred_time[min(which(1 - pred_surv[i, ] >= beta))]
#' For the (always satisfied in this pipeline) case of a nonincreasing survival
#' curve, the first index with cdf >= beta equals #(cdf < beta) + 1, computed in
#' a single findInterval() call per subject; a non-monotone row falls back to
#' the literal which() scan so the result is identical in all cases.
surv_quantile_grid <- function(pred_time, pred_surv, beta_grid) {
    if (!is.matrix(pred_surv)) {
        pred_surv <- matrix(pred_surv, nrow = 1)
    }
    m <- length(pred_time)
    tmax <- max(pred_time)
    n <- nrow(pred_surv)
    Q <- matrix(tmax, nrow = n, ncol = length(beta_grid))
    for (i in seq_len(n)) {
        cdf_i <- 1 - pred_surv[i, ]
        feas <- beta_grid <= cdf_i[m]
        if (!any(feas)) next
        if (!is.unsorted(cdf_i)) {
            idx <- findInterval(beta_grid[feas], cdf_i, left.open = TRUE) + 1L
        } else {
            idx <- vapply(beta_grid[feas], function(b) min(which(cdf_i >= b)), integer(1))
        }
        Q[i, feas] <- pred_time[idx]
    }
    Q
}








# ----------------- Build static CI from IPCW from dCP at tau = 0 ---------------------
build_static_from_dynamic <- function(results_bounds, tau_grid, TT.name = "TT") {
    
    R <- length(results_bounds)
    results_bounds_static  <- vector("list", R)
    results_summary_static <- vector("list", R)
    
    stopifnot(tau_grid[1] == 0)
    
    for (r in 1:R) {
        
        bounds_static_r <- vector("list", length(tau_grid))
        names(bounds_static_r) <- paste0("tau_", tau_grid)
        
        tab_static_r <- data.frame(
            rep = r,
            tau = tau_grid,
            ok = TRUE,
            coverage_cal_ipcw  = NA_real_,
            coverage_test_ipcw = NA_real_,
            coverage_cal_true  = NA_real_,
            coverage_test_true = NA_real_,
            mean_len_cal  = NA_real_,
            mean_len_test = NA_real_
        )
        
        b0_test <- results_bounds[[r]][[1]]$test
        b0_cal  <- results_bounds[[r]][[1]]$cal
        
        for (k in seq_along(tau_grid)) {
            
            bk_test <- results_bounds[[r]][[k]]$test
            bk_cal  <- results_bounds[[r]][[k]]$cal
            
            # if missing or wrong type, mark fail and store empty data.frames
            if (!is.data.frame(bk_test) || !is.data.frame(bk_cal)) {
                tab_static_r$ok[k] <- FALSE
                bounds_static_r[[k]] <- list(test = data.frame(), cal = data.frame())
                next
            }
            
            idx_test <- match(bk_test$id, b0_test$id)
            idx_cal  <- match(bk_cal$id,  b0_cal$id)
            
            if (anyNA(idx_test) || anyNA(idx_cal)) {
                tab_static_r$ok[k] <- FALSE
                next
            }
            
            lower_test <- b0_test$lower[idx_test]
            upper_test <- b0_test$upper[idx_test]
            lower_cal  <- b0_cal$lower[idx_cal]
            upper_cal  <- b0_cal$upper[idx_cal]
            
            bt <- bk_test
            bc <- bk_cal
            
            bt$lower <- lower_test
            bt$upper <- upper_test
            bc$lower <- lower_cal
            bc$upper <- upper_cal
            
            bt$covered_X <- as.integer(bt$lower <= bt$X & bt$X <= bt$upper)
            bc$covered_X <- as.integer(bc$lower <= bc$X & bc$X <= bc$upper)
            
            bt$covered_T <- if (!is.null(TT.name) && "T_true" %in% names(bt))
                as.integer(bt$lower <= bt$T_true & bt$T_true <= bt$upper) else rep(NA_integer_, nrow(bt))
            bc$covered_T <- if (!is.null(TT.name) && "T_true" %in% names(bc))
                as.integer(bc$lower <= bc$T_true & bc$T_true <= bc$upper) else rep(NA_integer_, nrow(bc))
            
            bounds_static_r[[k]] <- list(test = bt, cal = bc)
            
            tab_static_r$coverage_test_ipcw[k] <- if (nrow(bt) == 0 || !"w_ipcw" %in% names(bt)) NA_real_ else
                mean(bt$w_ipcw * bt$covered_X, na.rm = TRUE) / mean(bt$w_ipcw, na.rm = TRUE)
            tab_static_r$coverage_cal_ipcw[k]  <- if (nrow(bc) == 0 || !"w_ipcw" %in% names(bc)) NA_real_ else
                mean(bc$w_ipcw * bc$covered_X, na.rm = TRUE) / mean(bc$w_ipcw, na.rm = TRUE)
            
            tab_static_r$coverage_test_true[k] <- if (all(is.na(bt$covered_T))) NA_real_ else mean(bt$covered_T, na.rm = TRUE)
            tab_static_r$coverage_cal_true[k]  <- if (all(is.na(bc$covered_T))) NA_real_ else mean(bc$covered_T, na.rm = TRUE)
            
            tab_static_r$mean_len_test[k] <- if (nrow(bt) == 0) NA_real_ else mean(bt$upper - bt$lower, na.rm = TRUE)
            tab_static_r$mean_len_cal[k]  <- if (nrow(bc) == 0) NA_real_ else mean(bc$upper - bc$lower, na.rm = TRUE)
        }
        
        results_bounds_static[[r]]  <- bounds_static_r
        results_summary_static[[r]] <- tab_static_r
    }
    
    list(results_bounds_static  = results_bounds_static,
         results_summary_static = results_summary_static)
}
