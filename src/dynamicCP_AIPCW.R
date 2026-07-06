## Implementation for dynamic conformal prediction (CP) algorithm with sample splitting - general time-varying covariates setting
##  -- 2 splits: one for training the prediction model and estimate the censoring model (nuisance parameter); 
##               one for computing the miscoverage rate

.select_theta_index_from_curve <- function(theta_grid, curve,
                                          selector = c("first_crossing", "rightmost_nonnegative"),
                                          tol = 0) {
    selector <- match.arg(selector)
    if (length(theta_grid) != length(curve)) {
        stop(".select_theta_index_from_curve(): theta_grid and curve must have the same length.")
    }
    tol <- as.numeric(tol[[1L]])
    if (!is.finite(tol) || tol < 0) {
        stop(".select_theta_index_from_curve(): tol must be nonnegative.")
    }

    ok <- is.finite(curve)
    if (selector == "first_crossing" && length(curve) >= 2L) {
        for (jj in seq_len(length(curve) - 1L)) {
            if (ok[jj] && ok[jj + 1L] &&
                curve[jj] >= -tol && curve[jj + 1L] < -tol) {
                return(list(index = jj, reason = "first_downward_crossing"))
            }
        }
        idx_nonneg <- which(ok & curve >= -tol)
        if (length(idx_nonneg) > 0L) {
            return(list(index = max(idx_nonneg), reason = "no_downward_crossing_rightmost_nonnegative"))
        }
        return(list(index = NA_integer_, reason = "no_nonnegative_theta"))
    }

    idx_nonneg <- which(ok & curve >= -tol)
    if (length(idx_nonneg) > 0L) {
        return(list(index = max(idx_nonneg), reason = "rightmost_nonnegative"))
    }
    list(index = NA_integer_, reason = "no_nonnegative_theta")
}

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
#' @param model.pred.args Optional list of prediction-model-specific arguments.
#'   Used by \code{model.pred = "ddh"} to configure the R/Python bridge and
#'   neural survival backend.
#' @param model.C the censoring model
#' @param rsf_args.C A list of extra arguments passed to \code{randomForestSRC::rfsrc}
#'   when \code{model.C = "rsf"} in piecewise censoring estimation.
#' @param rsf_args.S A list of extra arguments passed to \code{randomForestSRC::rfsrc}
#'   when \code{model.S = "rsf"} for the event-time nuisance model.
#' @param rsf_args.Xi A list of extra arguments passed to \code{randomForestSRC::rfsrc}
#'   when \code{model.Xi = "rsf"} for Xi regression.
#' @param covname.pred.baseline names for the baseline covariates to be used in the prediction
#' @param covname.pred.timevarying names for the time-varying covariates to be used in the prediction
#' @param TT.name the variable name for the true event time
#' @param visit_times a vector of visit times $(t_1,...,t_K)$
#' @param localization whether localization is used for the AIPCW approach
#' @param seed seed for randomly split the data
#' @param train_ratio the proportion of data used for the training set
#' @param theta_grid a vector of theta values to be selected from
#' @param alpha the miscoverage tolerance level
#' @param pred_use_history A logical indicator for whether to use the full covariate history up to the landmark time \code{tau} in the prediction model. 
#'   If \code{FALSE} (default), the prediction model uses only the most recent covariate values at \code{tau}. 
#'   If \code{TRUE}, covariates measured at \code{visit_times[1:k_tau]} are incorporated, where \code{k_tau} is the index such that \code{visit_times[k_tau] == tau}.
#'
#' @param pred_history_mode A character string specifying how the covariate history is encoded when \code{pred_use_history = TRUE}. 
#'   \describe{
#'     \item{\code{"wide"}}{Include all available measurements up to \code{tau} as separate inputs, e.g., \code{L_t1, ..., L_t{k_tau}}.}
#'     \item{\code{"summary"}}{Use summary features of the history, including the mean, the last measurement, and a linear slope over visit index (or time index) up to \code{tau}.}
#'   }
#'
#' @param censoring_method A character string specifying how censoring survival probabilities are estimated for IPCW. 
#'   \describe{
#'     \item{\code{"global"}}{Fit a single censoring model over the full follow-up (optionally with time-varying covariates) and construct subject-level censoring survival probabilities by multiplying interval-specific contributions.}
#'     \item{\code{"piecewise"}}{Fit censoring models separately within each visit interval \code{[t_k, t_{k+1})} using the piecewise approach \code{G_k}, and compute the conditional censoring survival from \code{tau} to the observed endpoint as \code{H_{k_tau}(X)} (product over intervals).}
#'   }
#'
#' @param G_use_history A logical indicator for whether to use the full covariate history up to each interval start \code{t_k} when fitting the piecewise censoring models \code{G_k}. 
#'   If \code{FALSE} (default), \code{G_k} uses covariates measured at \code{t_k} only. 
#'   If \code{TRUE}, \code{G_k} uses covariates measured at \code{visit_times[1:k]} (encoded according to \code{G_history_mode}).
#'
#' @param G_history_mode A character string specifying how covariate history is encoded for the censoring models when \code{G_use_history = TRUE}. 
#'   \describe{
#'     \item{\code{"wide"}}{Include all measurements up to \code{t_k} as separate inputs, e.g., \code{L_t1, ..., L_tk}.}
#'     \item{\code{"summary"}}{Use summary features of the history up to \code{t_k}, including mean, last measurement, and linear slope over visit index (or time index).}
#'   }
#'
#' @param S_use_history A logical indicator for whether to use the full covariate history up to each landmark \code{t_k} when fitting the event-time models \code{S_k}. 
#'   If \code{FALSE} (default), \code{S_k} uses covariates measured at \code{t_k} only. 
#'   If \code{TRUE}, \code{S_k} uses covariates measured at \code{visit_times[1:k]} (encoded according to \code{S_history_mode}).
#'
#' @param S_history_mode A character string specifying how covariate history is encoded for the event-time models when \code{S_use_history = TRUE}. 
#'   \describe{
#'     \item{\code{"wide"}}{Include all measurements up to \code{t_k} as separate inputs, e.g., \code{L_t1, ..., L_tk}.}
#'     \item{\code{"summary"}}{Use summary features of the history up to \code{t_k}, including mean, last measurement, and linear slope over visit index (or time index).}
#'   }
#' @param S_fit_sample_mode A character string controlling how \code{S_k} is fit.
#'   \describe{
#'     \item{\code{"event_only_ipcw"}}{Current complete-case IPCW fit using only subjects with \code{\Delta=1}.}
#'     \item{\code{"all_at_risk"}}{Use all subjects at risk at \code{t_k} with ordinary right-censoring.}
#'     \item{\code{"hybrid_interval_ipcw"}}{Use all-at-risk fitting on \eqn{(t_k, t_{k+1}]} and an IPCW tail fit beyond \eqn{t_{k+1}} while keeping predictors fixed at \eqn{\bar L_k}.}
#'   }
#' @param S_small_event_threshold Nonnegative integer. If positive, and the number of
#'   events used to fit an \code{S_k} model falls below
#'   \code{max(S_small_event_threshold, S_small_event_min_events_per_coef * p)}, the
#'   requested \code{model.S} can be replaced by a simpler fallback model.
#' @param S_small_event_min_events_per_coef Minimum effective events per fitted
#'   covariate coefficient used in the adaptive \code{S_k} fallback threshold.
#' @param S_small_event_fallback_model Fallback model used when
#'   the \code{S_k} small-event rule is triggered. Supported values are
#'   \code{"cox"} and \code{"km"}.
#' @param S_small_event_metric Character string controlling how the small-event
#'   safeguard is evaluated. One of \code{"count"}, \code{"ess"}, or
#'   \code{"min_count_ess"}.
#' @param S_small_event_guard_models Character vector giving the \code{model.S}
#'   values subject to the small-event safeguard. By default this targets the
#'   more flexible survival learners.
#' @param Gtau_mode A character string specifying how \code{\tilde G(\tau|\bar L_\nu)} is handled.
#'   \describe{
#'     \item{\code{"one"}}{Set \code{\tilde G(\tau|\bar L_\nu)=1} for all subjects (default).}
#'     \item{\code{"estimated"}}{Use subject-specific \code{G(\tau|\bar L_\nu)} estimated from \code{Gfits}.}
#'     \item{\code{"tilted"}}{Use subject-specific exponentially tilted \code{\tilde G_\delta(\tau|\bar L_\nu)} based on weighted censoring mass from \code{G}.}
#'   }
#' @param Gtau_delta Numeric tilt parameter \code{\delta} used when \code{Gtau_mode="tilted"}.
#'   \code{\delta=0} recovers \code{\tilde G = G}.
#' @param localization A logical indicator for whether to use localization when computing the augmentation term
#' @param AIPCW_theta_selector Rule used to select the AIPCW theta from the
#'   AIPCW estimating curve. \code{"first_crossing"} scans theta from left to
#'   right and selects the grid point immediately before the first downward
#'   zero crossing; \code{"rightmost_nonnegative"} uses the previous behavior.
#' @param AIPCW_theta_tol Nonnegative tolerance used when checking whether the
#'   AIPCW estimating curve is nonnegative.
#' @param return_on_aipcw_fail Logical. If \code{TRUE}, return AIPCW selector diagnostics
#'   instead of stopping when no theta satisfies the AIPCW estimating equation.
dynamicCP_AIPCW_split <- function(dat, dat_test = NULL, start.name, stop.name, event.name, event.C.name, id.name, tau,
                            model.pred = c("cox", "rsf", "hal", "xgb_cox", "xgb_aft", "ddh", "ddh_py"),
                            model.pred.args = list(),
                            covname.pred.baseline = NULL, covname.pred.timevarying = NULL, # pred_tvtime = FALSE,
                            model.C = "km", covname.C.baseline = NULL, covname.C.timevarying = NULL, rsf_args.C = list(),
                            G_cox_min_censor_events = 10L,
                            G_cox_min_censor_events_per_coef = 5L,
                            G_rsf_min_censor_events = 10L,
                            G_rsf_min_censor_events_per_coef = 5L,
                            rsf_args.S = list(ntree = 1000L, nodesize = 30L, nsplit = 10L, forest = TRUE),
                            rsf_args.Xi = list(),
                            Xi_fallback_to_mean = TRUE,
                            Xi_min_n = 20L,
                            Xi_min_ess = 20,
                            Xi_min_obs_per_coef = 5,
                            Xi_min_class_n = 5L,
                            Xi_min_class_ess = 5,
                            model.S = "rsf", 
                            model.Xi = c("lm", "rsf", "xgb_reg", "xgb_multiclass", "hal_reg"),
                            TT.name = NULL, 
                            visit_times = NULL, 
                            seed = NULL, train_ratio = 0.5,
                            theta_grid, 
                            alpha = 0.1, trim.C = 0.05,
                            
                            pred_use_history = FALSE,
                            pred_history_mode = c("wide", "summary"),
                            
                            censoring_method = c("global", "piecewise"),  
                             G_use_history = FALSE,                         # history for Gk
                             G_history_mode = c("wide", "summary"),
                            
                             S_use_history = FALSE,                         # history for Sk
                             S_history_mode = c("wide", "summary"),
                             S_fit_sample_mode = c("event_only_ipcw", "all_at_risk", "hybrid_interval_ipcw"),
                             S_small_event_threshold = 10L,
                             S_small_event_min_events_per_coef = 5L,
                             S_small_event_fallback_model = "km",
                             S_small_event_metric = c("min_count_ess", "count", "ess"),
                             S_small_event_guard_models = c("cox", "rsf", "xgb_cox", "xgb_aft"),
                             xgb_nrounds = 200L,
                             xgb_use_cv = FALSE,
                             xgb_cv_nfold = 5L,
                            xgb_cv_early_stopping_rounds = 20L,
                            xgb_cv_param_grid = NULL,
                            xgb_cv_seed = NULL,
                            xgb_cv_verbose = 0,
                            
                             localization = TRUE,
                             Gtau_mode = c("one", "estimated", "tilted"),
                             Gtau_delta = 0,
                             ipcw_only = FALSE,
                             AIPCW_theta_selector = c("first_crossing", "rightmost_nonnegative"),
                             AIPCW_theta_tol = 0,
                            return_Gfits = FALSE,
                            return_on_aipcw_fail = FALSE){
    
    if (is.null(visit_times) || length(visit_times) == 0) stop("visit_times must be provided.")
    k_tau <- match(tau, visit_times)
    if (is.na(k_tau)) stop("tau must be one of visit_times.")
    
    
    # Check that the data set is regularly spaced
    if( sum(!(unique(dat[, start.name]) %in% visit_times))>0 | sum(!(unique(dat_test[, start.name]) %in% visit_times))>0){
        stop("The data sets should be regularly spaced at the time points in visit_times")
    }
    
    pred_history_mode <- match.arg(pred_history_mode)
    censoring_method  <- match.arg(censoring_method)
    G_history_mode    <- match.arg(G_history_mode)
    S_history_mode    <- match.arg(S_history_mode)
    S_fit_sample_mode <- match.arg(S_fit_sample_mode)
    model.pred        <- match.arg(model.pred)
    model.Xi          <- match.arg(model.Xi)
    Gtau_mode         <- match.arg(Gtau_mode)
    AIPCW_theta_selector <- match.arg(AIPCW_theta_selector)
    Gtau_delta        <- as.numeric(Gtau_delta)
    AIPCW_theta_tol   <- as.numeric(AIPCW_theta_tol[[1L]])
    if (!is.finite(AIPCW_theta_tol) || AIPCW_theta_tol < 0) {
        stop("AIPCW_theta_tol must be nonnegative.")
    }
    if (!is.list(rsf_args.C)) stop("rsf_args.C must be a list.")
    if (!is.list(rsf_args.S)) stop("rsf_args.S must be a list.")
    if (!is.list(rsf_args.Xi)) stop("rsf_args.Xi must be a list.")
    if (!is.list(model.pred.args)) stop("model.pred.args must be a list.")
    xgb_nrounds <- as.integer(xgb_nrounds)
    if (!is.finite(xgb_nrounds) || xgb_nrounds < 1L) {
        stop("xgb_nrounds must be >= 1.")
    }
    
    # theta grid used for C_{tau,theta}: enforce scalar domain in [0, 0.5]
    theta_grid <- as.numeric(theta_grid)
    if (length(theta_grid) == 0L || any(!is.finite(theta_grid))) {
        stop("theta_grid must be a non-empty numeric vector with finite values.")
    }
    theta_grid <- sort(unique(theta_grid))
    if (any(theta_grid < 0 | theta_grid > 0.5)) {
        stop("theta_grid must lie in [0, 0.5].")
    }
    if (length(Gtau_delta) != 1L || !is.finite(Gtau_delta)) {
        stop("Gtau_delta must be a finite scalar numeric value.")
    }
    if (Gtau_mode %in% c("estimated", "tilted") && censoring_method != "piecewise") {
        stop("Gtau_mode='", Gtau_mode, "' currently requires censoring_method='piecewise'.")
    }
    if (censoring_method == "global") {
        stop(
            "censoring_method='global' is currently not implemented in dynamicCP_AIPCW_split(). ",
            "The current AIPCW/IPCW-new pipeline depends on piecewise G_k objects (Gfits). ",
            "Please use censoring_method='piecewise'."
        )
    }
    
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
    
    # dat_tr_tau <- prepare_data_tau(dat_tr, start.name, stop.name, event.name, id.name, tau)
    # dat_cal_tau <- prepare_data_tau(dat_cal, start.name, stop.name, event.name, id.name, tau)
    # dat_test_tau <- prepare_data_tau(dat_test, start.name, stop.name, event.name, id.name, tau)
    dat_tr_tau <- make_tau_df(dat_tr, tau, visit_times, start.name, stop.name, event.name, id.name,
                              cov_tv = covname.pred.timevarying,
                              use_history = pred_use_history,
                              history_mode = pred_history_mode)
    
    dat_cal_tau <- make_tau_df(dat_cal, tau, visit_times, start.name, stop.name, event.name, id.name,
                               cov_tv = covname.pred.timevarying,
                               use_history = pred_use_history,
                               history_mode = pred_history_mode)
    
    dat_test_tau <- make_tau_df(dat_test, tau, visit_times, start.name, stop.name, event.name, id.name,
                                cov_tv = covname.pred.timevarying,
                                use_history = pred_use_history,
                                history_mode = pred_history_mode)
    
    # Train the prediction model 
    if(model.pred == "cox"){
        
        cov_base <- covname.pred.baseline
        cov_tv   <- covname.pred.timevarying
        
        if (!pred_use_history) {
            cov_pred_use <- c(cov_base, cov_tv)
        } else {
            cov_pred_use <- c(cov_base, expand_hist_covars(cov_tv, k_tau, pred_history_mode))
        }
        
        cov_pred_use <- keep_nonconstant_covars(dat_tr_tau, cov_pred_use)
        
        
        rhs <- if (length(cov_pred_use) == 0) "1" else paste(cov_pred_use, collapse = "+")
        formula.pred <- as.formula(paste0("Surv(", start.name, ",", stop.name, ",", event.name, ") ~ ", rhs))
        
        fit.pred <- coxph(formula.pred, data = dat_tr_tau)
        
        # prediction on the calibration data and test data
        # (se.fit/conf.type disabled: only $time/$surv are used; $surv is identical
        # and survfit() is several times faster without the SE/CI computation)
        sf_tr <- survfit(fit.pred, newdata = dat_tr_tau, se.fit = FALSE, conf.type = "none")
        pred_surv_tr <- t(sf_tr$surv)
        if (!is.matrix(pred_surv_tr)) {
            pred_surv_tr <- matrix(pred_surv_tr, nrow = nrow(dat_tr_tau), byrow = TRUE)
        }

        sf_cal <- survfit(fit.pred, newdata = dat_cal_tau, se.fit = FALSE, conf.type = "none")
        pred_time <- sf_cal$time
        pred_surv_cal <- t(sf_cal$surv)  # row - individual, col - times points

        sf_test <- survfit(fit.pred, newdata = dat_test_tau, se.fit = FALSE, conf.type = "none")
        pred_surv_test <- t(sf_test$surv)  # row - individual, col - times points
        
    } else if (model.pred == "rsf") {
        
        cov_base <- covname.pred.baseline
        cov_tv   <- covname.pred.timevarying
        
        if (!pred_use_history) {
            cov_pred_use <- c(cov_base, cov_tv)
        } else {
            cov_pred_use <- c(cov_base, expand_hist_covars(cov_tv, k_tau, pred_history_mode))
        }
        
        cov_pred_use <- keep_nonconstant_covars(dat_tr_tau, cov_pred_use)
        
        
        make_rsf_df <- function(dat_tau) {
            # remaining time after tau (relative time scale for landmarking)
            dat_tau$.time_rsf <- dat_tau[[stop.name]] - tau
            dat_tau
        }
        
        dat_tr_rsf <- make_rsf_df(dat_tr_tau)
        dat_cal_rsf  <- make_rsf_df(dat_cal_tau)
        dat_test_rsf <- make_rsf_df(dat_test_tau)
        
        if (length(cov_pred_use) == 0) {
            # KM on remaining time after tau, then shift time grid back to original scale
            km_fit <- survfit(Surv(.time_rsf, dat_tr_rsf[[event.name]]) ~ 1, data = dat_tr_rsf)
            
            pred_time_rel <- km_fit$time                 # time since tau
            pred_time     <- tau + pred_time_rel         # ORIGINAL time scale
            S_km <- km_fit$surv
            
            pred_surv_tr   <- matrix(rep(S_km, nrow(dat_tr_rsf)),  nrow = nrow(dat_tr_rsf),  byrow = TRUE)
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
            
            pr_tr <- predict(
                fit.pred,
                newdata = dat_tr_rsf,
                time.interest = pred_time_rel
            )
            pred_surv_tr <- pr_tr$survival
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
        
        if (isTRUE(pred_use_history)) {
            warning(
                "model.pred='ddh' uses the full long-format covariate trajectory ",
                "up to tau through the Python backend; pred_history_mode is ignored."
            )
        }
        
        cov_pred_use <- c(covname.pred.baseline, covname.pred.timevarying)
        cov_pred_use <- keep_nonconstant_covars(dat_tr_tau, cov_pred_use)
        
        ddh_pred <- predict_fn(
            dat_tr = dat_tr,
            dat_cal = dat_cal,
            dat_test = dat_test,
            dat_tr_tau = dat_tr_tau,
            dat_cal_tau = dat_cal_tau,
            dat_test_tau = dat_test_tau,
            start.name = start.name,
            stop.name = stop.name,
            event.name = event.name,
            id.name = id.name,
            tau = tau,
            cov_pred_use = cov_pred_use,
            model.pred.args = model.pred.args
        )
        
        pred_time <- ddh_pred$pred_time
        pred_surv_tr <- ddh_pred$pred_surv_tr
        pred_surv_cal <- ddh_pred$pred_surv_cal
        pred_surv_test <- ddh_pred$pred_surv_test
    } else if (model.pred == "hal") {
        
        cov_base <- covname.pred.baseline
        cov_tv   <- covname.pred.timevarying
        
        if (!pred_use_history) {
            cov_pred_use <- c(cov_base, cov_tv)
        } else {
            cov_pred_use <- c(cov_base, expand_hist_covars(cov_tv, k_tau, pred_history_mode))
        }
        cov_pred_use <- keep_nonconstant_covars(dat_tr_tau, cov_pred_use)
        
        if (length(cov_pred_use) == 0) {
            # With no usable covariates, HAL Cox reduces to a marginal survival model; use KM fallback.
            km_fit <- survfit(Surv(dat_tr_tau[[stop.name]] - tau, dat_tr_tau[[event.name]]) ~ 1, data = dat_tr_tau)
            pred_time_rel <- km_fit$time
            pred_time     <- tau + pred_time_rel
            S_km <- km_fit$surv
            
            pred_surv_tr   <- matrix(rep(S_km, nrow(dat_tr_tau)),  nrow = nrow(dat_tr_tau),  byrow = TRUE)
            pred_surv_cal  <- matrix(rep(S_km, nrow(dat_cal_tau)),  nrow = nrow(dat_cal_tau),  byrow = TRUE)
            pred_surv_test <- matrix(rep(S_km, nrow(dat_test_tau)), nrow = nrow(dat_test_tau), byrow = TRUE)
            
        } else {
            Y_hal <- survival::Surv(dat_tr_tau[[stop.name]] - tau, dat_tr_tau[[event.name]])
            hal_bundle_pred <- .fit_hal_cox_bundle(
                dat = dat_tr_tau,
                cov_use = cov_pred_use,
                Y_surv = Y_hal
            )
            
            pred_time_rel <- hal_bundle_pred$bh_time
            pred_time     <- tau + pred_time_rel
            
            build_hal_surv_mat <- function(newdata_tau) {
                lp_new <- .predict_hal_cox_lp(hal_bundle_pred, newdata_tau)
                n_new <- length(lp_new)
                m_grid <- length(pred_time_rel)
                out <- matrix(NA_real_, nrow = n_new, ncol = m_grid)
                for (ii in seq_len(n_new)) {
                    out[ii, ] <- .predict_surv_hal_cox(
                        hal_bundle = hal_bundle_pred,
                        lp = rep(lp_new[ii], m_grid),
                        t_stop = pred_time_rel,
                        t_start = rep(0, m_grid)
                    )
                }
                out
            }
            
            pred_surv_tr   <- build_hal_surv_mat(dat_tr_tau)
            pred_surv_cal  <- build_hal_surv_mat(dat_cal_tau)
            pred_surv_test <- build_hal_surv_mat(dat_test_tau)
        }
        
    } else if (model.pred == "grf") {
        
        if (!requireNamespace("grf", quietly = TRUE)) {
            stop("Package 'grf' is required for model.pred='grf'.")
        }
        
        cov_base <- covname.pred.baseline
        cov_tv   <- covname.pred.timevarying
        
        if (!pred_use_history) {
            cov_pred_use <- c(cov_base, cov_tv)
        } else {
            cov_pred_use <- c(cov_base, expand_hist_covars(cov_tv, k_tau, pred_history_mode))
        }
        cov_pred_use <- keep_nonconstant_covars(dat_tr_tau, cov_pred_use)
        
        if (length(cov_pred_use) == 0) {
            km_fit <- survfit(Surv(dat_tr_tau[[stop.name]] - tau, dat_tr_tau[[event.name]]) ~ 1, data = dat_tr_tau)
            pred_time_rel <- km_fit$time
            pred_time     <- tau + pred_time_rel
            S_km <- km_fit$surv
            
            pred_surv_tr   <- matrix(rep(S_km, nrow(dat_tr_tau)),  nrow = nrow(dat_tr_tau),  byrow = TRUE)
            pred_surv_cal  <- matrix(rep(S_km, nrow(dat_cal_tau)),  nrow = nrow(dat_cal_tau),  byrow = TRUE)
            pred_surv_test <- matrix(rep(S_km, nrow(dat_test_tau)), nrow = nrow(dat_test_tau), byrow = TRUE)
        } else {
            dm_tr <- .fit_design_matrix(dat_tr_tau, cov_use = cov_pred_use)
            if (is.null(dm_tr$X) || nrow(dm_tr$X) == 0) {
                stop("model.pred='grf': no complete-case rows left after design matrix construction.")
            }
            dat_tr_ml <- dat_tr_tau[dm_tr$keep, , drop = FALSE]
            
            dm_cal <- .predict_design_matrix(dat_cal_tau, cov_use = cov_pred_use, design_info = dm_tr$design_info)
            dm_test <- .predict_design_matrix(dat_test_tau, cov_use = cov_pred_use, design_info = dm_tr$design_info)
            if (any(!dm_cal$keep) || any(!dm_test$keep)) {
                stop("model.pred='grf': missing values in covariates are not supported for prediction.")
            }
            
            fit.pred <- grf::survival_forest(
                X = dm_tr$X,
                Y = as.numeric(dat_tr_ml[[stop.name]] - tau),
                D = as.numeric(dat_tr_ml[[event.name]])
            )
            
            pr_cal <- predict(fit.pred, newdata = dm_cal$X)
            pred_time_rel <- pr_cal$failure.times
            pred_time <- tau + pred_time_rel
            pred_surv_cal <- pr_cal$predictions
            
            pr_test <- predict(fit.pred, newdata = dm_test$X, failure.times = pred_time_rel)
            pred_surv_test <- pr_test$predictions
            
            pr_tr <- predict(fit.pred, newdata = dm_tr$X, failure.times = pred_time_rel)
            pred_surv_tr <- pr_tr$predictions
        }
        
    } else if (model.pred %in% c("xgb_cox", "xgb_aft")) {
        
        if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("Package 'xgboost' is required for model.pred='", model.pred, "'.")
        }
        
        cov_base <- covname.pred.baseline
        cov_tv   <- covname.pred.timevarying
        
        if (!pred_use_history) {
            cov_pred_use <- c(cov_base, cov_tv)
        } else {
            cov_pred_use <- c(cov_base, expand_hist_covars(cov_tv, k_tau, pred_history_mode))
        }
        cov_pred_use <- keep_nonconstant_covars(dat_tr_tau, cov_pred_use)
        
        if (length(cov_pred_use) == 0) {
            km_fit <- survfit(Surv(dat_tr_tau[[stop.name]] - tau, dat_tr_tau[[event.name]]) ~ 1, data = dat_tr_tau)
            pred_time_rel <- km_fit$time
            pred_time     <- tau + pred_time_rel
            S_km <- km_fit$surv
            
            pred_surv_tr   <- matrix(rep(S_km, nrow(dat_tr_tau)),  nrow = nrow(dat_tr_tau),  byrow = TRUE)
            pred_surv_cal  <- matrix(rep(S_km, nrow(dat_cal_tau)),  nrow = nrow(dat_cal_tau),  byrow = TRUE)
            pred_surv_test <- matrix(rep(S_km, nrow(dat_test_tau)), nrow = nrow(dat_test_tau), byrow = TRUE)
        } else {
            dm_tr <- .fit_design_matrix(dat_tr_tau, cov_use = cov_pred_use)
            if (is.null(dm_tr$X) || nrow(dm_tr$X) == 0) {
                stop("model.pred='", model.pred, "': no complete-case rows left after design matrix construction.")
            }
            dat_tr_ml <- dat_tr_tau[dm_tr$keep, , drop = FALSE]
            
            dm_cal <- .predict_design_matrix(dat_cal_tau, cov_use = cov_pred_use, design_info = dm_tr$design_info)
            dm_test <- .predict_design_matrix(dat_test_tau, cov_use = cov_pred_use, design_info = dm_tr$design_info)
            if (any(!dm_cal$keep) || any(!dm_test$keep)) {
                stop("model.pred='", model.pred, "': missing values in covariates are not supported for prediction.")
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
            
            if (model.pred == "xgb_cox") {
                time_rel_tr <- as.numeric(dat_tr_ml[[stop.name]] - tau)
                event_tr <- as.integer(dat_tr_ml[[event.name]])
                # xgboost Cox uses sign(label) for censoring indicator.
                y_cox <- ifelse(event_tr == 1L, time_rel_tr, -time_rel_tr)
                dtrain <- xgboost::xgb.DMatrix(data = dm_tr$X, label = y_cox)
                
                params <- c(list(objective = "survival:cox", eval_metric = "cox-nloglik"), xgb_params_common)
                xgb_fit <- .xgb_train_optional_cv(
                    dtrain = dtrain,
                    params = params,
                    nrounds = xgb_nrounds,
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
                
                lp_fit <- as.numeric(stats::predict(fit_xgb, newdata = dm_tr$X, outputmargin = TRUE))
                lp_fit[!is.finite(lp_fit)] <- 0
                lp_fit <- pmax(pmin(lp_fit, 30), -30)
                
                fit_anchor <- survival::coxph(
                    survival::Surv(time_rel_tr, event_tr) ~ offset(lp_fit),
                    ties = "breslow",
                    x = TRUE
                )
                bh <- survival::basehaz(fit_anchor, centered = FALSE)
                pred_time_rel <- bh$time
                pred_time <- tau + pred_time_rel
                S0_grid <- exp(-bh$hazard)
                
                lp_tr <- as.numeric(stats::predict(fit_xgb, newdata = dm_tr$X, outputmargin = TRUE))
                lp_tr[!is.finite(lp_tr)] <- 0
                lp_tr <- pmax(pmin(lp_tr, 30), -30)
                pred_surv_tr <- outer(exp(lp_tr), S0_grid, FUN = function(r, s0) s0^r)
                
                lp_cal <- as.numeric(stats::predict(fit_xgb, newdata = dm_cal$X, outputmargin = TRUE))
                lp_cal[!is.finite(lp_cal)] <- 0
                lp_cal <- pmax(pmin(lp_cal, 30), -30)
                pred_surv_cal <- outer(exp(lp_cal), S0_grid, FUN = function(r, s0) s0^r)
                
                lp_test <- as.numeric(stats::predict(fit_xgb, newdata = dm_test$X, outputmargin = TRUE))
                lp_test[!is.finite(lp_test)] <- 0
                lp_test <- pmax(pmin(lp_test, 30), -30)
                pred_surv_test <- outer(exp(lp_test), S0_grid, FUN = function(r, s0) s0^r)
                
            } else {
                time_rel_tr <- as.numeric(dat_tr_ml[[stop.name]] - tau)
                event_tr <- as.integer(dat_tr_ml[[event.name]])
                
                dtrain <- xgboost::xgb.DMatrix(data = dm_tr$X)
                xgboost::setinfo(dtrain, "label_lower_bound", time_rel_tr)
                xgboost::setinfo(
                    dtrain,
                    "label_upper_bound",
                    ifelse(event_tr == 1L, time_rel_tr, Inf)
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
                    nrounds = xgb_nrounds,
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
                
                # Use observed relative times as evaluation grid for quantile inversion.
                pred_time_rel <- sort(unique(time_rel_tr[is.finite(time_rel_tr) & time_rel_tr >= 0]))
                pred_time <- tau + pred_time_rel
                
                build_aft_surv <- function(mu_vec) {
                    tiny <- .Machine$double.eps
                    t_pos <- pmax(pred_time_rel, tiny)
                    z <- outer(mu_vec, log(t_pos), FUN = function(mu, lt) (lt - mu) / aft_scale)
                    Fz <- switch(
                        aft_dist,
                        normal = stats::pnorm(z),
                        logistic = stats::plogis(z),
                        extreme = exp(-exp(-z))
                    )
                    S <- 1 - Fz
                    S[, pred_time_rel <= 0] <- 1
                    S[S < 0] <- 0
                    S[S > 1] <- 1
                    S
                }
                
                mu_tr <- as.numeric(stats::predict(fit_xgb, newdata = dm_tr$X))
                mu_cal <- as.numeric(stats::predict(fit_xgb, newdata = dm_cal$X))
                mu_test <- as.numeric(stats::predict(fit_xgb, newdata = dm_test$X))
                
                pred_surv_tr <- build_aft_surv(mu_tr)
                pred_surv_cal <- build_aft_surv(mu_cal)
                pred_surv_test <- build_aft_surv(mu_test)
            }
        }
        
    } else {
        stop("This prediction model is not implemented in the function")
    }
    
    
    # Fit the censoring model on training data, compute the IPCW weights for the calibration data and test data
    denom_col <- NULL  # name of censoring survival column used for IPCW
    Gfits <- NULL
    
    if (censoring_method == "global") {
        
        if (model.C == "km") {
            
            dat_tr_0 <- prepare_data_tau(dat_tr, start.name, stop.name, event.C.name, id.name, tau = 0)
            formula_C <- as.formula(paste0("Surv(", stop.name, ", ", event.C.name, ") ~ 1"))
            kmfit_C <- survfit(formula_C, data = dat_tr_0)
            km_C <- summary(kmfit_C)
            Sc_stepfun <- stepfun(km_C$time, c(1, km_C$surv))
            
            dat_cal$surv.C.subj <- ave(dat_cal[, stop.name], dat_cal[, id.name],
                                       FUN = function(x) Sc_stepfun(max(x)))
            dat_test$surv.C.subj <- ave(dat_test[, stop.name], dat_test[, id.name],
                                        FUN = function(x) Sc_stepfun(max(x)))
            
        } else if (model.C == "cox") {
            
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
            
        } else if (model.C == "hal") {
            
            cov_C_use <- c(covname.C.baseline, covname.C.timevarying)
            cov_C_use <- keep_nonconstant_covars(dat_tr, cov_C_use)
            
            if (length(cov_C_use) == 0) {
                # HAL with no usable covariates collapses to a marginal model; use KM fallback.
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
                Y_surv_C <- survival::Surv(dat_tr[[start.name]], dat_tr[[stop.name]], dat_tr[[event.C.name]])
                hal_bundle_C <- .fit_hal_cox_bundle(
                    dat = dat_tr,
                    cov_use = cov_C_use,
                    Y_surv = Y_surv_C
                )
                
                lp_cal_C <- .predict_hal_cox_lp(hal_bundle_C, dat_cal)
                lp_test_C <- .predict_hal_cox_lp(hal_bundle_C, dat_test)
                
                dat_cal$surv.C <- .predict_surv_hal_cox(
                    hal_bundle = hal_bundle_C,
                    lp = lp_cal_C,
                    t_start = dat_cal[[start.name]],
                    t_stop = dat_cal[[stop.name]]
                )
                dat_test$surv.C <- .predict_surv_hal_cox(
                    hal_bundle = hal_bundle_C,
                    lp = lp_test_C,
                    t_start = dat_test[[start.name]],
                    t_stop = dat_test[[stop.name]]
                )
                
                dat_cal$surv.C.subj <- ave(dat_cal[, "surv.C"], dat_cal[, id.name], 
                                           FUN = function(x) prod(x))
                dat_test$surv.C.subj <- ave(dat_test[, "surv.C"], dat_test[, id.name], 
                                            FUN = function(x) prod(x))
            }
            
        } else {
            stop("For censoring_method='global', model.C must be 'km', 'cox', or 'hal'.")
        }
        
        # attach subject-level surv.C.subj onto landmark rows at tau
        survC_cal_subj <- tapply(dat_cal$surv.C.subj, dat_cal[[id.name]], function(x) x[1])
        survC_test_subj <- tapply(dat_test$surv.C.subj, dat_test[[id.name]], function(x) x[1])
        
        dat_cal_tau$surv.C.subj <- as.numeric(survC_cal_subj[as.character(dat_cal_tau[[id.name]])])
        dat_test_tau$surv.C.subj <- as.numeric(survC_test_subj[as.character(dat_test_tau[[id.name]])])
        
        denom_col <- "surv.C.subj"
        
        
    } else if (censoring_method == "piecewise") {
        
        # We use H_1(X), which equals G(0 -> X), provided visit_times[1] = 0
        if (length(visit_times) == 0 || abs(visit_times[1] - 0) > 1e-12) {
            stop("For censoring_method='piecewise' with IPCW denominator G(0->X), require visit_times[1] = 0.")
        }
        
        # Fit piecewise censoring models G_k on training data
        Gfits <- fit_Gk_list(
            dat_tr, visit_times,
            start.name, stop.name, event.name, event.C.name, id.name,
            model.C = model.C,  # "km" / "cox" / "rsf" / "hal"
            covname.C.baseline = covname.C.baseline,
            covname.C.timevarying = covname.C.timevarying,
            use_history = G_use_history,
            history_mode = G_history_mode,
            cox_min_censor_events = G_cox_min_censor_events,
            cox_min_censor_events_per_coef = G_cox_min_censor_events_per_coef,
            rsf_min_censor_events = G_rsf_min_censor_events,
            rsf_min_censor_events_per_coef = G_rsf_min_censor_events_per_coef,
            rsf_args = rsf_args.C,
            xgb_nrounds = xgb_nrounds,
            xgb_use_cv = xgb_use_cv,
            xgb_cv_nfold = xgb_cv_nfold,
            xgb_cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
            xgb_cv_param_grid = xgb_cv_param_grid,
            xgb_cv_seed = xgb_cv_seed,
            xgb_cv_verbose = xgb_cv_verbose
        )
        
        # Compute full censoring survival G(0->X) on calibration and test sets
        # via H_1(X) (since visit_times[1] = 0)
        G0X_cal_df <- compute_Hk_at_tk(
            dat_long   = dat_cal,
            visit_times = visit_times,
            k          = 1,
            start.name = start.name,
            stop.name  = stop.name,
            event.name = event.name,
            event.C.name = event.C.name,
            id.name    = id.name,
            Gfits      = Gfits,
            trim.C     = 1e-7
        )
        
        G0X_test_df <- compute_Hk_at_tk(
            dat_long   = dat_test,
            visit_times = visit_times,
            k          = 1,
            start.name = start.name,
            stop.name  = stop.name,
            event.name = event.name,
            event.C.name = event.C.name,
            id.name    = id.name,
            Gfits      = Gfits,
            trim.C     = 1e-7
        )
        
        # Attach subject-level G(0->X) to tau-landmark rows
        # compute_Hk_at_tk() returns columns: id, X, Hk
        idx_cal_match  <- match(as.character(dat_cal_tau[[id.name]]),  as.character(G0X_cal_df$id))
        idx_test_match <- match(as.character(dat_test_tau[[id.name]]), as.character(G0X_test_df$id))
        
        dat_cal_tau$surv.C.subj  <- G0X_cal_df$Hk[idx_cal_match]
        dat_test_tau$surv.C.subj <- G0X_test_df$Hk[idx_test_match]
        
        # Sanity checks
        if (anyNA(dat_cal_tau$surv.C.subj)) {
            miss_ids <- unique(dat_cal_tau[[id.name]][is.na(dat_cal_tau$surv.C.subj)])
            stop("piecewise censoring: failed to compute G(0->X) for some calibration ids. Example ids: ",
                 paste(utils::head(miss_ids, 10), collapse = ", "))
        }
        if (anyNA(dat_test_tau$surv.C.subj)) {
            miss_ids <- unique(dat_test_tau[[id.name]][is.na(dat_test_tau$surv.C.subj)])
            stop("piecewise censoring: failed to compute G(0->X) for some test ids. Example ids: ",
                 paste(utils::head(miss_ids, 10), collapse = ", "))
        }
        
        # Final IPCW denominator trimming (user-specified)
        dat_cal_tau$surv.C.subj  <- pmax(dat_cal_tau$surv.C.subj,  trim.C)
        dat_test_tau$surv.C.subj <- pmax(dat_test_tau$surv.C.subj, trim.C)
        
        # # Optional: also store on long data (useful for later estimating-function code/debugging)
        # cal_map  <- stats::setNames(G0X_cal_df$Hk,  as.character(G0X_cal_df$id))
        # test_map <- stats::setNames(G0X_test_df$Hk, as.character(G0X_test_df$id))
        # 
        # dat_cal$surv.C.subj  <- unname(cal_map[as.character(dat_cal[[id.name]])])
        # dat_test$surv.C.subj <- unname(test_map[as.character(dat_test[[id.name]])])
        
        denom_col <- "surv.C.subj"
        
    } else {
        stop(paste0("Unknown censoring_method: ", censoring_method))
    }

    if (isTRUE(ipcw_only)) {
        if (!identical(Gtau_mode, "one")) {
            stop("ipcw_only=TRUE currently supports only Gtau_mode='one'.")
        }

        X_cal <- dat_cal_tau[[stop.name]]
        Delta_cal <- dat_cal_tau[[event.name]]
        Gtau_cal <- rep(1.0, length(X_cal))

        hat_W_IPCW <- rep(NA_real_, length(theta_grid))
        for (jj in seq_along(theta_grid)) {
            theta_j <- theta_grid[jj]

            lo_j <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = theta_j)
            hi_j <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = 1 - theta_j)

            inside_j <- as.numeric(X_cal > lo_j & X_cal < hi_j)
            U_vec <- (Delta_cal / pmax(dat_cal_tau[[denom_col]], trim.C)) *
                Gtau_cal *
                as.numeric(X_cal > tau) *
                (inside_j - (1 - alpha))

            hat_W_IPCW[jj] <- sum(U_vec, na.rm = TRUE)
        }

        idx_nonneg <- which(hat_W_IPCW >= 0)
        if (length(idx_nonneg) == 0L) {
            stop("No theta in theta_grid satisfies sum_i U_IPCW,i(theta) >= 0.")
        }
        theta_index_lower_IPCW <- max(idx_nonneg)
        theta_hat_IPCW <- theta_grid[theta_index_lower_IPCW]

        lower_IPCW_cal <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = theta_hat_IPCW)
        lower_IPCW_test <- surv_quantile_vec_from_pred(pred_time, pred_surv_test, beta = theta_hat_IPCW)
        upper_IPCW_cal <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = 1 - theta_hat_IPCW)
        upper_IPCW_test <- surv_quantile_vec_from_pred(pred_time, pred_surv_test, beta = 1 - theta_hat_IPCW)

        lower_model_cal <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = alpha / 2)
        lower_model_test <- surv_quantile_vec_from_pred(pred_time, pred_surv_test, beta = alpha / 2)
        upper_model_cal <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = 1 - alpha / 2)
        upper_model_test <- surv_quantile_vec_from_pred(pred_time, pred_surv_test, beta = 1 - alpha / 2)

        lower_IPCW_DR_cal <- pmin(lower_IPCW_cal, lower_model_cal)
        lower_IPCW_DR_test <- pmin(lower_IPCW_test, lower_model_test)
        upper_IPCW_DR_cal <- pmax(upper_IPCW_cal, upper_model_cal)
        upper_IPCW_DR_test <- pmax(upper_IPCW_test, upper_model_test)

        w_cal <- dat_cal_tau[[event.name]] / pmax(dat_cal_tau[[denom_col]], trim.C)
        w_test <- dat_test_tau[[event.name]] / pmax(dat_test_tau[[denom_col]], trim.C)

        coverage_IPCW_cal <- mean(w_cal * (lower_IPCW_cal <= dat_cal_tau[[stop.name]] &
                                               dat_cal_tau[[stop.name]] <= upper_IPCW_cal)) / mean(w_cal)
        coverage_IPCW_test <- mean(w_test * (lower_IPCW_test <= dat_test_tau[[stop.name]] &
                                                 dat_test_tau[[stop.name]] <= upper_IPCW_test)) / mean(w_test)

        if (!is.null(TT.name)) {
            coverage_IPCW_cal_T <- mean(lower_IPCW_cal <= dat_cal_tau[[TT.name]] &
                                            dat_cal_tau[[TT.name]] <= upper_IPCW_cal)
            coverage_IPCW_test_T <- mean(lower_IPCW_test <= dat_test_tau[[TT.name]] &
                                             dat_test_tau[[TT.name]] <= upper_IPCW_test)
        } else {
            coverage_IPCW_cal_T <- NA_real_
            coverage_IPCW_test_T <- NA_real_
        }

        return(list(
            lower_IPCW_cal = lower_IPCW_cal,
            upper_IPCW_cal = upper_IPCW_cal,
            lower_IPCW_test = lower_IPCW_test,
            upper_IPCW_test = upper_IPCW_test,
            lower_model_cal = lower_model_cal,
            upper_model_cal = upper_model_cal,
            lower_model_test = lower_model_test,
            upper_model_test = upper_model_test,
            lower_IPCW_DR_cal = lower_IPCW_DR_cal,
            upper_IPCW_DR_cal = upper_IPCW_DR_cal,
            lower_IPCW_DR_test = lower_IPCW_DR_test,
            upper_IPCW_DR_test = upper_IPCW_DR_test,
            coverage_IPCW_cal = coverage_IPCW_cal,
            coverage_IPCW_test = coverage_IPCW_test,
            coverage_IPCW_cal_T = coverage_IPCW_cal_T,
            coverage_IPCW_test_T = coverage_IPCW_test_T,
            coverage_cal = coverage_IPCW_cal,
            coverage_test = coverage_IPCW_test,
            coverage_cal_T = coverage_IPCW_cal_T,
            coverage_test_T = coverage_IPCW_test_T,
            theta_hat_IPCW = theta_hat_IPCW,
            theta_index_IPCW = theta_index_lower_IPCW,
            hat_W_IPCW = hat_W_IPCW,
            theta_grid_used = theta_grid,
            id_tr = id_tr,
            id_cal_tau = dat_cal_tau[[id.name]],
            id_test_tau = dat_test_tau[[id.name]],
            survC_cal = dat_cal_tau[[denom_col]],
            survC_test = dat_test_tau[[denom_col]],
            Gfits = if (isTRUE(return_Gfits)) Gfits else NULL,
            ipcw_only = TRUE
        ))
    }
    
    ## AIPCW scaffolding: build integration grid u_j
    ## u_j = union(all observed censoring times in calibration, jump times of G)
    summ_cal <- subject_summary(dat_cal, id.name, stop.name, event.name, event.C.name)
    censor_times_cal <- as.numeric(summ_cal$X[summ_cal$DeltaC == 1])
    censor_times_cal <- sort(unique(censor_times_cal[is.finite(censor_times_cal)]))
    
    G_jump_times <- get_G_jump_times_from_Gfits(Gfits)
    G_jump_times <- sort(unique(as.numeric(G_jump_times[is.finite(G_jump_times)])))
    
    time_vec_AIPCW <- sort(unique(c(censor_times_cal, G_jump_times)))
    
    ## AIPCW scaffolding: all-subject calibration indexing (not only X > tau)
    ids_cal_all <- as.character(summ_cal[[id.name]])
    X_cal_all <- as.numeric(summ_cal$X)
    Delta_cal_all <- as.integer(summ_cal$Delta)
    DeltaC_cal_all <- as.integer(summ_cal$DeltaC)
    n_cal_all <- length(ids_cal_all)
    
    # map tau-risk subset rows into all-subject index (for later linking with pred_surv_cal rows)
    idx_tau_in_all <- match(as.character(dat_cal_tau[[id.name]]), ids_cal_all)
    if (anyNA(idx_tau_in_all)) {
        stop("Failed to align dat_cal_tau ids to all calibration ids.")
    }
    
    # Discrete integral increment terms over all calibration subjects:
    # I(X=u_j, DeltaC=1)/G(u_j) + I(X>=u_j) * dG(u_j)/G(u_j)^2
    m_u <- length(time_vec_AIPCW)
    G_cal_u_mat_all <- matrix(numeric(0), nrow = n_cal_all, ncol = 0)
    dG_cal_u_mat_all <- matrix(numeric(0), nrow = n_cal_all, ncol = 0)
    aipcw_increment_all <- matrix(numeric(0), nrow = n_cal_all, ncol = 0)
    
    if (m_u > 0L) {
        G_cal_obj <- predict_G_from_Gk_matrix(
            dat_new = dat_cal,
            time_vec = time_vec_AIPCW,
            start.name = start.name, stop.name = stop.name, event.name = event.name,
            event.C.name = event.C.name, id.name = id.name,
            Gfits = Gfits,
            trim.G = NULL,
            return_subject_info = TRUE
        )
        
        idx_G_all <- match(ids_cal_all, as.character(G_cal_obj$ids))
        if (anyNA(idx_G_all)) {
            stop("Failed to align G(u_j) rows to all calibration ids.")
        }
        G_cal_u_mat_all <- G_cal_obj$G[idx_G_all, , drop = FALSE]
        G_denom_all <- pmax(G_cal_u_mat_all, trim.C)
        # Left-limit denominator G(u-) for jump-based increment term
        G_prev_u_mat_all <- G_cal_u_mat_all
        G_prev_u_mat_all[, 1] <- 1
        if (m_u > 1L) {
            G_prev_u_mat_all[, 2:m_u] <- G_cal_u_mat_all[, 1:(m_u - 1L), drop = FALSE]
        }
        G_denom_prev_all <- pmax(G_prev_u_mat_all, trim.C)
        
        dG_cal_u_mat_all <- G_cal_u_mat_all
        dG_cal_u_mat_all[, 1] <- G_cal_u_mat_all[, 1] - 1
        if (m_u > 1L) {
            dG_cal_u_mat_all[, 2:m_u] <- G_cal_u_mat_all[, 2:m_u, drop = FALSE] -
                G_cal_u_mat_all[, 1:(m_u - 1L), drop = FALSE]
        }
        
        I_X_eq_u_cens_all <- outer(X_cal_all, time_vec_AIPCW, FUN = "==") &
            matrix(DeltaC_cal_all == 1L, nrow = n_cal_all, ncol = m_u)
        I_X_ge_u_all <- outer(X_cal_all, time_vec_AIPCW, FUN = ">=")
        
        term1_all <- I_X_eq_u_cens_all / G_denom_all
        term2_all <- I_X_ge_u_all * (dG_cal_u_mat_all / (G_denom_prev_all * G_denom_all))
        
        term1_all[!is.finite(term1_all)] <- 0
        term2_all[!is.finite(term2_all)] <- 0
        aipcw_increment_all <- term1_all + term2_all
    }
    
    # dat_cal_tau <- prepare_data_tau(dat_cal, start.name, stop.name, event.name, id.name, tau)
    # dat_test_tau <- prepare_data_tau(dat_test, start.name, stop.name, event.name, id.name, tau)
    
    
    ## Fit the event time model for S_k on training data, 
    
    if(censoring_method == "piecewise"){
        # compute H_k(X) on training (for fitting S_k)
        Hk_list_tr <- compute_Hk_list(
            dat_tr, visit_times,
            start.name, stop.name, event.name, event.C.name, id.name,
            Gfits, trim.C = trim.C
        )
        
        # fit S_k models on training
        Sfits <- fit_Sk_list(
            dat_tr, visit_times,
            start.name, stop.name, event.name, event.C.name, id.name,
            covname.S.baseline = covname.pred.baseline,     # or separate names
            covname.S.timevarying = covname.pred.timevarying,
            Hk_list_tr = Hk_list_tr,
            model.S = model.S,
            trim.H = trim.C,
            fit_sample_mode = S_fit_sample_mode,
            small_event_threshold = S_small_event_threshold,
            small_event_min_events_per_coef = S_small_event_min_events_per_coef,
            small_event_fallback_model = S_small_event_fallback_model,
            small_event_metric = S_small_event_metric,
            small_event_guard_models = S_small_event_guard_models,
            rsf_args = rsf_args.S,
            xgb_nrounds = xgb_nrounds,
            xgb_use_cv = xgb_use_cv,
            xgb_cv_nfold = xgb_cv_nfold,
            xgb_cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
            xgb_cv_param_grid = xgb_cv_param_grid,
            xgb_cv_seed = xgb_cv_seed,
            xgb_cv_verbose = xgb_cv_verbose
        )
        
    } else if (censoring_method == "global"){
        
        stop("censoring_method == 'global' is not implemented yet.")
        
    } else {
        stop(paste0("Unknown censoring_method: ", censoring_method))
    }
    
    
    
    
    
    ## On the calibration data -----------------------------------------------------
    
    # Select theta_hat_IPCW = max{theta in Theta: sum_i U_IPCW,i(theta) >= 0}
    summ_tr <- subject_summary(dat_tr, id.name, stop.name, event.name, event.C.name)
    ids_tr_all <- as.character(summ_tr[[id.name]])
    idx_tr_tau_in_all <- match(as.character(dat_tr_tau[[id.name]]), ids_tr_all)
    if (anyNA(idx_tr_tau_in_all)) {
        stop("Failed to align dat_tr_tau ids to all training ids.")
    }
    
    n_cal_tau <- nrow(dat_cal_tau)
    n_test_tau <- nrow(dat_test_tau)
    
    # c_max from data support used by fitted G: maximum jump time
    cmax_censor_tr <- if (length(G_jump_times) > 0L) max(G_jump_times) else NA_real_
    if (Gtau_mode == "tilted" && !is.finite(cmax_censor_tr)) {
        stop("Gtau_mode='tilted' requires non-empty G jump times in training fits.")
    }
    
    get_Gtau_estimated <- function(dat_long, ids_target) {
        G_tau_obj <- predict_G_from_Gk_matrix(
            dat_new = dat_long,
            time_vec = tau,
            start.name = start.name, stop.name = stop.name, event.name = event.name,
            event.C.name = event.C.name, id.name = id.name,
            Gfits = Gfits,
            trim.G = NULL,
            return_subject_info = TRUE
        )
        idx_gt <- match(as.character(ids_target), as.character(G_tau_obj$ids))
        if (anyNA(idx_gt)) {
            miss_ids <- ids_target[is.na(idx_gt)]
            stop("Failed to align G(tau|.) to requested ids. Example ids: ",
                 paste(utils::head(miss_ids, 10), collapse = ", "))
        }
        as.numeric(G_tau_obj$G[idx_gt, 1])
    }
    
    get_Gtau_tilted <- function(dat_long, ids_target) {
        g_est <- get_Gtau_estimated(dat_long, ids_target)
        
        # exact special case requested: delta = 0 returns estimated G(tau | .)
        if (abs(Gtau_delta) <= 1e-12) return(g_est)
        
        tilt_grid <- sort(unique(c(tau, G_jump_times[G_jump_times <= cmax_censor_tr])))
        tilt_grid <- as.numeric(tilt_grid[is.finite(tilt_grid)])
        if (length(tilt_grid) == 0L) return(g_est)
        
        G_obj <- predict_G_from_Gk_matrix(
            dat_new = dat_long,
            time_vec = tilt_grid,
            start.name = start.name, stop.name = stop.name, event.name = event.name,
            event.C.name = event.C.name, id.name = id.name,
            Gfits = Gfits,
            trim.G = NULL,
            return_subject_info = TRUE
        )
        idx_gt <- match(as.character(ids_target), as.character(G_obj$ids))
        if (anyNA(idx_gt)) {
            miss_ids <- ids_target[is.na(idx_gt)]
            stop("Failed to align tilted G rows to requested ids. Example ids: ",
                 paste(utils::head(miss_ids, 10), collapse = ", "))
        }
        G_mat <- G_obj$G[idx_gt, , drop = FALSE]
        
        out <- rep(NA_real_, nrow(G_mat))
        for (ii in seq_len(nrow(G_mat))) {
            gi <- as.numeric(G_mat[ii, ])
            fin <- which(is.finite(gi))
            if (length(fin) == 0L) {
                out[ii] <- g_est[ii]
                next
            }
            last <- max(fin)
            gi <- gi[seq_len(last)]
            tt <- tilt_grid[seq_len(last)]
            
            # dF_C(t_j | .) = G(t_{j-} | .) - G(t_j | .), with G(0-) = 1
            if (length(gi) == 1L) {
                dF <- 1 - gi[1]
            } else {
                dF <- c(1 - gi[1], gi[-length(gi)] - gi[-1])
            }
            dF[!is.finite(dF)] <- 0
            dF <- pmax(dF, 0)
            
            eta <- Gtau_delta * tt
            eta <- eta - max(eta)
            w <- exp(eta)
            
            denom <- sum(w * dF)
            if (!is.finite(denom) || denom <= 0) {
                out[ii] <- g_est[ii]
                next
            }
            num <- sum(w[tt > tau] * dF[tt > tau])
            out[ii] <- num / denom
        }
        # Numerical guard; theoretically out is already in [0,1].
        out <- pmin(pmax(out, 0), 1)
        out
    }
    
    get_Gtau_at_tau <- function(dat_long, ids_target) {
        summ_loc <- subject_summary(dat_long, id.name, stop.name, event.name, event.C.name)
        idx_loc <- match(as.character(ids_target), as.character(summ_loc[[id.name]]))
        if (anyNA(idx_loc)) {
            miss_ids <- ids_target[is.na(idx_loc)]
            stop("get_Gtau_at_tau(): ids_target contains ids not found in dat_long. Example ids: ",
                 paste(utils::head(miss_ids, 10), collapse = ", "))
        }
        X_target <- as.numeric(summ_loc$X[idx_loc])
        
        out <- if (Gtau_mode == "one") {
            rep(1.0, length(ids_target))
        } else if (Gtau_mode == "estimated") {
            get_Gtau_estimated(dat_long, ids_target)
        } else if (Gtau_mode == "tilted") {
            get_Gtau_tilted(dat_long, ids_target)
        } else {
            stop("Unknown Gtau_mode: ", Gtau_mode)
        }
        
        out[X_target <= tau] <- NA_real_
        out
    }
    
    Gtau_tr_tau  <- get_Gtau_at_tau(dat_tr, dat_tr_tau[[id.name]])
    Gtau_cal     <- get_Gtau_at_tau(dat_cal, dat_cal_tau[[id.name]])
    Gtau_cal_all <- get_Gtau_at_tau(dat_cal, ids_cal_all)
    
    X_cal <- dat_cal_tau[[stop.name]]
    Delta_cal <- dat_cal_tau[[event.name]]
    
    hat_W_IPCW <- rep(NA_real_, length(theta_grid))
    for (jj in seq_along(theta_grid)) {
        theta_j <- theta_grid[jj]
        
        lo_j <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = theta_j)
        hi_j <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = 1 - theta_j)
        
        # Open interval indicator for C_{tau,theta}
        inside_j <- as.numeric(X_cal > lo_j & X_cal < hi_j)
        
        U_vec <- (Delta_cal / pmax(dat_cal_tau[["surv.C.subj"]], trim.C)) *
            Gtau_cal *
            as.numeric(X_cal > tau) *
            (inside_j - (1 - alpha))
        
        hat_W_IPCW[jj] <- sum(U_vec, na.rm = TRUE)
    }
    
    idx_nonneg <- which(hat_W_IPCW >= 0)
    if (length(idx_nonneg) == 0L) {
        stop("No theta in theta_grid satisfies sum_i U_IPCW,i(theta) >= 0.")
    }
    theta_index_lower_IPCW <- max(idx_nonneg)
    theta_hat_IPCW <- theta_grid[theta_index_lower_IPCW]
    
    
    # Compute the IPCW lower bound
    lower_IPCW_cal <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = theta_hat_IPCW)
    lower_IPCW_test <- surv_quantile_vec_from_pred(pred_time, pred_surv_test, beta = theta_hat_IPCW)
    
    # Compute the IPCW upper bound
    upper_IPCW_cal <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = 1 - theta_hat_IPCW)
    upper_IPCW_test <- surv_quantile_vec_from_pred(pred_time, pred_surv_test, beta = 1 - theta_hat_IPCW)

    lower_model_cal <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = alpha / 2)
    lower_model_test <- surv_quantile_vec_from_pred(pred_time, pred_surv_test, beta = alpha / 2)
    upper_model_cal <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = 1 - alpha / 2)
    upper_model_test <- surv_quantile_vec_from_pred(pred_time, pred_surv_test, beta = 1 - alpha / 2)

    lower_IPCW_DR_cal <- pmin(lower_IPCW_cal, lower_model_cal)
    lower_IPCW_DR_test <- pmin(lower_IPCW_test, lower_model_test)
    upper_IPCW_DR_cal <- pmax(upper_IPCW_cal, upper_model_cal)
    upper_IPCW_DR_test <- pmax(upper_IPCW_test, upper_model_test)
    
    
    # Estimate empirical coverage (IPCW)
    w_cal  <- dat_cal_tau[[event.name]] / pmax(dat_cal_tau[["surv.C.subj"]], trim.C)
    w_test <- dat_test_tau[[event.name]] / pmax(dat_test_tau[["surv.C.subj"]], trim.C)
    
    coverage_IPCW_cal  <- mean(w_cal  * (lower_IPCW_cal  <= dat_cal_tau[[stop.name]]  & dat_cal_tau[[stop.name]]  <= upper_IPCW_cal)) / mean(w_cal)
    coverage_IPCW_test <- mean(w_test * (lower_IPCW_test <= dat_test_tau[[stop.name]] & dat_test_tau[[stop.name]] <= upper_IPCW_test)) / mean(w_test)
    
    if(!is.null(TT.name)){
        coverage_IPCW_cal_T  <- mean((lower_IPCW_cal  <= dat_cal_tau[[TT.name]]  & dat_cal_tau[[TT.name]]  <= upper_IPCW_cal))
        coverage_IPCW_test_T <- mean((lower_IPCW_test <= dat_test_tau[[TT.name]] & dat_test_tau[[TT.name]] <= upper_IPCW_test))
    }else{
        coverage_IPCW_cal_T  <- NA
        coverage_IPCW_test_T <- NA
    }
    
    ## AIPCW selector on the augmented estimating curve.
    build_Ctau_df <- function(ids_all, idx_tau_rows, q_lo_tau, q_hi_tau, g_tau_tau = NULL) {
        out <- data.frame(
            id = ids_all,
            q_lower = NA_real_,
            q_upper = NA_real_
        )
        out$q_lower[idx_tau_rows] <- q_lo_tau
        out$q_upper[idx_tau_rows] <- q_hi_tau
        if (!is.null(g_tau_tau)) {
            if (length(g_tau_tau) != length(idx_tau_rows)) {
                stop("build_Ctau_df(): g_tau_tau length must match number of tau rows.")
            }
            out$g_tau <- NA_real_
            out$g_tau[idx_tau_rows] <- as.numeric(g_tau_tau)
        }
        out
    }
    
    compute_aug_sum_at_theta <- function(theta_aug) {
        # empty integration grid => zero augmentation
        if (m_u == 0L) {
            return(list(
                aug_sum = 0,
                aug_vec = rep(0, n_cal_all),
                h_obj = NULL
            ))
        }
        
        # Build C_tau(theta) on training ids (for fitting xi_k on training data)
        q_lo_tr <- surv_quantile_vec_from_pred(pred_time, pred_surv_tr, beta = theta_aug)
        q_hi_tr <- surv_quantile_vec_from_pred(pred_time, pred_surv_tr, beta = 1 - theta_aug)
        Ctau_tr_df <- build_Ctau_df(ids_tr_all, idx_tr_tau_in_all, q_lo_tr, q_hi_tr,
                                    g_tau_tau = Gtau_tr_tau)
        
        Xi_fits_theta <- NULL
        if (k_tau > 1L) {
            Xi_fits_theta <- vector("list", k_tau - 1L)
            for (kk in seq_len(k_tau - 1L)) {
                Xi_fits_theta[[kk]] <- fit_Yk_condmean_at_tk(
                    dat_long = dat_tr,
                    visit_times = visit_times,
                    k = kk,
                    tau = tau,
                    start.name = start.name,
                    stop.name = stop.name,
                    event.name = event.name,
                    event.C.name = event.C.name,
                    id.name = id.name,
                    Hk_df = Hk_list_tr[[kk]],
                    Ctau_df = Ctau_tr_df,
                    q_lower.name = "q_lower",
                    q_upper.name = "q_upper",
                    alpha = alpha,
                    covname.baseline = covname.pred.baseline,
                    covname.timevarying = covname.pred.timevarying,
                    use_history = pred_use_history,
                    history_mode = pred_history_mode,
                    model = model.Xi,
                    rsf_args = rsf_args.Xi,
                    xgb_nrounds = xgb_nrounds,
                    xgb_use_cv = xgb_use_cv,
                    xgb_cv_nfold = xgb_cv_nfold,
                    xgb_cv_early_stopping_rounds = xgb_cv_early_stopping_rounds,
                    xgb_cv_param_grid = xgb_cv_param_grid,
                    xgb_cv_seed = xgb_cv_seed,
                    xgb_cv_verbose = xgb_cv_verbose,
                    Xi_fallback_to_mean = Xi_fallback_to_mean,
                    Xi_min_n = Xi_min_n,
                    Xi_min_ess = Xi_min_ess,
                    Xi_min_obs_per_coef = Xi_min_obs_per_coef,
                    Xi_min_class_n = Xi_min_class_n,
                    Xi_min_class_ess = Xi_min_class_ess
                )
            }
        }
        
        # h_tau on calibration (all subjects)
        h_obj <- build_h_tau_matrix_cal(
            dat_cal = dat_cal,
            dat_cal_tau = dat_cal_tau,
            pred_time = pred_time,
            pred_surv_cal = pred_surv_cal,
            visit_times = visit_times,
            tau = tau,
            theta = theta_aug,
            alpha = alpha,
            start.name = start.name,
            stop.name = stop.name,
            event.name = event.name,
            event.C.name = event.C.name,
            id.name = id.name,
            Gfits = Gfits,
            Sfits = Sfits,
            Xi_fits = Xi_fits_theta,
            time_vec = time_vec_AIPCW,
            Gtau_vec = Gtau_cal_all,
            row_data = "all",
            mask_after_X = TRUE,
            after_X_value = NA_real_
        )
        
        idx_h_all <- match(ids_cal_all, as.character(h_obj$ids))
        if (anyNA(idx_h_all)) {
            stop("Failed to align h_tau rows to all calibration ids.")
        }
        h_all <- h_obj$h[idx_h_all, , drop = FALSE]
        if (ncol(h_all) != ncol(aipcw_increment_all)) {
            stop("Mismatch between h_tau columns and augmentation increment columns.")
        }
        
        # For computation: NA placeholders correspond to out-of-support u > X and should contribute 0
        h_use <- h_all
        h_use[is.na(h_use)] <- 0
        
        aug_vec <- rowSums(h_use * aipcw_increment_all)
        list(
            aug_sum = sum(aug_vec),
            aug_vec = aug_vec,
            h_obj = h_obj
        )
    } # End of the function compute_aug_sum_at_theta()
    
    
    hat_W_AIPCW <- rep(NA_real_, length(theta_grid))
    aug_sum_grid <- rep(NA_real_, length(theta_grid))
    h_ref_AIPCW <- NULL
    
    if (isTRUE(localization)) {
        aug_ref <- compute_aug_sum_at_theta(theta_hat_IPCW)
        aug_sum_grid[] <- aug_ref$aug_sum
        h_ref_AIPCW <- aug_ref$h_obj
        hat_W_AIPCW <- hat_W_IPCW + aug_ref$aug_sum
    } else {
        for (jj in seq_along(theta_grid)) {
            aug_j <- compute_aug_sum_at_theta(theta_grid[jj])
            aug_sum_grid[jj] <- aug_j$aug_sum
            if (is.null(h_ref_AIPCW)) h_ref_AIPCW <- aug_j$h_obj
            hat_W_AIPCW[jj] <- hat_W_IPCW[jj] + aug_j$aug_sum
        }
    }
    
    AIPCW_selection <- .select_theta_index_from_curve(
        theta_grid = theta_grid,
        curve = hat_W_AIPCW,
        selector = AIPCW_theta_selector,
        tol = AIPCW_theta_tol
    )
    theta_index_AIPCW <- AIPCW_selection$index
    if (is.na(theta_index_AIPCW)) {
        if (isTRUE(return_on_aipcw_fail)) {
            return(list(
                error = TRUE,
                msg = paste0("No theta selected by AIPCW selector '", AIPCW_theta_selector, "'."),
                lower_IPCW_cal = lower_IPCW_cal,
                upper_IPCW_cal = upper_IPCW_cal,
                lower_IPCW_test = lower_IPCW_test,
                upper_IPCW_test = upper_IPCW_test,
                lower_model_cal = lower_model_cal,
                upper_model_cal = upper_model_cal,
                lower_model_test = lower_model_test,
                upper_model_test = upper_model_test,
                lower_IPCW_DR_cal = lower_IPCW_DR_cal,
                upper_IPCW_DR_cal = upper_IPCW_DR_cal,
                lower_IPCW_DR_test = lower_IPCW_DR_test,
                upper_IPCW_DR_test = upper_IPCW_DR_test,
                coverage_IPCW_cal = coverage_IPCW_cal,
                coverage_IPCW_test = coverage_IPCW_test,
                coverage_IPCW_cal_T = coverage_IPCW_cal_T,
                coverage_IPCW_test_T = coverage_IPCW_test_T,
                theta_grid_used = theta_grid,
                theta_hat_IPCW = theta_hat_IPCW,
                theta_index_IPCW = theta_index_lower_IPCW,
                hat_W_IPCW = hat_W_IPCW,
                aug_sum_grid = aug_sum_grid,
                hat_W_AIPCW = hat_W_AIPCW,
                AIPCW_theta_selector = AIPCW_theta_selector,
                AIPCW_theta_selector_reason = AIPCW_selection$reason,
                AIPCW_theta_tol = AIPCW_theta_tol,
                aipcw_increment_all = aipcw_increment_all,
                G_cal_u_mat_all = G_cal_u_mat_all,
                dG_cal_u_mat_all = dG_cal_u_mat_all,
                time_vec_AIPCW = time_vec_AIPCW,
                censor_times_cal = censor_times_cal,
                G_jump_times = G_jump_times,
                ids_cal_all = ids_cal_all,
                X_cal_all = X_cal_all,
                Delta_cal_all = Delta_cal_all,
                DeltaC_cal_all = DeltaC_cal_all,
                idx_tau_in_all = idx_tau_in_all,
                id_tr = id_tr,
                id_cal_tau = dat_cal_tau[[id.name]],
                id_test_tau = dat_test_tau[[id.name]],
                h_ref_AIPCW = h_ref_AIPCW,
                Gtau_mode = Gtau_mode,
                Gtau_delta = Gtau_delta,
                cmax_censor_tr = cmax_censor_tr,
                model.Xi = model.Xi,
                survC_cal = dat_cal_tau[["surv.C.subj"]],
                survC_test = dat_test_tau[["surv.C.subj"]]
            ))
        }
        stop("No theta selected by AIPCW selector '", AIPCW_theta_selector, "'.")
    }
    theta_hat_AIPCW <- theta_grid[theta_index_AIPCW]
    
    lower_AIPCW_cal <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = theta_hat_AIPCW)
    lower_AIPCW_test <- surv_quantile_vec_from_pred(pred_time, pred_surv_test, beta = theta_hat_AIPCW)
    upper_AIPCW_cal <- surv_quantile_vec_from_pred(pred_time, pred_surv_cal, beta = 1 - theta_hat_AIPCW)
    upper_AIPCW_test <- surv_quantile_vec_from_pred(pred_time, pred_surv_test, beta = 1 - theta_hat_AIPCW)
    
    coverage_AIPCW_cal  <- mean(w_cal  * (lower_AIPCW_cal  <= dat_cal_tau[[stop.name]]  & dat_cal_tau[[stop.name]]  <= upper_AIPCW_cal)) / mean(w_cal)
    coverage_AIPCW_test <- mean(w_test * (lower_AIPCW_test <= dat_test_tau[[stop.name]] & dat_test_tau[[stop.name]] <= upper_AIPCW_test)) / mean(w_test)
    
    if(!is.null(TT.name)){
        coverage_AIPCW_cal_T  <- mean((lower_AIPCW_cal  <= dat_cal_tau[[TT.name]]  & dat_cal_tau[[TT.name]]  <= upper_AIPCW_cal))
        coverage_AIPCW_test_T <- mean((lower_AIPCW_test <= dat_test_tau[[TT.name]] & dat_test_tau[[TT.name]] <= upper_AIPCW_test))
    } else {
        coverage_AIPCW_cal_T  <- NA
        coverage_AIPCW_test_T <- NA
    }
    
    
    return(list(lower_IPCW_cal = lower_IPCW_cal, upper_IPCW_cal = upper_IPCW_cal,
                lower_IPCW_test = lower_IPCW_test, upper_IPCW_test = upper_IPCW_test,
                lower_model_cal = lower_model_cal, upper_model_cal = upper_model_cal,
                lower_model_test = lower_model_test, upper_model_test = upper_model_test,
                lower_IPCW_DR_cal = lower_IPCW_DR_cal, upper_IPCW_DR_cal = upper_IPCW_DR_cal,
                lower_IPCW_DR_test = lower_IPCW_DR_test, upper_IPCW_DR_test = upper_IPCW_DR_test,
                coverage_IPCW_cal = coverage_IPCW_cal, coverage_IPCW_test = coverage_IPCW_test,
                coverage_IPCW_cal_T = coverage_IPCW_cal_T, coverage_IPCW_test_T = coverage_IPCW_test_T,
                # Backward-compatible aliases
                coverage_cal = coverage_IPCW_cal, coverage_test = coverage_IPCW_test,
                coverage_cal_T = coverage_IPCW_cal_T, coverage_test_T = coverage_IPCW_test_T,
                # IPCW selector diagnostics
                theta_hat_IPCW = theta_hat_IPCW,
                theta_index_IPCW = theta_index_lower_IPCW,
                hat_W_IPCW = hat_W_IPCW,
                theta_grid_used = theta_grid,
                # AIPCW integration-grid diagnostics
                lower_AIPCW_cal = lower_AIPCW_cal, upper_AIPCW_cal = upper_AIPCW_cal,
                lower_AIPCW_test = lower_AIPCW_test, upper_AIPCW_test = upper_AIPCW_test,
                coverage_AIPCW_cal = coverage_AIPCW_cal, coverage_AIPCW_test = coverage_AIPCW_test,
                coverage_AIPCW_cal_T = coverage_AIPCW_cal_T, coverage_AIPCW_test_T = coverage_AIPCW_test_T,
                theta_hat_AIPCW = theta_hat_AIPCW,
                theta_index_AIPCW = theta_index_AIPCW,
                AIPCW_theta_selector = AIPCW_theta_selector,
                AIPCW_theta_selector_reason = AIPCW_selection$reason,
                AIPCW_theta_tol = AIPCW_theta_tol,
                hat_W_AIPCW = hat_W_AIPCW, aug_sum_grid = aug_sum_grid,
                time_vec_AIPCW = time_vec_AIPCW,
                censor_times_cal = censor_times_cal,
                G_jump_times = G_jump_times,
                ids_cal_all = ids_cal_all,
                X_cal_all = X_cal_all,
                Delta_cal_all = Delta_cal_all,
                DeltaC_cal_all = DeltaC_cal_all,
                idx_tau_in_all = idx_tau_in_all,
                id_tr = id_tr,
                id_cal_tau = dat_cal_tau[[id.name]],
                id_test_tau = dat_test_tau[[id.name]],
                h_ref_AIPCW = h_ref_AIPCW,
                Gtau_mode = Gtau_mode,
                Gtau_delta = Gtau_delta,
                cmax_censor_tr = cmax_censor_tr,
                model.Xi = model.Xi,
                survC_cal = dat_cal_tau[["surv.C.subj"]], survC_test = dat_test_tau[["surv.C.subj"]]))

}














## Functions from Farina et al. (2025) ----------------------------------------------------------------------------------

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



#' Compute the estimated survival probability at a given time
#'
#' Retrieves the estimated survival probability for an individual `i` at time `t`,
#' based on the predicted survival curve.
#'
#' @param i Index of the individual in the calibration dataset
#' @param t Time at which survival probability is evaluated
#' @param pred_time Vector of time points for which survival probabilities were estimated
#' @param pred_surv Matrix of survival probabilities (rows correspond to individuals, columns to time points)
#' @return Estimated survival probability at time `t` for individual `i`
#'
surv_curve <- function(
        i, t, pred_time, pred_surv
) {
    # If the requested time t is before the first prediction time, return survival probability 1
    if (t < min(pred_time)) {
        return(1)
    } else {
        # Find the largest time point in pred_time that is less than or equal to t
        index <- max(which(pred_time <= t))
        return(pred_surv[i, index])  # Return the corresponding survival probability
    }
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
