# Bridge utilities for Python Dynamic-DeepHit style models.
# This file is only used when model.pred == "ddh" (or legacy "ddh_py").

`%||%` <- function(a, b) if (!is.null(a)) a else b

.ddh_bridge_state <- new.env(parent = emptyenv())

load_ddh_bridge <- function(py_module_path = "python/ddh_bridge.py",
                            force_reload = FALSE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package 'reticulate' is required for model.pred='ddh'.")
    }

    py_module_path <- normalizePath(py_module_path, mustWork = TRUE)
    cached <- !is.null(.ddh_bridge_state$py_module_path) &&
        identical(.ddh_bridge_state$py_module_path, py_module_path) &&
        !is.null(.ddh_bridge_state$bridge) &&
        !isTRUE(force_reload)
    if (cached) return(.ddh_bridge_state$bridge)

    bridge_env <- new.env(parent = baseenv())
    reticulate::source_python(py_module_path, envir = bridge_env)
    .ddh_bridge_state$bridge <- bridge_env
    .ddh_bridge_state$py_module_path <- py_module_path
    bridge_env
}

.ddh_numeric_prediction_data <- function(dat_list, covars, train_name = "dat_tr") {
    covars <- covars[!is.na(covars) & nzchar(covars)]
    dat_list <- lapply(dat_list, as.data.frame)
    if (length(covars) == 0) {
        return(list(dat_list = dat_list, covars = character(0)))
    }

    train_dat <- dat_list[[train_name]]
    new_covars <- character(0)

    for (v in covars) {
        if (!v %in% names(train_dat)) next
        x_train <- train_dat[[v]]
        is_numeric_like <- is.numeric(x_train) || is.integer(x_train) || is.logical(x_train)

        if (is_numeric_like) {
            x_num <- suppressWarnings(as.numeric(x_train))
            med <- suppressWarnings(stats::median(x_num[is.finite(x_num)], na.rm = TRUE))
            if (!is.finite(med)) med <- 0
            for (nm in names(dat_list)) {
                x <- suppressWarnings(as.numeric(dat_list[[nm]][[v]]))
                x[!is.finite(x)] <- med
                dat_list[[nm]][[v]] <- x
            }
            new_covars <- c(new_covars, v)
            next
        }

        x_chr <- as.character(x_train)
        has_missing <- any(is.na(x_chr))
        x_chr[is.na(x_chr)] <- "__MISSING__"
        if (is.factor(x_train)) {
            lev <- levels(droplevels(x_train))
            if (has_missing) lev <- c(lev, "__MISSING__")
        } else {
            lev <- sort(unique(x_chr))
        }
        if (length(lev) <= 1L) next

        lev_use <- lev[-1L]
        dummy_names <- make.names(paste0(v, "_ddh_", lev_use), unique = TRUE)
        for (nm in names(dat_list)) {
            x <- as.character(dat_list[[nm]][[v]])
            x[is.na(x)] <- "__MISSING__"
            for (j in seq_along(lev_use)) {
                dat_list[[nm]][[dummy_names[j]]] <- as.numeric(x == lev_use[j])
            }
        }
        new_covars <- c(new_covars, dummy_names)
    }

    list(dat_list = dat_list, covars = new_covars)
}

predict_surv_ddh <- function(dat_tr, dat_cal, dat_test,
                                dat_tr_tau, dat_cal_tau, dat_test_tau,
                                start.name, stop.name, event.name, id.name, tau,
                                cov_pred_use, model.pred.args = list()) {
    py_path <- model.pred.args$ddh_bridge_path %||% "python/ddh_bridge.py"
    py_reload <- isTRUE(model.pred.args$py_force_reload)

    bridge <- load_ddh_bridge(py_module_path = py_path, force_reload = py_reload)

    # remove bridge-loading keys from the payload passed to python
    payload <- model.pred.args
    payload$ddh_bridge_path <- NULL
    payload$py_force_reload <- NULL

    if (is.null(bridge$fit_predict_ddh)) {
        stop("Function fit_predict_ddh() not found in python bridge: ", py_path)
    }

    numeric_dat <- .ddh_numeric_prediction_data(
        dat_list = list(
            dat_tr = dat_tr,
            dat_cal = dat_cal,
            dat_test = dat_test,
            dat_tr_tau = dat_tr_tau,
            dat_cal_tau = dat_cal_tau,
            dat_test_tau = dat_test_tau
        ),
        covars = cov_pred_use,
        train_name = "dat_tr"
    )
    dat_tr <- numeric_dat$dat_list$dat_tr
    dat_cal <- numeric_dat$dat_list$dat_cal
    dat_test <- numeric_dat$dat_list$dat_test
    dat_tr_tau <- numeric_dat$dat_list$dat_tr_tau
    dat_cal_tau <- numeric_dat$dat_list$dat_cal_tau
    dat_test_tau <- numeric_dat$dat_list$dat_test_tau
    cov_pred_use <- numeric_dat$covars

    out <- bridge$fit_predict_ddh(
        dat_tr = dat_tr,
        dat_cal = dat_cal,
        dat_test = dat_test,
        dat_tr_tau = dat_tr_tau,
        dat_cal_tau = dat_cal_tau,
        dat_test_tau = dat_test_tau,
        tau = as.numeric(tau),
        start_col = start.name,
        stop_col = stop.name,
        event_col = event.name,
        id_col = id.name,
        covars = as.list(cov_pred_use),
        args = payload
    )

    pred_time <- as.numeric(out$pred_time %||% out$time_grid)
    pred_surv_tr <- as.matrix(out$pred_surv_tr %||% out$surv_tr)
    pred_surv_cal <- as.matrix(out$pred_surv_cal %||% out$surv_cal)
    pred_surv_test <- as.matrix(out$pred_surv_test %||% out$surv_test)

    if (length(pred_time) == 0) stop("ddh bridge returned empty pred_time.")
    if (ncol(pred_surv_tr) != length(pred_time)) {
        stop("ddh: ncol(pred_surv_tr) != length(pred_time).")
    }
    if (ncol(pred_surv_cal) != length(pred_time)) {
        stop("ddh: ncol(pred_surv_cal) != length(pred_time).")
    }
    if (ncol(pred_surv_test) != length(pred_time)) {
        stop("ddh: ncol(pred_surv_test) != length(pred_time).")
    }

    # Optional id-based row alignment if python returns ids.
    id_tr_py <- out$id_tr
    id_cal_py <- out$id_cal
    id_test_py <- out$id_test
    if (!is.null(id_tr_py)) {
        idx <- match(dat_tr_tau[[id.name]], as.vector(id_tr_py))
        if (anyNA(idx)) stop("ddh: id mismatch for training predictions.")
        pred_surv_tr <- pred_surv_tr[idx, , drop = FALSE]
    } else if (nrow(pred_surv_tr) != nrow(dat_tr_tau)) {
        stop("ddh: row mismatch in training predictions and id_tr not provided.")
    }
    if (!is.null(id_cal_py)) {
        idx <- match(dat_cal_tau[[id.name]], as.vector(id_cal_py))
        if (anyNA(idx)) stop("ddh: id mismatch for calibration predictions.")
        pred_surv_cal <- pred_surv_cal[idx, , drop = FALSE]
    } else if (nrow(pred_surv_cal) != nrow(dat_cal_tau)) {
        stop("ddh: row mismatch in calibration predictions and id_cal not provided.")
    }
    if (!is.null(id_test_py)) {
        idx <- match(dat_test_tau[[id.name]], as.vector(id_test_py))
        if (anyNA(idx)) stop("ddh: id mismatch for test predictions.")
        pred_surv_test <- pred_surv_test[idx, , drop = FALSE]
    } else if (nrow(pred_surv_test) != nrow(dat_test_tau)) {
        stop("ddh: row mismatch in test predictions and id_test not provided.")
    }

    # Keep a sorted increasing time grid.
    ord <- order(pred_time)
    pred_time <- pred_time[ord]
    pred_surv_tr <- pred_surv_tr[, ord, drop = FALSE]
    pred_surv_cal <- pred_surv_cal[, ord, drop = FALSE]
    pred_surv_test <- pred_surv_test[, ord, drop = FALSE]

    list(
        pred_time = pred_time,
        pred_surv_tr = pred_surv_tr,
        pred_surv_cal = pred_surv_cal,
        pred_surv_test = pred_surv_test
    )
}

# Backward compatibility alias
predict_surv_ddh_py <- predict_surv_ddh
