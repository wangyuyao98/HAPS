
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


# -------- helper: only keep bariables that are not constant --------
# keep_nonconstant_covars <- function(df, vars) {
#     vars <- vars[!is.na(vars) & nzchar(vars)]
#     if (length(vars) == 0) return(character(0))
#     
#     keep <- vapply(vars, function(v) {
#         if (!v %in% names(df)) return(FALSE)
#         x <- df[[v]]
#         if (is.factor(x)) {
#             nlevels(droplevels(x)) >= 2
#         } else {
#             ux <- unique(x[is.finite(x)])
#             length(ux) >= 2
#         }
#     }, logical(1))
#     
#     vars[keep]
# }

keep_nonconstant_covars <- function(df, vars, min_unique = 2) {
    vars <- vars[!is.na(vars) & nzchar(vars)]
    if (length(vars) == 0) return(character(0))
    
    keep <- vapply(vars, function(v) {
        if (!v %in% names(df)) return(FALSE)
        
        x <- df[[v]]
        
        # drop missing values (and for numeric drop non-finite)
        if (is.factor(x)) {
            x2 <- x[!is.na(x)]
            if (length(x2) == 0) return(FALSE)
            nlevels(droplevels(x2)) >= min_unique
            
        } else if (is.character(x)) {
            x2 <- x[!is.na(x)]
            if (length(x2) == 0) return(FALSE)
            length(unique(x2)) >= min_unique
            
        } else if (inherits(x, c("Date", "POSIXct", "POSIXlt"))) {
            x2 <- x[!is.na(x)]
            if (length(x2) == 0) return(FALSE)
            length(unique(x2)) >= min_unique
            
        } else if (is.logical(x)) {
            x2 <- x[!is.na(x)]
            if (length(x2) == 0) return(FALSE)
            length(unique(x2)) >= min_unique
            
        } else {
            # numeric / integer and other atomic vectors
            x2 <- x[!is.na(x)]
            if (length(x2) == 0) return(FALSE)
            if (is.numeric(x2)) x2 <- x2[is.finite(x2)]
            if (length(x2) == 0) return(FALSE)
            length(unique(x2)) >= min_unique
        }
    }, logical(1))
    
    vars[keep]
}



# -------- helper: long -> subject (for preparing data set to apply bootCP) --------
long_to_subject_df <- function(dat_long, id.name, start.name, stop.name, event.name,
                               covars, TT.name = NULL) {
    dat_long <- dat_long[order(dat_long[[id.name]], dat_long[[start.name]]), ]
    
    # baseline row per subject
    base <- dat_long %>%
        dplyr::group_by(.data[[id.name]]) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
    
    # subject-level observed time + event indicator
    ydelta <- dat_long %>%
        dplyr::group_by(.data[[id.name]]) %>%
        dplyr::summarise(
            y = max(.data[[stop.name]]),
            delta = max(.data[[event.name]]),
            .groups = "drop"
        )
    
    # columns to carry from the baseline row
    cols <- unique(c(id.name, covars, TT.name))
    cols <- cols[!is.na(cols) & nzchar(cols)]  # drops NULLs/NA/""
    
    # optional: strict check (recommended)
    missing_cols <- setdiff(cols, names(base))
    if (length(missing_cols) > 0) {
        stop("long_to_subject_df: these columns are missing in dat_long/base: ",
             paste(missing_cols, collapse = ", "))
    }
    
    out <- dplyr::left_join(ydelta, base[, cols, drop = FALSE], by = id.name)
    as.data.frame(out)
}

## ---- helper: wrapper for one tau (same signature as your simulation wrapper) ----
run_one_tau <- function(dat_long, dat_test, tau,
                        start.name, stop.name, event.name, event.C.name, id.name,
                        covname.pred.baseline, covname.pred.timevarying, # pred_tvtime,
                        model.pred, model.C, covname.C.baseline, covname.C.timevarying,
                        theta_grid, alpha, trim.C,
                        TT.name,
                        seed, train_ratio) {
    
    out <- tryCatch({
        dynamicCP_split(
            dat = dat_long, dat_test = dat_test,
            start.name = start.name, stop.name = stop.name,
            event.name = event.name, event.C.name = event.C.name, id.name = id.name,
            tau = tau,
            model.pred = model.pred,
            covname.pred.baseline = covname.pred.baseline,
            covname.pred.timevarying = covname.pred.timevarying,
            # pred_tvtime = pred_tvtime,
            model.C = model.C,
            covname.C.baseline = covname.C.baseline,
            covname.C.timevarying = covname.C.timevarying,
            TT.name = TT.name,
            seed = seed, train_ratio = train_ratio,
            theta_grid = theta_grid,
            alpha = alpha, trim.C = trim.C
        )
    }, error = function(e) {
        list(error = TRUE, msg = conditionMessage(e))
    })
    
    out
}


## ---- helper: Function for adding baseline measurements of the time-varying covariates as baseline covariates ---
add_baseline_covariates <- function(dat, id.name, tv_covs, suffix = "0") {
    dat <- dat[order(dat[[id.name]], dat[["start"]]), ]   # ensure ordered by time within id
    for (v in tv_covs) {
        base_name <- paste0(v, suffix)
        dat[[base_name]] <- ave(dat[[v]], dat[[id.name]], FUN = function(x) x[1])
    }
    dat
}



# --------------------- helpers for plots ----------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ============================================================
# Helper: pool bounds over reps by tau
# results_bounds: list over reps; each rep is list over tau-index; each tau has $test data.frame
# ============================================================
pool_bounds_by_tau <- function(results_bounds, tau_list,
                               which_set = c("test", "cal"),
                               what = c("lower", "upper", "len")) {
    which_set <- match.arg(which_set)
    what <- match.arg(what)
    
    out <- vector("list", length(tau_list))
    names(out) <- paste0("tau_", tau_list)
    
    for (k in seq_along(tau_list)) {
        vecs <- lapply(results_bounds, function(brep) {
            obj <- brep[[k]][[which_set]]
            if (is.null(obj) || nrow(obj) == 0) return(numeric(0))
            if (what == "lower") return(obj$lower)
            if (what == "upper") return(obj$upper)
            obj$upper - obj$lower
        })
        out[[k]] <- unlist(vecs, use.names = FALSE)
    }
    out
}

multi_pos <- function(tau_list, m, offset = 0.18) {
    shifts <- seq(-offset, offset, length.out = m)
    rep(seq_along(tau_list), each = m) + rep(shifts, times = length(tau_list))
}

interleave_k_by_tau <- function(lists, tau_list, names) {
    K <- length(lists)
    out <- vector("list", K * length(tau_list))
    nm  <- character(K * length(tau_list))
    idx <- 1
    for (k in seq_along(tau_list)) {
        for (j in seq_len(K)) {
            out[[idx]] <- lists[[j]][[k]]
            nm[idx] <- sprintf("%s_tau=%.1f", names[j], tau_list[k])
            idx <- idx + 1
        }
    }
    names(out) <- nm
    out
}

build_cov_list_k <- function(res_list, tau_list, names, col = "coverage_test_true") {
    K <- length(res_list)
    out <- vector("list", K * length(tau_list))
    nm  <- character(K * length(tau_list))
    idx <- 1
    for (k in seq_along(tau_list)) {
        t <- tau_list[k]
        for (j in seq_len(K)) {
            out[[idx]] <- res_list[[j]][[col]][res_list[[j]]$tau == t]
            nm[idx] <- sprintf("%s_tau=%.1f", names[j], t)
            idx <- idx + 1
        }
    }
    names(out) <- nm
    out
}

drop_empty_groups <- function(vals_list, pos, cols) {
    keep <- sapply(vals_list, function(v) any(is.finite(v)))
    list(vals = vals_list[keep], pos = pos[keep], cols = cols[keep], keep = keep)
}

# ============================================================
# Helper: bootCP expand to all taus (TEST ONLY)
# Assumes boot_obj$bootCP$results_bounds_bootCP[[r]]$tau_0$test exists with T_true
# ============================================================
build_bootCP_static_testonly <- function(results_bounds_tau0, tau_grid,
                                         TT.col = "T_true",
                                         X.col  = "X") {
    R <- length(results_bounds_tau0)
    out_bounds  <- vector("list", R)
    out_summary <- vector("list", R)
    
    stopifnot(tau_grid[1] == 0)
    
    for (r in seq_len(R)) {
        b0wrap <- results_bounds_tau0[[r]]
        if (is.null(b0wrap) || is.null(b0wrap$tau_0) || is.null(b0wrap$tau_0$test)) {
            out_bounds[[r]] <- vector("list", length(tau_grid))
            out_summary[[r]] <- data.frame(rep=r, tau=tau_grid, ok=FALSE,
                                           coverage_test_true=NA_real_,
                                           mean_len_test=NA_real_,
                                           msg="Missing tau_0/test")
            next
        }
        
        df0 <- b0wrap$tau_0$test
        if (!(TT.col %in% names(df0))) {
            out_bounds[[r]] <- vector("list", length(tau_grid))
            out_summary[[r]] <- data.frame(rep=r, tau=tau_grid, ok=FALSE,
                                           coverage_test_true=NA_real_,
                                           mean_len_test=NA_real_,
                                           msg=paste0("Missing ", TT.col))
            next
        }
        
        tab_r <- data.frame(
            rep = r, tau = tau_grid, ok = TRUE, msg = NA_character_,
            coverage_test_true = NA_real_,
            mean_len_test = NA_real_
        )
        
        bounds_r <- vector("list", length(tau_grid))
        names(bounds_r) <- paste0("tau_", tau_grid)
        
        for (k in seq_along(tau_grid)) {
            tau <- tau_grid[k]
            
            keep <- which(df0[[TT.col]] > tau)
            dfk <- df0[keep, , drop = FALSE]
            # Extreme case: nobody in the test set survives past tau. The bare
            # `dfk$tau <- tau` below would fail on a 0-row data frame with a
            # cryptic "replacement has 1 row, data has 0"; error informatively.
            if (nrow(dfk) == 0) {
                stop("build_bootCP_static_testonly(): no test subjects with ",
                     TT.col, " > tau = ", tau, " (replication ", r,
                     "); cannot compute static bounds at this tau.")
            }
            dfk$tau <- tau
            
            if (X.col %in% names(dfk)) {
                dfk$covered_X <- as.integer(dfk$lower <= dfk[[X.col]] & dfk[[X.col]] <= dfk$upper)
            } else {
                dfk$covered_X <- NA_integer_
            }
            dfk$covered_T <- as.integer(dfk$lower <= dfk[[TT.col]] & dfk[[TT.col]] <= dfk$upper)
            
            bounds_r[[k]] <- list(test = dfk)
            
            if (nrow(dfk) == 0) {
                tab_r$coverage_test_true[k] <- NA_real_
                tab_r$mean_len_test[k] <- NA_real_
            } else {
                tab_r$coverage_test_true[k] <- mean(dfk$covered_T, na.rm = TRUE)
                tab_r$mean_len_test[k] <- mean(dfk$upper - dfk$lower, na.rm = TRUE)
            }
        }
        
        out_bounds[[r]]  <- bounds_r
        out_summary[[r]] <- tab_r
    }
    
    list(results_bounds_static = out_bounds,
         results_summary_static = out_summary)
}


# bootCP expand to all taus (TEST ONLY) for application:
# subset survivors at tau using observed X, and compute IPCW coverage on X
build_bootCP_static_testonly_ipcw <- function(results_bounds_tau0, tau_grid,
                                              X.col = "X", w.col = "w_ipcw") {
    R <- length(results_bounds_tau0)
    out_bounds  <- vector("list", R)
    out_summary <- vector("list", R)
    
    stopifnot(tau_grid[1] == 0)
    
    for (r in seq_len(R)) {
        b0wrap <- results_bounds_tau0[[r]]
        if (is.null(b0wrap) || is.null(b0wrap$tau_0) || is.null(b0wrap$tau_0$test)) {
            out_bounds[[r]] <- vector("list", length(tau_grid))
            out_summary[[r]] <- data.frame(rep=r, tau=tau_grid, ok=FALSE,
                                           coverage_test_ipcw=NA_real_,
                                           mean_len_test=NA_real_,
                                           msg="Missing tau_0/test")
            next
        }
        
        df0 <- b0wrap$tau_0$test
        if (!(X.col %in% names(df0)) || !(w.col %in% names(df0))) {
            out_bounds[[r]] <- vector("list", length(tau_grid))
            out_summary[[r]] <- data.frame(rep=r, tau=tau_grid, ok=FALSE,
                                           coverage_test_ipcw=NA_real_,
                                           mean_len_test=NA_real_,
                                           msg=paste0("Missing ", X.col, " or ", w.col))
            next
        }
        
        tab_r <- data.frame(
            rep = r, tau = tau_grid, ok = TRUE, msg = NA_character_,
            coverage_test_ipcw = NA_real_,
            mean_len_test = NA_real_
        )
        
        bounds_r <- vector("list", length(tau_grid))
        names(bounds_r) <- paste0("tau_", tau_grid)
        
        for (k in seq_along(tau_grid)) {
            tau <- tau_grid[k]
            
            # survivors at tau: observed time X > tau
            keep <- which(df0[[X.col]] > tau)
            dfk <- df0[keep, , drop = FALSE]
            # Extreme case: no observed survivors past tau (see the analogous
            # guard in build_bootCP_static_testonly); error informatively.
            if (nrow(dfk) == 0) {
                stop("build_bootCP_static_testonly_ipcw(): no test subjects with ",
                     X.col, " > tau = ", tau, " (replication ", r,
                     "); cannot compute static bounds at this tau.")
            }
            dfk$tau <- tau
            
            dfk$covered_X <- as.integer(dfk$lower <= dfk[[X.col]] & dfk[[X.col]] <= dfk$upper)
            
            bounds_r[[k]] <- list(test = dfk)
            
            if (nrow(dfk) == 0) {
                tab_r$coverage_test_ipcw[k] <- NA_real_
                tab_r$mean_len_test[k] <- NA_real_
            } else {
                w <- dfk[[w.col]]
                tab_r$coverage_test_ipcw[k] <- mean(w * dfk$covered_X, na.rm = TRUE) / mean(w, na.rm = TRUE)
                tab_r$mean_len_test[k] <- mean(dfk$upper - dfk$lower, na.rm = TRUE)
            }
        }
        
        out_bounds[[r]]  <- bounds_r
        out_summary[[r]] <- tab_r
    }
    
    list(results_bounds_static = out_bounds,
         results_summary_static = out_summary)
}



# ============================================================
# Plot function: 4-panel boxplots
# ============================================================

make_4panel_boxplots <- function(methods, tau_list, alpha, out_pdf,
                                 method_order = names(methods),
                                 method_label = NULL,
                                 boxfill_cols, boxborder_cols,
                                 offset = 0.18, boxwex = 0.05,
                                 ylim_coverage = NULL, ylim_length = NULL,
                                 ylim_lower = NULL, ylim_upper = NULL,
                                 alpha_line_col = "grey60", ref_tau_col = "grey60",
                                 align = c("h", "v"),
                                 width = 9, height = 8,
                                 plot_cex = list(),
                                 coverage_col = "coverage_test_true",
                                 coverage_ticks = NULL) {
    
    # ---- font size defaults ----
    cex_legend <- plot_cex$legend %||% 1.0
    cex_axis   <- plot_cex$axis   %||% 1.0
    cex_lab    <- plot_cex$lab    %||% 1.1
    cex_main   <- plot_cex$main   %||% 1.1
    
    align <- match.arg(align)
    
    stopifnot(length(method_order) == length(boxfill_cols),
              length(method_order) == length(boxborder_cols))
    
    K <- length(method_order)
    
    res_list <- lapply(method_order, function(nm) methods[[nm]]$res_ok)
    bnd_list <- lapply(method_order, function(nm) methods[[nm]]$bounds)
    
    # # positions and colors
    # pos <- multi_pos(tau_list, m = K, offset = offset)
    # cols <- rep(boxfill_cols, times = length(tau_list))
    # border_cols_full <- rep(boxborder_cols, times = length(tau_list))
    # positions and colors (IMPORTANT: use index positions, not tau numeric values)
    tau_centers <- seq_along(tau_list)
    pos <- multi_pos(tau_centers, m = K, offset = offset)
    
    cols <- rep(boxfill_cols, times = length(tau_list))
    border_cols_full <- rep(boxborder_cols, times = length(tau_list))
    
    
    # Coverage list
    cov_list <- build_cov_list_k(res_list, tau_list, method_order, col = coverage_col)
    
    # Bounds pooled
    lowers <- lapply(bnd_list, pool_bounds_by_tau, tau_list = tau_list, which_set="test", what="lower")
    uppers <- lapply(bnd_list, pool_bounds_by_tau, tau_list = tau_list, which_set="test", what="upper")
    lens   <- lapply(bnd_list, pool_bounds_by_tau, tau_list = tau_list, which_set="test", what="len")
    
    lower_list <- interleave_k_by_tau(lowers, tau_list, method_order)
    upper_list <- interleave_k_by_tau(uppers, tau_list, method_order)
    len_list   <- interleave_k_by_tau(lens,   tau_list, method_order)
    
    # drop empty groups
    cov_obj   <- drop_empty_groups(cov_list,   pos, cols)
    upper_obj <- drop_empty_groups(upper_list, pos, cols)
    len_obj   <- drop_empty_groups(len_list,   pos, cols)
    lower_obj <- drop_empty_groups(lower_list, pos, cols)
    
    # tau_centers <- seq_along(tau_list)
    
    # helper to map "keep" to borders
    border_cov   <- border_cols_full[cov_obj$keep]
    border_upper <- border_cols_full[upper_obj$keep]
    border_len   <- border_cols_full[len_obj$keep]
    border_lower <- border_cols_full[lower_obj$keep]
    
    pdf(out_pdf, width = width, height = height)
    
    if (align == "h") {
        # Original horizontal layout: 2x2, tau on y-axis (horizontal=TRUE)
        par(mfrow = c(2, 2),
            mar = c(2.5, 2.0, 2.3, 1.2),
            oma = c(0, 2.5, 0, 0),
            cex.axis = cex_axis, cex.lab = cex_lab)
        
        # 1) Coverage
        boxplot(cov_obj$vals,
                at = -cov_obj$pos,
                col = cov_obj$cols,
                border = border_cov,
                horizontal = TRUE,
                yaxt = "n",
                xaxt = if (is.null(coverage_ticks)) "s" else "n",
                main = "Coverage", cex.main = cex_main,
                ylim = ylim_coverage,
                boxwex = boxwex,
                outline = FALSE)
        abline(v = 1 - alpha, lty = 2, col = alpha_line_col)
        if (!is.null(coverage_ticks)) {
            axis(1, at = coverage_ticks, labels = coverage_ticks, cex.axis = cex_axis)
        }
        axis(2, at = -tau_centers, labels = tau_list, las = 1, cex.axis = cex_axis)
        
        legend_text <- method_order
        if (!is.null(method_label)) {
            legend_text <- unname(method_label[method_order])
            legend_text[is.na(legend_text)] <- method_order[is.na(legend_text)]
        }
        legend("topleft",
               legend = legend_text,
               fill = boxfill_cols,
               border = boxborder_cols,
               bty = "n",
               cex = cex_legend)
        
        # 3) Length
        boxplot(len_obj$vals,
                at = -len_obj$pos,
                col = len_obj$cols,
                border = border_len,
                horizontal = TRUE,
                yaxt = "n",
                main = "Length", cex.main = cex_main,
                ylim = ylim_length,
                boxwex = boxwex,
                outline = FALSE)
        
        # 4) Lower
        boxplot(lower_obj$vals,
                at = -lower_obj$pos,
                col = lower_obj$cols,
                border = border_lower,
                horizontal = TRUE,
                yaxt = "n",
                main = "Lower", cex.main = cex_main,
                ylim = ylim_lower,
                boxwex = boxwex,
                outline = FALSE)
        abline(v = tau_list, lty = 2, col = ref_tau_col)
        axis(2, at = -tau_centers, labels = tau_list, las = 1, cex.axis = cex_axis)
        
        
        # 2) Upper
        boxplot(upper_obj$vals,
                at = -upper_obj$pos,
                col = upper_obj$cols,
                border = border_upper,
                horizontal = TRUE,
                yaxt = "n",
                main = "Upper", cex.main = cex_main,
                ylim = ylim_upper,
                boxwex = boxwex,
                outline = FALSE)
        
        
        # shared y-label
        mtext(expression("Prediction time " * tau),
              side = 2, outer = TRUE, line = 1.2, cex = cex_lab)
        
    } else {
        # Vertical layout with requested panel arrangement:
        #   [ Coverage | Upper ]
        #   [ Length   | Lower ]
        par(mfrow = c(2, 2),
            mar = c(3.0, 3.5, 2.3, 1.2),
            oma = c(2.0, 0.5, 0, 0),
            cex.axis = cex_axis, cex.lab = cex_lab)
        
        # top-left: Coverage
        boxplot(cov_obj$vals,
                at = cov_obj$pos,
                col = cov_obj$cols,
                border = border_cov,
                xaxt = "n",
                yaxt = if (is.null(coverage_ticks)) "s" else "n",
                main = "Coverage", cex.main = cex_main,
                ylim = ylim_coverage,
                boxwex = boxwex,
                outline = FALSE)
        abline(h = 1 - alpha, lty = 2, col = alpha_line_col)
        axis(1, at = tau_centers, labels = tau_list, cex.axis = cex_axis)
        if (!is.null(coverage_ticks)) {
            axis(2, at = coverage_ticks, labels = coverage_ticks, cex.axis = cex_axis)
        }
        
        legend_text <- method_order
        if (!is.null(method_label)) {
            legend_text <- unname(method_label[method_order])
            legend_text[is.na(legend_text)] <- method_order[is.na(legend_text)]
        }
        legend("topleft",
               legend = legend_text,
               fill = boxfill_cols,
               border = boxborder_cols,
               bty = "n",
               cex = cex_legend)
        
        # top-right: Upper
        boxplot(upper_obj$vals,
                at = upper_obj$pos,
                col = upper_obj$cols,
                border = border_upper,
                xaxt = "n",
                main = "Upper", cex.main = cex_main,
                ylim = ylim_upper,
                boxwex = boxwex,
                outline = FALSE)
        axis(1, at = tau_centers, labels = tau_list, cex.axis = cex_axis)
        
        # bottom-left: Length
        boxplot(len_obj$vals,
                at = len_obj$pos,
                col = len_obj$cols,
                border = border_len,
                xaxt = "n",
                main = "Length", cex.main = cex_main, 
                ylim = ylim_length,
                boxwex = boxwex,
                outline = FALSE)
        axis(1, at = tau_centers, labels = tau_list, cex.axis = cex_axis)
        
        # bottom-right: Lower
        boxplot(lower_obj$vals,
                at = lower_obj$pos,
                col = lower_obj$cols,
                border = border_lower,
                xaxt = "n",
                main = "Lower", cex.main = cex_main,
                ylim = ylim_lower,
                boxwex = boxwex,
                outline = FALSE)
        # in vertical plot: reference is horizontal line at y=tau
        abline(h = tau_list, lty = 2, col = ref_tau_col)
        axis(1, at = tau_centers, labels = tau_list, cex.axis = cex_axis)
        
        # shared x-label
        mtext(expression("Prediction time " * tau),
              side = 1, outer = TRUE, line = 0.8, cex = cex_lab)
    }
    
    dev.off()
    cat("Saved:", out_pdf, "\n")
}


# ============================================================
# Combine (Coverage + Length) panels from Fig B (top) and Fig A (bottom)
# with horizontal boxplots
# ============================================================
make_horiz_cov_len_twofigs <- function(
        methods, tau_list, alpha,
        # Figure A spec
        figA_order, figA_fill, figA_edge,
        # Figure B spec (plotted on TOP)
        figB_order, figB_fill, figB_edge,
        out_pdf,
        method_label = NULL,
        offsetA = 0.18, offsetB = 0.25,
        boxwex = 0.05,
        ylim_coverage = NULL, ylim_length = NULL,
        alpha_line_col = "grey60",
        plot_cex = list(),
        coverage_col = "coverage_test_true",
        coverage_ticks = NULL,
        width = 10, height = 7
) {
    # ---- font size defaults ----
    cex_legend <- plot_cex$legend %||% 1.0
    cex_axis   <- plot_cex$axis   %||% 1.0
    cex_lab    <- plot_cex$lab    %||% 1.1
    cex_main   <- plot_cex$main   %||% 1.1
    
    tau_centers <- seq_along(tau_list)
    
    # ---- build plotting objects for one figure spec ----
    build_cov_len_objs <- function(method_order, boxfill_cols, boxborder_cols, offset) {
        K <- length(method_order)
        
        res_list <- lapply(method_order, function(nm) methods[[nm]]$res_ok)
        bnd_list <- lapply(method_order, function(nm) methods[[nm]]$bounds)
        
        # IMPORTANT: use index positions, not tau numeric values
        pos <- multi_pos(tau_centers, m = K, offset = offset)
        cols <- rep(boxfill_cols, times = length(tau_list))
        border_cols_full <- rep(boxborder_cols, times = length(tau_list))
        
        # Coverage
        cov_list <- build_cov_list_k(res_list, tau_list, method_order, col = coverage_col)
        
        # Length pooled
        lens <- lapply(bnd_list, pool_bounds_by_tau, tau_list = tau_list, which_set="test", what="len")
        len_list <- interleave_k_by_tau(lens, tau_list, method_order)
        
        # drop empty groups
        cov_obj <- drop_empty_groups(cov_list, pos, cols)
        len_obj <- drop_empty_groups(len_list, pos, cols)
        
        border_cov <- border_cols_full[cov_obj$keep]
        border_len <- border_cols_full[len_obj$keep]
        
        legend_text <- method_order
        if (!is.null(method_label)) {
            legend_text <- unname(method_label[method_order])
            legend_text[is.na(legend_text)] <- method_order[is.na(legend_text)]
        }
        
        list(
            method_order = method_order,
            boxfill_cols = boxfill_cols,
            boxborder_cols = boxborder_cols,
            cov_obj = cov_obj,
            len_obj = len_obj,
            border_cov = border_cov,
            border_len = border_len,
            legend_text = legend_text
        )
    }
    
    # Build TOP = B, BOTTOM = A
    B <- build_cov_len_objs(figB_order, figB_fill, figB_edge, offsetB)
    A <- build_cov_len_objs(figA_order, figA_fill, figA_edge, offsetA)
    
    pdf(out_pdf, width = width, height = height)
    
    par(mfrow = c(2, 2),
        mar = c(2.5, 2.0, 2.3, 1.2),
        oma = c(0, 2.5, 0, 0),
        cex.axis = cex_axis, cex.lab = cex_lab)
    
    # ---------------- TOP ROW: Figure B ----------------
    # (1) Coverage (B)  [TITLE ON]
    boxplot(B$cov_obj$vals,
            at = -B$cov_obj$pos,
            col = B$cov_obj$cols,
            border = B$border_cov,
            horizontal = TRUE,
            yaxt = "n",
            xaxt = if (is.null(coverage_ticks)) "s" else "n",
            main = "Coverage", cex.main = cex_main,
            ylim = ylim_coverage,
            boxwex = boxwex,
            outline = FALSE)
    abline(v = 1 - alpha, lty = 2, col = alpha_line_col)
    if (!is.null(coverage_ticks)) {
        axis(1, at = coverage_ticks, labels = coverage_ticks, cex.axis = cex_axis)
    }
    axis(2, at = -tau_centers, labels = tau_list, las = 1, cex.axis = cex_axis)
    legend("topleft",
           legend = B$legend_text,
           fill = B$boxfill_cols,
           border = B$boxborder_cols,
           bty = "n",
           cex = cex_legend)
    
    # y-axis label for TOP ROW (panel 1 only)
    mtext(expression("Prediction time " * tau),
          side = 2, line = 2.3, cex = cex_lab)
    
    # (2) Length (B)    [TITLE ON]
    boxplot(B$len_obj$vals,
            at = -B$len_obj$pos,
            col = B$len_obj$cols,
            border = B$border_len,
            horizontal = TRUE,
            yaxt = "n",
            main = "Length", cex.main = cex_main,
            ylim = ylim_length,
            boxwex = boxwex,
            outline = FALSE)
    
    # ---------------- BOTTOM ROW: Figure A ----------------
    # (3) Coverage (A)  [TITLE OFF]
    boxplot(A$cov_obj$vals,
            at = -A$cov_obj$pos,
            col = A$cov_obj$cols,
            border = A$border_cov,
            horizontal = TRUE,
            yaxt = "n",
            xaxt = if (is.null(coverage_ticks)) "s" else "n",
            main = "",  # <- no title on bottom row
            ylim = ylim_coverage,
            boxwex = boxwex,
            outline = FALSE)
    abline(v = 1 - alpha, lty = 2, col = alpha_line_col)
    if (!is.null(coverage_ticks)) {
        axis(1, at = coverage_ticks, labels = coverage_ticks, cex.axis = cex_axis)
    }
    axis(2, at = -tau_centers, labels = tau_list, las = 1, cex.axis = cex_axis)
    legend("topleft",
           legend = A$legend_text,
           fill = A$boxfill_cols,
           border = A$boxborder_cols,
           bty = "n",
           cex = cex_legend)
    
    # y-axis label for BOTTOM ROW (panel 3 only)
    mtext(expression("Prediction time " * tau),
          side = 2, line = 2.3, cex = cex_lab)
    
    # (4) Length (A)    [TITLE OFF]
    boxplot(A$len_obj$vals,
            at = -A$len_obj$pos,
            col = A$len_obj$cols,
            border = A$border_len,
            horizontal = TRUE,
            yaxt = "n",
            main = "",  # <- no title on bottom row
            ylim = ylim_length,
            boxwex = boxwex,
            outline = FALSE)
    
    # REMOVE the old shared outer label:
    # mtext(expression("Prediction time " * tau),
    #       side = 2, outer = TRUE, line = 1.2, cex = cex_lab)
    
    
    dev.off()
    cat("Saved:", out_pdf, "\n")
}


# ============================================================
# Make a 1x2 (Coverage | Length) plot for ONE figure spec
# with horizontal boxplots (tau on y-axis)
# ============================================================
make_horiz_cov_len_onefig <- function(
        methods, tau_list, alpha,
        out_pdf,
        method_order,
        boxfill_cols, boxborder_cols,
        method_label = NULL,
        offset = 0.22,
        boxwex = 0.05,
        ylim_coverage = NULL, ylim_length = NULL,
        alpha_line_col = "grey60",
        plot_cex = list(),
        coverage_col = "coverage_test_true",
        coverage_ticks = NULL,
        width = 10, height = 4.5,
        title_prefix = NULL
) {
    # ---- font size defaults ----
    cex_legend <- plot_cex$legend %||% 1.0
    cex_axis   <- plot_cex$axis   %||% 1.0
    cex_lab    <- plot_cex$lab    %||% 1.1
    cex_main   <- plot_cex$main   %||% 1.1
    
    K <- length(method_order)
    stopifnot(length(boxfill_cols) == K, length(boxborder_cols) == K)
    
    tau_centers <- seq_along(tau_list)
    
    res_list <- lapply(method_order, function(nm) methods[[nm]]$res_ok)
    bnd_list <- lapply(method_order, function(nm) methods[[nm]]$bounds)
    
    # positions & colors (use index positions, not numeric tau)
    pos <- multi_pos(tau_centers, m = K, offset = offset)
    cols <- rep(boxfill_cols, times = length(tau_list))
    border_cols_full <- rep(boxborder_cols, times = length(tau_list))
    
    # Coverage + Length lists
    cov_list <- build_cov_list_k(res_list, tau_list, method_order, col = coverage_col)
    lens     <- lapply(bnd_list, pool_bounds_by_tau, tau_list = tau_list, which_set="test", what="len")
    len_list <- interleave_k_by_tau(lens, tau_list, method_order)
    
    # drop empty groups
    cov_obj <- drop_empty_groups(cov_list, pos, cols)
    len_obj <- drop_empty_groups(len_list, pos, cols)
    
    border_cov <- border_cols_full[cov_obj$keep]
    border_len <- border_cols_full[len_obj$keep]
    
    legend_text <- method_order
    if (!is.null(method_label)) {
        legend_text <- unname(method_label[method_order])
        legend_text[is.na(legend_text)] <- method_order[is.na(legend_text)]
    }
    
    main_cov <- "Coverage"
    main_len <- "Length"
    if (!is.null(title_prefix)) {
        main_cov <- paste0(main_cov, " (", title_prefix, ")")
        main_len <- paste0(main_len, " (", title_prefix, ")")
    }
    
    pdf(out_pdf, width = width, height = height)
    
    par(mfrow = c(1, 2),
        mar = c(2.5, 2.0, 2.3, 1.2),
        oma = c(0, 2.5, 0, 0),
        cex.axis = cex_axis, cex.lab = cex_lab)
    
    # ---- Coverage ----
    boxplot(cov_obj$vals,
            at = -cov_obj$pos,
            col = cov_obj$cols,
            border = border_cov,
            horizontal = TRUE,
            yaxt = "n",
            xaxt = if (is.null(coverage_ticks)) "s" else "n",
            main = main_cov, cex.main = cex_main,
            ylim = ylim_coverage,
            boxwex = boxwex,
            outline = FALSE)
    abline(v = 1 - alpha, lty = 2, col = alpha_line_col)
    if (!is.null(coverage_ticks)) {
        axis(1, at = coverage_ticks, labels = coverage_ticks, cex.axis = cex_axis)
    }
    axis(2, at = -tau_centers, labels = tau_list, las = 1, cex.axis = cex_axis)
    
    legend("topleft",
           legend = legend_text,
           fill = boxfill_cols,
           border = boxborder_cols,
           bty = "n",
           cex = cex_legend)
    
    # ---- Length ----
    boxplot(len_obj$vals,
            at = -len_obj$pos,
            col = len_obj$cols,
            border = border_len,
            horizontal = TRUE,
            yaxt = "n",
            main = main_len, cex.main = cex_main,
            ylim = ylim_length,
            boxwex = boxwex,
            outline = FALSE)
    
    # shared y-label
    mtext(expression("Prediction time " * tau),
          side = 2, outer = TRUE, line = 1.2, cex = cex_lab)
    
    dev.off()
    cat("Saved:", out_pdf, "\n")
}







# -------------------- Subset the bounds at tau_list ------------------------------
subset_bounds_to_tau <- function(bounds, tau_list, tol = 1e-12) {
    # bounds: list over reps; each rep is list over taus (dynamic) or named list (static)
    out <- vector("list", length(bounds))
    for (r in seq_along(bounds)) {
        br <- bounds[[r]]
        
        # If already named like "tau_0", "tau_1", ...
        nms <- names(br)
        if (!is.null(nms) && all(grepl("^tau_", nms))) {
            # match by numeric tau parsed from names
            tau_in <- suppressWarnings(as.numeric(sub("^tau_", "", nms)))
            keep <- vapply(tau_list, function(t) any(abs(tau_in - t) < tol), logical(1))
            tau_keep <- tau_list[keep]
            
            # rebuild in tau_list order
            idx <- vapply(tau_keep, function(t) which.min(abs(tau_in - t)), integer(1))
            out[[r]] <- br[idx]
            names(out[[r]]) <- paste0("tau_", tau_keep)
            next
        }
        
        # Otherwise treat as positional list corresponding to tau_grid.
        # In that case we MUST have per-element tau stored inside (common in your structure):
        # bounds[[r]][[k]]$test$tau or bounds[[r]][[k]]$cal$tau or similar.
        tau_in <- rep(NA_real_, length(br))
        for (k in seq_along(br)) {
            bk <- br[[k]]
            tval <- NA_real_
            if (is.list(bk) && !is.null(bk$test) && nrow(bk$test) > 0 && "tau" %in% names(bk$test)) {
                tval <- unique(bk$test$tau)[1]
            } else if (is.list(bk) && !is.null(bk$cal) && nrow(bk$cal) > 0 && "tau" %in% names(bk$cal)) {
                tval <- unique(bk$cal$tau)[1]
            } else if (!is.null(attr(bk, "tau"))) {
                tval <- attr(bk, "tau")
            }
            tau_in[k] <- tval
        }
        
        if (anyNA(tau_in)) {
            stop("subset_bounds_to_tau: cannot infer tau values for rep r=", r,
                 ". Please ensure each bounds[[r]][[k]] carries tau in $test/$cal.")
        }
        
        # select closest matches
        idx <- vapply(tau_list, function(t) {
            j <- which.min(abs(tau_in - t))
            if (abs(tau_in[j] - t) > tol) NA_integer_ else j
        }, integer(1))
        
        if (anyNA(idx)) {
            miss <- tau_list[is.na(idx)]
            stop("subset_bounds_to_tau: these tau values not found in bounds: ",
                 paste(miss, collapse = ", "))
        }
        
        out[[r]] <- br[idx]
        names(out[[r]]) <- paste0("tau_", tau_list)
    }
    out
}
