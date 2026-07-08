## Evaluation-side helpers for Gtau_mode = "estimated" / "tilted" simulations.
##
## These compute the ORACLE (true-DGM-law) censoring survival weight
##   w_i = \tilde G_delta(tau | path_i) = P_delta(\tilde C_i > tau | path_i)
## used to evaluate coverage on the at-risk population {T > tau, \tilde C > tau}
## by reweighting the uncensored test subjects that survived to tau:
##   coverage = sum_i w_i 1{T_i in C_hat_i} / sum_i w_i    over {T_i > tau}.
##
## The censoring DGM (linWB1 / "linear_weibull", see gen_ICML_simu.R) is an
## interval-wise Weibull renewal process: on interval k = (t_k, t_{k+1}] the
## within-interval censoring waiting time is Weibull(shape a_{C,k}, scale
## sigma_{C,k}(L_k)) with sigma_{C,k} = exp(-(beta_C0[k] + beta_CL[k] L_k)/a_{C,k}).
## C crosses in interval k iff it did not cross in intervals 1..k-1.
##
## delta is the exponential tilt: the censoring density is reweighted by exp(delta t)
## and renormalized. delta = 0 recovers the untilted true law (closed form below).
## These mirror get_Gtau_tilted() inside dynamicCP_AIPCW_split(), but use the exact
## continuous true Weibull law instead of a fitted/discretized G.

# ---- shared parameter extraction --------------------------------------------
.gtau_linWB1_pars <- function(change_times, par) {
    ct <- sort(unique(c(0, change_times)))
    ct_ext <- c(ct, Inf)
    K <- length(ct_ext) - 1L          # number of intervals
    w <- diff(ct_ext)                 # interval widths (last = Inf)
    a  <- par$shape_C
    b0 <- par$beta_C0
    bL <- par$beta_CL
    if (length(a) != K || length(b0) != K || length(bL) != K) {
        stop(".gtau_linWB1_pars(): par$shape_C/beta_C0/beta_CL must have length K = ",
             K, " (one per interval).")
    }
    list(ct = ct, ct_ext = ct_ext, K = K, w = w, a = a, b0 = b0, bL = bL)
}

# subject x interval Weibull censoring scale sigma_{C,k}(L_{i,k})
.gtau_sigma_mat <- function(L, p) {
    # L: n x K covariate-path matrix; returns n x K scales
    n <- nrow(L)
    sig <- matrix(NA_real_, n, p$K)
    for (k in seq_len(p$K)) {
        lin <- p$b0[k] + p$bL[k] * L[, k]
        sig[, k] <- exp(-lin / p$a[k])
    }
    sig
}

#' True untilted censoring survival at tau: G(tau | path) = P(C > tau | path).
#'
#' tau must be one of the change times t_nu. Returns w_i = prod_{k < nu}
#' S_{C,k}(width_k | L_{i,k}) = exp(-sum_{k<nu} (w_k / sigma_{C,k})^{a_k}).
#' Needs only L_{.,1:(nu-1)} (available for survivors on standard uncensored data).
#'
#' @param L n x K covariate-path matrix (attr "L_full" from the generator).
#' @param tau prediction time (a change time).
#' @param change_times finite change-time vector (e.g. c(0,3,6)).
#' @param par DGM parameter list with shape_C, beta_C0, beta_CL.
#' @return numeric vector of length nrow(L) in [0, 1].
compute_Gtau_true_linWB1 <- function(L, tau, change_times, par) {
    L <- as.matrix(L)
    p <- .gtau_linWB1_pars(change_times, par)
    nu <- match(tau, p$ct)
    if (is.na(nu)) stop("compute_Gtau_true_linWB1(): tau must be one of the change times.")
    if (ncol(L) < p$K) {
        stop("compute_Gtau_true_linWB1(): L must have K = ", p$K, " columns (full path).")
    }
    if (nu == 1L) return(rep(1.0, nrow(L)))       # everyone survives censoring to t_1 = 0
    sig <- .gtau_sigma_mat(L, p)
    cumhaz <- rep(0.0, nrow(L))
    for (k in seq_len(nu - 1L)) {                 # intervals strictly before tau
        cumhaz <- cumhaz + (p$w[k] / sig[, k])^p$a[k]
    }
    exp(-cumhaz)
}

#' True EXPONENTIALLY TILTED censoring survival at tau under the linWB1 DGM.
#'
#'   w_i = [ \int_{t>tau} e^{delta t} dF_C(t|path_i) ] / [ \int_0^\infty e^{delta t} dF_C(t|path_i) ]
#'
#' where F_C is the interval-wise-Weibull renewal censoring CDF along the FULL
#' path. delta = 0 returns the closed form above (exactly). The tilt normalizer
#' integrates over ALL intervals, so the full n x K path is required.
#'
#' @param L n x K FULL covariate-path matrix.
#' @param tau prediction time (a change time t_nu).
#' @param delta tilt parameter (delta = 0 -> untilted true law).
#' @param change_times finite change-time vector.
#' @param par DGM parameter list (shape_C, beta_C0, beta_CL).
#' @param n_grid grid points per finite interval for trapezoidal integration.
#' @param tail_tol max acceptable un-captured Weibull tail mass on the last
#'   (infinite) interval; a larger cap is used until below this (warns if not met).
#' @return numeric vector of length nrow(L) in [0, 1].
compute_Gtau_true_tilted_linWB1 <- function(L, tau, delta, change_times, par,
                                            n_grid = 2000L, tail_tol = 1e-8) {
    if (abs(delta) <= 1e-12) {
        return(compute_Gtau_true_linWB1(L, tau, change_times, par))
    }
    L <- as.matrix(L)
    p <- .gtau_linWB1_pars(change_times, par)
    nu <- match(tau, p$ct)
    if (is.na(nu)) stop("compute_Gtau_true_tilted_linWB1(): tau must be one of the change times.")
    if (ncol(L) < p$K) {
        stop("compute_Gtau_true_tilted_linWB1(): L must have K = ", p$K,
             " columns (full path required for the tilt normalizer).")
    }
    n <- nrow(L)
    sig <- .gtau_sigma_mat(L, p)

    # prefix_{i,k} = P(no censoring crossing before interval k) = prod_{j<k} S_j(w_j)
    prefix <- matrix(1.0, n, p$K)
    if (p$K >= 2L) {
        for (k in 2:p$K) {
            prefix[, k] <- prefix[, k - 1L] *
                exp(-(p$w[k - 1L] / sig[, k - 1L])^p$a[k - 1L])
        }
    }

    # I_{i,k} = \int_0^{w_k} e^{delta s} f_k(s | L_{i,k}) ds  (relative time s)
    # Interval contribution to \int e^{delta t} dF is prefix_{i,k} e^{delta t_k} I_{i,k}.
    # Factor out e^{delta * t_ref} (t_ref = max change time considered) -> cancels in ratio.
    t_ref <- max(p$ct[nu], p$ct[p$K])
    I_k <- matrix(0.0, n, p$K)
    for (k in seq_len(p$K)) {
        upper <- if (is.finite(p$w[k])) p$w[k] else {
            # infinite last interval: integrate to a cap where the Weibull tail is
            # negligible. Set the cap well past the tail_tol quantile so the residual
            # mass sits comfortably below tail_tol (avoids a boundary-triggered warning).
            cap <- max(stats::qweibull(1 - tail_tol * 1e-3, shape = p$a[k], scale = sig[, k]))
            tail_mass <- max(stats::pweibull(cap, shape = p$a[k], scale = sig[, k],
                                             lower.tail = FALSE))
            if (tail_mass > tail_tol) {
                warning(sprintf(
                    "compute_Gtau_true_tilted_linWB1(): last-interval tail mass %.2e exceeds tail_tol at delta=%.3g.",
                    tail_mass, delta))
            }
            cap
        }
        s <- seq(0, upper, length.out = n_grid)
        ds <- s[2] - s[1]
        # integrand matrix n x n_grid: e^{delta s} * dweibull(s; a_k, sig_{.,k})
        # weight e^{delta s} shared across subjects; density varies by subject scale.
        ew <- exp(delta * s)
        # trapezoidal weights
        tw <- rep(ds, n_grid); tw[1] <- ds / 2; tw[n_grid] <- ds / 2
        ewt <- ew * tw
        for (i in seq_len(n)) {
            fi <- stats::dweibull(s, shape = p$a[k], scale = sig[i, k])
            I_k[i, k] <- sum(ewt * fi)
        }
        I_k[, k] <- I_k[, k] * exp(delta * (p$ct_ext[k] - t_ref))  # e^{delta t_k}, ref-shifted
    }

    contrib <- prefix * I_k                 # n x K, up to the shared e^{delta t_ref}
    denom <- rowSums(contrib)
    num   <- rowSums(contrib[, nu:p$K, drop = FALSE])
    w <- num / denom
    w[!is.finite(w)] <- compute_Gtau_true_linWB1(L, tau, change_times, par)[!is.finite(w)]
    pmin(pmax(w, 0), 1)
}

#' Weighted survivor-conditional coverage and length under a censoring weight.
#'
#' @param covered integer/logical 1{T in C_hat} on the survivor rows (T > tau).
#' @param len interval lengths on the same rows.
#' @param w censoring-survival weights on the same rows.
#' @return list(coverage, med_len, q25_len, q75_len, ess).
weighted_coverage_summary <- function(covered, len, w) {
    ok <- is.finite(w) & w >= 0 & is.finite(covered)
    covered <- covered[ok]; w <- w[ok]; len <- len[ok]
    sw <- sum(w)
    if (sw <= 0) {
        return(list(coverage = NA_real_, med_len = NA_real_,
                    q25_len = NA_real_, q75_len = NA_real_, ess = 0))
    }
    wq <- function(x, prob) {
        o <- order(x); x <- x[o]; ww <- w[o]
        cw <- cumsum(ww) / sw
        x[which(cw >= prob)[1]]
    }
    list(
        coverage = sum(w * covered) / sw,
        med_len  = wq(len, 0.5),
        q25_len  = wq(len, 0.25),
        q75_len  = wq(len, 0.75),
        ess      = sw^2 / sum(w^2)
    )
}
