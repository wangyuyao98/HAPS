## =====================================================================
## Comparison of evaluation strategies for the at-risk coverage estimand
##     beta_delta = P( T in C_hat | T > tau, C_tilde > tau ),
## for one FIXED interval rule, against a near-exact truth from a 50k
## reference population. Estimators per fresh test set (n = 500):
##   (a) est-subgroup : censored test data, at-risk {X>tau}, mean covered_T
##                      [the legacy Gtau_mode="estimated" evaluation]
##   (b) est-weighted : uncensored survivors {T>tau}, weights w0 = G_true(tau)
##                      [gtau_eval_method="weighted", delta = 0]
##   (c) tilt-weighted: survivors {T>tau}, weights w_delta (tilted true law)
##                      [gtau_eval_method="weighted", Gtau_mode="tilted"]
##   (d) tilt-ratio   : at-risk {X>tau}, importance ratio w_delta/w0
##   (e) tilt-sample  : physically draw C_tilde once, subset, plain mean
##
## Expected result (seeds fixed): all ~unbiased; SD ordering
## weighted < sampled < ratio-subgroup; this motivates the "weighted" default
## of gtau_eval_method in the drivers.
## Runtime: ~4-6 minutes. Run from the repository root:
##     Rscript validation/compare_gtau_evaluators.R
## =====================================================================
suppressMessages(library(survival))
if (!file.exists("src/gen_ICML_simu.R")) {
    stop("Run validation/compare_gtau_evaluators.R from the repository root.")
}
source("src/gen_ICML_simu.R"); source("src/helpers.R")
source("src/dynamicCP.R"); source("src/gtau_eval_helpers.R")

ct <- c(3, 6); tau <- 3; delta <- 0.15
par <- list(shape_T = c(4,4,5), shape_C = c(3,3,3), beta_T0 = c(-8,-8,-5),
            beta_TL = c(1,2,3), beta_C0 = c(-6,-6,-5), beta_CL = c(2,2,2))

## ---- fixed HAPS-style interval rule: landmark Cox at tau, quantile band ----
dat_tr <- simulate_dataset_long(1000, seed = 31, change_times = ct, no_censoring = FALSE,
                                par = par, rho = 0.3)
tr_tau <- prepare_data_tau(dat_tr, "start","stop","event","id", tau)
cox_fit <- coxph(Surv(start, stop, event) ~ L, data = tr_tau)
levels_lohi <- c(0.06, 0.94)
make_bounds <- function(dat) {
    dt <- prepare_data_tau(dat, "start","stop","event","id", tau)
    sf <- survfit(cox_fit, newdata = dt, se.fit = FALSE, conf.type = "none")
    Q <- surv_quantile_grid(sf$time, t(sf$surv), levels_lohi)
    data.frame(id = dt$id, lo = Q[,1], hi = Q[,2], TT = dt$TT, X = dt$stop)
}
cov_ind <- function(b) as.integer(b$lo <= b$TT & b$TT <= b$hi)

## ---- truth from a 50k reference population (uncensored; direct weighting) ---
ref_u <- simulate_dataset_long(50000, seed = 777, change_times = ct, no_censoring = TRUE,
                               par = par, rho = 0.3, return_L_full = TRUE)
Lref <- attr(ref_u, "L_full"); rownames(Lref) <- as.character(seq_len(50000))
bru <- make_bounds(ref_u)                     # rows: T>tau survivors
cvu <- cov_ind(bru)
ids_u <- as.character(bru$id)
w0_ref <- compute_Gtau_true_linWB1(Lref[ids_u,,drop=FALSE], tau, ct, par)
wd_ref <- compute_Gtau_true_tilted_linWB1(Lref[ids_u,,drop=FALSE], tau, delta, ct, par)
truth0 <- sum(w0_ref*cvu)/sum(w0_ref)
truthd <- sum(wd_ref*cvu)/sum(wd_ref)
cat(sprintf("TRUTH (50k ref): beta_0 = %.4f   beta_%.2f = %.4f   (unweighted surv: %.4f)\n",
            truth0, delta, truthd, mean(cvu)))

## ---- tilted C_tilde sampler for (e): inverse CDF from raw DGM pieces ------
sample_Ct <- function(Lrow) {
    tg <- seq(1e-4, 40, length.out = 1500)
    dens <- numeric(length(tg)); pre <- 1; ce <- c(0, ct, Inf)
    for (k in 1:3) {
        sig <- exp(-(par$beta_C0[k] + par$beta_CL[k]*Lrow[k])/par$shape_C[k])
        inb <- tg > ce[k] & tg <= ce[k+1]
        dens[inb] <- pre * dweibull(tg[inb]-ce[k], par$shape_C[k], sig)
        if (is.finite(ce[k+1])) pre <- pre * pweibull(ce[k+1]-ce[k], par$shape_C[k], sig, lower.tail=FALSE)
    }
    dens <- dens * exp(delta * tg)
    cdf <- cumsum(dens); cdf <- cdf/cdf[length(cdf)]
    tg[findInterval(runif(1), cdf) + 1L]
}

## ---- M fresh test sets of n=500 --------------------------------------------
M <- 300; n_te <- 500
est <- matrix(NA_real_, M, 5, dimnames = list(NULL, c("a_est_subgrp","b_est_wt","c_tilt_wt","d_tilt_ratio","e_tilt_sample")))
n_surv <- n_risk <- numeric(M)
set.seed(4242)
for (m in seq_len(M)) {
    sd_m <- 20000 + m
    dc <- simulate_dataset_long(n_te, seed = sd_m, change_times = ct, no_censoring = FALSE,
                                par = par, rho = 0.3, return_L_full = TRUE)
    Lm <- attr(dc, "L_full"); rownames(Lm) <- as.character(seq_len(n_te))
    du <- simulate_dataset_long(n_te, seed = sd_m, change_times = ct, no_censoring = TRUE,
                                par = par, rho = 0.3, return_L_full = TRUE)
    # same seed => same latent (T, path); dc additionally realizes C via X
    bu <- make_bounds(du); cu <- cov_ind(bu); idu <- as.character(bu$id)       # survivors T>tau
    bc <- make_bounds(dc); cc <- cov_ind(bc); idc <- as.character(bc$id)       # at-risk X>tau
    n_surv[m] <- length(idu); n_risk[m] <- length(idc)
    w0_u <- compute_Gtau_true_linWB1(Lm[idu,,drop=FALSE], tau, ct, par)
    wd_u <- compute_Gtau_true_tilted_linWB1(Lm[idu,,drop=FALSE], tau, delta, ct, par)
    w0_c <- compute_Gtau_true_linWB1(Lm[idc,,drop=FALSE], tau, ct, par)
    wd_c <- compute_Gtau_true_tilted_linWB1(Lm[idc,,drop=FALSE], tau, delta, ct, par)
    est[m,"a_est_subgrp"] <- mean(cc)
    est[m,"b_est_wt"]     <- sum(w0_u*cu)/sum(w0_u)
    est[m,"c_tilt_wt"]    <- sum(wd_u*cu)/sum(wd_u)
    r <- wd_c / pmax(w0_c, 1e-12)
    est[m,"d_tilt_ratio"] <- sum(r*cc)/sum(r)
    keep <- vapply(idu, function(i) sample_Ct(Lm[i,]) > tau, logical(1))
    est[m,"e_tilt_sample"] <- if (any(keep)) mean(cu[keep]) else NA_real_
}

## ---- report -----------------------------------------------------------------
tgt <- c(truth0, truth0, truthd, truthd, truthd)
cat(sprintf("\n%-14s %-22s %8s %8s %8s\n", "estimator", "target", "bias", "SD", "RMSE"))
for (j in 1:5) {
    e <- est[,j]; b <- mean(e, na.rm=TRUE) - tgt[j]; s <- sd(e, na.rm=TRUE)
    cat(sprintf("%-14s %-22s %+8.4f %8.4f %8.4f\n",
        colnames(est)[j],
        ifelse(j<=2, sprintf("beta_0=%.4f",truth0), sprintf("beta_%.2f=%.4f",delta,truthd)),
        b, s, sqrt(b^2+s^2)))
}
cat(sprintf("\navg survivors/test set (T>tau): %.0f | avg at-risk subgroup (X>tau): %.0f\n",
            mean(n_surv), mean(n_risk)))
