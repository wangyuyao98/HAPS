## =====================================================================
## Validation: the weighted evaluation of the at-risk estimand equals the
## physical-subgroup evaluation with C_tilde DRAWN from the tilted TRUE law.
##
## Estimand (for a fixed interval rule C_hat and a specified censoring
## mechanism C_tilde independent of T given the covariate path):
##     beta_delta = P( T in C_hat | T > tau, C_tilde > tau ).
## Weighted estimator (used by the pipeline): on uncensored test survivors
## {T > tau}, sum_i w_i 1{T_i in C_hat_i} / sum_i w_i with
## w_i = G_tilde_delta(tau | path_i).
##
## This script draws FULL C_tilde trajectories by inverse-CDF from the raw DGM
## pieces (prefix survival x interval Weibull density x exp(delta t),
## normalized) -- the sampler never touches w_i -- and shows the mean of the
## physical-subgroup coverages matches the weighted estimate.
##
## Expected output (seeds fixed): |weighted - sampled| ~ 1e-4, with ~95% of
## subgroup replicates within +/- 2 SD of the weighted value.
## Runtime: ~1-2 minutes. Run from the repository root:
##     Rscript validation/validate_gtau_weighted_evaluation.R
## =====================================================================
suppressMessages(library(survival))
if (!file.exists("src/gen_ICML_simu.R")) {
    stop("Run validation/validate_gtau_weighted_evaluation.R from the repository root.")
}
source("src/gen_ICML_simu.R"); source("src/helpers.R"); source("src/helpers_AIPCW.R")
source("src/dynamicCP.R"); source("src/dynamicCP_AIPCW.R"); source("src/gtau_eval_helpers.R")

ct <- c(3, 6); tau <- 3; delta <- 0.15
par <- list(shape_T = c(4,4,5), shape_C = c(3,3,3), beta_T0 = c(-8,-8,-5),
            beta_TL = c(1,2,3), beta_C0 = c(-6,-6,-5), beta_CL = c(2,2,2))

## 1) real pipeline: train (censored) + test (uncensored, full path), one HAPS-A fit
dat_tr <- simulate_dataset_long(1000, seed = 31, change_times = ct, no_censoring = FALSE,
                                par = par, rho = 0.3)
dat_te <- simulate_dataset_long(4000, seed = 32, change_times = ct, no_censoring = TRUE,
                                par = par, rho = 0.3, return_L_full = TRUE)
L <- attr(dat_te, "L_full"); rownames(L) <- as.character(seq_len(4000))
TT <- setNames(dat_te$TT[match(seq_len(4000), dat_te$id)], as.character(seq_len(4000)))

fit <- dynamicCP_AIPCW_split(
    dat = dat_tr, dat_test = dat_te, tau = tau,
    start.name = "start", stop.name = "stop", event.name = "event",
    event.C.name = "event.C", id.name = "id",
    covname.pred.timevarying = "L", model.pred = "cox",
    model.C = "cox", covname.C.timevarying = "L",
    model.S = "cox", model.Xi = "xgb_reg",
    theta_grid = seq(0, 0.5, 0.01), alpha = 0.1, trim.C = 0.05,
    TT.name = "TT", visit_times = c(0, ct), seed = 77, train_ratio = 0.5,
    localization = FALSE, censoring_method = "piecewise",
    Gtau_mode = "tilted", Gtau_delta = delta)

ids <- as.character(fit$id_test_tau)
TTs <- TT[ids]; keep <- TTs > tau
ids <- ids[keep]
covred <- as.integer(fit$lower_AIPCW_test[keep] <= TTs[keep] &
                     TTs[keep] <= fit$upper_AIPCW_test[keep])
Ls <- L[ids, , drop = FALSE]
n_s <- length(ids)
cat(sprintf("survivors at tau=%g: %d;  raw (unweighted) coverage on {T>tau}: %.4f\n",
            tau, n_s, mean(covred)))

## 2) weighted estimate (the implementation under validation)
w <- compute_Gtau_true_tilted_linWB1(Ls, tau, delta, ct, par)
w_est <- sum(w * covred) / sum(w)
cat(sprintf("weighted estimate  (sum w*cov / sum w):            %.4f\n", w_est))

## 3) physical sampling of full C_tilde from the tilted TRUE law (identity-free)
set.seed(99)
B <- 400
tgrid <- seq(1e-4, 40, length.out = 3000)
ct_ext <- c(0, ct, Inf)
num_b <- den_b <- numeric(B)
for (i in seq_len(n_s)) {
    dens <- numeric(length(tgrid)); prefix <- 1
    for (k in 1:3) {
        sig <- exp(-(par$beta_C0[k] + par$beta_CL[k] * Ls[i, k]) / par$shape_C[k])
        lo <- ct_ext[k]; hi <- ct_ext[k + 1]
        inb <- tgrid > lo & tgrid <= hi
        dens[inb] <- prefix * dweibull(tgrid[inb] - lo, par$shape_C[k], sig)
        if (is.finite(hi)) prefix <- prefix * pweibull(hi - lo, par$shape_C[k], sig, lower.tail = FALSE)
    }
    dens <- dens * exp(delta * tgrid)                    # exponential tilt
    cdf <- cumsum(dens); cdf <- cdf / cdf[length(cdf)]
    Ct <- tgrid[findInterval(runif(B), cdf) + 1L]        # B full C_tilde draws
    at_risk <- Ct > tau                                  # subgroup membership only
    num_b <- num_b + at_risk * covred[i]
    den_b <- den_b + at_risk
}
mc <- num_b / den_b
cat(sprintf("physical C_tilde sampling, B=%d replicates:        %.4f  (MC SE %.4f)\n",
            B, mean(mc), sd(mc) / sqrt(B)))
cat(sprintf("|weighted - sampled| = %.5f\n", abs(w_est - mean(mc))))
cat(sprintf("share of subgroup replicates within +/-2SD of weighted: %.2f\n",
            mean(abs(mc - w_est) <= 2 * sd(mc))))
ok <- abs(w_est - mean(mc)) < 4 * sd(mc) / sqrt(B)
cat(if (ok) "VALIDATION: PASS\n" else "VALIDATION: FAIL (difference exceeds 4x MC SE)\n")
