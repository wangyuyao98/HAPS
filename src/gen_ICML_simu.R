simulate_dataset_long_linear_weibull <- function(
        n, seed = NULL,
        change_times = c(0, 3, 6),   # finite change times
        tau_max = Inf,               # admin cap for observed data
        no_censoring = FALSE,
        par,                         # parameters in the DGM

        # ---- L process: stationary AR(1); iid when rho = 0 ----
        rho = 0,                     # Corr(L_k, L_{k-1})
        L_mean = 0,
        L_sd = 1,
        return_L_full = FALSE,       # opt-in: attach the full n x K covariate path
                                     # matrix as attr(., "L_full"). Default off keeps
                                     # the returned object byte-identical to before;
                                     # used only by the Gtau tilted-evaluation helpers
                                     # (the full path is needed for the tilt normalizer).
        T_max = Inf                  # opt-in bounded event-time support: when finite,
                                     # the event time is drawn from its conditional law
                                     # given T <= T_max (TRUNCATED generation via
                                     # inverse-CDF in the last interval -- no atom at
                                     # T_max, so T stays absolutely continuous and every
                                     # draw satisfies T < T_max strictly). Default Inf
                                     # keeps the generator byte-identical to before.
) {
    if (!is.null(seed)) set.seed(seed)
    
    stopifnot(is.numeric(change_times), all(change_times >= 0))
    stopifnot(all(diff(change_times) > 0))
    stopifnot(is.infinite(tau_max) || tau_max > max(change_times))
    stopifnot(is.finite(rho), abs(rho) <= 1)
    
    # Extend with Inf to define last interval
    change_times_ext <- sort(unique(c(0, change_times, Inf)))
    K <- length(change_times_ext) - 1
    w <- diff(change_times_ext)
    
    # Parameters must be length K
    shape_T <- par$shape_T
    shape_C <- par$shape_C
    beta_T0 <- par$beta_T0
    beta_TL <- par$beta_TL
    beta_C0 <- par$beta_C0
    beta_CL <- par$beta_CL
    
    stopifnot(length(shape_T) == K, length(shape_C) == K,
              length(beta_T0) == K, length(beta_TL) == K,
              length(beta_C0) == K, length(beta_CL) == K)

    stopifnot(is.infinite(T_max) || T_max > max(change_times))
    # Truncation binds only in the last (infinite) interval: T <= T_max is
    # equivalent to (last-interval waiting time) <= T_max - t_K for subjects
    # reaching it; earlier intervals have finite widths below T_max.
    trunc_last_T <- is.finite(T_max)
    w_last_cap <- if (trunc_last_T) T_max - change_times_ext[K] else Inf
    
    # ---- Stationary AR(1) generator for L path ----
    # z_k ~ AR(1) with Var(z_k)=1; L_k = L_mean + L_sd * z_k
    gen_L_path <- function(K) {
        z <- numeric(K)
        
        # stationary initial draw
        z[1] <- rnorm(1)
        
        if (K >= 2) {
            eps <- rnorm(K - 1)
            innov_sd <- sqrt(1 - rho^2)
            for (k in 2:K) {
                z[k] <- rho * z[k - 1] + innov_sd * eps[k - 1]
            }
        }
        
        L_mean + L_sd * z
    }
    
    out_list <- vector("list", n)
    TT <- CC <- X <- rep(NA_real_, n)
    Delta <- rep(NA_integer_, n)
    L_full <- if (isTRUE(return_L_full)) matrix(NA_real_, nrow = n, ncol = K) else NULL

    for (i in seq_len(n)) {

        # ---- Step 1: simulate latent TT and C ----
        TT_i <- 0
        CC_i <- 0
        T_done <- FALSE
        C_done <- FALSE

        L_path <- gen_L_path(K)
        if (isTRUE(return_L_full)) L_full[i, ] <- L_path
        
        for (k in seq_len(K)) {
            Lk <- L_path[k]
            len_k <- w[k]
            
            # Event time
            linpred_T <- beta_T0[k] + beta_TL[k] * Lk
            scale_T   <- exp(-linpred_T / shape_T[k])
            if (trunc_last_T && k == K && !T_done) {
                # truncated Weibull via inverse-CDF: W | W <= w_last_cap
                # (keeps T absolutely continuous; T < T_max strictly a.s.)
                p_cap <- pweibull(w_last_cap, shape = shape_T[k], scale = scale_T)
                T_raw <- qweibull(runif(1) * p_cap, shape = shape_T[k], scale = scale_T)
            } else {
                T_raw <- rweibull(1, shape = shape_T[k], scale = scale_T)
            }
            
            # Censoring time
            linpred_C <- beta_C0[k] + beta_CL[k] * Lk
            scale_C   <- exp(-linpred_C / shape_C[k])
            C_raw     <- rweibull(1, shape = shape_C[k], scale = scale_C)
            
            if (!T_done) {
                if (T_raw <= len_k) {
                    TT_i <- TT_i + T_raw
                    T_done <- TRUE
                } else {
                    TT_i <- TT_i + len_k
                }
            }
            
            if (!C_done) {
                if (C_raw <= len_k) {
                    CC_i <- CC_i + C_raw
                    C_done <- TRUE
                } else {
                    CC_i <- CC_i + len_k
                }
            }
            
            if (T_done && C_done) break
        }
        
        TT[i] <- TT_i
        CC[i] <- CC_i
        
        # ---- Step 2: observed endpoint ----
        CC_obs_i <- if (no_censoring) Inf else CC_i
        X_i <- min(TT_i, CC_obs_i, tau_max)
        X[i] <- X_i
        
        delta_i <- as.integer(TT_i <= CC_obs_i && TT_i <= tau_max)
        Delta[i] <- delta_i
        cens_i <- as.integer(!no_censoring && CC_i < TT_i && CC_i <= tau_max)
        
        # ---- Build long data ----
        rows <- list()
        for (k in seq_len(K)) {
            start_k <- change_times_ext[k]
            stop_k  <- change_times_ext[k + 1]
            
            if (start_k >= X_i) break
            
            end_k <- min(stop_k, X_i)
            
            rows[[length(rows) + 1L]] <- data.frame(
                id = i,
                k = k,
                start = start_k,
                stop = end_k,
                L = L_path[k],
                event = 0L,
                event.C = 0L,
                stringsAsFactors = FALSE
            )
            
            if (end_k >= X_i - 1e-12) break
        }
        
        dat_i <- do.call(rbind, rows)
        dat_i$event[nrow(dat_i)]   <- delta_i
        dat_i$event.C[nrow(dat_i)] <- cens_i
        
        out_list[[i]] <- dat_i
    }
    
    dat_long <- do.call(rbind, out_list)
    
    subj <- data.frame(
        id = seq_len(n),
        TT = TT, C = CC, X = X, delta = Delta,
        stringsAsFactors = FALSE
    )
    
    dat_long <- merge(dat_long, subj, by = "id", all.x = TRUE)
    dat_long <- dat_long[order(dat_long$id, dat_long$start), ]
    rownames(dat_long) <- NULL
    
    # Attributes
    attr(dat_long, "censor_rate") <- mean(subj$delta == 0)
    attr(dat_long, "change_times") <- change_times_ext
    attr(dat_long, "tau_max") <- tau_max
    attr(dat_long, "L_process") <- "stationary AR(1); iid when rho = 0"
    attr(dat_long, "rho") <- rho
    attr(dat_long, "dgm_name") <- "linear_weibull"
    attr(dat_long, "T_max") <- T_max
    if (isTRUE(return_L_full)) {
        rownames(L_full) <- as.character(seq_len(n))
        attr(dat_long, "L_full") <- L_full
    }

    dat_long
}


.simulate_dataset_long_mixC_uniform_cox <- function(
        n, seed = NULL,
        change_times = c(0, 3, 6),
        tau_max = Inf,
        no_censoring = FALSE,
        par,
        rho = 0,
        L_mean = 0,
        L_sd = 1
) {
    if (!is.null(seed)) set.seed(seed)
    
    stopifnot(is.numeric(change_times), all(change_times >= 0))
    stopifnot(all(diff(change_times) > 0))
    stopifnot(is.infinite(tau_max) || tau_max > max(change_times))
    stopifnot(is.finite(rho), abs(rho) <= 1)
    
    change_times_ext <- sort(unique(c(0, change_times, Inf)))
    K <- length(change_times_ext) - 1
    w <- diff(change_times_ext)
    
    shape_T <- par$shape_T
    shape_C <- par$shape_C
    beta_T0 <- par$beta_T0
    beta_TL <- par$beta_TL
    beta_TB <- if (!is.null(par$beta_TB)) par$beta_TB else rep(0, K)
    beta_C0 <- par$beta_C0
    beta_CL <- par$beta_CL
    beta_CL2 <- if (!is.null(par$beta_CL2)) par$beta_CL2 else rep(0, K)
    prob_cox <- if (!is.null(par$prob_cox)) as.numeric(par$prob_cox) else 0.5
    uniform_cmax <- if (!is.null(par$uniform_cmax)) as.numeric(par$uniform_cmax) else (max(change_times) + 2)
    
    stopifnot(length(shape_T) == K, length(shape_C) == K,
              length(beta_T0) == K, length(beta_TL) == K, length(beta_TB) == K,
              length(beta_C0) == K, length(beta_CL) == K,
              length(beta_CL2) == K)
    if (length(prob_cox) != 1L || !is.finite(prob_cox) || prob_cox < 0 || prob_cox > 1) {
        stop("For dgm_name='mixC_uniform_cox', par$prob_cox must be a scalar in [0, 1].")
    }
    if (length(uniform_cmax) != 1L || !is.finite(uniform_cmax) || uniform_cmax <= 0) {
        stop("For dgm_name='mixC_uniform_cox', par$uniform_cmax must be a positive finite scalar.")
    }
    
    gen_L_path <- function(K) {
        z <- numeric(K)
        z[1] <- rnorm(1)
        
        if (K >= 2) {
            eps <- rnorm(K - 1)
            innov_sd <- sqrt(1 - rho^2)
            for (k in 2:K) {
                z[k] <- rho * z[k - 1] + innov_sd * eps[k - 1]
            }
        }
        
        L_mean + L_sd * z
    }
    
    gen_event_time <- function(L_path, B_i) {
        TT_i <- 0
        T_done <- FALSE
        
        for (k in seq_len(K)) {
            Lk <- L_path[k]
            len_k <- w[k]
            linpred_T <- beta_T0[k] + beta_TL[k] * Lk + beta_TB[k] * B_i
            scale_T <- exp(-linpred_T / shape_T[k])
            T_raw <- rweibull(1, shape = shape_T[k], scale = scale_T)
            
            if (T_raw <= len_k) {
                TT_i <- TT_i + T_raw
                T_done <- TRUE
                break
            } else {
                TT_i <- TT_i + len_k
            }
        }
        
        if (!T_done) TT_i else TT_i
    }
    
    gen_censor_time <- function(L_path, B_i) {
        if (B_i == 0L) {
            return(runif(1, min = 0, max = uniform_cmax))
        }
        
        CC_i <- 0
        C_done <- FALSE
        
        for (k in seq_len(K)) {
            Lk <- L_path[k]
            len_k <- w[k]
            linpred_C <- beta_C0[k] + beta_CL[k] * Lk + beta_CL2[k] * (Lk^2)
            scale_C <- exp(-linpred_C / shape_C[k])
            C_raw <- rweibull(1, shape = shape_C[k], scale = scale_C)
            
            if (C_raw <= len_k) {
                CC_i <- CC_i + C_raw
                C_done <- TRUE
                break
            } else {
                CC_i <- CC_i + len_k
            }
        }
        
        if (!C_done) CC_i else CC_i
    }
    
    out_list <- vector("list", n)
    TT <- CC <- X <- rep(NA_real_, n)
    Delta <- B <- rep(NA_integer_, n)
    
    for (i in seq_len(n)) {
        B_i <- rbinom(1, size = 1, prob = prob_cox)
        L_path <- gen_L_path(K)
        TT_i <- gen_event_time(L_path, B_i)
        CC_i <- gen_censor_time(L_path, B_i)
        
        TT[i] <- TT_i
        CC[i] <- CC_i
        B[i] <- B_i
        
        CC_obs_i <- if (no_censoring) Inf else CC_i
        X_i <- min(TT_i, CC_obs_i, tau_max)
        X[i] <- X_i
        
        delta_i <- as.integer(TT_i <= CC_obs_i && TT_i <= tau_max)
        Delta[i] <- delta_i
        cens_i <- as.integer(!no_censoring && CC_i < TT_i && CC_i <= tau_max)
        
        rows <- list()
        for (k in seq_len(K)) {
            start_k <- change_times_ext[k]
            stop_k  <- change_times_ext[k + 1]
            
            if (start_k >= X_i) break
            
            end_k <- min(stop_k, X_i)
            rows[[length(rows) + 1L]] <- data.frame(
                id = i,
                k = k,
                start = start_k,
                stop = end_k,
                B = B_i,
                L = L_path[k],
                event = 0L,
                event.C = 0L,
                stringsAsFactors = FALSE
            )
            
            if (end_k >= X_i - 1e-12) break
        }
        
        dat_i <- do.call(rbind, rows)
        dat_i$event[nrow(dat_i)] <- delta_i
        dat_i$event.C[nrow(dat_i)] <- cens_i
        out_list[[i]] <- dat_i
    }
    
    dat_long <- do.call(rbind, out_list)
    
    subj <- data.frame(
        id = seq_len(n),
        B = B,
        TT = TT,
        C = CC,
        X = X,
        delta = Delta,
        stringsAsFactors = FALSE
    )
    
    dat_long <- merge(dat_long, subj, by = c("id", "B"), all.x = TRUE)
    dat_long <- dat_long[order(dat_long$id, dat_long$start), ]
    rownames(dat_long) <- NULL
    
    attr(dat_long, "censor_rate") <- mean(subj$delta == 0)
    attr(dat_long, "change_times") <- change_times_ext
    attr(dat_long, "tau_max") <- tau_max
    attr(dat_long, "L_process") <- "stationary AR(1); iid when rho = 0"
    attr(dat_long, "rho") <- rho
    attr(dat_long, "dgm_name") <- "mixC_uniform_cox"
    attr(dat_long, "censoring_process") <- "mixture: B=0 uniform censoring, B=1 nonlinear Cox-like Weibull censoring"
    attr(dat_long, "prob_cox") <- prob_cox
    attr(dat_long, "uniform_cmax") <- uniform_cmax
    
    dat_long
}


simulate_dataset_long <- function(
        n, seed = NULL,
        change_times = c(0, 3, 6),   # finite change times
        tau_max = Inf,               # admin cap for observed data
        no_censoring = FALSE,
        par,                         # parameters in the DGM
        
        # ---- L process: stationary AR(1); iid when rho = 0 ----
        rho = 0,                     # Corr(L_k, L_{k-1})
        L_mean = 0,
        L_sd = 1,
        dgm_name = c("linear_weibull", "mixC_uniform_cox"),
        return_L_full = FALSE,
        T_max = Inf
) {
    dgm_name <- match.arg(dgm_name)
    if (isTRUE(return_L_full) && dgm_name != "linear_weibull") {
        stop("simulate_dataset_long(): return_L_full=TRUE is only supported for ",
             "dgm_name='linear_weibull'.")
    }
    if (is.finite(T_max) && dgm_name != "linear_weibull") {
        stop("simulate_dataset_long(): finite T_max is only supported for ",
             "dgm_name='linear_weibull'.")
    }

    switch(
        dgm_name,
        linear_weibull = simulate_dataset_long_linear_weibull(
            n = n,
            seed = seed,
            change_times = change_times,
            tau_max = tau_max,
            no_censoring = no_censoring,
            par = par,
            rho = rho,
            L_mean = L_mean,
            L_sd = L_sd,
            return_L_full = return_L_full,
            T_max = T_max
        ),
        mixC_uniform_cox = .simulate_dataset_long_mixC_uniform_cox(
            n = n,
            seed = seed,
            change_times = change_times,
            tau_max = tau_max,
            no_censoring = no_censoring,
            par = par,
            rho = rho,
            L_mean = L_mean,
            L_sd = L_sd
        )
    )
}
