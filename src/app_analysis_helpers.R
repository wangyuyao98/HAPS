`%||%` <- function(x, y) if (is.null(x)) y else x

app_dataset_config <- function(dataset_name, alpha = 0.2) {
    dataset_name <- match.arg(dataset_name, c("pbcseq", "colon"))
    alpha <- as.numeric(alpha)

    if (identical(dataset_name, "pbcseq")) {
        list(
            dataset_name = "pbcseq",
            out_setup = "app_pbcseq",
            n_all = 312L,
            covname.pred.baseline = c("age"),
            covname.pred.timevarying = c("log_bili", "log_protime", "log_albumin", "edema"),
            covname.C.baseline = c("age"),
            covname.C.timevarying = c("log_bili", "log_protime", "log_albumin", "edema"),
            tau_list_plot = seq(0, 5, by = 1),
            ylim_coverage = c(0.25, 1),
            ylim_length = c(0, 14.3),
            ylim_lower = c(0, 14.3),
            ylim_upper = if (alpha == 0.1) c(10, 14.3) else c(7, 14.3)
        )
    } else {
        list(
            dataset_name = "colon",
            out_setup = "app_colon",
            n_all = 929L,
            covname.pred.baseline = c("rx", "age", "sex", "node4"),
            covname.pred.timevarying = c("recur"),
            covname.C.baseline = c("rx", "age", "sex", "node4"),
            covname.C.timevarying = c("recur"),
            tau_list_plot = seq(0, 5, by = 1),
            ylim_coverage = c(0.5, 1),
            ylim_length = c(0, 8.5),
            ylim_lower = c(0, 9),
            ylim_upper = c(8.5, 9.2)
        )
    }
}

prepare_app_data <- function(dataset_name) {
    dataset_name <- match.arg(dataset_name, c("pbcseq", "colon"))

    if (identical(dataset_name, "pbcseq")) {
        data(pbc, package = "survival")
        dat0 <- survival::pbcseq
        dat0$futime <- dat0$futime / 365.25
        dat0$day <- dat0$day / 365.25

        first <- with(dat0, c(TRUE, diff(id) != 0))
        last <- c(first[-1], TRUE)

        dat0$start <- with(dat0, ifelse(first, 0, day))
        dat0$stop <- with(dat0, ifelse(last, futime, c(day[-1], 0)))
        dat0$event <- with(dat0, ifelse(last & status == 2, 1, 0))
        dat0$event.C <- 0L
        for (i in unique(dat0$id)) {
            idx <- which(dat0$id == i)
            if (dat0$event[max(idx)] == 0) dat0$event.C[max(idx)] <- 1L
        }

        dat0 <- as.data.frame(dat0)
        dat0$log_bili <- log(dat0$bili)
        dat0$log_protime <- log(dat0$protime)
        dat0$log_albumin <- log(dat0$albumin)
        return(dat0)
    }

    data(cancer, package = "survival")
    dat_raw <- survival::colon
    dat_raw$time <- dat_raw$time / 365.25

    by_id <- dat_raw |>
        dplyr::group_by(id) |>
        dplyr::summarise(
            t_rec = suppressWarnings(min(time[etype == 1 & status == 1], na.rm = TRUE)),
            has_rec = any(etype == 1 & status == 1),
            t_death = time[etype == 2][1],
            event = status[etype == 2][1],
            rx = rx[etype == 2][1],
            age = age[etype == 2][1],
            sex = sex[etype == 2][1],
            nodes = nodes[etype == 2][1],
            node4 = node4[etype == 2][1],
            surg = surg[etype == 2][1],
            .groups = "drop"
        ) |>
        dplyr::mutate(t_rec = ifelse(is.infinite(t_rec), NA_real_, t_rec))

    colon_death <- by_id |>
        dplyr::rowwise() |>
        dplyr::do({
            s <- .
            if (!isTRUE(s$has_rec) || is.na(s$t_rec) || s$t_rec >= s$t_death) {
                data.frame(
                    id = s$id, start = 0, stop = s$t_death, recur = 0, event = s$event,
                    rx = s$rx, age = s$age, sex = s$sex, nodes = s$nodes,
                    node4 = s$node4, surg = s$surg
                )
            } else {
                rbind(
                    data.frame(
                        id = s$id, start = 0, stop = s$t_rec, recur = 0, event = 0,
                        rx = s$rx, age = s$age, sex = s$sex, nodes = s$nodes,
                        node4 = s$node4, surg = s$surg
                    ),
                    data.frame(
                        id = s$id, start = s$t_rec, stop = s$t_death,
                        recur = 1, event = s$event,
                        rx = s$rx, age = s$age, sex = s$sex, nodes = s$nodes,
                        node4 = s$node4, surg = s$surg
                    )
                )
            }
        }) |>
        dplyr::ungroup() |>
        as.data.frame()

    colon_death$event.C <- 0L
    for (i in unique(colon_death$id)) {
        idx <- which(colon_death$id == i)
        if (colon_death$event[max(idx)] == 0) colon_death$event.C[max(idx)] <- 1L
    }

    colon_death
}

app_base_config <- function(cfg, alpha = 0.2) {
    list(
        setup = cfg$out_setup,
        out_setup = cfg$out_setup,
        dataset_name = cfg$dataset_name,
        n_all = cfg$n_all,
        covname.pred.baseline = cfg$covname.pred.baseline,
        covname.pred.timevarying = cfg$covname.pred.timevarying,
        covname.C.baseline = cfg$covname.C.baseline,
        covname.C.timevarying = cfg$covname.C.timevarying,
        start.name = "start",
        stop.name = "stop",
        event.name = "event",
        event.C.name = "event.C",
        id.name = "id",
        train_ratio = 0.5,
        test_ratio = 0.2,
        tau_grid = cfg$tau_list_plot,
        theta_grid = seq(0, 1, by = 0.01),
        alpha = alpha,
        trim.C = 0.05,
        TT.name = NULL
    )
}

app_make_split_seeds <- function(R, seed = 123L) {
    R <- as.integer(R)
    if (!is.finite(R) || R < 1L) stop("R must be a positive integer.")

    R_pool <- max(R, 200L)
    set.seed(seed)
    seeds_test <- sample(10^7, R_pool, replace = FALSE)
    seeds_split <- sample(10^7, R_pool, replace = FALSE)
    list(
        seed = seed,
        seeds_test = seeds_test[seq_len(R)],
        seeds_split = seeds_split[seq_len(R)]
    )
}

app_base_run <- function(cfg, alpha = 0.2, R = 200L, seed = 123L) {
    list(
        path = "generated_from_survival_package_data",
        config = app_base_config(cfg, alpha = alpha),
        seeds = app_make_split_seeds(R, seed = seed),
        R_full = as.integer(R)
    )
}

app_out_dir <- function(cfg) {
    out_dir <- file.path("results", cfg$out_setup)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_dir
}

empty_app_bounds <- function() {
    list(test = data.frame(), cal = data.frame())
}

subject_summary <- function(dat, id.name, stop.name, event.name, event.C.name) {
    ids <- unique(dat[[id.name]])
    out <- lapply(ids, function(id) {
        di <- dat[dat[[id.name]] == id, , drop = FALSE]
        data.frame(
            id_tmp = id,
            X = max(di[[stop.name]], na.rm = TRUE),
            Delta = as.integer(any(di[[event.name]] == 1L, na.rm = TRUE)),
            DeltaC = as.integer(any(di[[event.C.name]] == 1L, na.rm = TRUE))
        )
    })
    out <- do.call(rbind, out)
    names(out)[names(out) == "id_tmp"] <- id.name
    out
}

subject_ids_at_risk_app <- function(dat, tau, cfg) {
    summ <- subject_summary(dat, cfg$id.name, cfg$stop.name, cfg$event.name, cfg$event.C.name)
    as.character(summ[[cfg$id.name]][summ$X > tau])
}

make_landmark_residual_data_app <- function(dat, tau, cfg) {
    dat_tau <- prepare_data_tau(dat, cfg$start.name, cfg$stop.name, cfg$event.name, cfg$id.name, tau)
    if (nrow(dat_tau) == 0L) return(dat_tau)

    summ <- subject_summary(dat, cfg$id.name, cfg$stop.name, cfg$event.name, cfg$event.C.name)
    idx <- match(as.character(dat_tau[[cfg$id.name]]), as.character(summ[[cfg$id.name]]))
    if (anyNA(idx)) stop("Failed to align landmark rows to subject summary.")

    dat_tau[[cfg$start.name]] <- 0
    dat_tau[[cfg$stop.name]] <- pmax(summ$X[idx] - tau, 0)
    dat_tau[[cfg$event.name]] <- summ$Delta[idx]
    dat_tau[[cfg$event.C.name]] <- summ$DeltaC[idx]
    dat_tau
}

make_app_bounds_from_fit <- function(dat_tau, lower, upper, survC, rep_id, tau, cfg) {
    w_ipcw <- dat_tau[[cfg$event.name]] / pmax(survC, cfg$trim.C)
    covered_X <- as.integer(lower <= dat_tau[[cfg$stop.name]] & dat_tau[[cfg$stop.name]] <= upper)
    data.frame(
        rep = rep_id,
        tau = tau,
        id = dat_tau[[cfg$id.name]],
        lower = lower,
        upper = upper,
        X = dat_tau[[cfg$stop.name]],
        w_ipcw = w_ipcw,
        covered_X = covered_X
    )
}

summarize_app_bounds <- function(bt, bc) {
    len_test <- bt$upper - bt$lower
    len_cal <- bc$upper - bc$lower
    data.frame(
        coverage_cal_ipcw = mean(bc$w_ipcw * bc$covered_X, na.rm = TRUE) / mean(bc$w_ipcw, na.rm = TRUE),
        coverage_test_ipcw = mean(bt$w_ipcw * bt$covered_X, na.rm = TRUE) / mean(bt$w_ipcw, na.rm = TRUE),
        mean_len_cal = mean(len_cal, na.rm = TRUE),
        mean_len_test = mean(len_test, na.rm = TRUE),
        q25_len_test = as.numeric(stats::quantile(len_test, 0.25, na.rm = TRUE)),
        q50_len_test = as.numeric(stats::quantile(len_test, 0.50, na.rm = TRUE)),
        q75_len_test = as.numeric(stats::quantile(len_test, 0.75, na.rm = TRUE))
    )
}
