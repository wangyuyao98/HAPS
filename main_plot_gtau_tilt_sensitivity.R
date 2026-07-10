rm(list = ls())

## Plot the Gtau tilted-censoring sensitivity study produced by
## main_simu_gtau_tilt_sensitivity.R. For each tau, plots mean coverage (and
## median length) vs the evaluation tilt delta_eval, with one curve per
## calibration arm:
##   - "one"                 (paper method; ignores the shift)
##   - "estimated"           (tilt family, delta_cal = 0)
##   - "matched" tilt        (delta_cal == delta_eval; the ideal oracle target)
##   - fixed-delta_cal tilts (mismatch robustness), one curve per delta_cal.
## Usage: Rscript main_plot_gtau_tilt_sensitivity.R <R> <n> <n_test> [alpha] [setup] [results_root] [models]
## setup in {linWB1 (default), linWB2} selects <results_root>/<setup>/gtau_tilt/.
## results_root defaults to "results"; pass e.g. results/osg/gtau_grid/collected
## to plot OSG-collected files without moving them. models (default "cox") is
## the model-case tag: non-default cases read/write "_<models>"-tagged files.

if (!file.exists("src/gen_ICML_simu.R")) {
    stop("Please run main_plot_gtau_tilt_sensitivity.R from the repository root directory.")
}
suppressMessages({ library(ggplot2); library(dplyr) })

args  <- commandArgs(trailingOnly = TRUE)
R      <- if (length(args) >= 1L) as.integer(args[[1L]]) else 200L
n      <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1000L
n_test <- if (length(args) >= 3L) as.integer(args[[3L]]) else 1000L
alpha  <- if (length(args) >= 4L) as.numeric(args[[4L]]) else 0.1
setup  <- if (length(args) >= 5L && nzchar(args[[5L]])) args[[5L]] else "linWB1"
if (!setup %in% c("linWB1", "linWB2")) stop("setup must be 'linWB1' or 'linWB2'.")
results_root <- if (length(args) >= 6L && nzchar(args[[6L]])) args[[6L]] else "results"
models <- if (length(args) >= 7L && nzchar(args[[7L]])) args[[7L]] else "cox"
mtag   <- if (identical(models, "cox")) "" else paste0("_", models)

folder  <- file.path(results_root, setup, "gtau_tilt")
infile  <- file.path(folder, sprintf("gtau_tilt_sensitivity_R%d_n%d_ntest%d_alpha%s%s.rds",
                                     R, n, n_test, format(alpha), mtag))
if (!file.exists(infile)) stop("Missing results file: ", infile)
obj <- readRDS(infile)
res_all <- obj$results
res <- res_all[res_all$ok, ]
delta_grid <- obj$config$delta_grid
tau_grid   <- obj$config$tau_grid

## Infeasibility rate per (tau, arm): calibrations with no feasible theta are
## excluded from the coverage/length means, so report the exclusion rate
## alongside the figures (CSV + plot caption).
cal_rows <- res_all[abs(res_all$delta_eval) < 1e-9, ]
cal_rows$arm_lab <- ifelse(cal_rows$method == "tilted",
                           sprintf("tilted(dc=%+.2f)", cal_rows$delta_cal), cal_rows$method)
infeas <- aggregate(list(infeasible_rate = !cal_rows$ok),
                    by = list(tau = cal_rows$tau, arm = cal_rows$arm_lab), FUN = mean)
max_infeas <- max(infeas$infeasible_rate)

# Build an arm label per row.
res$arm <- with(res, ifelse(
    method == "one", "one (paper)",
    ifelse(abs(delta_cal - delta_eval) < 1e-9, "tilted (matched)",
           ifelse(method == "estimated", "estimated (delta_cal=0)",
                  sprintf("tilted (delta_cal=%+.2f)", delta_cal)))))
# Keep: one, matched, estimated(=delta_cal 0 line across delta_eval), and the two
# canonical mismatch calibrations for a readable figure.
keep_fixed <- c(sprintf("tilted (delta_cal=%+.2f)", -0.15),
                sprintf("tilted (delta_cal=%+.2f)",  0.15),
                "estimated (delta_cal=0)")
res <- res[res$arm %in% c("one (paper)", "tilted (matched)", keep_fixed), ]

agg <- res %>%
    group_by(tau, arm, delta_eval) %>%
    summarise(coverage = mean(coverage, na.rm = TRUE),
              med_len  = mean(med_len,  na.rm = TRUE),
              n_rep    = sum(is.finite(coverage)), .groups = "drop")
agg$tau_lab <- factor(paste0("tau == ", agg$tau), levels = paste0("tau == ", tau_grid))

if (!dir.exists(file.path(folder, "plots"))) dir.create(file.path(folder, "plots"))

p_cov <- ggplot(agg, aes(delta_eval, coverage, colour = arm, shape = arm)) +
    geom_hline(yintercept = 1 - alpha, linetype = "dashed", colour = "grey40") +
    geom_line() + geom_point() +
    facet_wrap(~ tau_lab, labeller = label_parsed) +
    labs(x = expression(delta[eval]~"(evaluation censoring tilt)"),
         y = "Mean survivor-conditional coverage",
         colour = "Calibration arm", shape = "Calibration arm",
         title = sprintf("Censoring-shift sensitivity (%s, R=%d, n=%d)", setup, R, n),
         caption = sprintf(
             "Infeasible calibrations (no feasible theta) excluded from means; max rate %.1f%% (see %s).",
             100 * max_infeas, "gtau_tilt_infeasibility_*.csv")) +
    theme_bw() + theme(legend.position = "bottom")

p_len <- ggplot(agg, aes(delta_eval, med_len, colour = arm, shape = arm)) +
    geom_line() + geom_point() +
    facet_wrap(~ tau_lab, labeller = label_parsed) +
    labs(x = expression(delta[eval]), y = "Mean of median interval length",
         colour = "Calibration arm", shape = "Calibration arm") +
    theme_bw() + theme(legend.position = "bottom")

f_cov <- file.path(folder, "plots", sprintf("gtau_tilt_coverage_R%d_n%d%s.pdf", R, n, mtag))
f_len <- file.path(folder, "plots", sprintf("gtau_tilt_length_R%d_n%d%s.pdf", R, n, mtag))
f_inf <- file.path(folder, "plots", sprintf("gtau_tilt_infeasibility_R%d_n%d%s.csv", R, n, mtag))
ggsave(f_cov, p_cov, width = 10, height = 4)
ggsave(f_len, p_len, width = 10, height = 4)
write.csv(infeas[order(infeas$tau, infeas$arm), ], f_inf, row.names = FALSE)
cat("Saved:\n ", f_cov, "\n ", f_len, "\n ", f_inf, "\n")
