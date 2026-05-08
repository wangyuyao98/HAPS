rm(list = ls())

# Run this script from the repo root.
if (!file.exists("src/helpers.R")) {
    stop("Please run main_plot_mainpaper_simu_covlen.R from the repository root directory.")
}

source("src/helpers.R")
source("src/helpers_AIPCW.R")
source("src/dynamicCP_AIPCW.R")

## Main-paper simulation figures:
## coverage + length only, one PDF per selected simulation setting.

alpha <- 0.1
rho <- 0.3
n <- 1000L
n_test <- 500L
R <- 200L

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && nzchar(args[[1L]])) R <- as.integer(args[[1L]])
if (length(args) >= 2L && nzchar(args[[2L]])) n <- as.integer(args[[2L]])
if (length(args) >= 3L && nzchar(args[[3L]])) n_test <- as.integer(args[[3L]])
if (!is.finite(R) || R < 1L) stop("R must be a positive integer.")
if (!is.finite(n) || n < 10L) stop("n must be at least 10.")
if (!is.finite(n_test) || n_test < 10L) stop("n_test must be at least 10.")

tau_list <- c(0, 3, 6)
main_plot_width <- 8
main_two_panel_height <- 3.3
main_stacked_height <- 6
main_plot_cex <- list(legend = 1.4, axis = 1.4, lab = 1.3, main = 1.4)
coverage_ticks <- c(0.7, 0.8, 0.9, 1)

method_label <- c(
    DynaCoPS = "HAPS",
    DynaCoPS_A = "HAPS-A",
    DynaCoPS_DR = "HAPS-DR",
    IPCW0 = "IPCW",
    IPCW_trunc = "trunc",
    LM_IPCW = "LM"
)

boxfill_cols <- c(
    DynaCoPS = adjustcolor("#D73027", 0.78),
    DynaCoPS_A = adjustcolor("#1A9850", 0.78),
    DynaCoPS_DR = adjustcolor("#2C7BB6", 0.78),
    IPCW0 = adjustcolor("#737373", 0.72),
    IPCW_trunc = adjustcolor("#BDBDBD", 0.72),
    LM_IPCW = adjustcolor("#F6E8B1", 0.85)
)

boxborder_cols <- c(
    DynaCoPS = "#8E1B16",
    DynaCoPS_A = "#0E6330",
    DynaCoPS_DR = "#174B78",
    IPCW0 = "#4D4D4D",
    IPCW_trunc = "#7A7A7A",
    LM_IPCW = "#A67C52"
)

load_method_bundle <- function(path) {
    if (!file.exists(path)) {
        stop("Missing result file: ", path)
    }
    obj <- readRDS(path)
    list(
        res_ok = obj$results$res_ok,
        bounds = obj$results$results_bounds,
        config = obj$config
    )
}

load_static_bundle <- function(path, source = c("static_ready", "dynamic_tau0")) {
    source <- match.arg(source)
    if (!file.exists(path)) {
        stop("Missing result file: ", path)
    }

    obj <- readRDS(path)
    if (identical(source, "static_ready")) {
        return(list(
            res_ok = obj$results$res_ok,
            bounds = obj$results$results_bounds,
            config = obj$config
        ))
    }

    stat <- build_static_from_dynamic(
        results_bounds = obj$results$results_bounds,
        tau_grid = obj$config$tau_grid,
        TT.name = obj$config$TT.name
    )
    res_static <- do.call(rbind, stat$results_summary_static)
    list(
        res_ok = res_static[res_static$ok %in% TRUE, , drop = FALSE],
        bounds = stat$results_bounds_static,
        config = obj$config
    )
}

format_dynamic_file <- function(setup, method, model.pred, model.C, model.S, model.Xi, loc) {
    file.path(
        "results", setup,
        paste0(
            method, "_", model.pred, "_", model.C,
            "_S", model.S, "_Xi", model.Xi,
            "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
            "_alpha", alpha, "_loc", loc, "_Gtauone_d0.rds"
        )
    )
}

format_linwb_dr_file <- function(tag) {
    file.path(
        "results", "linWB1",
        paste0(
            "derived_DynaCoPS_DR_", tag,
            "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
            "_alpha", alpha, ".rds"
        )
    )
}

format_static_file <- function(setup, method, tag) {
    file.path(
        "results", setup,
        paste0(
            method, "_", tag,
            "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
            "_alpha", alpha, ".rds"
        )
    )
}

format_ddh_tuning_file <- function(method) {
    file.path(
        "results", "mixC2_DGM1", "ddh_tuning",
        paste0(
            method,
            "_ddh_rsf_Sxgb_cox_Xixgb_multiclass_tunewide64_epochs50",
            "_rho", rho, "_R", R, "_n", n, "_ntest", n_test,
            "_alpha", alpha, "_loc0_Gtauone_d0.rds"
        )
    )
}

plot_main_covlen <- function(setting) {
    out_dir <- file.path("results", setting$setup, "plots")
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    method_order <- c("DynaCoPS", "DynaCoPS_DR", "DynaCoPS_A")
    methods <- list(
        DynaCoPS = load_method_bundle(setting$dyna_file),
        DynaCoPS_DR = load_method_bundle(setting$dr_file),
        DynaCoPS_A = load_method_bundle(setting$aipcw_file)
    )

    out_pdf <- file.path(
        out_dir,
        paste0("h_mainpaper_covlen_versions_", setting$tag, ".pdf")
    )

    make_horiz_cov_len_onefig(
        methods = methods,
        tau_list = tau_list,
        alpha = alpha,
        out_pdf = out_pdf,
        method_order = method_order,
        boxfill_cols = unname(boxfill_cols[method_order]),
        boxborder_cols = unname(boxborder_cols[method_order]),
        method_label = method_label,
        offset = 0.18,
        boxwex = 0.12,
        ylim_coverage = c(0.7, 1),
        ylim_length = NULL,
        alpha_line_col = "grey60",
        plot_cex = main_plot_cex,
        coverage_col = "coverage_test_true",
        coverage_ticks = coverage_ticks,
        width = main_plot_width,
        height = main_two_panel_height
    )

    cat("\n", setting$label, "\n", sep = "")
    for (m in method_order) {
        res_ok <- methods[[m]]$res_ok
        cat(
            method_label[[m]], ": coverage=",
            round(mean(res_ok$coverage_test_true, na.rm = TRUE), 4),
            ", length=",
            round(mean(res_ok$mean_len_test, na.rm = TRUE), 4),
            ", n_ok=", nrow(res_ok), "\n",
            sep = ""
        )
    }
    invisible(out_pdf)
}

plot_stacked_main_covlen <- function(setting) {
    out_dir <- file.path("results", setting$setup, "plots")
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    methods <- list(
        DynaCoPS = load_method_bundle(setting$dyna_file),
        DynaCoPS_DR = load_method_bundle(setting$dr_file),
        DynaCoPS_A = load_method_bundle(setting$aipcw_file),
        IPCW0 = load_static_bundle(setting$static_file, setting$static_source),
        IPCW_trunc = load_static_bundle(setting$ipcw_trunc_file, "static_ready"),
        LM_IPCW = load_static_bundle(setting$lm_ipcw_file, "static_ready")
    )

    fig_top_order <- c("DynaCoPS", "IPCW0", "IPCW_trunc", "LM_IPCW")
    fig_bottom_order <- c("DynaCoPS", "DynaCoPS_DR", "DynaCoPS_A")

    out_pdf <- file.path(
        out_dir,
        paste0("h_mainpaper_covlen_stacked_", setting$tag, ".pdf")
    )

    make_horiz_cov_len_twofigs(
        methods = methods,
        tau_list = tau_list,
        alpha = alpha,
        figA_order = fig_bottom_order,
        figA_fill = unname(boxfill_cols[fig_bottom_order]),
        figA_edge = unname(boxborder_cols[fig_bottom_order]),
        figB_order = fig_top_order,
        figB_fill = unname(boxfill_cols[fig_top_order]),
        figB_edge = unname(boxborder_cols[fig_top_order]),
        out_pdf = out_pdf,
        method_label = method_label,
        offsetA = 0.18,
        offsetB = 0.28,
        boxwex = 0.12,
        ylim_coverage = c(0.7, 1),
        ylim_length = NULL,
        alpha_line_col = "grey60",
        plot_cex = main_plot_cex,
        coverage_col = "coverage_test_true",
        coverage_ticks = coverage_ticks,
        width = main_plot_width,
        height = main_stacked_height
    )

    cat("\n", setting$label, " stacked static + versions\n", sep = "")
    invisible(out_pdf)
}

settings <- list(
    linWB1_cox_cox = list(
        setup = "linWB1",
        label = "linWB1: Cox prediction with Cox nuisance models",
        tag = "cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone",
        dyna_file = format_dynamic_file("linWB1", "dCP_IPCW", "cox", "cox", "cox", "xgb_multiclass", 0),
        dr_file = format_linwb_dr_file("cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone"),
        aipcw_file = format_dynamic_file("linWB1", "dCP_AIPCW", "cox", "cox", "cox", "xgb_multiclass", 0)
    ),
    mixC2_DGM1_ddh_flex = list(
        setup = "mixC2_DGM1",
        label = "mixC2_DGM1: tuned DDH prediction with flexible nuisance models",
        tag = "ddh_rsf_Sxgb_cox_Xixgb_multiclass_tunewide64_epochs50_loc0_Gtauone",
        dyna_file = format_ddh_tuning_file("dCP_IPCW_new"),
        dr_file = format_ddh_tuning_file("derived_DynaCoPS_DR"),
        aipcw_file = format_ddh_tuning_file("dCP_AIPCW")
    )
)

invisible(lapply(settings, plot_main_covlen))

stacked_settings <- list(
    linWB1_cox_cox = c(
        settings$linWB1_cox_cox,
        list(
            static_source = "static_ready",
            static_file = format_static_file("linWB1", "static_IPCW0", "cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone"),
            ipcw_trunc_file = format_static_file("linWB1", "static_IPCW_trunc", "cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone"),
            lm_ipcw_file = format_static_file("linWB1", "static_LM_IPCW", "cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone")
        )
    ),
    mixC2_DGM1_ddh_flex_static_rsf = list(
        setup = "mixC2_DGM1",
        label = "mixC2_DGM1: tuned DDH prediction with flexible nuisance models and RSF static methods",
        tag = "ddh_rsf_Sxgb_cox_Xixgb_multiclass_tunewide64_epochs50_static_rsf_rsf_globalG_loc0_Gtauone",
        dyna_file = format_ddh_tuning_file("dCP_IPCW_new"),
        dr_file = format_ddh_tuning_file("derived_DynaCoPS_DR"),
        aipcw_file = format_ddh_tuning_file("dCP_AIPCW"),
        static_source = "static_ready",
        static_file = format_static_file("mixC2_DGM1", "static_IPCW0_globalG", "rsf_rsf_Srsf_Xixgb_reg_loc0_Gtauone"),
        ipcw_trunc_file = format_static_file("mixC2_DGM1", "static_IPCW_trunc_globalG", "rsf_rsf_Srsf_Xixgb_reg_loc0_Gtauone"),
        lm_ipcw_file = format_static_file("mixC2_DGM1", "static_LM_IPCW", "rsf_rsf_Srsf_Xixgb_reg_loc0_Gtauone")
    )
)

invisible(lapply(stacked_settings, plot_stacked_main_covlen))
