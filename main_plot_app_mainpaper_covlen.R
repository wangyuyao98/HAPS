rm(list = ls())

# Run this script from the repository root.
if (!file.exists("src/helpers.R")) {
    stop("Please run main_plot_app_mainpaper_covlen.R from the repository root directory.")
}

source("src/helpers.R")
source("src/app_analysis_helpers.R")
source("src/ddh_tuning_presets.R")

args <- commandArgs(trailingOnly = TRUE)

alpha <- if (length(args) >= 1L && nzchar(args[[1L]])) as.numeric(args[[1L]]) else 0.2
R <- if (length(args) >= 2L && nzchar(args[[2L]])) as.integer(args[[2L]]) else 200L
ddh_arg <- if (length(args) >= 3L && nzchar(args[[3L]])) args[[3L]] else "wide64_epochs50"
time_feature_tag <- if (length(args) >= 4L && nzchar(args[[4L]])) args[[4L]] else "time"
G_model <- if (length(args) >= 5L && nzchar(args[[5L]])) args[[5L]] else "cox"

ddh_tuning_preset <- NULL
ddh_epochs <- suppressWarnings(as.integer(ddh_arg))
if (is.na(ddh_epochs)) {
    ddh_tuning_preset <- ddh_arg
    ddh_epochs <- as.integer(ddh_get_tuning_preset(ddh_tuning_preset)$args$epochs)
}

if (!is.finite(alpha)) stop("alpha must be numeric.")
if (!is.finite(R) || R < 1L) stop("R must be a positive integer.")
time_feature_tag <- match.arg(time_feature_tag, c("time", "notime"))
G_model <- match.arg(G_model, c("cox", "rsf"))

ddh_plot_tag <- if (is.null(ddh_tuning_preset)) {
    paste0("epochs", ddh_epochs)
} else {
    paste0("tune", ddh_tuning_tag(ddh_tuning_preset))
}

main_plot_width <- 8
main_two_panel_height <- 3.5
full_tau_height <- 5.8
main_plot_cex <- list(legend = 1.4, axis = 1.4, lab = 1.3, main = 1.4)
coverage_ticks <- NULL
coverage_xlim <- c(0.2, 1)
app_method_offset <- 0.26
app_boxwex <- 0.09

method_label <- c(
    DynaCoPS = "HAPS",
    DynaCoPS_DR = "HAPS-DR",
    IPCW = "IPCW",
    IPCW_trunc = "trunc",
    LM_IPCW = "LM"
)

boxfill_cols <- c(
    DynaCoPS = grDevices::adjustcolor("#D73027", 0.78),
    DynaCoPS_DR = grDevices::adjustcolor("#2C7BB6", 0.78),
    IPCW = grDevices::adjustcolor("#737373", 0.72),
    IPCW_trunc = grDevices::adjustcolor("#BDBDBD", 0.72),
    LM_IPCW = grDevices::adjustcolor("#F6E8B1", 0.85)
)

boxborder_cols <- c(
    DynaCoPS = "#8E1B16",
    DynaCoPS_DR = "#174B78",
    IPCW = "#4D4D4D",
    IPCW_trunc = "#7A7A7A",
    LM_IPCW = "#A67C52"
)

load_app_bundle <- function(path) {
    if (!file.exists(path)) stop("Missing application result file: ", path)
    obj <- readRDS(path)
    list(
        res_ok = obj$results$res_ok,
        bounds = obj$results$results_bounds,
        config = obj$config,
        path = path
    )
}

app_file <- function(dataset_name, file) {
    file.path("results", paste0("app_", dataset_name), file)
}

ddh_tune_part <- function() {
    if (is.null(ddh_tuning_preset)) return("")
    paste0("_tune", ddh_tuning_tag(ddh_tuning_preset))
}

ddh_file <- function(dataset_name, method, G_model) {
    cfg <- app_dataset_config(dataset_name, alpha)
    app_file(
        dataset_name,
        paste0(
            method, "_ddh_", G_model, "_", time_feature_tag,
            ddh_tune_part(),
            "_R", R, "_n", cfg$n_all,
            "_alpha", alpha, "_epochs", ddh_epochs, ".rds"
        )
    )
}

static_file <- function(dataset_name, method) {
    cfg <- app_dataset_config(dataset_name, alpha)
    if (identical(G_model, "rsf")) {
        new_file <- app_file(
            dataset_name,
            paste0(method, "_rsfpred_rsfG_R", R, "_n", cfg$n_all, "_alpha", alpha, ".rds")
        )
        old_file <- app_file(
            dataset_name,
            paste0(method, "_rsfG_R", R, "_n", cfg$n_all, "_alpha", alpha, ".rds")
        )
        if (file.exists(new_file) || !file.exists(old_file)) new_file else old_file
    } else {
        new_file <- app_file(
            dataset_name,
            paste0(method, "_coxpred_coxG_R", R, "_n", cfg$n_all, "_alpha", alpha, ".rds")
        )
        old_file <- app_file(
            dataset_name,
            paste0(method, "_coxcox_R", R, "_n", cfg$n_all, "_alpha", alpha, ".rds")
        )
        if (file.exists(new_file) || !file.exists(old_file)) new_file else old_file
    }
}

build_methods <- function(dataset_name) {
    list(
        DynaCoPS = load_app_bundle(ddh_file(dataset_name, "dCP_IPCW", G_model)),
        DynaCoPS_DR = load_app_bundle(ddh_file(dataset_name, "derived_DynaCoPS_DR", G_model)),
        IPCW = load_app_bundle(static_file(dataset_name, "static_IPCW")),
        IPCW_trunc = load_app_bundle(static_file(dataset_name, "static_IPCW_trunc")),
        LM_IPCW = load_app_bundle(static_file(dataset_name, "static_LM_IPCW"))
    )
}

subset_method_tau <- function(obj, tau_list) {
    obj$bounds <- subset_bounds_to_tau(obj$bounds, tau_list)
    obj$res_ok <- obj$res_ok[obj$res_ok$tau %in% tau_list, , drop = FALSE]
    obj
}

plot_one_dataset <- function(dataset_name, tau_list, suffix, height) {
    cfg <- app_dataset_config(dataset_name, alpha)
    methods <- lapply(build_methods(dataset_name), subset_method_tau, tau_list = tau_list)
    method_order <- c("DynaCoPS", "DynaCoPS_DR", "IPCW", "IPCW_trunc", "LM_IPCW")
    out_dir <- file.path("results", paste0("app_", dataset_name), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    g_suffix <- if (identical(G_model, "cox")) "" else paste0("_", G_model, "G")

    out_pdf <- file.path(
        out_dir,
        paste0(
            dataset_name,
            "_h_mainpaper_covlen_HAPS_vs_static",
            g_suffix,
            "_",
            suffix,
            "_", ddh_plot_tag,
            "_alpha", alpha,
            ".pdf"
        )
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
        offset = app_method_offset,
        boxwex = app_boxwex,
        ylim_coverage = coverage_xlim,
        ylim_length = cfg$ylim_length,
        alpha_line_col = "grey60",
        plot_cex = main_plot_cex,
        coverage_col = "coverage_test_ipcw",
        coverage_ticks = coverage_ticks,
        width = main_plot_width,
        height = height
    )
    invisible(out_pdf)
}

datasets <- c("pbcseq", "colon")
generated <- c(
    unlist(lapply(datasets, plot_one_dataset,
                  tau_list = c(0, 3, 5),
                  suffix = "tau0_3_5",
                  height = main_two_panel_height)),
    unlist(lapply(datasets, plot_one_dataset,
                  tau_list = 0:5,
                  suffix = "tau0_5",
                  height = full_tau_height))
)

cat("Generated main-paper-style application coverage/length plots:\n")
cat(paste0("  ", generated, collapse = "\n"), "\n")
