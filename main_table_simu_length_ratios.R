rm(list = ls())

# Run this script from the repository root.
if (!file.exists("src/helpers.R")) {
    stop("Please run main_table_simu_length_ratios.R from the repository root directory.")
}

source("src/helpers.R")
source("src/helpers_AIPCW.R")
source("src/dynamicCP_AIPCW.R")

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
main_tau_list <- c(3, 6)
out_dir <- file.path("results", "mainpaper_tables")
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

method_label <- c(
    DynaCoPS = "HAPS",
    DynaCoPS_A = "HAPS-A",
    DynaCoPS_DR = "HAPS-DR",
    IPCW0 = "IPCW",
    IPCW_trunc = "trunc",
    LM_IPCW = "LM"
)

dynamic_methods <- c("DynaCoPS", "DynaCoPS_A", "DynaCoPS_DR")
static_methods <- c("IPCW0", "IPCW_trunc", "LM_IPCW")

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

path_exists <- function(paths) {
    vapply(paths, file.exists, logical(1))
}

missing_paths <- function(setting) {
    paths <- c(
        setting$dyna_file,
        setting$dr_file,
        setting$aipcw_file,
        setting$static_file,
        setting$ipcw_trunc_file,
        setting$lm_ipcw_file
    )
    paths[!path_exists(paths)]
}

median_len_by_tau <- function(bundle, tau_values = tau_list) {
    res <- bundle$res_ok
    out <- stats::aggregate(
        mean_len_test ~ tau,
        data = res[res$tau %in% tau_values & is.finite(res$mean_len_test), , drop = FALSE],
        FUN = stats::median,
        na.rm = TRUE
    )
    stats::setNames(out$mean_len_test, as.character(out$tau))
}

mean_coverage_by_tau <- function(bundle, tau_values = tau_list) {
    res <- bundle$res_ok
    out <- stats::aggregate(
        coverage_test_true ~ tau,
        data = res[res$tau %in% tau_values & is.finite(res$coverage_test_true), , drop = FALSE],
        FUN = mean,
        na.rm = TRUE
    )
    stats::setNames(out$coverage_test_true, as.character(out$tau))
}

load_setting_methods <- function(setting) {
    list(
        DynaCoPS = load_method_bundle(setting$dyna_file),
        DynaCoPS_DR = load_method_bundle(setting$dr_file),
        DynaCoPS_A = load_method_bundle(setting$aipcw_file),
        IPCW0 = load_static_bundle(setting$static_file, setting$static_source),
        IPCW_trunc = load_static_bundle(setting$ipcw_trunc_file, "static_ready"),
        LM_IPCW = load_static_bundle(setting$lm_ipcw_file, "static_ready")
    )
}

build_summary_for_setting <- function(setting_name, setting) {
    miss <- missing_paths(setting)
    if (length(miss) > 0L) {
        warning(
            "Skipping ", setting_name, " because these files are missing:\n  ",
            paste(miss, collapse = "\n  "),
            call. = FALSE
        )
        return(NULL)
    }

    methods <- load_setting_methods(setting)
    len_by_method <- lapply(methods, median_len_by_tau)
    cov_by_method <- lapply(methods, mean_coverage_by_tau)

    length_rows <- do.call(rbind, lapply(names(methods), function(method) {
        data.frame(
            setting_key = setting_name,
            setting_label = setting$label,
            method = method,
            method_label = unname(method_label[[method]]),
            tau = tau_list,
            median_mean_len_test = as.numeric(len_by_method[[method]][as.character(tau_list)]),
            mean_coverage_test_true = as.numeric(cov_by_method[[method]][as.character(tau_list)]),
            stringsAsFactors = FALSE
        )
    }))

    ratio_rows <- do.call(rbind, lapply(dynamic_methods, function(dyn) {
        do.call(rbind, lapply(static_methods, function(sta) {
            data.frame(
                setting_key = setting_name,
                setting_label = setting$label,
                dynamic_method = dyn,
                dynamic_label = unname(method_label[[dyn]]),
                static_method = sta,
                static_label = unname(method_label[[sta]]),
                tau = tau_list,
                dynamic_median_mean_len_test = as.numeric(len_by_method[[dyn]][as.character(tau_list)]),
                static_median_mean_len_test = as.numeric(len_by_method[[sta]][as.character(tau_list)]),
                length_ratio = as.numeric(len_by_method[[dyn]][as.character(tau_list)] /
                    len_by_method[[sta]][as.character(tau_list)]),
                dynamic_mean_coverage_test_true = as.numeric(cov_by_method[[dyn]][as.character(tau_list)]),
                static_mean_coverage_test_true = as.numeric(cov_by_method[[sta]][as.character(tau_list)]),
                stringsAsFactors = FALSE
            )
        }))
    }))

    list(lengths = length_rows, ratios = ratio_rows)
}

fmt_num <- function(x, digits = 2) {
    ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}

escape_latex <- function(x) {
    x <- gsub("\\\\", "\\\\textbackslash{}", x)
    x <- gsub("&", "\\\\&", x)
    x <- gsub("%", "\\\\%", x)
    x <- gsub("_", "\\\\_", x)
    x
}

write_main_latex <- function(df, path) {
    lines <- c(
        "\\begin{table}[!ht]",
        "\\centering",
        "\\caption{Ratio of median test-set prediction interval length for HAPS relative to landmark IPCW. Values below one indicate shorter HAPS intervals.}",
        "\\label{tab:main_simu_len_ratio_haps_ipcwlm}",
        "{\\fontsize{8}{9}\\selectfont",
        "\\begin{tabular}{lcc}",
        "\\toprule",
        paste0("Setting & $\\tau=", main_tau_list[1], "$ & $\\tau=", main_tau_list[2], "$ \\\\"),
        "\\midrule"
    )

    for (i in seq_len(nrow(df))) {
        lines <- c(
            lines,
            paste0(
                escape_latex(df$setting_label[i]), " & ",
                fmt_num(df[[paste0("tau_", main_tau_list[1])]][i]), " & ",
                fmt_num(df[[paste0("tau_", main_tau_list[2])]][i]), " \\\\"
            )
        )
    }

    lines <- c(
        lines,
        "\\bottomrule",
        "\\end{tabular}",
        "}",
        "\\end{table}"
    )
    writeLines(lines, path)
}

write_appendix_latex <- function(df, path) {
    static_header_short <- c(
        IPCW0 = "IPCW",
        IPCW_trunc = "trunc",
        LM_IPCW = "LM"
    )
    lines <- c(
        "\\begin{table}[!ht]",
        "\\centering",
        "\\caption{Ratio of median test-set prediction interval length for each HAPS version relative to static comparison methods. Values below one indicate shorter dynamic intervals. Here trunc denotes IPCW-trunc and LM denotes landmark IPCW.}",
        "\\label{tab:appendix_simu_len_ratio_haps_static}",
        "{\\fontsize{7}{8}\\selectfont",
        "\\begin{tabular}{llccccccccc}",
        "\\toprule",
        "& & \\multicolumn{3}{c}{$\\tau=0$} & \\multicolumn{3}{c}{$\\tau=3$} & \\multicolumn{3}{c}{$\\tau=6$} \\\\",
        "\\cmidrule(lr){3-5} \\cmidrule(lr){6-8} \\cmidrule(lr){9-11}",
        paste0(
            "Setting & Dynamic method & ",
            paste(rep(unname(static_header_short[static_methods]), times = length(tau_list)), collapse = " & "),
            " \\\\"
        ),
        "\\midrule"
    )

    setting_levels <- unique(df$setting_label)
    for (s in setting_levels) {
        dfs <- df[df$setting_label == s, , drop = FALSE]
        dyn_levels <- unique(dfs$dynamic_label)
        for (j in seq_along(dyn_levels)) {
            d <- dyn_levels[j]
            row_df <- dfs[dfs$dynamic_label == d, , drop = FALSE]
            first_col <- if (j == 1L) escape_latex(s) else ""
            vals <- c()
            for (tau in tau_list) {
                for (sta in unname(method_label[static_methods])) {
                    idx <- row_df$tau == tau & row_df$static_label == sta
                    vals <- c(vals, fmt_num(row_df$length_ratio[idx][1]))
                }
            }
            lines <- c(
                lines,
                paste0(first_col, " & ", escape_latex(d), " & ", paste(vals, collapse = " & "), " \\\\")
            )
        }
        if (!identical(s, tail(setting_levels, 1))) {
            lines <- c(lines, "\\addlinespace")
        }
    }

    lines <- c(
        lines,
        "\\bottomrule",
        "\\end{tabular}",
        "}",
        "\\end{table}"
    )
    writeLines(lines, path)
}

settings <- list(
    linWB1_cox_cox = list(
        setup = "linWB1",
        label = "linWB1",
        tag = "cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone",
        dyna_file = format_dynamic_file("linWB1", "dCP_IPCW", "cox", "cox", "cox", "xgb_multiclass", 0),
        dr_file = format_linwb_dr_file("cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone"),
        aipcw_file = format_dynamic_file("linWB1", "dCP_AIPCW", "cox", "cox", "cox", "xgb_multiclass", 0),
        static_source = "static_ready",
        static_file = format_static_file("linWB1", "static_IPCW0", "cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone"),
        ipcw_trunc_file = format_static_file("linWB1", "static_IPCW_trunc", "cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone"),
        lm_ipcw_file = format_static_file("linWB1", "static_LM_IPCW", "cox_cox_Scox_Xixgb_multiclass_loc0_Gtauone")
    ),
    mixC2_DGM1_ddh_flex_static_rsf_globalG = list(
        setup = "mixC2_DGM1",
        label = "mixC2 DGM1",
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

summaries <- lapply(names(settings), function(nm) build_summary_for_setting(nm, settings[[nm]]))
summaries <- summaries[!vapply(summaries, is.null, logical(1))]
if (length(summaries) == 0L) {
    stop("No complete settings found; length-ratio tables were not generated.")
}

length_summary <- do.call(rbind, lapply(summaries, `[[`, "lengths"))
ratio_summary <- do.call(rbind, lapply(summaries, `[[`, "ratios"))

main_df_long <- ratio_summary[
    ratio_summary$dynamic_method == "DynaCoPS" &
        ratio_summary$static_method == "LM_IPCW" &
        ratio_summary$tau %in% main_tau_list,
    ,
    drop = FALSE
]

main_settings <- unique(main_df_long[, c("setting_key", "setting_label")])
main_table <- main_settings
for (tau in main_tau_list) {
    vals <- main_df_long[main_df_long$tau == tau, c("setting_key", "length_ratio")]
    names(vals)[2] <- paste0("tau_", tau)
    main_table <- merge(main_table, vals, by = "setting_key", all.x = TRUE, sort = FALSE)
}

write.csv(length_summary, file.path(out_dir, "simu_length_coverage_summary.csv"), row.names = FALSE)
write.csv(ratio_summary, file.path(out_dir, "simu_length_ratio_haps_appendix_full.csv"), row.names = FALSE)
write.csv(main_table, file.path(out_dir, "simu_length_ratio_main_haps_vs_ipcwlm.csv"), row.names = FALSE)

write_main_latex(
    main_table,
    file.path(out_dir, "simu_length_ratio_main_haps_vs_ipcwlm.tex")
)
write_appendix_latex(
    ratio_summary,
    file.path(out_dir, "simu_length_ratio_haps_appendix_full.tex")
)

cat("Saved length-ratio tables to:\n")
cat("  ", file.path(out_dir, "simu_length_coverage_summary.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "simu_length_ratio_haps_appendix_full.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "simu_length_ratio_haps_appendix_full.tex"), "\n", sep = "")
cat("  ", file.path(out_dir, "simu_length_ratio_main_haps_vs_ipcwlm.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "simu_length_ratio_main_haps_vs_ipcwlm.tex"), "\n", sep = "")
