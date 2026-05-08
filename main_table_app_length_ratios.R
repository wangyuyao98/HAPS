rm(list = ls())

# Ratio of median test-set interval length for application results.
# Run this script from the repository root after the application result files exist.

if (!file.exists("src/app_analysis_helpers.R")) {
    stop("Please run main_table_app_length_ratios.R from the repository root directory.")
}

source("src/app_analysis_helpers.R")
source("src/ddh_tuning_presets.R")

args <- commandArgs(trailingOnly = TRUE)

alpha <- if (length(args) >= 1L && nzchar(args[[1L]])) as.numeric(args[[1L]]) else 0.2
R <- if (length(args) >= 2L && nzchar(args[[2L]])) as.integer(args[[2L]]) else 200L
ddh_arg <- if (length(args) >= 3L && nzchar(args[[3L]])) args[[3L]] else "wide64_epochs50"
time_feature_tag <- if (length(args) >= 4L && nzchar(args[[4L]])) args[[4L]] else "time"
G_model <- if (length(args) >= 5L && nzchar(args[[5L]])) args[[5L]] else "cox"

if (!is.finite(alpha)) stop("alpha must be numeric.")
if (!is.finite(R) || R < 1L) stop("R must be a positive integer.")
time_feature_tag <- match.arg(time_feature_tag, c("time", "notime"))
G_model <- match.arg(G_model, c("cox", "rsf"))

ddh_tuning_preset <- NULL
ddh_epochs <- suppressWarnings(as.integer(ddh_arg))
if (is.na(ddh_epochs)) {
    ddh_tuning_preset <- ddh_arg
    ddh_epochs <- as.integer(ddh_get_tuning_preset(ddh_tuning_preset)$args$epochs)
}

out_dir <- file.path("results", "mainpaper_tables")
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

tau_keep <- c(3, 5)

ddh_tune_part <- function() {
    if (is.null(ddh_tuning_preset)) return("")
    paste0("_tune", ddh_tuning_tag(ddh_tuning_preset))
}

app_file <- function(dataset_name, file) {
    file.path("results", paste0("app_", dataset_name), file)
}

ddh_file <- function(dataset_name, method) {
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
    static_tag <- if (identical(G_model, "rsf")) "rsfpred_rsfG" else "coxpred_coxG"
    app_file(
        dataset_name,
        paste0(method, "_", static_tag, "_R", R, "_n", cfg$n_all, "_alpha", alpha, ".rds")
    )
}

median_len_by_tau <- function(path, tau_keep = c(3, 5)) {
    if (!file.exists(path)) {
        stop("Missing result file: ", path)
    }
    obj <- readRDS(path)
    res <- obj$results$res_ok
    out <- aggregate(
        mean_len_test ~ tau,
        data = res[res$tau %in% tau_keep & is.finite(res$mean_len_test), , drop = FALSE],
        FUN = median,
        na.rm = TRUE
    )
    stats::setNames(out$mean_len_test, as.character(out$tau))
}

dataset_labels <- c(pbcseq = "PBC", colon = "Colon")
result_files <- setNames(vector("list", length(dataset_labels)), names(dataset_labels))
for (dataset in names(dataset_labels)) {
    result_files[[dataset]] <- list(
        label = unname(dataset_labels[[dataset]]),
        HAPS = ddh_file(dataset, "dCP_IPCW"),
        HAPS_DR = ddh_file(dataset, "derived_DynaCoPS_DR"),
        IPCW = static_file(dataset, "static_IPCW"),
        trunc = static_file(dataset, "static_IPCW_trunc"),
        LM = static_file(dataset, "static_LM_IPCW")
    )
}

dynamic_methods <- c(HAPS = "HAPS", HAPS_DR = "HAPS-DR")
static_methods <- c("IPCW", "trunc", "LM")

ratio_rows <- lapply(names(result_files), function(dataset) {
    files <- result_files[[dataset]]
    med <- lapply(files[setdiff(names(files), "label")], median_len_by_tau, tau_keep = tau_keep)

    do.call(rbind, lapply(names(dynamic_methods), function(dynamic_method) {
        do.call(rbind, lapply(tau_keep, function(tau) {
            tau_chr <- as.character(tau)
            data.frame(
                dataset = dataset,
                dataset_label = files$label,
                dynamic_method = dynamic_methods[[dynamic_method]],
                tau = tau,
                IPCW = med[[dynamic_method]][tau_chr] / med[["IPCW"]][tau_chr],
                trunc = med[[dynamic_method]][tau_chr] / med[["trunc"]][tau_chr],
                LM = med[[dynamic_method]][tau_chr] / med[["LM"]][tau_chr],
                row.names = NULL
            )
        }))
    }))
})

ratio_table <- do.call(rbind, ratio_rows)
ratio_table[, static_methods] <- round(ratio_table[, static_methods], 3)

suffix <- paste0(
    if (identical(G_model, "cox")) "coxG" else "rsfG",
    "_R", R,
    "_alpha", alpha
)
csv_file <- file.path(out_dir, paste0("app_length_ratio_table_", suffix, ".csv"))
write.csv(ratio_table, csv_file, row.names = FALSE)

fmt <- function(x) sprintf("%.3f", as.numeric(x))

latex_lines <- c(
    "\\begin{table}[!ht]",
    "\\centering",
    "\\caption{Ratio of median test-set prediction interval length for HAPS and HAPS-DR relative to baseline methods in the two data applications. Values below one indicate shorter intervals.}",
    "\\label{tab:app_length_ratio}",
    "{\\fontsize{8}{9}\\selectfont",
    "\\begin{tabular}{llcccccc}",
    "\\toprule",
    "&  & \\multicolumn{3}{c}{$\\tau=3$} & \\multicolumn{3}{c}{$\\tau=5$} \\\\",
    "\\cmidrule(lr){3-5} \\cmidrule(lr){6-8}",
    "Data & Method & IPCW & trunc & LM & IPCW & trunc & LM \\\\",
    "\\midrule"
)

for (dataset in names(result_files)) {
    d <- ratio_table[ratio_table$dataset == dataset, , drop = FALSE]
    for (i in seq_along(dynamic_methods)) {
        method_label <- unname(dynamic_methods[[i]])
        row <- d[d$dynamic_method == method_label, , drop = FALSE]
        row <- row[match(tau_keep, row$tau), , drop = FALSE]
        first_col <- if (i == 1L) result_files[[dataset]]$label else ""
        latex_lines <- c(
            latex_lines,
            paste0(
                first_col, " & ", method_label, " & ",
                paste(c(
                    fmt(row$IPCW[1]), fmt(row$trunc[1]), fmt(row$LM[1]),
                    fmt(row$IPCW[2]), fmt(row$trunc[2]), fmt(row$LM[2])
                ), collapse = " & "),
                " \\\\"
            )
        )
    }
    if (!identical(dataset, tail(names(result_files), 1))) {
        latex_lines <- c(latex_lines, "\\midrule")
    }
}

latex_lines <- c(
    latex_lines,
    "\\bottomrule",
    "\\end{tabular}",
    "}",
    "\\end{table}"
)

tex_file <- file.path(out_dir, paste0("app_length_ratio_table_", suffix, ".tex"))
writeLines(latex_lines, tex_file)

cat("Saved application length-ratio CSV to:\n  ", csv_file, "\n", sep = "")
cat("Saved application length-ratio LaTeX table to:\n  ", tex_file, "\n", sep = "")
print(ratio_table)
