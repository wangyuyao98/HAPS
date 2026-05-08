ddh_tuning_presets <- function() {
    list(
        base20 = list(
            description = "Current DDH baseline used in previous runs.",
            args = list(
                epochs = 20L,
                hidden_size = 32L,
                num_layers = 1L,
                dropout = 0.1,
                lr = 1e-3,
                weight_decay = 1e-5,
                batch_size = 64L,
                num_time_bins = 40L
            )
        ),
        epochs50 = list(
            description = "Same architecture as base20 with longer training.",
            args = list(
                epochs = 50L,
                hidden_size = 32L,
                num_layers = 1L,
                dropout = 0.1,
                lr = 1e-3,
                weight_decay = 1e-5,
                batch_size = 64L,
                num_time_bins = 40L
            )
        ),
        wide64_epochs50 = list(
            description = "Larger GRU hidden size with longer training.",
            args = list(
                epochs = 50L,
                hidden_size = 64L,
                num_layers = 1L,
                dropout = 0.1,
                lr = 1e-3,
                weight_decay = 1e-5,
                batch_size = 64L,
                num_time_bins = 40L
            )
        ),
        wide64_epochs100 = list(
            description = "Larger GRU hidden size with substantially longer training.",
            args = list(
                epochs = 100L,
                hidden_size = 64L,
                num_layers = 1L,
                dropout = 0.1,
                lr = 1e-3,
                weight_decay = 1e-5,
                batch_size = 64L,
                num_time_bins = 40L
            )
        ),
        wide64_lowreg50 = list(
            description = "Larger GRU with weaker regularization to diagnose underfitting.",
            args = list(
                epochs = 50L,
                hidden_size = 64L,
                num_layers = 1L,
                dropout = 0,
                lr = 1e-3,
                weight_decay = 1e-6,
                batch_size = 64L,
                num_time_bins = 40L
            )
        ),
        lstmish100_layers2_100 = list(
            description = "Closer in capacity to the official Dynamic-DeepHit example, but still using the local GRU backend.",
            args = list(
                epochs = 100L,
                hidden_size = 100L,
                num_layers = 2L,
                dropout = 0.1,
                lr = 1e-4,
                weight_decay = 1e-5,
                batch_size = 32L,
                num_time_bins = 40L
            )
        ),
        wide64_bins60_50 = list(
            description = "Larger hidden size plus a finer residual-time grid.",
            args = list(
                epochs = 50L,
                hidden_size = 64L,
                num_layers = 1L,
                dropout = 0.1,
                lr = 1e-3,
                weight_decay = 1e-5,
                batch_size = 64L,
                num_time_bins = 60L
            )
        )
    )
}

ddh_tuning_preset_names <- function() {
    names(ddh_tuning_presets())
}

ddh_tuning_preset_table <- function() {
    presets <- ddh_tuning_presets()
    do.call(
        rbind,
        lapply(names(presets), function(nm) {
            x <- presets[[nm]]
            data.frame(
                preset = nm,
                description = x$description,
                epochs = x$args$epochs,
                hidden_size = x$args$hidden_size,
                num_layers = x$args$num_layers,
                dropout = x$args$dropout,
                lr = x$args$lr,
                weight_decay = x$args$weight_decay,
                batch_size = x$args$batch_size,
                num_time_bins = x$args$num_time_bins,
                stringsAsFactors = FALSE
            )
        })
    )
}

ddh_get_tuning_preset <- function(name) {
    presets <- ddh_tuning_presets()
    if (!name %in% names(presets)) {
        stop(
            "Unknown DDH tuning preset: ", name,
            ". Available presets: ", paste(names(presets), collapse = ", ")
        )
    }
    presets[[name]]
}

ddh_apply_tuning_preset <- function(base_args, preset_name) {
    preset <- ddh_get_tuning_preset(preset_name)
    utils::modifyList(base_args, preset$args)
}

ddh_tuning_tag <- function(preset_name) {
    gsub("[^A-Za-z0-9]+", "_", preset_name)
}
