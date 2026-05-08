mixC_dgm_registry <- function() {
    mixC_variants <- list(
        DGM1 = list(
            shape_T = c(4, 4, 5),
            shape_C = c(3, 3, 3),
            beta_T0 = c(-8, -8, -5),
            beta_TL = c(1, 2, 3),
            beta_TB = c(0.7, 0.9, 1.1),
            beta_C0 = c(-8.0, -7.5, -7.0),
            beta_CL = c(0.5, 0.6, 0.7),
            beta_CL2 = c(0.4, 0.5, 0.6),
            prob_cox = 0.5,
            uniform_cmax = 12
        ),
        DGM2 = list(
            shape_T = c(4, 4, 5),
            shape_C = c(2, 3, 5),
            beta_T0 = c(-8, -8, -5),
            beta_TL = c(1, 2, 3),
            beta_TB = c(0.7, 0.9, 1.1),
            beta_C0 = c(-8.0, -7.3, -6.7),
            beta_CL = c(0.7, 0.9, 1.1),
            beta_CL2 = c(0.6, 0.8, 1.0),
            prob_cox = 0.5,
            uniform_cmax = 9
        ),
        DGM3 = list(
            shape_T = c(4, 4, 5),
            shape_C = c(1.5, 3, 5),
            beta_T0 = c(-8, -8, -5),
            beta_TL = c(1, 2, 3),
            beta_TB = c(0.7, 0.9, 1.1),
            beta_C0 = c(-7.9, -7.2, -6.6),
            beta_CL = c(0.8, 1.0, 1.2),
            beta_CL2 = c(0.7, 0.9, 1.1),
            prob_cox = 0.5,
            uniform_cmax = 9
        ),
        DGM4 = list(
            shape_T = c(4, 4, 5),
            shape_C = c(1, 1, 1),
            beta_T0 = c(-8, -8, -5),
            beta_TL = c(1, 2, 3),
            beta_TB = c(0.7, 0.9, 1.1),
            beta_C0 = c(-2.8, -2.6, -2.4),
            beta_CL = c(-0.8, -1.0, -1.2),
            beta_CL2 = c(0, 0, 0),
            prob_cox = 0.4,
            uniform_cmax = 6.5
        )
    )

    mixC2_variants <- list(
        DGM1 = utils::modifyList(mixC_variants$DGM1, list(
            uniform_cmax = 15
        )),
        DGM2 = utils::modifyList(mixC_variants$DGM2, list(
            uniform_cmax = 15
        )),
        DGM3 = utils::modifyList(mixC_variants$DGM4, list(
            uniform_cmax = 15
        )),
        DGM4 = utils::modifyList(mixC_variants$DGM4, list(
            beta_C0 = c(-0.8, -0.6, -0.4),
            uniform_cmax = 15
        ))
    )

    list(
        mixC = list(
            alias_map = c(
                base = "DGM1",
                shape_inc_moderate = "DGM2",
                shape_inc_stronger = "DGM3",
                scenario34_tempered = "DGM4"
            ),
            variants = mixC_variants
        ),
        mixC2 = list(
            alias_map = c(
                base = "DGM1",
                shape_inc_moderate = "DGM2",
                scenario34_tempered = "DGM3",
                scenario34_tempered_highcens = "DGM4"
            ),
            variants = mixC2_variants
        )
    )
}

mixC_dgm_family_names <- function() {
    names(mixC_dgm_registry())
}

mixC_dgm_variant_map <- function(dgm_family = "mixC") {
    registry <- mixC_dgm_registry()
    if (!dgm_family %in% names(registry)) {
        stop(
            "Unknown dgm_family: ", dgm_family,
            ". Choose one of: ", paste(names(registry), collapse = ", ")
        )
    }
    registry[[dgm_family]]$variants
}

mixC_dgm_variant_alias_map <- function(dgm_family = "mixC") {
    registry <- mixC_dgm_registry()
    if (!dgm_family %in% names(registry)) {
        stop(
            "Unknown dgm_family: ", dgm_family,
            ". Choose one of: ", paste(names(registry), collapse = ", ")
        )
    }
    registry[[dgm_family]]$alias_map
}

resolve_mixC_dgm_variant <- function(dgm_variant, dgm_family = "mixC") {
    alias_map <- mixC_dgm_variant_alias_map(dgm_family = dgm_family)
    if (dgm_variant %in% names(alias_map)) {
        dgm_variant <- unname(alias_map[[dgm_variant]])
    }
    dgm_variant
}

get_mixC_dgm_par <- function(dgm_variant, dgm_family = "mixC") {
    dgm_variant <- resolve_mixC_dgm_variant(dgm_variant, dgm_family = dgm_family)
    dgm_map <- mixC_dgm_variant_map(dgm_family = dgm_family)
    if (!dgm_variant %in% names(dgm_map)) {
        stop(
            "Unknown dgm_variant: ", dgm_variant,
            " for dgm_family=", dgm_family,
            ". Choose one of: ", paste(names(dgm_map), collapse = ", ")
        )
    }
    dgm_map[[dgm_variant]]
}
