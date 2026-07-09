## Registry of linear-Weibull ("linWB") DGM presets.
##
## linWB1 -- the original paper DGM1 (UNCHANGED; kept for comparability with the
##   published results and the existing sensitivity runs under results/linWB1/).
##   Note: under this DGM the censoring positivity assumption holds only up to
##   trimming -- G(T|path) < 0.05 for ~4-6% of subjects (see
##   docs/dgm_positivity_notes.md) -- which is absorbed by trim.C = 0.05.
##
## linWB2 -- positivity-respecting variant for the general-Gtau paper:
##   * event-time law IDENTICAL to linWB1, but generated from its conditional
##     law given T <= T_max = 20 (truncated inverse-CDF; no atom; removes 0.31%
##     tail mass), so the support is bounded and known;
##   * censoring retuned (last-interval shape 3 -> 2; beta_C0 = (-5,-5,-5.8);
##     beta_CL = 0.3) so that G(T_max | path) >= ~0.10 across essentially the
##     entire covariate-path distribution (empirical 0.01%-quantile 0.104 over
##     100k paths) with censoring rate ~31%;
##   * the widest candidate prediction set can then be taken as (tau, T_max]
##     -- a KNOWN constant, so the candidate family is fully fixed given the
##     training data and Algorithm 1 is feasible by construction.
get_linWB_dgm_par <- function(setup = c("linWB1", "linWB2")) {
    setup <- match.arg(setup)
    switch(
        setup,
        linWB1 = list(
            par = list(
                shape_T = c(4, 4, 5),
                shape_C = c(3, 3, 3),
                beta_T0 = c(-8, -8, -5),
                beta_TL = c(1, 2, 3),
                beta_C0 = c(-6, -6, -5),
                beta_CL = c(2, 2, 2)
            ),
            T_max = Inf,
            change_times = c(3, 6),
            rho = 0.3,
            description = "Paper DGM1 (unchanged); positivity via trim.C only."
        ),
        linWB2 = list(
            par = list(
                shape_T = c(4, 4, 5),
                shape_C = c(3, 3, 2),
                beta_T0 = c(-8, -8, -5),
                beta_TL = c(1, 2, 3),
                beta_C0 = c(-5, -5, -5.8),
                beta_CL = c(0.3, 0.3, 0.3)
            ),
            T_max = 20,
            change_times = c(3, 6),
            rho = 0.3,
            description = paste0(
                "Positivity-respecting variant: same event law truncated at T_max=20; ",
                "censoring retuned so G(T_max|path) >= ~0.10 (eta) with ~31% censoring; ",
                "widest candidate set (tau, T_max].")
        )
    )
}
