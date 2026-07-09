# General G_tau options: research directions and planned extensions

Working notes for the follow-up paper on generalized `Gtau_mode` options
(`"one"`, `"estimated"`, `"tilted"`) in HAPS-A. Status: **planning document** —
the DGM2 extension and the studies below are scoped but intentionally not yet
implemented (author decision 2026-07-08). See `IMPLEMENTATION_NOTES.md` §5 for
what IS implemented (evaluation options for DGM1/linWB1).

## The estimand

For a specified residual at-risk mechanism C̃ (independent of T given the
covariate path), target

    P( T ∈ Ĉ(Z̄_τ) | T > τ, C̃ > τ ).

The three modes instantiate C̃ as: `"one"` → C̃ ≡ ∞ (all survivors; the original
HAPS estimand), `"estimated"` → the observed censoring law, `"tilted"` → the
exp(δt)-tilted law. The tilted family indexes a neighborhood of at-risk
mechanisms, i.e. a sensitivity analysis over "who is still under follow-up at
the decision time"; δ>0 shifts censoring later (less dropout by τ), δ<0 earlier.

## What is needed to run estimated/tilted + sensitivity under DGM2 (`mixC2_DGM1`)

DGM2 censoring is a subgroup mixture (`.simulate_dataset_long_mixC_uniform_cox`
in `src/gen_ICML_simu.R`): B~Bern(0.5); B=0 → C~Unif(0, cmax=15); B=1 →
interval-wise Weibull renewal with scale exp(−(β_C0+β_CL·L+β_CL2·L²)/a).
The oracle-weight machinery is currently DGM1-specific, so:

1. `return_L_full` flag for the mixC generator (+ allow it in the
   `simulate_dataset_long` dispatcher; B is already emitted). Default off ⇒
   existing outputs unchanged.
2. mixC true-law weight helpers in `src/gtau_eval_helpers.R`, branching on B:
   - B=0 (Uniform): closed forms G(τ)=(cmax−τ)/cmax and
     G̃_δ(τ)=(e^{δ·cmax}−e^{δτ})/(e^{δ·cmax}−1).
   - B=1: reuse the interval-integral logic with the quadratic linear predictor;
     δ=0 must equal the closed form (unit test, as for DGM1).
3. A DGM dispatch wrapper (e.g. `compute_Gtau_true_tilted(dgm_name, ...)`) so
   drivers select the right law.
4. Generalize `main_simu_gtau_tilt_sensitivity.R` to a DGM argument. Default
   prediction model for the DGM2 study: fast Cox/RSF (the study's object is the
   censoring/G_tau law, not the predictor); DDH behind a flag for a
   paper-faithful Setup B run.
5. (optional) tilted evaluation wiring in `main_simu_dynamic_DDH_tuning.R`
   (it currently passes Gtau_mode to calibration only).

Note: DGM2's Cox censoring model is misspecified for the mixture — DGM2 is the
natural setting where the G_tau choice and censoring-model robustness can matter.

## Candidate simulation studies

| Study | DGM / lever | knob | metric | expected story |
|---|---|---|---|---|
| A validity | DGM1 & DGM2 | G_tau ∈ {one,estimated,tilted} × τ | coverage on {T>τ,C̃>τ} | all ≈ nominal |
| B sensitivity | DGM1 (robust) vs DGM2 (mixture) | δ_cal × δ_eval | coverage & width vs δ | DGM1 flat; DGM2 shows gaps |
| C when-it-matters | DGM2, amplify B=0/B=1 contrast | one vs estimated | coverage gap of `one` | estimated fixes a real gap |
| D double robustness | DGM2, misspecify G and/or S,ξ | method × nuisance-correctness | coverage | AIPCW robust where IPCW isn't |
| E transportability | source vs shifted follow-up | tilted δ = shift | coverage under transport | tilt restores validity |
| F informativeness | any, fixed mode | τ | q25/median/q75 width | narrows with τ |

Details:
- **A. Validity across target populations.** Establish nominal coverage for each
  G_tau mode across DGM1 (correct Cox) and DGM2 (mixture). Figure: coverage vs τ,
  faceted by mode × DGM.
- **B. Sensitivity to the at-risk mechanism (headline).** Tilted δ as a
  principled sensitivity/transportability parameter; matched vs mismatched
  δ_cal. The existing sensitivity driver already produces the δ_cal×δ_eval sweep
  for DGM1 (preliminary R=40 run: mild sensitivity — a robustness story — even
  under ±25–34pp at-risk shifts).
- **C. When does accounting for the at-risk mechanism matter?** Construct a DGM
  where `one` visibly miscovers on the true at-risk population and
  estimated/tilted corrects it — lever: informative censoring selecting a
  subpopulation with different (T|path) tail behavior that the predictor misses
  (DGM2's B-subgroups are a ready, already-implemented lever).
- **D. Double robustness under the general estimand.** HAPS vs HAPS-A vs
  HAPS-DR crossing correct/misspecified G and (S, ξ); does the AIPCW
  augmentation retain double robustness once g_tau ≠ 1 enters the estimating
  equation?
- **E. Transportability framing.** Tilted G_tau ↔ target population with a
  specified follow-up regime; connect to covariate-shift conformal prediction.
- **F. Informativeness vs τ.** Interval width vs τ WITHIN a fixed G_tau target
  (never across modes — different modes target different populations, so
  cross-mode width comparisons are not an efficiency statement). Requires
  restoring q25/median/q75 and an unweighted per-τ width summary among T>τ
  survivors. Preliminary (R=40, `one` arm): median width 9.90 (τ=0) → 7.77
  (τ=3) → 2.24 (τ=6).

## Related decisions on record

- File naming: keep existing `main_*.R` names; add plain-English header
  docblocks + a README script map instead of renaming (renames would break
  README/derive/plot references and the arXiv-linked history).
- Evaluation: `gtau_eval_method = "weighted"` default; see
  `IMPLEMENTATION_NOTES.md` §5 and `validation/`.
