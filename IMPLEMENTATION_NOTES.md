# Implementation notes

This file records places where the code's *default* behavior differs from the
algorithm as stated in the paper, plus inactive/optional code paths that are
intentionally kept for future development. It is meant to help anyone comparing
the code against the paper (arXiv:2605.06581) reconcile the two. Line numbers are
approximate and refer to the `code-improvements` branch.

These items were reviewed and **left as-is on purpose** (the code is what produced
the reported results); nothing here is a pending code change.

## Repository history

- `v1.0-original` (tag, commit `e8eeef0`) — the original code as of the
  arXiv:2605.06581v1 submission. Use this for exact reproduction of the paper.
- `v1.1-code-improvements` (tag, commit `1dee016`) — after the code review +
  speedup pass: robustness fixes and result-preserving speedups; the full
  default-path smoke suite verified bit-identical to `v1.0-original`. Full diff:
  https://github.com/wangyuyao98/HAPS/compare/v1.0-original...v1.1-code-improvements
- General-G_tau development (completed `Gtau_mode = "estimated"`/`"tilted"`,
  weighted at-risk evaluation, censoring-shift sensitivity study, `linWB2` DGM,
  `support_upper` candidate bound) — branch `gtau-modes`, reviewed and merged via
  GitHub pull request (see the repository's pull-request history).

## 1. Deviations from the paper's stated algorithm (in the shipped pipeline)

These are deliberate implementation choices that differ from the paper text and
that would change reported numbers if altered. Documented here for transparency.

### 1.1 IPCW censoring-weight trimming `trim.C = 0.05`
- Where: `src/dynamicCP.R` — `dynamicCP_split(..., trim.C = 0.05)` (default, arg list
  line ~30). The floor is applied to the estimated censoring survival `Ĝ` in both
  the θ-calibration estimating function and the reported IPCW coverage
  (`pmax(surv.C.subj, trim.C)`, lines ~536, ~585-586). The AIPCW path applies an
  analogous floor via `trim.G`/`trim.S`/`trim.Gtau`.
- Paper: Algorithm 1 (Step 4) and the coverage estimator eq. (24) use `Ĝ`
  directly, with no trimming constant mentioned.
- Effect: caps IPCW weights at `1/0.05 = 20`; changes θ selection and reported
  coverage only when `Ĝ(X_i) < 0.05` (heavy-censoring tails). A standard,
  variance-reducing stabilization, but not described in the paper.

### 1.2 Last-interval `S_K` fit uses event-only IPCW, not a plain right-censored fit
- Where: `src/helpers_AIPCW.R` — `fit_Sk_list(..., fit_sample_mode = "event_only_ipcw")`
  (default, choices at line ~1541; the per-interval fit at line ~1612 subsets to
  `Δ == 1` with weights `1/H_k(X)` for every k, including `k = K`). The default is
  set in `dynamicCP_AIPCW_split(S_fit_sample_mode = "event_only_ipcw")` (line ~175)
  and is not overridden by any driver.
- Paper: Appendix C.2.2 ("Estimation of S_k") specifies the IPCW-weighted
  event-only fit for `k = 1, ..., K-1`, and for the **last interval** a *standard
  right-censored survival regression* among subjects with `X > t_K` (no IPCW
  weighting, since no future covariates remain).
- Effect: `S_K` (which drives `h_τ` for all `u` in the last interval — the only
  interval when `τ = t_K`, e.g. τ=6 with K=3) is a different estimator than the
  paper's, and it makes `S_K` depend on the censoring model. The code already
  contains the paper's behavior under `fit_sample_mode = "hybrid_interval_ipcw"`
  (which fits `all_at_risk` at `k = K`); it is simply not the default.

### 1.3 AIPCW θ selector defaults to `first_crossing`, not the paper's `inf` rule
- Where: `src/dynamicCP_AIPCW.R` — `AIPCW_theta_selector = "first_crossing"`
  (default, choices at line ~193). The alternative `"rightmost_nonnegative"`
  implements `θ̂ = inf{θ : Σ U_AIPCW ≥ 0}` under the code's θ→(1−θ) grid
  parametrization.
- Paper: Section 4 says HAPS-A replaces `U` by `U_AIPCW` in Step 4 of Algorithm 1,
  i.e. the `inf` rule.
- Effect: the two agree when the AIPCW estimating curve crosses zero once; they can
  differ when it is non-monotone (the augmentation term need not be monotone in θ),
  in which case `first_crossing` yields a wider/more conservative interval.

### 1.4 `ξ_k` silently falls back to a constant weighted mean under small-sample guards
- Where: `src/helpers_AIPCW.R` — `fit_Yk_condmean_at_tk(..., Xi_fallback_to_mean = TRUE)`
  with `Xi_min_n`, `Xi_min_ess`, `Xi_min_obs_per_coef`, `Xi_min_class_n`,
  `Xi_min_class_ess` guards. Defaults set in `dynamicCP_AIPCW_split` (line ~152).
- Paper: Appendix F.3.2 specifies an XGBoost `multi:softprob` model for `ξ_k`.
- Effect: when the uncensored landmark fitting sample is too small (plausible in the
  small-n appendix runs, e.g. n=300 at τ=6), `ξ_k` loses covariate dependence and
  becomes a constant. Consider reporting the fallback frequency when using small n.

## 2. Inactive / optional code paths kept for future development

These are intentionally retained (not dead weight to remove) — they are hooks for
alternative nuisance-estimation options under development. Listed so reviewers know
they are currently unreached.

- `model.C = "xgb_piecewise"` in `src/dynamicCP.R` (branch at line ~406) depends on
  `fit_Gk_list_xgb_piecewise()` / `compute_Hk_at_tk_xgb_piecewise()`, which are not
  defined in the repository; selecting it currently errors. (A piecewise-XGB
  analogue of eq. (10) / Appendix C.2.1.)
- `grf` backends: `fit_Gk_list` / `fit_Sk_list` / `fit_Yk_condmean_at_tk` contain
  `grf` branches that are excluded from their `match.arg` choices, so they are
  unreachable; `main_simu_dynamic_DDH_tuning.R` (lines ~95-96) accepts `grf` as a
  CLI value, so selecting it would error until `grf` is added to those choices.
- `censoring_method = "global"` and `model.pred = "grf"` blocks in
  `src/dynamicCP_AIPCW.R` are similarly gated off.
- Duplicated helpers exist in more than one sourced file (`surv_quantile`,
  `prepare_data_tau`, `build_static_from_dynamic`, and an older `keep_nonconstant_covars`
  in `dynamicCP.R`). The active copy depends on `source()` order; they currently
  agree except the old `keep_nonconstant_covars`, which drops character-typed
  covariates (all covariates in the shipped pipelines are numeric or factor, so this
  is inert). Consolidate before relying on character covariates.

## 3. DGM naming vs the paper

- Paper DGM1 = code `linWB1` (`linear_weibull`); paper DGM2 = code `mixC2_DGM1`;
  paper DGM3 = code `mixC2_DGM3` (see `src/mixC_dgm_registry.R`). The registry also
  ships some unused near-miss presets (e.g. a `mixC2/DGM2` entry that matches no
  paper DGM, and the `mixC` family); the parameter sets actually selected by the
  drivers reproduce the paper's reported censoring/survival marginals.
- A hard-coded `mixC` block in `main_simu_dynamic_linWB1_cox_multiclass_loc0.R`
  (line ~164) uses `uniform_cmax = 12`, which does not match paper DGM2's
  `Unif(0, 15)`. It is unreachable in shipped runs because `dgm_name` is fixed to
  `"linear_weibull"` in that script; only relevant if that block is ever activated.

## 4. Optional speed flag: `fast_rsf_predict` (default OFF)

`predict.rfsrc()` consumes the R RNG stream on every call, so the number of RSF
prediction calls affects the seeds of all downstream randomized fits in the same
replication. The AIPCW pipeline can hoist RSF predictions (one call per interval /
per τ instead of one per evaluation time / per θ); the predicted values are
identical, and in the DGM2 DDH smoke run this cut the driver from ~100s to ~41s
(~2.4x), but results become bitwise-different (statistically equivalent) from the
published `.rds` files.

- Default `FALSE` — bit-identical to the original implementation (use for exact
  reproduction of the paper's results).
- Enable per call via `dynamicCP_AIPCW_split(..., fast_rsf_predict = TRUE)` or
  globally via the environment variable `HAPS_FAST_RSF_PREDICT=TRUE` (e.g.
  `HAPS_FAST_RSF_PREDICT=TRUE Rscript main_simu_dynamic_DDH_tuning.R ...`) for new
  experiments where exact bitwise reproduction is not required.

## 5. Gtau evaluation options (`gtau_eval_method`, default `"weighted"`)

For `Gtau_mode = "estimated"` / `"tilted"`, the simulation drivers evaluate coverage on
the at-risk population {T>τ, C̃>τ} (C̃ ⊥ T | covariate path; C̃ ~ the true censoring law,
tilted by exp(δt) under `"tilted"`). Because C̃ enters only through 1(C̃>τ), conditioning on
the subgroup is equivalent to weighting uncensored survivors by w = G̃_δ(τ|path) — the
weighted estimator is the Rao–Blackwellization of physically constructing the subgroup.

`gtau_eval_method` in `main_simu_dynamic_linWB1_cox_multiclass_loc0.R` (persisted in the
result config and honored by `main_derive_HAPS_DR_linWB1.R`):

- `"weighted"` (default): oracle-weighted coverage on uncensored T>τ survivors →
  `coverage_test_gtau` (δ=0 weights for `"estimated"`). Chosen as default for consistency
  across the two modes and lowest variance — in the evaluator comparison
  (`validation/compare_gtau_evaluators.R`; M=300 test sets of n=500) all strategies were
  unbiased with SDs: weighted 0.0153 < physical-C̃-sampled 0.0159 < ratio-from-subgroup
  0.0170 (tilted target), and weighted 0.0155 < subgroup 0.0162 (estimated target).
- `"subgroup"`: physical at-risk construction.
  * estimated: legacy behavior — censored test data, subgroup X>τ, coverage in
    `coverage_test_true` (pre-existing outputs reproduce exactly).
  * tilted: membership indicators I(C̃>τ) ~ Bernoulli(w) (the indicator's exact law) drawn
    on an isolated deterministic RNG stream (calibration outputs bit-identical to a
    weighted run); subgroup mean → `coverage_test_gtau`, weighted counterpart →
    `coverage_test_gtau_wt`; bounds carry `at_risk_gtau` alongside `w_gtau`.

The identity behind the weighted evaluator is validated end-to-end (full C̃ trajectories
drawn from the tilted true law by inverse-CDF, never using the weights) in
`validation/validate_gtau_weighted_evaluation.R`. True-law weights are implemented for the
linWB1 DGM only (`src/gtau_eval_helpers.R`); other DGMs currently require
estimated+subgroup.

## 6. DDH prediction-model caveats

- CPU threads: `src/ddh_backend.py` now pins `torch.set_num_threads` (default 1,
  override via `args["torch_num_threads"]`) for cross-run/cross-machine
  reproducibility of the neural model.
- Extrapolation beyond the training horizon: test residual times past the last
  training time-bin are clipped to that bin (`ddh_backend.py`, bin assignment),
  so survival beyond the training range reuses the last estimated hazard.
- Dropout placement: a single-layer GRU cannot take PyTorch inter-layer dropout, so
  the specified 0.1 dropout is applied in the prediction head instead.
