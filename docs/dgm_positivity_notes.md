# Why linWB2: censoring positivity, bounded support, and the widest candidate set

Notes for collaborators explaining the design of the `linWB2` DGM (used in the
general-G_tau paper's simulations) and why `linWB1` (the original paper's DGM1)
is kept unchanged alongside it. All numbers below are Monte-Carlo estimates on
100k draws; scripts to reproduce them live in the session notes and can be
re-derived from `src/gtau_eval_helpers.R` + `src/linWB_dgm_registry.R`.

## 1. The problem that started this

At small n and late prediction times (n=300, τ=6), a few percent of simulation
replications ended with Algorithm 1 **infeasible**: no θ on the grid satisfied
the calibration inequality. Replaying every failure showed a single mechanism:
the candidate intervals are quantile bands read off the **training** prediction
grid, so the widest candidate is (first training event time after τ, last
training event time after τ) — and some uncensored **calibration** subjects fall
outside that range (e.g., one failed replication had 13 of 78 calibration events
beyond the training maximum). No candidate set can ever cover them, and with
IPCW/G̃ weighting a few such subjects push the estimating function negative on
the whole grid.

Fix requirement: the candidate family must contain a "large enough" widest set
(the feasibility premise of the coverage theorem), and — per the theory — the
family must be **fixed given the training data only**. A data-dependent widest
set (e.g., (τ, max pooled event time]) restores feasibility but makes one member
depend on the calibration data; workable with a conformal-max argument, but not
as clean.

## 2. The cleanest resolution: bounded, known support

If T has bounded support with known upper bound T_max, the widest candidate can
be **(τ, T_max] — a known constant**: the family is fully training-fixed and
Algorithm 1 is feasible by construction. We implement this by **truncated
generation** (draw T from its law conditional on T ≤ T_max via inverse-CDF in
the last interval): no point mass at T_max, T stays absolutely continuous, and
every draw satisfies T < T_max strictly. (Naive min-capping T∧T_max would create
an atom at T_max — ties, broken absolute-continuity, and open-interval issues.)

Under linWB1's event law, T_max = 20 truncates only **0.31%** of the event mass
(0.36% among τ=6 survivors), versus 1.95% at T_max = 15.

## 3. The catch: censoring positivity constrains more than T_max

Assumption 2.3(ii) needs G(T | path) = P(C > T | path) ≥ η on the support. Under
**linWB1's censoring law** (β_CL = 2: censoring hazard scales as e^{2L}) this
fails badly regardless of any T_max choice:

- G(T_max | path) across the path distribution: at T_max = 12 the **median**
  path has G = 0.11 and 42% of paths are below 0.05; at T_max = 20, 84% are.
- The operative quantity G(T | path) at the subject's own event time: minimum
  ≈ 1e-86; **~4.4% of subjects below 0.05** — essentially independent of T_max
  (it is mid-range events on high-L paths, not the far tail of T).

Consequence worth stating plainly: **linWB1 has always relied on the weight
trimming `trim.C = 0.05`** (documented as deviation §1.1 in
IMPLEMENTATION_NOTES.md) to absorb a real positivity violation affecting ~4–6%
of subjects. That is a fine benchmark to keep — it matches the published paper
— but it is not a setting where the censoring assumptions genuinely hold.

## 4. linWB2: retuned censoring so positivity actually holds

Keeping the event law identical (so T-marginals match linWB1 up to the 0.31%
truncation), we retuned only the censoring:

| | linWB1 | **linWB2** |
|---|---|---|
| censoring shapes a_C | (3, 3, 3) | (3, 3, **2**) |
| β_C0 | (−6, −6, −5) | (**−5, −5, −5.8**) |
| β_CL | (2, 2, 2) | (**0.3, 0.3, 0.3**) |
| T support | unbounded | **(0, 20)**, truncated generation |
| censoring rate | ~36% | **~31%** |
| G(T_max\|path), 0.01% path quantile | ~0 | **0.104** |
| P(G(T_max\|path) < 0.05) | 0.84 | **0 / 100k** |
| widest candidate set | training-grid extremes | **(τ, 20], known constant** |

The two levers: the flatter last-interval shape (t² instead of t³) slows the
late hazard growth, and β_CL = 0.3 bounds the across-path hazard ratio (e^{0.3L}
instead of e^{2L}) — this is precisely what makes a uniform-in-path η possible.
The trade-off, stated openly: linWB2's censoring depends on covariates more
mildly than linWB1's. For studies that *want* violent covariate–censoring
dependence (e.g., stress settings where baseline methods visibly fail), design
that as its own DGM; linWB2 is the assumptions-hold benchmark.

Residual caveat: L is Gaussian (unbounded), so "η ≈ 0.10" is an all-but-
negligible-mass statement (0.01%-quantile over 100k paths), not a literal a.s.
bound; truncating L at ±4 (mass 6e-5) would make it exact if ever needed.

## 5. Code structure

- `linWB1` is untouched; all its results (paper + sensitivity runs under
  `results/linWB1/`) remain exactly reproducible.
- `linWB2` is an additive preset in `src/linWB_dgm_registry.R`; the generator
  gains an opt-in `T_max` (default `Inf` = byte-identical legacy behavior); the
  calibration gains an opt-in known-support widest set (τ, T_max]; outputs go
  to `results/linWB2/`.
