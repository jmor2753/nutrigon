# nutrigon 0.1.0

Initial release.

## Core workflow

* `nutrigon_fit()` fits a response surface and locates the nutritional
  optimum (closed-form where possible, with a Monte Carlo fallback near
  the boundary of the design). A `model_type` argument (`"lm"`, `"glm"`,
  or `"lmer"`) selects the underlying fitter:
  - `"lm"` (default): ordinary least squares.
  - `"glm"`: any family, including quasi-families (`quasipoisson()`,
    `quasibinomial()`) for overdispersed counts and proportions.
    Curvature, elongation, and the optimum are all located on the
    **link scale**, since that is the scale on which the fitted surface
    is quadratic.
  - `"lmer"`: a linear mixed-effects model via \pkg{lme4}/\pkg{lmerTest}
    (both Suggests -- install with
    `install.packages(c("lme4", "lmerTest"))` if needed), locating the
    population-level optimum. Requires `re_formula`, e.g.
    `re_formula = "+ (1 | block)"`, and `re_group_cols` naming the
    original column(s) that formula refers to (e.g.
    `re_group_cols = "block"`), since fitting internally reduces the
    data to just the x/y/z columns beforehand. Any valid \pkg{lme4}
    random-effects syntax is supported, including crossed and nested
    terms and random slopes (see `?nutrigon_fit` for the `.x`/`.y`
    naming convention random slopes must use).
* `nutrigon_compare()` computes the angular separation (theta) between
  two fitted optima, with a bootstrap percentile confidence interval.
  Supports both paired (same sampling units) and unpaired designs, and
  the two surfaces can use different `model_type`s (e.g. a binomial GLM
  for viability compared against an `lm` fit for a continuous trait).
* `nutrigon_coefs()` is an accessor for the fitted link-scale quadratic
  coefficients, working identically across all three model types.
* `nutrigon_elongation()` estimates the curvature elongation ratio of a
  fitted surface from its Hessian eigenvalues (link scale for `"glm"`
  fits).
* `nutrigon_curvature()` computes local Gaussian and mean curvature (via
  the first and second fundamental forms of the fitted surface treated
  as a 2-manifold) at every point of a grid restricted to the design
  hull, plus a surface-area ratio comparing the fitted surface's area to
  a flat reference surface over the same design. Available for all three
  model types via the same coefficient-extraction layer used by
  `nutrigon_elongation()`. The surface-area ratio requires \pkg{pracma}
  (Suggests).
* `nutrigon_check()` provides a per-result reliability diagnostic,
  comparing a result's risk-factor profile against a simulation-based
  validation grid of configurations with known ground-truth bootstrap
  coverage. Defaults to the validation grid shipped with the package
  (`nutrigon_validation_grid`, `seed = 2025`) for immediate use without
  running a fresh simulation.
* `nutrigon_build_grid()` builds a validation grid from scratch, for
  users working with a non-standard design.

## Plotting helpers

* `nutrigon_grid()`, `nutrigon_hull()`, and `nutrigon_region_polygon()`
  prepare a fitted surface for visualisation with ggplot2.

## Example data

* `lee_drosophila` -- Lee et al. (2008) Drosophila lifespan/fecundity
  data, used in `vignette("getting-started")`.
* `zanco_survival`, `zanco_wingsize`, and `zanco_tradeoff_matrix` --
  Zanco et al. (2025) viability and wing size data, and a pre-computed
  12x12 pairwise angular trade-off matrix, used in
  `vignette("zanco-tradeoff-matrix")`.
* `nutrigon_validation_grid` -- a pre-built simulation validation grid
  (`seed = 2025`, 1000 configurations) for `nutrigon_check()`.

## Notes

* `model` (a component of the object returned by `nutrigon_fit()`) is a
  `nutrigon_lm_fit`/`nutrigon_glm_fit`/`nutrigon_lmer_fit` wrapper, not a
  bare `lm`/`glm`/`merMod` object. The raw fitted model is available as
  `fit$model$fit` (e.g. for `summary()`, residual diagnostics, or
  `AIC()`); use `nutrigon_coefs(fit)` rather than calling `coef()` on
  `fit$model` directly.
* `nutrigon_compare()` and `nutrigon_build_grid()` implement their
  optional `seed` argument via `withr::local_seed()`; neither writes to
  `.GlobalEnv`, and results remain fully reproducible given the same
  seed.
