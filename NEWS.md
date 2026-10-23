# nutrigon 0.1.0

First release.

## Fitting

* `nutrigon_fit()` fits a quadratic response surface and locates the
  nutritional optimum with a hybrid estimator: the closed-form stationary
  point when it is admissible (inside the design hull, correct curvature
  sign), otherwise a constrained Monte Carlo search summarised by the point
  with the highest (or lowest) fitted value. `model_type` selects ordinary
  least squares (`"lm"`), a generalized linear model with any family
  including quasi-families (`"glm"`; optimum located on the link scale), or
  a linear mixed-effects model (`"lmer"`; population-level optimum;
  requires lme4 and lmerTest).
* `nutrigon_coefs()`, `nutrigon_elongation()`, `nutrigon_snr()` and
  `nutrigon_curvature()` return the fitted quadratic coefficients, the
  curvature elongation ratio, the signal-to-noise ratio, and local
  Gaussian/mean curvature with a surface-area ratio.

## Inference

* `nutrigon_compare()` estimates the angular separation between two optima
  with a percentile bootstrap interval that refits both surfaces in every
  replicate. Resampling is by case (default) or stratified within design
  points, paired or unpaired. Returns the locations of the optima, each
  surface's angular support, the bootstrap distribution and the Monte
  Carlo fallback rate.
* `nutrigon_matrix()` builds the N x N angular trade-off matrix with joint
  bootstrap uncertainty (the same replicate of each surface is used in
  every pair) and the bootstrap stability of a hierarchical clustering.
* `nutrigon_equivalence()` tests whether two optima are aligned within a
  pre-specified angle.

## Coverage check

* `nutrigon_build_grid()` simulates a validation grid, either on the
  user's own design and at the signal-to-noise ratio of the fitted
  surfaces (recommended) or on a generic rail design.
* `nutrigon_check()` estimates the coverage expected for one result from
  the nearest grid configurations, with a Wilson interval, excluding
  near-zero angles and warning when the grid comes from a different
  design.

## Plotting helpers

* `nutrigon_grid()`, `nutrigon_hull()` and `nutrigon_region_polygon()`
  prepare fitted surfaces for ggplot2.

## Data

* `lee_drosophila`, `zanco_survival`, `zanco_wingsize`,
  `zanco_tradeoff_matrix` and `nutrigon_validation_grid`.
