# nutrigon 0.1.0

Initial release.

## Core workflow

* `nutrigon_fit()` fits a quadratic response surface and locates the
  nutritional optimum (closed-form where possible, with a Monte Carlo
  fallback near the boundary of the design).
* `nutrigon_compare()` computes the angular separation (theta) between
  two fitted optima, with a bootstrap percentile confidence interval.
  Supports both paired (same sampling units) and unpaired designs.
* `nutrigon_elongation()` estimates the curvature elongation ratio of a
  fitted surface from its Hessian eigenvalues.
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
