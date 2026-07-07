# CRAN submission comments

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.5.0
* (add further environments here before submission, e.g. via
  `rhub::check_for_cran()`, GitHub Actions `r-lib/actions/check-r-package`
  on windows-latest / macOS-latest / ubuntu-latest, and win-builder)

## R CMD check results

0 errors | 0 warnings | 0 notes

(Update this section with the actual output of your final
`R CMD check --as-cran` run before submitting -- see below.)

## Downstream dependencies

This is a new submission; there are no reverse dependencies.

## Additional notes for CRAN maintainers

* Two exported functions (`nutrigon_compare()`, `nutrigon_build_grid()`)
  accept an optional `seed` argument for reproducibility. Both save and
  restore the caller's `.Random.seed` on exit, so calling them does not
  alter the caller's global RNG state.
* The package ships five small pre-built data objects (`lee_drosophila`,
  `zanco_survival`, `zanco_wingsize`, `nutrigon_validation_grid`,
  `zanco_tradeoff_matrix`), documented individually in `R/data.R`. The
  latter two are pre-computed simulation/bootstrap results included so
  that `nutrigon_check()` and the second vignette do not need to re-run
  expensive simulations (several minutes) at package load or vignette
  build time.
* All examples run in well under 5 seconds; none are wrapped in
  `\donttest{}` or `\dontrun{}`.
