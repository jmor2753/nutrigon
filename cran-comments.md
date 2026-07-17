# CRAN submission comments

## Resubmission

This is a resubmission of nutrigon 0.1.0, addressing two issues raised by
CRAN on the previous submission:

1. **DESCRIPTION references.** The `Description` field now formats both
   references as `authors (year) "Title" <doi:...>`, per
   <https://contributor.r-project.org/cran-cookbook/description_issues.html#references>.
   Lee et al. (2008), previously mentioned without a link, is now linked
   via its DOI.
2. **Writing to .GlobalEnv.** `nutrigon_compare()` and
   `nutrigon_build_grid()` previously implemented seed reproducibility
   with a hand-rolled save/restore of `.Random.seed` directly in
   `.GlobalEnv` (a strict save-then-restore-on-exit that left no net
   change, but is still flagged by CRAN's automated check on the
   package's source text -- see
   <https://contributor.r-project.org/cran-cookbook/code_issues.html#writing-to-the-.globalenv>).
   Both functions now use `withr::local_seed()` instead (new `Imports`
   dependency on `withr`). Behaviour is unchanged: results remain fully
   reproducible given the same `seed`, verified by direct comparison
   against the previous implementation's output.

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
  accept an optional `seed` argument for reproducibility, implemented via
  `withr::local_seed()`. Neither function writes to `.GlobalEnv`.
* The package ships five small pre-built data objects (`lee_drosophila`,
  `zanco_survival`, `zanco_wingsize`, `nutrigon_validation_grid`,
  `zanco_tradeoff_matrix`), documented individually in `R/data.R`. The
  latter two are pre-computed simulation/bootstrap results included so
  that `nutrigon_check()` and the second vignette do not need to re-run
  expensive simulations (several minutes) at package load or vignette
  build time.
* All examples run in well under 5 seconds; none are wrapped in
  `\donttest{}` or `\dontrun{}`, with the exception of a small number of
  `model_type = "glm"`/`"lmer"` examples that depend on Suggests
  packages not being guaranteed available on all check platforms.
