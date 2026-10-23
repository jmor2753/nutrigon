## Resubmission

This is a resubmission of nutrigon 0.1.0. Changes in response to the
previous review:

* The `Description` field now formats both references as
  `authors (year) "Title" <doi:...>`, and Lee et al. (2008) is linked via
  its DOI.
* The package no longer modifies `.Random.seed` in `.GlobalEnv`. Optional
  `seed` arguments in `nutrigon_compare()`, `nutrigon_matrix()` and
  `nutrigon_build_grid()` are applied locally with `withr::local_seed()` /
  `withr::with_seed()` (`withr` added to Imports).

Since the previous submission the package also gained `nutrigon_matrix()`,
`nutrigon_equivalence()` and `nutrigon_snr()` (see NEWS.md).

## Test environments

* local: macOS Tahoe 26.6.2 (aarch64-apple-darwin23), R 4.6.1
* local: Ubuntu 24.04 (x86_64-pc-linux-gnu), R 4.3.3
* win-builder: R-devel            <- confirm after devtools::check_win_devel()

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Juliano Morimoto <juliano.morimoto@gmail.com>'
  New submission

## Notes for the CRAN team

* Long-running examples (bootstrap simulations) are wrapped in
  `\donttest{}`; the package passes `R CMD check --as-cran`, which runs
  them with `--run-donttest`.
* Suggested packages are used conditionally: `lme4`/`lmerTest` (mixed
  models) and `pracma` (surface-area ratio) are checked with
  `requireNamespace()` before use, and plotting code in the vignettes
  (`ggplot2`, `ggdendro`, `ggrepel`) is only evaluated when those packages
  are installed.
