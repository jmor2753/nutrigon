# nutrigon <a href="https://github.com/jmor2753/nutrigon"><img src="man/figures/logo.png" align="right" height="139" alt="nutrigon logo" /></a>

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
<!-- badges: end -->

## Overview

**nutrigon** fits quadratic response surfaces to trait data collected
under the Geometric Framework for Nutrition (GFN) and computes the
angular separation (theta) between the nutritional optima of two traits,
with percentile bootstrap confidence intervals that refit the surfaces in
every replicate. Theta measures the difference in nutrient balance
(direction) between optima; it does not measure their radial distance
(total intake) or the fitness cost of a trade-off, so report it with the
locations of the optima.

The package also provides a design-matched coverage check
(`nutrigon_build_grid()` + `nutrigon_check()`), which simulates validation
configurations on your own design and at the signal-to-noise ratio of your
fitted surfaces, and an equivalence test for alignment of optima
(`nutrigon_equivalence()`), since the unsigned interval cannot by itself
show that two optima coincide.

## Installation

```r
# from a local source tarball
install.packages("path/to/nutrigon", repos = NULL, type = "source")
```

## Core workflow

```r
library(nutrigon)
data(lee_drosophila)

fit1 <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                      "lifespan", design_id_col = "treatment")
fit2 <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                      "dailyeggs", design_id_col = "treatment")

result <- nutrigon_compare(
  lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
  lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
  paired = TRUE, B = 500, seed = 1        # case resampling (validated)
)
print(result)

# Design-matched coverage check (slow: build once, then saveRDS())
grid <- nutrigon_build_grid(
  design = lee_drosophila, x_col = "protein_eaten", y_col = "carb_eaten",
  design_id_col = "treatment",
  snr = c(nutrigon_snr(fit1), nutrigon_snr(fit2)), seed = 99
)
nutrigon_check(result,
               elong1 = nutrigon_elongation(fit1),
               elong2 = nutrigon_elongation(fit2),
               validation_grid = grid)

# Are two optima aligned within 10 degrees?
nutrigon_equivalence(result, delta = 10)
```

See `vignette("getting-started")` for the full walkthrough, and
`vignette("zanco-tradeoff-matrix")` for the N x N trade-off matrix,
clustering, and overlay-plot example built on Zanco et al. (2025).

## Example data

- `lee_drosophila` -- Lee et al. (2008) Drosophila lifespan/fecundity data.
- `zanco_survival`, `zanco_wingsize` -- Zanco et al. (2025) viability and
  wing size data.
- `nutrigon_validation_grid` -- a generic pre-built validation grid
  (seed = 2025, 5-rail design) for trying `nutrigon_check()`; build a
  design-matched grid for applied results.
- `zanco_tradeoff_matrix` -- the pre-computed 12x12 pairwise angular
  trade-off matrix for the Zanco example.

## Citation

> Morimoto, J. (in review). Bootstrap inference for angular separation
> between nutritional optima in the Geometric Framework for Nutrition.
> *Methods in Ecology and Evolution*.

And the original Nutrigonometry paper:

> Morimoto, J. et al. (2023). Nutrigonometry I: using right-angle
> triangles to quantify nutritional trade-offs in performance landscapes.
> *The American Naturalist*, 201(5), 725-740.

## License

MIT

## Author

Juliano Morimoto ([ORCID: 0000-0003-3561-1920](https://orcid.org/0000-0003-3561-1920))
