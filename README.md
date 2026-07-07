# nutrigon <a href="https://github.com/jmor2753/nutrigon"><img src="man/figures/logo.png" align="right" height="139" alt="nutrigon logo" /></a>

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
<!-- badges: end -->

## Overview

**nutrigon** fits quadratic response surfaces to trait data collected
under the Geometric Framework for Nutrition (GFN) and computes the
angular separation (theta) between the nutritional optima of two traits,
with bootstrap confidence intervals. It also provides a simulation-based
per-result reliability diagnostic that flags whether a given result
resembles configurations that achieved nominal coverage or dropped below nominal in extensive
simulation (see Morimoto, in review, *Environmental and Ecological
Statistics*).

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
  design_id_col1 = "treatment", design_id_col2 = "treatment",
  paired = TRUE, B = 500, seed = 1
)
print(result)

check <- nutrigon_check(result,
                         elong1 = nutrigon_elongation(fit1),
                         elong2 = nutrigon_elongation(fit2))
check$flag
```

See `vignette("getting-started")` for the full walkthrough, and
`vignette("zanco-tradeoff-matrix")` for the N x N trade-off matrix,
clustering, and overlay-plot example built on Zanco et al. (2025).

## Example data

- `lee_drosophila` -- Lee et al. (2008) Drosophila lifespan/fecundity data.
- `zanco_survival`, `zanco_wingsize` -- Zanco et al. (2025) viability and
  wing size data.
- `nutrigon_validation_grid` -- a pre-built simulation validation grid
  (seed = 2025) for `nutrigon_check()`.
- `zanco_tradeoff_matrix` -- the pre-computed 12x12 pairwise angular
  trade-off matrix for the Zanco example.

## Citation

> Morimoto, J. (in review). Bootstrap inference for angular separation
> between constrained optima of nutritional ecology response surfaces.
> *Environmental and Ecological Statistics*.

And the original Nutrigonometry paper:

> Morimoto, J. et al. (2023). Nutrigonometry I: using right-angle
> triangles to quantify nutritional trade-offs in performance landscapes.
> *The American Naturalist*, 201(5), 725-740.

## License

MIT

## Author

Juliano Morimoto ([ORCID: 0000-0003-3561-1920](https://orcid.org/0000-0003-3561-1920))

## License

MIT

## Author

Juliano Morimoto ([ORCID: 0000-0003-3561-1920](https://orcid.org/0000-0003-3561-1920))
