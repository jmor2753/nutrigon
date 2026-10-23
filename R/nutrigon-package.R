#' nutrigon: Angular Separation Between Nutritional Optima
#'
#' Fits response surfaces (ordinary least squares, generalized linear
#' models with any family, or linear mixed-effects models) to trait data
#' collected under the Geometric Framework for Nutrition (GFN) and computes
#' the angular separation (theta) between the nutritional optima of two
#' traits, with bootstrap confidence intervals.
#'
#' @section Core workflow:
#' \describe{
#'   \item{[nutrigon_fit()]}{Fit a quadratic response surface (`"lm"`,
#'     `"glm"`, or `"lmer"`, via `model_type`) and locate the nutritional
#'     optimum.}
#'   \item{[nutrigon_compare()]}{Compare two surfaces (optionally with
#'     different `model_type`s); returns the angular separation and its
#'     bootstrap CI (the main result).}
#'   \item{[nutrigon_coefs()]}{Extract the fitted link-scale quadratic
#'     coefficients from a fit, regardless of model type.}
#'   \item{[nutrigon_elongation()]}{Estimate the curvature elongation ratio
#'     of a fitted surface.}
#'   \item{[nutrigon_curvature()]}{Local Gaussian and mean curvature, and a
#'     surface-area ratio, for a fitted surface.}
#'   \item{[nutrigon_check()]}{Design-matched coverage check for one
#'     result against a validation grid.}
#'   \item{[nutrigon_equivalence()]}{Equivalence test for alignment of two
#'     optima (small angles).}
#'   \item{[nutrigon_matrix()]}{Angular trade-off matrix for N surfaces,
#'     with joint bootstrap uncertainty and clustering stability.}
#' }
#'
#' @section Utilities:
#' \describe{
#'   \item{[nutrigon_build_grid()]}{Build a validation grid by simulation,
#'     on your own design and at your surfaces' signal-to-noise ratio.}
#'   \item{[nutrigon_snr()]}{Signal-to-noise ratio of a fitted surface.}
#'   \item{[nutrigon_grid()], [nutrigon_hull()], [nutrigon_region_polygon()]}{
#'     Plot helpers (require \pkg{ggplot2}).}
#' }
#'
#' @section Example data:
#' `lee_drosophila`, from Lee et al. (2008), and the pre-built validation
#' grid `nutrigon_validation_grid` (`seed = 2025`).
#'
#' @section Optional dependencies:
#' `model_type = "lmer"` in [nutrigon_fit()]/[nutrigon_compare()] requires
#' the \pkg{lme4} and \pkg{lmerTest} packages; the surface-area ratio
#' returned by [nutrigon_curvature()] requires \pkg{pracma}. None of these
#' are required for `"lm"`/`"glm"` fits or for local curvature alone, so
#' they are Suggests, not hard dependencies -- install them with
#' `install.packages(c("lme4", "lmerTest", "pracma"))` if needed.
#'
#' @references
#' Lee, K.P., Simpson, S.J., Clissold, F.J., Brooks, R., Ballard, J.W.O.,
#' Taylor, P.W., Soran, N. & Raubenheimer, D. (2008). Lifespan and
#' reproduction in Drosophila: new insights from nutritional geometry.
#' \emph{PNAS}, 105(7), 2498-2503.
#'
#' Morimoto, J. et al. (2023). Nutrigonometry I: using right-angle triangles
#' to quantify nutritional trade-offs in performance landscapes. \emph{The
#' American Naturalist}, 201(5), 725-740.
#'
#' Morimoto, J. (in review). Bootstrap inference for angular separation
#' between nutritional optima in the Geometric Framework for Nutrition.
#' \emph{Methods in Ecology and Evolution}.
#'
#' @keywords internal
"_PACKAGE"
