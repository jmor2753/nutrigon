#' nutrigon: Angular Separation Between Nutritional Optima
#'
#' Fits quadratic response surfaces to trait data collected under the
#' Geometric Framework for Nutrition (GFN) and computes the angular
#' separation (theta) between the nutritional optima of two traits, with
#' bootstrap confidence intervals.
#'
#' @section Core workflow:
#' \describe{
#'   \item{[nutrigon_fit()]}{Fit a quadratic response surface and locate the
#'     nutritional optimum.}
#'   \item{[nutrigon_compare()]}{Compare two surfaces; returns the angular
#'     separation and its bootstrap CI (the main result).}
#'   \item{[nutrigon_elongation()]}{Estimate the curvature elongation ratio
#'     of a fitted surface.}
#'   \item{[nutrigon_check()]}{Per-result reliability diagnostic against a
#'     pre-built validation grid.}
#' }
#'
#' @section Utilities:
#' \describe{
#'   \item{[nutrigon_build_grid()]}{Build a validation grid from
#'     simulation.}
#'   \item{[nutrigon_grid()], [nutrigon_hull()], [nutrigon_region_polygon()]}{
#'     Plot helpers (require \pkg{ggplot2}).}
#' }
#'
#' @section Example data:
#' `lee_drosophila`, from Lee et al. (2008), and the pre-built validation
#' grid `nutrigon_validation_grid` (`seed = 2025`).
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
#' between constrained optima of nutritional ecology response surfaces.
#' \emph{Environmental and Ecological Statistics}.
#'
#' @keywords internal
"_PACKAGE"
