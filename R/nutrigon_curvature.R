#' Local Gaussian and mean curvature, and surface-area ratio, of a fitted
#' surface
#'
#' Computes the local Gaussian curvature and mean curvature of a fitted
#' quadratic surface at every point of a grid restricted to the convex hull
#' of the design, together with a surface-area ratio (the fitted surface's
#' area relative to a flat surface spanning the same design). This treats
#' the fitted response surface as a 2-manifold embedded in 3-space (the
#' first and second fundamental forms of differential geometry), which is a
#' different notion of "curvature" from the Hessian-eigenvalue elongation
#' ratio returned by [nutrigon_elongation()] -- elongation describes the
#' *shape* of the quadratic (how stretched the peak/valley is along its
#' principal axes), while Gaussian/mean curvature describe the *local
#' bending* of the surface itself at each point.
#'
#' @details For `model_type = "glm"` fits (see [nutrigon_fit()]), curvature
#'   is computed on the **link scale**, consistent with how the optimum
#'   itself is located and how [nutrigon_elongation()] computes elongation.
#'   The fitted quadratic's coefficients are read via the same internal
#'   dispatch layer used throughout the package, so curvature is available
#'   identically for `"lm"`, `"glm"`, and `"lmer"` fits.
#'
#'   Gaussian curvature (`Gauss_curv`) is the product of the two principal
#'   curvatures at a point; mean curvature (`Mean_curv`) is their average.
#'   Both are computed in closed form from the first and second fundamental
#'   forms of the fitted quadratic surface, following standard differential
#'   geometry (see e.g. do Carmo, *Differential Geometry of Curves and
#'   Surfaces*).
#'
#'   The surface-area ratio (`surface_area_ratio`, a single number returned
#'   alongside the grid, not per-point) is the fitted surface's area over
#'   the design hull, divided by the area of a flat surface spanning the
#'   same hull. A ratio of 1 indicates a flat response; larger values
#'   indicate a more strongly curved (bowl-, dome-, or saddle-shaped)
#'   surface. This is computed by numerically integrating the surface-area
#'   element (via \pkg{pracma}) over the bounding box of the design, for
#'   both the fitted surface and a flat reference surface fitted to a
#'   constant response over the same design.
#'
#' @param fit A `"nutrigon_fit"` object (from [nutrigon_fit()]).
#' @param n Grid resolution per axis (default 50, matching the resolution
#'   used in the original curvature analysis this function is based on).
#'   Only points inside the convex hull of the design are returned.
#' @param quiet Suppress the link-scale reminder message for GLM fits
#'   (default `FALSE`).
#'
#' @return A list with components:
#'   \describe{
#'     \item{grid}{A data frame with columns named after the original
#'       `x_col` and `y_col`, plus `Gauss_curv` and `Mean_curv`, one row per
#'       grid point inside the design hull.}
#'     \item{surface_area_ratio}{A single number: the fitted surface's area
#'       divided by a flat reference surface's area over the same design
#'       bounding box (see Details). Requires \pkg{pracma}.}
#'   }
#'
#' @examples
#' data(lee_drosophila)
#' fit_lifespan <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                               "lifespan", design_id_col = "treatment")
#' curv_lifespan <- nutrigon_curvature(fit_lifespan, n = 30)
#' head(curv_lifespan$grid)
#' curv_lifespan$surface_area_ratio
#'
#' @seealso [nutrigon_fit()], [nutrigon_elongation()], [nutrigon_grid()]
#' @export
nutrigon_curvature <- function(fit, n = 50, quiet = FALSE) {
  if (!inherits(fit, "nutrigon_fit"))
    stop("fit must be a \"nutrigon_fit\" object (from nutrigon_fit()).")

  cf <- .get_link_coefs(fit$model)
  .check_quadratic_coefs(cf)

  x_col <- fit$x_col; y_col <- fit$y_col
  d  <- fit$data
  xs <- seq(min(d$.x), max(d$.x), length.out = n)
  ys <- seq(min(d$.y), max(d$.y), length.out = n)
  gd <- expand.grid(.x = xs, .y = ys)
  inside <- .point_in_hull(gd$.x, gd$.y, fit$hull$x, fit$hull$y)
  gd <- gd[inside, ]

  curv <- .curvature_at(gd$.x, gd$.y, cf)
  gd$Gauss_curv <- curv$Gauss_curv
  gd$Mean_curv  <- curv$Mean_curv

  names(gd)[1:2] <- c(x_col, y_col)

  if (fit$model_type == "glm" && !quiet) {
    message("nutrigon_curvature(): curvature was computed on the LINK ",
            "scale (family: ", .glm_family_name(fit$model), "), not the ",
            "original response scale. Pass quiet = TRUE to suppress this ",
            "message.")
  }

  surface_area_ratio <- .surface_area_ratio(d$.x, d$.y, cf)

  list(grid = gd, surface_area_ratio = surface_area_ratio)
}

# ============================================================================
# Internal: curvature and surface-area machinery.
#
# Gaussian/mean curvature and a surface-area ratio for a quadratic
# z = f(x, y) treated as a 2-manifold in 3-space (first/second fundamental
# forms). Coefficients come from .get_link_coefs(), so the same code works
# for lm/glm/lmer fits, on the scale each model estimates a quadratic on.
#
# Coefficient naming: a/b/c/d/e for the quadratic
#   z = a.x^2 + b.y^2 + c.x + d.y + e.xy  (+ intercept, irrelevant to
# curvature). These map onto the package's canonical coefficient names as:
#   a <- cf["I(.x^2)"]   b <- cf["I(.y^2)"]   e <- cf[".x:.y"]
#   c <- cf[".x"]        d <- cf[".y"]
# ============================================================================

# 3D cross product (vectors padded to length 3 if shorter).
.cross3 <- function(x, y) {
  create3D <- function(v) utils::head(c(v, rep(0, 3)), 3)
  x <- create3D(x); y <- create3D(y)
  j <- function(i) (i - 1) %% 3 + 1
  vapply(1:3, function(i) x[j(i + 1)] * y[j(i + 2)] - x[j(i + 2)] * y[j(i + 1)],
         numeric(1))
}

.norm_vec <- function(x) sqrt(sum(x^2))

# Gaussian and mean curvature at a set of (x, y) points, for a quadratic
# surface with link-scale coefficients cf (named as in .get_link_coefs()).
# Vectorised over x, y.
.curvature_at <- function(x, y, cf) {
  a <- unname(cf["I(.x^2)"]); b <- unname(cf["I(.y^2)"])
  e <- unname(cf[".x:.y"]);   cc <- unname(cf[".x"]); dd <- unname(cf[".y"])

  n <- length(x)
  gauss <- numeric(n); mean_c <- numeric(n)

  for (i in seq_len(n)) {
    partial_x_z <- 2 * a * x[i] + cc + e * y[i]
    partial_y_z <- 2 * b * y[i] + dd + e * x[i]

    partial_x <- c(1, 0, partial_x_z)
    partial_y <- c(0, 1, partial_y_z)

    Efirst <- sum(partial_x * partial_x)
    Ffirst <- sum(partial_x * partial_y)
    Gfirst <- sum(partial_y * partial_y)

    cr <- .cross3(partial_x, partial_y)
    normal_unit <- -cr / .norm_vec(cr)

    Lsecond <- sum(c(0, 0, -2 * a) * normal_unit)
    Msecond <- sum(c(0, 0, -e)     * normal_unit)
    Nsecond <- sum(c(0, 0, -2 * b) * normal_unit)

    first_fundForm  <- matrix(c(Efirst, Ffirst, Ffirst, Gfirst), 2, 2)
    second_fundForm <- matrix(c(Lsecond, Msecond, Msecond, Nsecond), 2, 2)

    detI  <- Efirst * Gfirst - Ffirst^2
    gauss[i]  <- det(second_fundForm) / detI
    mean_c[i] <- sum(diag(second_fundForm %*% solve(first_fundForm))) / 2
  }

  list(Gauss_curv = gauss, Mean_curv = mean_c)
}

# Surface-area element for numerical integration (pracma::quad2d).
.surface_area_integrand <- function(cf, x, y) {
  a <- unname(cf["I(.x^2)"]); b <- unname(cf["I(.y^2)"])
  e <- unname(cf[".x:.y"]);   cc <- unname(cf[".x"]); dd <- unname(cf[".y"])
  sqrt(1 + (2 * a * x + cc + e * y)^2 + (2 * b * y + dd + e * x)^2)
}

# Ratio of fitted-surface area to a flat reference surface's area, over the
# bounding box of the design. Requires pracma (Suggests).
.surface_area_ratio <- function(x, y, cf) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    warning("surface_area_ratio requires the pracma package; returning NA. ",
            "Install it with install.packages(\"pracma\").", call. = FALSE)
    return(NA_real_)
  }
  flat_cf <- cf
  flat_cf[c("I(.x^2)", "I(.y^2)", ".x:.y", ".x", ".y")] <- 0

  area_fit  <- pracma::quad2d(function(x, y) .surface_area_integrand(cf, x, y),
                               min(x), max(x), min(y), max(y))
  area_flat <- pracma::quad2d(function(x, y) .surface_area_integrand(flat_cf, x, y),
                               min(x), max(x), min(y), max(y))
  area_fit / area_flat
}
