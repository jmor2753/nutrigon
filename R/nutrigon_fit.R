#' Fit a quadratic response surface and locate the nutritional optimum
#'
#' Fits a full quadratic response surface (`z ~ x + y + x^2 + y^2 + x:y`) to
#' trait data collected under the Geometric Framework for Nutrition, and
#' locates the nutritional optimum (peak or valley) of the fitted surface.
#' The optimum is found in closed form when the fitted surface has the
#' correct curvature and the stationary point falls inside the convex hull
#' of the design; otherwise a Monte Carlo fallback samples points inside the
#' hull and takes the top (or bottom) quantile region.
#'
#' @param data A data frame containing nutrient intake and trait columns.
#' @param x_col Name of the protein (x-axis nutrient) column.
#' @param y_col Name of the carbohydrate (y-axis nutrient) column.
#' @param z_col Name of the trait (response) column.
#' @param design_id_col Name of a column identifying distinct design points
#'   (e.g. treatment, or rail x concentration combination). Used for
#'   within-design-point bootstrap resampling in [nutrigon_compare()]. If
#'   `NULL`, plain case-level resampling is used there instead.
#' @param peak `TRUE` (default) to locate the maximum; `FALSE` for the
#'   minimum.
#' @param n_cloud Number of candidate points for the Monte Carlo fallback.
#' @param n_reps Number of MC draws (the final region is the union of all
#'   draws).
#' @param region_quantile Top (or bottom) quantile defining the peak region.
#'   Default: 0.95 for peaks, 0.05 for valleys.
#'
#' @return An object of class `"nutrigon_fit"`, a list with components:
#'   \describe{
#'     \item{optimum}{Data frame with `x`, `y`, `theta` (angle from origin,
#'       degrees), and `used_fallback` (logical).}
#'     \item{model}{The `lm()` object for the fitted quadratic.}
#'     \item{data}{The cleaned input data.}
#'     \item{hull}{Design hull vertices (for plotting).}
#'     \item{region}{The top/bottom quantile region (for plotting).}
#'     \item{elongation}{Curvature elongation ratio (eigenvalue ratio of the
#'       Hessian).}
#'   }
#'
#' @examples
#' data(lee_drosophila)
#' fit_lifespan <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                               "lifespan", design_id_col = "treatment")
#' print(fit_lifespan)
#'
#' @seealso [nutrigon_compare()], [nutrigon_elongation()]
#' @export
nutrigon_fit <- function(data, x_col, y_col, z_col,
                          design_id_col = NULL,
                          peak          = TRUE,
                          n_cloud       = 500,
                          n_reps        = 5,
                          region_quantile = NULL) {

  fit <- .fit_quadratic(data, x_col, y_col, z_col)
  opt <- .locate_optimum(fit, peak, n_cloud, n_reps, region_quantile)

  if (is.null(region_quantile))
    region_quantile <- if (peak) 0.95 else 0.05

  # Build region for plotting
  cloud <- .sample_in_hull(fit$hull$x, fit$hull$y, n = n_cloud * n_reps)
  names(cloud) <- c(".x", ".y")
  cloud$pred   <- stats::predict(fit$model, newdata = cloud)
  cutoff <- stats::quantile(cloud$pred, region_quantile)
  region <- if (peak) cloud[cloud$pred >= cutoff, ] else
                       cloud[cloud$pred <= cutoff, ]
  region_out <- data.frame(x = region$.x, y = region$.y)
  names(region_out) <- c(x_col, y_col)

  # Elongation from Hessian eigenvalues
  cf  <- stats::coef(fit$model)
  H   <- matrix(c(2 * cf["I(.x^2)"],     cf[".x:.y"],
                      cf[".x:.y"],    2 * cf["I(.y^2)"]),
                nrow = 2)
  eig <- abs(eigen(H, only.values = TRUE)$values)
  elong <- if (min(eig) > 0) max(eig) / min(eig) else NA_real_

  optimum <- data.frame(x = opt$x, y = opt$y,
                         theta = .theta_deg(opt$x, opt$y),
                         used_fallback = opt$used_fallback)
  names(optimum)[1:2] <- c(x_col, y_col)

  structure(
    list(optimum     = optimum,
         model       = fit$model,
         data        = fit$data,
         hull        = fit$hull,
         region      = region_out,
         elongation  = elong,
         design_id_col = design_id_col,
         x_col       = x_col,
         y_col       = y_col,
         z_col       = z_col,
         peak        = peak),
    class = "nutrigon_fit"
  )
}

#' @export
print.nutrigon_fit <- function(x, ...) {
  opt <- x$optimum
  cat(sprintf("Nutrigonometry surface fit\n"))
  cat(sprintf("  Trait:   %s\n", x$z_col))
  cat(sprintf("  Optimum: %s = %.2f, %s = %.2f\n",
              x$x_col, opt[[x$x_col]], x$y_col, opt[[x$y_col]]))
  cat(sprintf("  Angle from origin (theta): %.2f degrees\n", opt$theta))
  cat(sprintf("  Curvature elongation:      %.2f\n", x$elongation))
  cat(sprintf("  Used MC fallback:          %s\n", opt$used_fallback))
  invisible(x)
}
