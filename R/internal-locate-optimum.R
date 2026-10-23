# ============================================================================
# internal-locate-optimum.R
#
# Hybrid optimum estimator: closed-form stationary point when admissible
# (inside the design hull, correct curvature sign); otherwise a Monte Carlo
# search over the hull, summarised by the point with the highest (or
# lowest) fitted value ("max" rule). Coefficients and predictions are read
# through the .get_link_coefs() / .predict_link() dispatch layer
# (internal-fit.R), so the same algorithm serves lm, glm and lmer fits.
# ============================================================================

.locate_optimum <- function(fit, peak = TRUE, n_cloud = 500, n_reps = 5,
                             region_quantile = NULL) {
  if (is.null(region_quantile))
    region_quantile <- if (peak) 0.95 else 0.05

  cf <- .get_link_coefs(fit$model_obj)
  .check_quadratic_coefs(cf)

  a  <- unname(cf["I(.x^2)"])
  b  <- unname(cf["I(.y^2)"])
  e  <- unname(cf[".x:.y"])
  cc <- unname(cf[".x"])
  d  <- unname(cf[".y"])
  det <- 4 * a * b - e^2

  correct_sign <- if (peak) a < 0 else a > 0
  valid_det    <- det > 0

  if (valid_det && correct_sign) {
    xs <- (e * d - 2 * b * cc) / det
    ys <- (e * cc - 2 * a * d) / det
    if (.point_in_hull(xs, ys, fit$hull$x, fit$hull$y)) {
      return(list(x = xs, y = ys, used_fallback = FALSE))
    }
  }

  # Monte Carlo fallback: sample cloud inside hull, evaluate on the LINK
  # scale (see internal-fit.R design note), take top/bottom region.
  pool <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    cloud <- .sample_in_hull(fit$hull$x, fit$hull$y, n = n_cloud)
    names(cloud) <- c(".x", ".y")
    cloud$pred <- .predict_link(fit$model_obj, newdata = cloud)
    cutoff <- stats::quantile(cloud$pred, region_quantile)
    region <- if (peak) cloud[cloud$pred >= cutoff, ] else
                         cloud[cloud$pred <= cutoff, ]
    pool[[i]] <- region
  }
  region <- do.call(rbind, pool)
  if (nrow(region) == 0L)
    stop("MC fallback returned empty region -- too few cloud points?")

  opt <- region[which.max(if (peak) region$pred else -region$pred), ]
  list(x = unname(opt$.x), y = unname(opt$.y), used_fallback = TRUE)
}
