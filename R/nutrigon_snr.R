#' Signal-to-noise ratio of a fitted surface
#'
#' The standard deviation of the fitted values over the observed data
#' divided by the residual standard deviation, equivalent to
#' `sqrt(R^2 / (1 - R^2))` up to degrees-of-freedom corrections. Pass the
#' SNRs of your two surfaces to `nutrigon_build_grid(snr = ...)` so that the
#' validation grid reflects how flat your surfaces are relative to their
#' noise; this matters most for small angles.
#'
#' Defined for Gaussian fits (`model_type = "lm"`, and population-level
#' predictions for `"lmer"`). For `"glm"` fits there is no single residual
#' standard deviation on the link scale, so an error is returned; the
#' validation grids assume Gaussian surfaces.
#'
#' @param fit A `"nutrigon_fit"` object from [nutrigon_fit()].
#' @return A single positive number.
#'
#' @examples
#' data(lee_drosophila)
#' fit <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten", "lifespan")
#' nutrigon_snr(fit)
#'
#' @seealso [nutrigon_build_grid()]
#' @export
nutrigon_snr <- function(fit) {
  if (!inherits(fit, "nutrigon_fit"))
    stop("fit must be a \"nutrigon_fit\" object (from nutrigon_fit()).")
  m <- fit$model$fit
  if (fit$model_type == "glm")
    stop("nutrigon_snr() is defined for Gaussian fits (\"lm\"/\"lmer\") only.")
  fitted_vals <- if (fit$model_type == "lmer")
    stats::predict(m, re.form = NA) else stats::fitted(m)
  s <- stats::sigma(m)
  if (!is.finite(s) || s <= 0)
    stop("Residual SD is zero or undefined; SNR cannot be computed.")
  stats::sd(fitted_vals) / s
}
