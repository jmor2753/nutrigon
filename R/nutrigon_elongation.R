#' Curvature elongation ratio of a fitted quadratic surface
#'
#' Returns the ratio of the larger to smaller absolute eigenvalue of the
#' fitted quadratic's Hessian. A value near 1 indicates roughly symmetric
#' curvature; large values (e.g. > 10) indicate a narrow, ridge-like peak.
#' This is used as one of the risk factors in [nutrigon_check()].
#'
#' @param fit A `"nutrigon_fit"` object (from [nutrigon_fit()]), or an `lm`
#'   object from a quadratic fit with terms `.x`, `.y`, `I(.x^2)`,
#'   `I(.y^2)`, `.x:.y`.
#'
#' @return A single positive number (the elongation ratio, >= 1), or `NA` if
#'   the Hessian is singular.
#'
#' @examples
#' data(lee_drosophila)
#' fit_lifespan <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                               "lifespan", design_id_col = "treatment")
#' nutrigon_elongation(fit_lifespan)
#'
#' @seealso [nutrigon_fit()], [nutrigon_check()]
#' @export
nutrigon_elongation <- function(fit) {
  mod <- if (inherits(fit, "nutrigon_fit")) fit$model else fit
  cf  <- stats::coef(mod)
  H   <- matrix(c(2 * cf["I(.x^2)"],     cf[".x:.y"],
                      cf[".x:.y"],    2 * cf["I(.y^2)"]),
                nrow = 2)
  eig <- abs(eigen(H, only.values = TRUE)$values)
  if (min(eig) == 0) return(NA_real_)
  max(eig) / min(eig)
}
