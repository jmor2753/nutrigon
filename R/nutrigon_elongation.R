# ============================================================================
# elongation-generalized.R
#
# Updated nutrigon_elongation() -- the real one, as pasted. Only ONE line
# changes from the original: `cf <- stats::coef(mod)` becomes
# `cf <- .get_link_coefs(mod)`, routed through the same dispatch layer as
# everywhere else. This is deliberately still accepting either a
# "nutrigon_fit" object OR a raw model object directly, matching the
# original function's contract -- but "raw model object" now means a
# nutrigon_*_fit wrapper (nutrigon_lm_fit/nutrigon_glm_fit/nutrigon_lmer_fit)
# rather than a bare lm object, since that's what fit$model now is.
#
# BREAKING CHANGE relative to the previous version: passing a bare lm()
# object built OUTSIDE nutrigon_fit() (as the old docs example/param text
# technically allowed: "or an lm object from a quadratic fit...") no longer
# works directly, because .get_link_coefs() dispatches on the wrapper
# class, not on "lm". If you need to elongation-check a bare lm object you
# built by hand, wrap it first: structure(list(fit = your_lm),
# class = "nutrigon_lm_fit"). This is documented below; flagging here in
# case it affects any existing user code or vignette examples.
# ============================================================================

#' Curvature elongation ratio of a fitted quadratic surface
#'
#' Returns the ratio of the larger to smaller absolute eigenvalue of the
#' fitted quadratic's Hessian (on the link scale for GLM fits -- see
#' Details). A value near 1 indicates roughly symmetric curvature; large
#' values (e.g. > 10) indicate a narrow, ridge-like peak. This is used as
#' one of the risk factors in [nutrigon_check()].
#'
#' @details For `model_type = "glm"` fits, the Hessian is computed on the
#'   **link scale**, consistent with how the optimum itself is located (see
#'   `?nutrigon_fit`). This means elongation values are directly comparable
#'   between `"lm"` and `"glm"` fits only insofar as both describe curvature
#'   on the scale each model actually estimates a quadratic on -- they are
#'   not automatically on a common footing with an `"lm"` fit of the same
#'   trait's raw values.
#'
#' @param fit A `"nutrigon_fit"` object (from [nutrigon_fit()]); or a
#'   `nutrigon_*_fit` wrapper object (i.e. `some_fit$model`) for `"lm"`,
#'   `"glm"`, or `"lmer"` fits.
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
  cf  <- .get_link_coefs(mod)
  .check_quadratic_coefs(cf)
  H   <- matrix(c(2 * cf["I(.x^2)"],     cf[".x:.y"],
                      cf[".x:.y"],    2 * cf["I(.y^2)"]),
                nrow = 2)
  eig <- abs(eigen(H, only.values = TRUE)$values)
  if (min(eig) == 0) return(NA_real_)
  max(eig) / min(eig)
}
