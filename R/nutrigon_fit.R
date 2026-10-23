# ============================================================================
# nutrigon_fit.R
#
# fit$model is a nutrigon_*_fit wrapper (not a raw lm/glm/merMod object);
# the raw fitted model is fit$model$fit. Use nutrigon_coefs(fit) for the
# link-scale quadratic coefficients.
# ============================================================================

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
#' Three model families are supported via `model_type`:
#' \describe{
#'   \item{`"lm"`}{(default) Ordinary least squares. Appropriate for
#'     approximately Gaussian, homoscedastic trait responses.}
#'   \item{`"glm"`}{Generalized linear model, any `family` (including
#'     quasi-families such as `quasipoisson()`, `quasibinomial()` for
#'     overdispersed counts/proportions). Pass `family = ` via `...`. The
#'     closed-form and Monte Carlo optimum-finding both operate on the
#'     **link scale** (see Details), which is the scale on which the
#'     surface is quadratic.}
#'   \item{`"lmer"`}{Linear mixed-effects model via \pkg{lme4}. Requires
#'     `re_formula` specifying the random-effects structure, e.g.
#'     `re_formula = "+ (1 | block)"`. The optimum located is the
#'     population-level optimum (random effects integrated out via
#'     `re.form = NA`), not a group-specific optimum.}
#' }
#'
#' @details
#' For `model_type = "glm"`, the fitted linear predictor
#' \eqn{\eta(x,y) = \beta_0 + \beta_x x + \beta_y y + \beta_{xx} x^2 +
#' \beta_{yy} y^2 + \beta_{xy} xy} is quadratic in \eqn{(x,y)}, but the
#' predicted trait value on the response scale generally is not (e.g. for a
#' logit or log link). Because the link function is monotonic, the point
#' that maximises \eqn{\eta} also maximises the response-scale prediction,
#' so locating the optimum on the link scale and reporting it in the
#' original nutrient units loses nothing, and lets the same closed-form
#' stationary-point algebra be reused unchanged across `"lm"`, `"glm"`, and
#' `"lmer"`.
#'
#' @param data A data frame containing nutrient intake and trait columns.
#' @param x_col Name of the protein (x-axis nutrient) column.
#' @param y_col Name of the carbohydrate (y-axis nutrient) column.
#' @param z_col Name of the trait (response) column.
#' @param design_id_col Name of a column identifying distinct design points
#'   (e.g. treatment, or rail x concentration combination). Stored with the
#'   fit for reference; design-point identifiers are used by
#'   [nutrigon_compare()] only with `resample = "stratified"` and by
#'   [nutrigon_build_grid()] for design-matched grids.
#' @param model_type One of `"lm"` (default), `"glm"`, or `"lmer"`. See
#'   Details.
#' @param re_formula For `model_type = "lmer"` only: the random-effects
#'   portion of the model formula, as a string starting with `"+"`, e.g.
#'   `"+ (1 | block)"`. Required when `model_type = "lmer"`; ignored
#'   otherwise. Any valid \pkg{lme4} random-effects syntax is supported --
#'   multiple/crossed terms (`"+ (1 | block) + (1 | batch)"`), nesting
#'   (`"+ (1 | block/subblock)"`), interaction grouping factors
#'   (`"+ (1 | block:batch)"`), and random slopes -- since this string is
#'   passed to \pkg{lme4} essentially unmodified. **For random slopes**,
#'   refer to the nutrient columns as `.x` and `.y` (not your original
#'   `x_col`/`y_col` names), e.g. `"+ (1 + .x | block)"`: internally the
#'   fitting data frame renames your columns to `.x`/`.y`/`.z` before
#'   fitting, so `"+ (1 + x | block)"` will fail with an `object 'x' not
#'   found` error if your column was actually named `x_col = "protein"`
#'   or similar. Random *intercepts* (`"+ (1 | block)"`) are unaffected,
#'   since they don't reference the nutrient columns by name.
#' @param re_group_cols For `model_type = "lmer"` only: character vector of
#'   column name(s) in `data` that `re_formula` refers to as grouping
#'   factors (e.g. `"block"` for `re_formula = "+ (1 | block)"`; both
#'   `c("block", "batch")` for crossed terms or `c("block", "subblock")`
#'   for nested terms). Required whenever `re_formula` is supplied.
#'   Needed because internally the data is reduced to just the x/y/z
#'   columns before fitting; any column `re_formula` needs must be named
#'   here so it is carried through.
#' @param peak `TRUE` (default) to locate the maximum; `FALSE` for the
#'   minimum.
#' @param n_cloud Number of candidate points for the Monte Carlo fallback.
#' @param n_reps Number of MC draws (the final region is the union of all
#'   draws).
#' @param region_quantile Top (or bottom) quantile defining the peak region.
#'   Default: 0.95 for peaks, 0.05 for valleys.
#' @param ... Forwarded to the underlying fitter: e.g. `family = binomial()`
#'   or `family = quasipoisson()` for `model_type = "glm"`; `weights = `,
#'   `control = `, etc. for `model_type = "glm"` or `"lmer"`. Ignored for
#'   `model_type = "lm"`.
#'
#' @return An object of class `"nutrigon_fit"`, a list with components:
#'   \describe{
#'     \item{optimum}{Data frame with `x`, `y`, `theta` (angle from origin,
#'       degrees), and `used_fallback` (logical).}
#'     \item{model}{A `nutrigon_*_fit` wrapper object. The raw fitted
#'       `lm`/`glm`/`merMod` object is available as `model$fit` (e.g. for
#'       `summary()`, residual diagnostics, or `AIC()`).}
#'     \item{model_type}{Character: `"lm"`, `"glm"`, or `"lmer"`.}
#'     \item{data}{The cleaned input data.}
#'     \item{hull}{Design hull vertices (for plotting).}
#'     \item{region}{The top/bottom quantile region (for plotting), on the
#'       link scale for `model_type = "glm"`.}
#'     \item{elongation}{Curvature elongation ratio (eigenvalue ratio of the
#'       Hessian), computed on the link scale.}
#'   }
#'
#' @examples
#' data(lee_drosophila)
#'
#' # Default: ordinary least squares
#' fit_lifespan <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                               "lifespan", design_id_col = "treatment")
#' print(fit_lifespan)
#'
#' \donttest{
#' # GLM example: a count-like trait fitted with a quasi-Poisson family
#' fit_eggs_glm <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                               "dailyeggs", design_id_col = "treatment",
#'                               model_type = "glm", family = quasipoisson())
#'
#' # Mixed-model example (requires lme4 and lmerTest; install with
#' # install.packages(c("lme4", "lmerTest")) if not already available).
#' # re_group_cols must name the original column(s) re_formula refers to,
#' # since the fitter only sees the x/y/z columns by default.
#' fit_lifespan_lmer <- nutrigon_fit(
#'   lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
#'   model_type = "lmer", re_formula = "+ (1 | treatment)",
#'   re_group_cols = "treatment"
#' )
#' }
#'
#' @seealso [nutrigon_compare()], [nutrigon_elongation()], [nutrigon_coefs()]
#' @export
nutrigon_fit <- function(data, x_col, y_col, z_col,
                          design_id_col   = NULL,
                          model_type      = c("lm", "glm", "lmer"),
                          re_formula      = NULL,
                          re_group_cols   = NULL,
                          peak            = TRUE,
                          n_cloud         = 500,
                          n_reps          = 5,
                          region_quantile = NULL,
                          ...) {
  model_type <- match.arg(model_type)

  fit <- .fit_surface(data, x_col, y_col, z_col,
                       model_type = model_type, re_formula = re_formula,
                       re_group_cols = re_group_cols, ...)
  opt <- .locate_optimum(fit, peak, n_cloud, n_reps, region_quantile)

  if (is.null(region_quantile))
    region_quantile <- if (peak) 0.95 else 0.05

  # Build region for plotting (on the link scale for glm; identity for lm)
  cloud <- .sample_in_hull(fit$hull$x, fit$hull$y, n = n_cloud * n_reps)
  names(cloud) <- c(".x", ".y")
  cloud$pred   <- .predict_link(fit$model_obj, newdata = cloud)
  cutoff <- stats::quantile(cloud$pred, region_quantile)
  region <- if (peak) cloud[cloud$pred >= cutoff, ] else
                       cloud[cloud$pred <= cutoff, ]
  region_out <- data.frame(x = region$.x, y = region$.y)
  names(region_out) <- c(x_col, y_col)

  # Elongation from Hessian eigenvalues -- now via the dispatch layer
  # instead of stats::coef(fit$model) directly, so it works for lm/glm/lmer
  # alike. Computed on the link scale, consistent with .locate_optimum().
  cf <- .get_link_coefs(fit$model_obj)
  .check_quadratic_coefs(cf)
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
    list(optimum       = optimum,
         model         = fit$model_obj,
         model_type    = model_type,
         data          = fit$data,
         hull          = fit$hull,
         region        = region_out,
         elongation    = elong,
         design_id_col = design_id_col,
         x_col         = x_col,
         y_col         = y_col,
         z_col         = z_col,
         peak          = peak),
    class = "nutrigon_fit"
  )
}

#' Extract fitted quadratic coefficients from a `nutrigon_fit` object
#'
#' Convenience accessor for the fitted link-scale quadratic coefficients
#' (`.x`, `.y`, `I(.x^2)`, `I(.y^2)`, `.x:.y`), regardless of whether the
#' underlying model is `"lm"`, `"glm"`, or `"lmer"`. Use this instead of
#' calling `coef()` on `fit$model` directly, since for `"lmer"` fits
#' `fit$model` is a wrapper object, not the raw `merMod`, and for `"glm"`
#' fits the coefficients are on the link scale by design (see
#' [nutrigon_fit()]).
#'
#' @param fit An object of class `"nutrigon_fit"`.
#'
#' @return A named numeric vector.
#'
#' @examples
#' data(lee_drosophila)
#' fit_lifespan <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                               "lifespan", design_id_col = "treatment")
#' nutrigon_coefs(fit_lifespan)
#'
#' @seealso [nutrigon_fit()]
#' @export
nutrigon_coefs <- function(fit) {
  if (!inherits(fit, "nutrigon_fit"))
    stop("fit must be a \"nutrigon_fit\" object (from nutrigon_fit()).")
  .get_link_coefs(fit$model)
}

#' @export
print.nutrigon_fit <- function(x, ...) {
  opt <- x$optimum
  cat(sprintf("Nutrigonometry surface fit (%s)\n", x$model_type))
  cat(sprintf("  Trait:   %s\n", x$z_col))
  cat(sprintf("  Optimum: %s = %.2f, %s = %.2f\n",
              x$x_col, opt[[x$x_col]], x$y_col, opt[[x$y_col]]))
  cat(sprintf("  Angle from origin (theta): %.2f degrees\n", opt$theta))
  cat(sprintf("  Curvature elongation:      %.2f\n", x$elongation))
  cat(sprintf("  Used MC fallback:          %s\n", opt$used_fallback))
  if (x$model_type == "glm")
    cat("  Note: optimum located on the link scale (see ?nutrigon_fit)\n")
  invisible(x)
}
