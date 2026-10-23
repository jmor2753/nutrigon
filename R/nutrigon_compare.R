#' Compare two nutritional optima and bootstrap the angular separation
#'
#' Fits two quadratic response surfaces, locates each nutritional optimum
#' with the hybrid estimator (closed-form stationary point when admissible,
#' otherwise a constrained Monte Carlo search with the `max` rule), and
#' computes the angular separation theta between the two optima together
#' with a percentile bootstrap confidence interval. Each bootstrap replicate
#' refits both surfaces, so uncertainty in the surface coefficients is
#' propagated into the location of each optimum and into theta.
#'
#' @section Resampling:
#' `resample = "case"` (default) draws observations (rows) with replacement
#' from the whole dataset. This is the scheme whose coverage was validated
#' by simulation. With `paired = TRUE` the same rows are used for both
#' traits, which preserves the covariance between the two fitted surfaces
#' when both traits are measured on the same individuals; with
#' `paired = FALSE` each trait is resampled independently.
#' `resample = "stratified"` resamples within each design point (requires
#' `design_id_col1`/`design_id_col2`), keeping the number of observations
#' per design point fixed. It is useful for balanced designs with several
#' replicates per design point, but its coverage has not been separately
#' validated.
#'
#' @section What theta measures:
#' theta is the difference between the polar angles of the two optima,
#' i.e. the difference in nutrient balance (e.g. protein:carbohydrate
#' direction). It does not measure the radial distance of the optima from
#' the origin (total intake) or the fitness cost of feeding between them.
#' Report it with the locations of the optima (`point1`, `point2`).
#'
#' @param data1,data2 Data frames for traits 1 and 2. Can be the same object
#'   if both traits were measured on the same individuals.
#' @param x1_col,y1_col,z1_col Column names for surface 1.
#' @param x2_col,y2_col,z2_col Column names for surface 2.
#' @param design_id_col1,design_id_col2 Design-point ID columns. Required
#'   for `resample = "stratified"`; otherwise optional.
#' @param paired `TRUE` (default) if both traits were measured on the same
#'   individuals (rows of `data1` and `data2` must correspond); `FALSE` for
#'   unpaired data.
#' @param resample `"case"` (default, validated) or `"stratified"`. See
#'   the Resampling section.
#' @param peak1,peak2 `TRUE` to locate maxima; `FALSE` for minima.
#' @param model_type1,model_type2 One of `"lm"` (default), `"glm"`, or
#'   `"lmer"`, for surfaces 1 and 2 respectively. See [nutrigon_fit()].
#' @param re_formula1,re_formula2 Random-effects formula strings, required
#'   only when the corresponding `model_type` is `"lmer"`.
#' @param re_group_cols1,re_group_cols2 Grouping column name(s) referred to
#'   by the corresponding `re_formula`.
#' @param family1,family2 Forwarded as `family = ` to [stats::glm()] when
#'   the corresponding `model_type` is `"glm"`.
#' @param extra_args1,extra_args2 Named lists of further arguments for the
#'   underlying fitters (e.g. `list(weights = ...)`).
#' @param B Number of bootstrap resamples (default 500).
#' @param conf_level Confidence level of the percentile interval (default
#'   0.95).
#' @param n_cloud Monte Carlo fallback cloud size.
#' @param n_reps_cloud Number of Monte Carlo draws per estimate.
#' @param seed Optional random seed (set locally; the global RNG state is
#'   restored on exit).
#'
#' @return An object of class `"nutrigon_comparison"`, a list with:
#'   \describe{
#'     \item{theta_hat}{Point estimate of angular separation (degrees),
#'       `|theta1 - theta2|`.}
#'     \item{ci}{Percentile bootstrap CI for theta (unsigned), named e.g.
#'       `c(lower_95, upper_95)`.}
#'     \item{theta1, theta2}{Polar angles of the two optima (degrees).}
#'     \item{point1, point2}{Locations of the two optima, `c(x, y)`, in the
#'       original nutrient units.}
#'     \item{used_fallback1, used_fallback2}{Whether the Monte Carlo
#'       fallback located each full-data optimum.}
#'     \item{angular_support}{List with `surface1` and `surface2`, each
#'       `c(theta_min, theta_max)`: the angular range of that surface's
#'       design hull (degrees). Used by [nutrigon_check()].}
#'     \item{boot_diffs}{The `B` bootstrap values of theta (`NA` for failed
#'       replicates).}
#'     \item{boot_signed}{The `B` bootstrap values of `theta1 - theta2`.}
#'     \item{fallback_rate}{Proportion of successful bootstrap replicates in
#'       which at least one optimum required the Monte Carlo fallback.}
#'     \item{n_failed}{Number of bootstrap replicates that failed (e.g. a
#'       rank-deficient refit) and were excluded.}
#'     \item{n_fallback}{Deprecated alias of `n_failed`, kept for
#'       compatibility with version 0.1.0.}
#'     \item{B, conf_level, paired, resample}{Settings used.}
#'   }
#'
#' @examples
#' data(lee_drosophila)
#' result <- nutrigon_compare(
#'   lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
#'   lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
#'   design_id_col1 = "treatment", design_id_col2 = "treatment",
#'   paired = TRUE, B = 50, seed = 1
#' )
#' print(result)
#'
#' \donttest{
#' # Different model types for the two surfaces
#' result_glm <- nutrigon_compare(
#'   lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
#'   lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
#'   paired = TRUE, B = 50, seed = 1,
#'   model_type2 = "glm", family2 = quasipoisson()
#' )
#' }
#'
#' @seealso [nutrigon_fit()], [nutrigon_check()], [nutrigon_equivalence()]
#' @export
nutrigon_compare <- function(data1, x1_col, y1_col, z1_col,
                             data2, x2_col, y2_col, z2_col,
                             design_id_col1 = NULL,
                             design_id_col2 = NULL,
                             paired         = TRUE,
                             resample       = c("case", "stratified"),
                             peak1          = TRUE,
                             peak2          = TRUE,
                             model_type1    = c("lm", "glm", "lmer"),
                             model_type2    = c("lm", "glm", "lmer"),
                             re_formula1    = NULL,
                             re_formula2    = NULL,
                             re_group_cols1 = NULL,
                             re_group_cols2 = NULL,
                             family1        = NULL,
                             family2        = NULL,
                             extra_args1    = list(),
                             extra_args2    = list(),
                             B              = 500,
                             conf_level     = 0.95,
                             n_cloud        = 500,
                             n_reps_cloud   = 5,
                             seed           = NULL) {

  model_type1 <- match.arg(model_type1)
  model_type2 <- match.arg(model_type2)
  resample    <- match.arg(resample)

  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be strictly between 0 and 1.")
  if (paired && nrow(data1) != nrow(data2))
    stop("paired = TRUE requires data1 and data2 to have the same rows ",
         "(the same individuals). Got ", nrow(data1), " and ", nrow(data2),
         " rows; use paired = FALSE for traits measured on different units.")
  if (resample == "stratified") {
    if (is.null(design_id_col1) || is.null(design_id_col2))
      stop("resample = \"stratified\" requires design_id_col1 and design_id_col2.")
    for (cc in list(list(data1, design_id_col1), list(data2, design_id_col2)))
      if (!cc[[2]] %in% names(cc[[1]]))
        stop("Design-point column \"", cc[[2]], "\" not found in data.")
  }

  if (!is.null(family1)) extra_args1$family <- family1
  if (!is.null(family2)) extra_args2$family <- family2

  fit_one <- function(data, x_col, y_col, z_col, model_type, re_formula,
                      re_group_cols, extra) {
    do.call(.fit_surface,
            c(list(data = data, x_col = x_col, y_col = y_col, z_col = z_col,
                   model_type = model_type, re_formula = re_formula,
                   re_group_cols = re_group_cols),
              extra))
  }

  if (!is.null(seed)) withr::local_seed(seed)

  # ---- point estimates ---------------------------------------------------
  fit1 <- fit_one(data1, x1_col, y1_col, z1_col, model_type1, re_formula1,
                  re_group_cols1, extra_args1)
  fit2 <- fit_one(data2, x2_col, y2_col, z2_col, model_type2, re_formula2,
                  re_group_cols2, extra_args2)
  opt1 <- .locate_optimum(fit1, peak1, n_cloud, n_reps_cloud)
  opt2 <- .locate_optimum(fit2, peak2, n_cloud, n_reps_cloud)

  theta1    <- .theta_deg(opt1$x, opt1$y)
  theta2    <- .theta_deg(opt2$x, opt2$y)
  theta_hat <- abs(theta1 - theta2)

  support1 <- .angular_support_xy(fit1$data$.x, fit1$data$.y)
  support2 <- .angular_support_xy(fit2$data$.x, fit2$data$.y)

  # ---- resampling --------------------------------------------------------
  strat_index <- function(groups) {
    idx <- split(seq_along(groups), groups)
    unlist(lapply(idx, function(ii) ii[sample.int(length(ii), length(ii),
                                                   replace = TRUE)]),
           use.names = FALSE)
  }
  draw_rows <- function(data, id_col) {
    if (resample == "stratified") strat_index(data[[id_col]]) else
      sample.int(nrow(data), nrow(data), replace = TRUE)
  }

  boot_signed   <- rep(NA_real_, B)
  boot_fallback <- rep(NA, B)

  for (b in seq_len(B)) {
    res_b <- tryCatch({
      if (paired) {
        idx <- draw_rows(data1, design_id_col1)
        b1 <- data1[idx, , drop = FALSE]; b2 <- data2[idx, , drop = FALSE]
      } else {
        b1 <- data1[draw_rows(data1, design_id_col1), , drop = FALSE]
        b2 <- data2[draw_rows(data2, design_id_col2), , drop = FALSE]
      }
      f1 <- fit_one(b1, x1_col, y1_col, z1_col, model_type1, re_formula1,
                    re_group_cols1, extra_args1)
      f2 <- fit_one(b2, x2_col, y2_col, z2_col, model_type2, re_formula2,
                    re_group_cols2, extra_args2)
      o1 <- .locate_optimum(f1, peak1, n_cloud, n_reps_cloud)
      o2 <- .locate_optimum(f2, peak2, n_cloud, n_reps_cloud)
      c(.theta_deg(o1$x, o1$y) - .theta_deg(o2$x, o2$y),
        o1$used_fallback || o2$used_fallback)
    }, error = function(e) c(NA_real_, NA))
    boot_signed[b]   <- res_b[1]
    boot_fallback[b] <- as.logical(res_b[2])
  }

  boot_diffs <- abs(boot_signed)
  valid      <- boot_diffs[!is.na(boot_diffs)]
  n_failed   <- sum(is.na(boot_diffs))
  alpha      <- 1 - conf_level
  ci <- if (length(valid) >= 10L)
    stats::quantile(valid, c(alpha / 2, 1 - alpha / 2), names = FALSE) else
      c(NA_real_, NA_real_)
  pct <- format(100 * conf_level, trim = TRUE)
  names(ci) <- paste0(c("lower_", "upper_"), pct)
  if (n_failed > 0.2 * B)
    warning(n_failed, " of ", B, " bootstrap replicates failed; the ",
            "interval may be unreliable.", call. = FALSE)

  structure(
    list(theta_hat       = theta_hat,
         ci              = ci,
         theta1          = theta1,
         theta2          = theta2,
         point1          = c(x = opt1$x, y = opt1$y),
         point2          = c(x = opt2$x, y = opt2$y),
         used_fallback1  = opt1$used_fallback,
         used_fallback2  = opt2$used_fallback,
         angular_support = list(surface1 = support1, surface2 = support2),
         boot_diffs      = boot_diffs,
         boot_signed     = boot_signed,
         fallback_rate   = mean(boot_fallback, na.rm = TRUE),
         n_failed        = n_failed,
         n_fallback      = n_failed,
         B               = B,
         conf_level      = conf_level,
         paired          = paired,
         resample        = resample,
         z1_col          = z1_col,
         z2_col          = z2_col,
         model_type1     = model_type1,
         model_type2     = model_type2),
    class = "nutrigon_comparison"
  )
}

#' @export
print.nutrigon_comparison <- function(x, ...) {
  pct <- if (is.null(x$conf_level)) "95" else format(100 * x$conf_level, trim = TRUE)
  cat("Angular separation between nutritional optima\n")
  cat(sprintf("  Trait 1 (%s, %s): theta = %.2f degrees\n",
              x$z1_col, x$model_type1, x$theta1))
  cat(sprintf("  Trait 2 (%s, %s): theta = %.2f degrees\n",
              x$z2_col, x$model_type2, x$theta2))
  cat(sprintf("  Angular separation: %.2f degrees\n", x$theta_hat))
  cat(sprintf("  %s%% CI: [%.2f, %.2f] degrees\n", pct, x$ci[1], x$ci[2]))
  if (!is.null(x$resample))
    cat(sprintf("  Bootstrap: B = %d (%s, %s resampling), %d failed\n",
                x$B, if (isTRUE(x$paired)) "paired" else "unpaired",
                x$resample, x$n_failed))
  if (!is.null(x$fallback_rate) && is.finite(x$fallback_rate))
    cat(sprintf("  Monte Carlo fallback used in %.0f%% of replicates\n",
                100 * x$fallback_rate))
  if (is.finite(x$theta_hat) && x$theta_hat < 5)
    cat("  Note: angle < 5 degrees; the unsigned interval cannot show that\n",
        "  the optima are aligned. See ?nutrigon_equivalence.\n", sep = "")
  invisible(x)
}
