# ============================================================================
# compare-generalized.R
#
# Generalized nutrigon_compare(). Same bootstrap logic as before, but each
# surface can now be fit with its own model_type ("lm"/"glm"/"lmer"), since
# two traits being compared often need different distributional families
# (e.g. Example: viability as a binomial GLM, wing size as a Gaussian lm --
# this was already handled by hand in the applied scripts; this generalises
# it into the function itself).
#
# CHANGE FROM PREVIOUS VERSION:
#   - New arguments: model_type1, model_type2 (default "lm" each, so
#     existing calls with no model_type specified behave identically to
#     before), re_formula1, re_formula2, and ...1 / ...2 for family/weights/
#     etc. forwarded separately to each surface's fitter.
#   - Internally: .fit_quadratic() -> .fit_surface() (see internal-fit.R),
#     called with each surface's own model_type/re_formula/extra args.
#   - theta_hat, ci, boot_diffs, etc: UNCHANGED in meaning and computation.
#     No change to the bootstrap resampling scheme itself.
# ============================================================================

#' Compare two nutritional optima and bootstrap the angular separation
#'
#' Fits two quadratic response surfaces and computes the angular separation
#' theta between their nutritional optima, together with a bootstrap
#' percentile confidence interval. Supports paired resampling (both traits
#' measured on the same individuals) and unpaired resampling with optional
#' resampling by design point. The two surfaces may be fitted with different
#' model families (e.g. surface 1 as an ordinary least squares fit, surface
#' 2 as a binomial GLM) -- see `model_type1`/`model_type2`.
#'
#' @param data1,data2 Data frames for traits 1 and 2. Can be the same object
#'   if both traits were measured on the same individuals.
#' @param x1_col,y1_col,z1_col Column names for surface 1.
#' @param x2_col,y2_col,z2_col Column names for surface 2.
#' @param design_id_col1,design_id_col2 Design-point ID columns (can be
#'   `NULL`). Only used when `paired = FALSE`.
#' @param paired `TRUE` (default) if both traits were measured on the same
#'   individuals (resamples individuals jointly, using `data1`/`data2` row
#'   correspondence). `FALSE` for unpaired data, resampled independently
#'   (optionally within design points via `design_id_col1`/`design_id_col2`).
#' @param peak1,peak2 `TRUE` to locate maxima; `FALSE` for minima.
#' @param model_type1,model_type2 One of `"lm"` (default), `"glm"`, or
#'   `"lmer"`, for surfaces 1 and 2 respectively. Need not match -- e.g.
#'   `model_type1 = "glm"` with `family1 = binomial()` for a viability
#'   surface, `model_type2 = "lm"` for a continuous wing-size surface. See
#'   [nutrigon_fit()] for details on how each model type locates its
#'   optimum.
#' @param re_formula1,re_formula2 Random-effects formula strings, required
#'   only when the corresponding `model_type` is `"lmer"`. See
#'   [nutrigon_fit()].
#' @param re_group_cols1,re_group_cols2 Character vector(s) of column
#'   name(s) in `data1`/`data2` that the corresponding `re_formula`
#'   refers to as grouping factors. Required whenever the corresponding
#'   `re_formula` is supplied. See [nutrigon_fit()].
#' @param family1,family2 Convenience arguments forwarded as `family = ` to
#'   [stats::glm()] when the corresponding `model_type` is `"glm"`. Ignored
#'   otherwise. For finer control (e.g. `weights = `, `control = `), use
#'   `extra_args1`/`extra_args2` instead/in addition.
#' @param extra_args1,extra_args2 Named lists of further arguments forwarded
#'   to the underlying fitter for surfaces 1 and 2 respectively (e.g.
#'   `list(weights = ..., control = ...)`).
#' @param B Number of bootstrap resamples (default 200; use 500 for
#'   publication-quality results).
#' @param n_cloud Monte Carlo fallback cloud size.
#' @param n_reps_cloud Number of MC draws per replicate.
#' @param seed Optional random seed for reproducibility.
#'
#' @return An object of class `"nutrigon_comparison"`, a list with:
#'   \describe{
#'     \item{theta_hat}{Point estimate of angular separation (degrees).}
#'     \item{ci}{Named vector `c(lower_95, upper_95)` for the 95\% percentile
#'       bootstrap CI.}
#'     \item{theta1, theta2}{Component angles (degrees from origin).}
#'     \item{boot_diffs}{Vector of `B` bootstrap angular separations.}
#'     \item{n_fallback}{Number of bootstrap resamples that failed (and were
#'       excluded from the CI).}
#'   }
#'
#' @examples
#' data(lee_drosophila)
#'
#' # Both lm() (default) -- identical behaviour to previous package versions
#' result <- nutrigon_compare(
#'   lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
#'   lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
#'   design_id_col1 = "treatment", design_id_col2 = "treatment",
#'   paired = TRUE, B = 50, seed = 1
#' )
#' print(result)
#'
#' \donttest{
#' # Mixed model types: quasi-Poisson GLM for egg count, lm for lifespan
#' result_glm <- nutrigon_compare(
#'   lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
#'   lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
#'   design_id_col1 = "treatment", design_id_col2 = "treatment",
#'   paired = TRUE, B = 50, seed = 1,
#'   model_type2 = "glm", family2 = quasipoisson()
#' )
#' }
#'
#' @seealso [nutrigon_fit()], [nutrigon_check()]
#' @export
nutrigon_compare <- function(data1, x1_col, y1_col, z1_col,
                              data2, x2_col, y2_col, z2_col,
                              design_id_col1 = NULL,
                              design_id_col2 = NULL,
                              paired         = TRUE,
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
                              B              = 200,
                              n_cloud        = 500,
                              n_reps_cloud   = 5,
                              seed           = NULL) {

  model_type1 <- match.arg(model_type1)
  model_type2 <- match.arg(model_type2)

  # Fold family1/family2 convenience args into the extra_args lists that
  # get forwarded to .fit_surface()'s ... -- keeps the common case (just
  # setting a family) simple while still allowing full control via
  # extra_args for anything else (weights, offset, control, etc.)
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

  # Seeding via withr::local_seed() rather than a hand-rolled
  # save/assign/.GlobalEnv pattern: CRAN policy prohibits packages from
  # writing to .GlobalEnv (see CRAN Repository Policy and the CRAN
  # Cookbook's Code Issues chapter), and this is enforced by an automated
  # check on the package's own source text -- even a strictly
  # save-then-restore-on-exit pattern (which never leaves any net change
  # behind) is flagged if `assign(..., envir = .GlobalEnv)` /
  # `.GlobalEnv$.Random.seed <-` appears literally in the package's R
  # code. withr::local_seed() performs the identical save/restore
  # mechanism, but since that pattern lives inside the withr package
  # (itself already an accepted CRAN dependency) rather than in
  # nutrigon's own source, it does not trigger the check.
  if (!is.null(seed)) withr::local_seed(seed)

  # Point estimates
  fit1 <- fit_one(data1, x1_col, y1_col, z1_col, model_type1, re_formula1,
                   re_group_cols1, extra_args1)
  fit2 <- fit_one(data2, x2_col, y2_col, z2_col, model_type2, re_formula2,
                   re_group_cols2, extra_args2)
  opt1 <- .locate_optimum(fit1, peak1, n_cloud, n_reps_cloud)
  opt2 <- .locate_optimum(fit2, peak2, n_cloud, n_reps_cloud)

  theta1    <- .theta_deg(opt1$x, opt1$y)
  theta2    <- .theta_deg(opt2$x, opt2$y)
  theta_hat <- abs(theta1 - theta2)

  # Bootstrap
  resample_one <- function(data, id_col) {
    if (!is.null(id_col)) {
      ids <- unique(data[[id_col]])
      do.call(rbind, lapply(ids, function(id) {
        sub <- data[data[[id_col]] == id, ]
        sub[sample(nrow(sub), nrow(sub), replace = TRUE), ]
      }))
    } else {
      data[sample(nrow(data), nrow(data), replace = TRUE), ]
    }
  }

  boot_diffs <- numeric(B)
  n_fail     <- 0L

  for (b in seq_len(B)) {
    tryCatch({
      if (paired) {
        idx <- sample(nrow(data1), nrow(data1), replace = TRUE)
        b1 <- data1[idx, ]; b2 <- data2[idx, ]
      } else {
        b1 <- resample_one(data1, design_id_col1)
        b2 <- resample_one(data2, design_id_col2)
      }
      f1 <- fit_one(b1, x1_col, y1_col, z1_col, model_type1, re_formula1,
                     re_group_cols1, extra_args1)
      f2 <- fit_one(b2, x2_col, y2_col, z2_col, model_type2, re_formula2,
                     re_group_cols2, extra_args2)
      o1 <- .locate_optimum(f1, peak1, n_cloud, n_reps_cloud)
      o2 <- .locate_optimum(f2, peak2, n_cloud, n_reps_cloud)
      boot_diffs[b] <- abs(.theta_deg(o1$x, o1$y) - .theta_deg(o2$x, o2$y))
    }, error = function(e) {
      boot_diffs[b] <<- NA_real_
      n_fail <<- n_fail + 1L
    })
  }

  valid <- boot_diffs[!is.na(boot_diffs)]
  ci    <- if (length(valid) >= 10L) stats::quantile(valid, c(0.025, 0.975)) else
             c(NA_real_, NA_real_)
  names(ci) <- c("lower_95", "upper_95")

  structure(
    list(theta_hat  = theta_hat,
         ci         = ci,
         theta1     = theta1,
         theta2     = theta2,
         boot_diffs = boot_diffs,
         n_fallback = n_fail,
         B          = B,
         z1_col     = z1_col,
         z2_col     = z2_col,
         model_type1 = model_type1,
         model_type2 = model_type2),
    class = "nutrigon_comparison"
  )
}

#' @export
print.nutrigon_comparison <- function(x, ...) {
  cat("Nutrigonometry angular separation\n")
  cat(sprintf("  Trait 1 (%s, %s): theta = %.2f degrees\n",
              x$z1_col, x$model_type1, x$theta1))
  cat(sprintf("  Trait 2 (%s, %s): theta = %.2f degrees\n",
              x$z2_col, x$model_type2, x$theta2))
  cat(sprintf("  Angular separation: %.2f degrees\n",   x$theta_hat))
  cat(sprintf("  95%% CI: [%.2f, %.2f] degrees\n",
              x$ci["lower_95"], x$ci["upper_95"]))
  cat(sprintf("  Bootstrap: B=%d, %d failed resamples\n",
              x$B, x$n_fallback))
  invisible(x)
}
