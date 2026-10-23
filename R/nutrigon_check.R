#' Design-matched coverage check for one result
#'
#' Estimates the coverage expected for the interval of a specific result
#' by averaging the simulated coverage of the `k` validation-grid
#' configurations closest to it in three risk factors: maximum curvature
#' elongation, minimum distance of either optimum from the edge of its
#' design's angular support, and the angle difference. Distances are
#' standardised by the grid's own standard deviations.
#'
#' @section Recommended use:
#' Build the grid on your own design and at the signal-to-noise ratio of
#' your fitted surfaces with [nutrigon_build_grid()] (`design = your data`,
#' `snr = c(nutrigon_snr(fit1), nutrigon_snr(fit2))`). The grid shipped
#' with the package (`validation_grid = "default"`) was simulated on a
#' generic 5-rail design with fixed radial positions and noise; a warning
#' is issued when its angular support differs from that of your result.
#'
#' @section Near-zero angles:
#' Near zero, the unsigned percentile interval cannot contain the true
#' value by construction. Grid configurations with an angle difference
#' below `min_angle` are therefore excluded from the neighbour search, and
#' results with `theta_hat < min_angle` are not given a coverage estimate.
#' To assess whether two optima are aligned, use [nutrigon_equivalence()].
#'
#' @section Limitations:
#' The check is conditional on the grid's generating model (quadratic
#' surfaces with Gaussian noise), does not account for model
#' misspecification, and is only as fine-grained as the grid near the
#' result. It is a diagnostic and does not change the interval.
#'
#' @param result A `"nutrigon_comparison"` object from [nutrigon_compare()].
#' @param elong1,elong2 Elongation ratios of the two surfaces (from
#'   [nutrigon_elongation()]).
#' @param k Number of nearest grid configurations (default 10).
#' @param validation_grid A grid from [nutrigon_build_grid()];
#'   `"default"` (the shipped generic grid); or `NULL` to build a small
#'   generic grid on the fly (slow; for testing only).
#' @param min_angle Near-zero threshold in degrees (default 5).
#' @param conf_level Confidence level of the Wilson interval for the
#'   estimated coverage (default 0.95).
#' @param theta_min,theta_max Angular support to use for the edge distance
#'   only if `result` does not carry its own (objects from version 0.1.0).
#'
#' @return An object of class `"nutrigon_check"`, a list with:
#'   \describe{
#'     \item{flag}{Plain-language summary.}
#'     \item{estimated_coverage}{Mean simulated coverage of the `k` nearest
#'       eligible configurations (`NA` for a near-zero result).}
#'     \item{coverage_ci}{Wilson interval for `estimated_coverage`, pooling
#'       the simulated datasets of the neighbours (`NA` if the grid has no
#'       `n_valid_sims` column).}
#'     \item{n_sims_pooled}{Number of simulated datasets behind
#'       `coverage_ci`.}
#'     \item{risk_factors}{The result's own risk factors.}
#'     \item{n_neighbours}{Number of neighbours used.}
#'     \item{n_grid_points_eligible}{Grid configurations with an angle
#'       difference of at least `min_angle`.}
#'     \item{neighbours, neighbour_distances}{The neighbouring grid rows and
#'       their standardised distances.}
#'   }
#'
#' @examples
#' data(lee_drosophila)
#' fit_lifespan <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                              "lifespan", design_id_col = "treatment")
#' fit_eggs <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                          "dailyeggs", design_id_col = "treatment")
#' result <- nutrigon_compare(
#'   lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
#'   lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
#'   paired = TRUE, B = 50, seed = 1
#' )
#' # Quick look with the shipped generic grid (warns: different design)
#' check <- suppressWarnings(
#'   nutrigon_check(result, elong1 = nutrigon_elongation(fit_lifespan),
#'                  elong2 = nutrigon_elongation(fit_eggs)))
#' check
#'
#' \donttest{
#' # Recommended: a grid matched to the design and signal-to-noise ratio
#' grid <- nutrigon_build_grid(
#'   n_grid_points = 10, n_sim_per_point = 5, B_boot = 20,   # use defaults in practice
#'   design = lee_drosophila, x_col = "protein_eaten", y_col = "carb_eaten",
#'   design_id_col = "treatment",
#'   snr = c(nutrigon_snr(fit_lifespan), nutrigon_snr(fit_eggs)), seed = 99)
#' nutrigon_check(result, nutrigon_elongation(fit_lifespan),
#'                nutrigon_elongation(fit_eggs), validation_grid = grid)
#' }
#'
#' @seealso [nutrigon_build_grid()], [nutrigon_compare()],
#'   [nutrigon_equivalence()]
#' @export
nutrigon_check <- function(result,
                           elong1,
                           elong2,
                           k               = 10,
                           validation_grid = "default",
                           min_angle       = 5,
                           conf_level      = 0.95,
                           theta_min       = 5.71,
                           theta_max       = 84.29) {

  if (!inherits(result, "nutrigon_comparison"))
    stop("result must be a \"nutrigon_comparison\" object from nutrigon_compare().")

  if (is.character(validation_grid) && identical(validation_grid, "default")) {
    validation_grid <- .load_default_validation_grid()
  } else if (is.null(validation_grid)) {
    message("Building a fresh validation grid (small, generic design)...")
    validation_grid <- nutrigon_build_grid(n_grid_points = 20,
                                           n_sim_per_point = 8, seed = 42,
                                           verbose = FALSE)
  }
  need <- c("max_elongation", "min_edge_distance", "angle_diff", "coverage")
  if (!is.data.frame(validation_grid) || !all(need %in% names(validation_grid)))
    stop("validation_grid must be a data frame from nutrigon_build_grid().")

  # ---- risk factors of this result -----------------------------------------
  th1 <- result$theta1; th2 <- result$theta2
  if (!is.null(result$angular_support)) {
    s1 <- result$angular_support$surface1; s2 <- result$angular_support$surface2
  } else {
    s1 <- s2 <- c(theta_min = theta_min, theta_max = theta_max)
    warning("result has no angular_support (created with nutrigon 0.1.0?); ",
            "using theta_min/theta_max. Re-run nutrigon_compare() to use ",
            "your design's own support.", call. = FALSE)
  }
  rf <- list(
    max_elongation    = max(elong1, elong2, na.rm = TRUE),
    min_edge_distance = min(th1 - s1[["theta_min"]], s1[["theta_max"]] - th1,
                            th2 - s2[["theta_min"]], s2[["theta_max"]] - th2),
    angle_diff        = abs(th1 - th2))

  # ---- does the grid describe the same design? -----------------------------
  gs <- attr(validation_grid, "angular_support")
  if (!is.null(gs)) {
    rs <- c(min(s1[["theta_min"]], s2[["theta_min"]]),
            max(s1[["theta_max"]], s2[["theta_max"]]))
    if (any(abs(unname(gs) - rs) > 1))
      warning(sprintf(paste0(
        "The validation grid was simulated on angular support [%.2f, %.2f] deg, ",
        "but this result's design spans [%.2f, %.2f] deg. Build a design-matched ",
        "grid with nutrigon_build_grid(design = your_data, ...)."),
        gs[1], gs[2], rs[1], rs[2]), call. = FALSE)
  }
  if (is.finite(rf$max_elongation) &&
      rf$max_elongation > max(validation_grid$max_elongation, na.rm = TRUE))
    warning(sprintf(paste0("Maximum elongation (%.1f) lies outside the range ",
                           "covered by the grid (up to %.1f); the surface is ",
                           "ridge-like and the estimate is an extrapolation."),
                    rf$max_elongation,
                    max(validation_grid$max_elongation, na.rm = TRUE)),
            call. = FALSE)

  eligible <- which(validation_grid$angle_diff >= min_angle &
                      !is.na(validation_grid$coverage))

  base <- list(risk_factors = rf,
               n_grid_points_total = nrow(validation_grid),
               n_grid_points_eligible = length(eligible))

  # ---- near-zero result ----------------------------------------------------
  if (rf$angle_diff < min_angle) {
    out <- c(list(
      flag = sprintf(paste0(
        "NEAR-ZERO ANGLE (%.2f deg < %.0f deg): coverage is not estimated. The ",
        "unsigned interval cannot show that the optima are aligned; use ",
        "nutrigon_equivalence() with a pre-specified smallest angle of interest."),
        rf$angle_diff, min_angle),
      estimated_coverage = NA_real_,
      coverage_ci = c(lower = NA_real_, upper = NA_real_),
      n_sims_pooled = 0L, n_neighbours = 0L,
      neighbours = validation_grid[0, ], neighbour_distances = numeric(0)), base)
    return(structure(out, class = "nutrigon_check"))
  }
  if (!length(eligible))
    stop("No grid configurations with angle difference >= min_angle.")

  g <- validation_grid[eligible, , drop = FALSE]
  vars <- c("max_elongation", "min_edge_distance", "angle_diff")
  sds <- vapply(vars, function(v) stats::sd(g[[v]]), numeric(1))
  sds[!is.finite(sds) | sds == 0] <- 1
  dists <- sqrt(Reduce(`+`, lapply(vars, function(v)
    ((rf[[v]] - g[[v]]) / sds[[v]])^2)))

  k_use <- min(k, nrow(g))
  nn    <- order(dists)[seq_len(k_use)]
  est_cov <- mean(g$coverage[nn])

  if ("n_valid_sims" %in% names(g)) {
    n_pool <- sum(g$n_valid_sims[nn])
    x_pool <- sum(round(g$coverage[nn] * g$n_valid_sims[nn]))
    cov_ci <- .wilson_ci(x_pool, n_pool, conf_level)
  } else {
    n_pool <- NA_integer_; cov_ci <- c(lower = NA_real_, upper = NA_real_)
  }

  flag <- if (est_cov >= 0.92) {
    sprintf("CLOSE TO NOMINAL: similar simulated configurations had ~%.0f%% coverage.", 100 * est_cov)
  } else if (est_cov >= 0.80) {
    sprintf("SOME CAUTION: similar simulated configurations had ~%.0f%% coverage (below nominal 95%%).", 100 * est_cov)
  } else {
    sprintf("REAL CAUTION: similar simulated configurations had only ~%.0f%% coverage; treat the interval as a rough guide.", 100 * est_cov)
  }

  out <- c(list(flag = flag,
                estimated_coverage = est_cov,
                coverage_ci = cov_ci,
                n_sims_pooled = n_pool,
                n_neighbours = k_use,
                neighbours = g[nn, , drop = FALSE],
                neighbour_distances = dists[nn]), base)
  structure(out, class = "nutrigon_check")
}

#' @export
print.nutrigon_check <- function(x, ...) {
  cat("Design-matched coverage check\n")
  cat("  ", x$flag, "\n", sep = "")
  if (is.finite(x$estimated_coverage)) {
    cat(sprintf("  Estimated coverage: %.1f%%", 100 * x$estimated_coverage))
    if (all(is.finite(x$coverage_ci)))
      cat(sprintf(" [%.1f, %.1f] (%d simulated datasets)",
                  100 * x$coverage_ci[1], 100 * x$coverage_ci[2],
                  as.integer(x$n_sims_pooled)))
    cat(sprintf("\n  Neighbours: %d of %d eligible grid configurations\n",
                x$n_neighbours, x$n_grid_points_eligible))
  }
  rf <- x$risk_factors
  cat(sprintf("  Risk factors: max elongation = %.2f, min edge distance = %.2f deg, angle = %.2f deg\n",
              rf$max_elongation, rf$min_edge_distance, rf$angle_diff))
  invisible(x)
}
