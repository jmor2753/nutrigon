#' Per-result reliability diagnostic
#'
#' Compares a specific result's risk-factor vector (maximum curvature
#' elongation, minimum angular edge distance, and angular separation) to a
#' pre-built validation grid of configurations with known simulation-based
#' ground-truth coverage. Reports whether the result resembles
#' configurations that validated well or poorly in simulation.
#'
#' This is a **diagnostic**, not a correction: it does not change the
#' confidence interval returned by [nutrigon_compare()].
#'
#' @param result A `"nutrigon_comparison"` object (from
#'   [nutrigon_compare()]).
#' @param elong1,elong2 Elongation ratios for the two surfaces (from
#'   [nutrigon_elongation()]).
#' @param k Number of nearest validation-grid neighbours (default 10).
#' @param validation_grid A pre-built grid, typically from
#'   [nutrigon_build_grid()]. Defaults to `"default"`, which loads and uses
#'   the grid shipped with the package, `nutrigon_validation_grid`
#'   (`seed = 2025`, 1000 configurations, `B_boot = 100`) -- this is fast
#'   (no simulation is run) and is the right choice for almost all users.
#'   Pass `NULL` to build a fresh grid from scratch instead (slow, several
#'   minutes; only needed if your design deviates from the standard 5-rail
#'   `[5.71, 84.29]` angular support the shipped grid was built under, or
#'   you want a custom `n`/`n_sim`/`B_boot`). Alternatively pass your own
#'   pre-built grid object (e.g. from a prior [nutrigon_build_grid()] call)
#'   to reuse it directly.
#' @param theta_min,theta_max Angular support of the design (degrees).
#'   Defaults to the standard 5-rail design used in validation,
#'   `[5.71, 84.29]`. Only used when a fresh grid is built (i.e. when
#'   `validation_grid = NULL`); ignored when using the shipped or a
#'   supplied grid.
#'
#' @return A list with:
#'   \describe{
#'     \item{flag}{Character description of estimated reliability.}
#'     \item{estimated_coverage}{Numeric, mean coverage among the `k`
#'       nearest validated neighbours.}
#'     \item{risk_factors}{List of the result's own risk-factor values.}
#'     \item{n_neighbours}{Number of neighbours actually used.}
#'     \item{neighbours}{The nearest-neighbour rows of the validation grid.}
#'   }
#'
#' @examples
#' data(lee_drosophila)
#' fit_lifespan <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                               "lifespan", design_id_col = "treatment")
#' fit_eggs <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                           "dailyeggs", design_id_col = "treatment")
#' result <- nutrigon_compare(
#'   lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
#'   lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
#'   design_id_col1 = "treatment", design_id_col2 = "treatment",
#'   paired = TRUE, B = 50, seed = 1
#' )
#' # Uses the shipped validation grid by default -- fast, no simulation.
#' check <- nutrigon_check(result,
#'                          elong1 = nutrigon_elongation(fit_lifespan),
#'                          elong2 = nutrigon_elongation(fit_eggs))
#' check$flag
#'
#' @seealso [nutrigon_compare()], [nutrigon_build_grid()],
#'   [nutrigon_elongation()]
#' @export
nutrigon_check <- function(result,
                            elong1,
                            elong2,
                            k                = 10,
                            validation_grid  = "default",
                            theta_min        = 5.71,
                            theta_max        = 84.29) {

  if (identical(validation_grid, "default")) {
    validation_grid <- .load_default_validation_grid()
  } else if (is.null(validation_grid)) {
    message("Building a fresh validation grid (n=20, ~1-2 min)...")
    validation_grid <- nutrigon_build_grid(n = 20, n_sim = 8, seed = 42,
                                            theta_min = theta_min,
                                            theta_max = theta_max)
  }

  # Risk factors for this result
  max_elong    <- max(elong1, elong2, na.rm = TRUE)
  angle_diff   <- result$theta_hat
  theta1       <- result$theta1; theta2 <- result$theta2
  min_edge     <- min(theta1 - theta_min, theta_max - theta1,
                       theta2 - theta_min, theta_max - theta2)

  rf <- c(max_elongation    = max_elong,
          min_edge_distance = min_edge,
          angle_diff        = angle_diff)

  # Standardise using grid's own SDs
  sds <- apply(validation_grid[, c("max_elongation",
                                    "min_edge_distance",
                                    "angle_diff")], 2, stats::sd)
  sds[sds == 0] <- 1  # guard against zero variance

  dists <- sqrt(
    ((rf["max_elongation"]    - validation_grid$max_elongation)    / sds["max_elongation"])^2    +
    ((rf["min_edge_distance"] - validation_grid$min_edge_distance) / sds["min_edge_distance"])^2 +
    ((rf["angle_diff"]        - validation_grid$angle_diff)        / sds["angle_diff"])^2
  )

  k_use <- min(k, nrow(validation_grid))
  nn    <- order(dists)[seq_len(k_use)]
  est_cov <- mean(validation_grid$coverage[nn], na.rm = TRUE)

  flag <- if (est_cov >= 0.92) {
    "RELIABLE: nearby validated configurations showed coverage near or above nominal 95%."
  } else if (est_cov >= 0.80) {
    "SOME CAUTION: nearby validated configurations showed coverage between 80-92%."
  } else {
    "REAL CAUTION: nearby validated configurations showed coverage below 80%."
  }

  list(flag               = flag,
       estimated_coverage = est_cov,
       risk_factors       = as.list(rf),
       n_neighbours       = k_use,
       neighbours         = validation_grid[nn, ])
}
