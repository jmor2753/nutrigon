#' Build a reusable validation grid for the reliability diagnostic
#'
#' Builds a grid of simulated configurations spanning a range of curvature
#' elongations and angular separations, runs real bootstrap simulations for
#' each, and records ground-truth empirical coverage. The result is used by
#' [nutrigon_check()] to assess whether a specific applied result resembles
#' configurations that validated well or poorly.
#'
#' This function is computationally intensive (it runs real bootstrap
#' simulations for every grid configuration). Build once and reuse for all
#' comparisons made under the same design. The package ships a pre-built
#' grid, `nutrigon_grid_2025` (`seed = 2025`), so most users will not need
#' to call this function directly -- see `data(nutrigon_grid_2025)`.
#'
#' @param n Number of random configurations (default 20 for a quick check;
#'   50-100 recommended for thorough diagnostics).
#' @param n_sim Simulated datasets per configuration (default 8; 20 for more
#'   reliable ground-truth coverage estimates).
#' @param B_boot Bootstrap resamples per simulated dataset (default 50).
#' @param seed Random seed.
#' @param theta_min,theta_max Angular support of your design (degrees).
#'
#' @return A data frame with columns `max_elongation`, `min_edge_distance`,
#'   `angle_diff`, and `coverage` (ground-truth empirical coverage).
#'
#' @examples
#' \donttest{
#' grid <- nutrigon_build_grid(n = 5, n_sim = 3, B_boot = 20, seed = 1)
#' head(grid)
#' }
#'
#' @seealso [nutrigon_check()]
#' @export
nutrigon_build_grid <- function(n         = 20,
                                 n_sim     = 8,
                                 B_boot    = 50,
                                 seed      = 42,
                                 theta_min = 5.71,
                                 theta_max = 84.29) {
  set.seed(seed)
  INTAKE_SD   <- 8
  RESPONSE_SD <- 2
  MAX_TOTAL   <- 200
  N_LEVELS    <- 3
  N_RAILS     <- 5
  rail_ratios <- exp(seq(log(0.1), log(10), length.out = N_RAILS))
  nominal     <- .build_nominal_design(rail_ratios, N_LEVELS, MAX_TOTAL)

  configs <- vector("list", n)
  for (i in seq_len(n)) {
    t1    <- stats::runif(1, theta_min, theta_max)
    t2    <- stats::runif(1, theta_min, theta_max)
    e1    <- stats::runif(1, 1, 40)
    e2    <- stats::runif(1, 1, 40)
    configs[[i]] <- list(theta1 = t1, theta2 = t2, elong1 = e1, elong2 = e2)
  }

  results <- vector("list", n)
  for (i in seq_len(n)) {
    cfg     <- configs[[i]]
    covered <- logical(n_sim)
    for (s in seq_len(n_sim)) {
      tryCatch({
        d <- .simulate_two_surfaces(nominal, cfg, INTAKE_SD, RESPONSE_SD)
        res <- nutrigon_compare(
          d, "x", "y", "z1", d, "x", "y", "z2",
          design_id_col1 = "id", design_id_col2 = "id",
          paired = TRUE, B = B_boot, n_cloud = 200, n_reps_cloud = 3
        )
        covered[s] <- res$ci["lower_95"] <= abs(cfg$theta1 - cfg$theta2) &&
                       abs(cfg$theta1 - cfg$theta2) <= res$ci["upper_95"]
      }, error = function(e) NULL)
    }
    max_elong  <- max(cfg$elong1, cfg$elong2)
    min_edge   <- min(cfg$theta1 - theta_min, theta_max - cfg$theta1,
                       cfg$theta2 - theta_min, theta_max - cfg$theta2)
    results[[i]] <- data.frame(
      max_elongation    = max_elong,
      min_edge_distance = min_edge,
      angle_diff        = abs(cfg$theta1 - cfg$theta2),
      coverage          = mean(covered, na.rm = TRUE)
    )
    if (i %% 5 == 0)
      message(sprintf("  Grid: %d / %d configurations done.", i, n))
  }
  do.call(rbind, results)
}
