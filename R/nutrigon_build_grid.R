#' Build a validation grid for the design-matched coverage check
#'
#' Simulates many two-surface configurations with known true optima, runs
#' the full bootstrap on repeated simulated datasets for each, and records
#' the empirical coverage of the interval. The grid is used by
#' [nutrigon_check()] to estimate the coverage expected for a specific
#' result.
#'
#' @section Two modes:
#' \describe{
#'   \item{Design-matched (recommended; `design = your data`)}{Configurations
#'     are simulated on your own experimental design: the observed design
#'     points (mean intake of each `design_id_col` group) with the observed
#'     number of replicates per point; the observed within-point standard
#'     deviations of intake as intake noise; true optima at angles drawn
#'     uniformly within the angular support of your design hull and at
#'     radial distances drawn uniformly between the `radius_quantiles` of
#'     the observed radial intakes (accepted only if inside the hull). If
#'     `snr` is supplied (see [nutrigon_snr()]), response noise is set so
#'     that each simulated surface has that signal-to-noise ratio.}
#'   \item{Generic (`design = NULL`)}{A generic rail design
#'     (`rail_ratios` x `n_levels`, one observation per design point,
#'     optima at radial distances 145 and 90, intake noise `intake_sd`,
#'     response noise `response_sd`), as in the original validation.}
#' }
#'
#' All configurations are drawn before any simulation, so grids built with
#' the same `seed` but different `snr` or `paired` settings contain
#' identical configurations and can be compared configuration by
#' configuration. Elongation ratios are drawn by `elongation_sampler`.
#'
#' This function is slow: it runs `n_grid_points * n_sim_per_point`
#' bootstraps of `B_boot` replicates each. Build once, save with
#' [saveRDS()], and reuse.
#'
#' @param n_grid_points Number of configurations (default 100).
#' @param n_sim_per_point Simulated datasets per configuration (default 20;
#'   coverage per configuration is estimated in steps of
#'   `1 / n_sim_per_point`).
#' @param B_boot Bootstrap resamples per simulated dataset (default 50).
#' @param design Optional data frame of your own data (switches on the
#'   design-matched mode).
#' @param x_col,y_col,design_id_col Column names in `design` (required with
#'   `design`).
#' @param snr Optional target signal-to-noise ratio(s), one or two positive
#'   numbers (surface 1, surface 2), e.g. from [nutrigon_snr()].
#' @param paired `TRUE` (default) to simulate both traits on the same
#'   individuals; `FALSE` to simulate them on independent individuals and
#'   resample them independently.
#' @param radius_quantiles Quantiles of the observed radial intake between
#'   which true optima are placed (design-matched mode).
#' @param elongation_sampler Function of no arguments returning one
#'   elongation ratio (default: uniform on 1-40).
#' @param rail_ratios,n_levels,intake_sd,response_sd Settings of the generic
#'   mode. `response_sd` is ignored when `snr` is supplied.
#' @param n_cloud,n_reps_cloud Monte Carlo fallback settings passed to
#'   [nutrigon_compare()].
#' @param seed Random seed (set locally).
#' @param verbose Print progress.
#' @param n,n_sim,theta_min,theta_max Deprecated (version 0.1.0 names).
#'   `n` and `n_sim` map onto `n_grid_points` and `n_sim_per_point`;
#'   `theta_min`/`theta_max` are ignored (the support is computed from the
#'   design).
#'
#' @return A data frame with one row per configuration: `max_elongation`,
#'   `min_edge_distance`, `angle_diff`, `coverage`, `n_valid_sims`,
#'   `theta1`, `theta2`, `radius1`, `radius2`. Attributes:
#'   `angular_support`, `grid_design` (`"user_design"` or
#'   `"default_rails"`), `snr`, `paired`.
#'
#' @examples
#' \donttest{
#' data(lee_drosophila)
#' # Tiny grid, for illustration only (use the defaults in practice)
#' fit1 <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten", "lifespan")
#' fit2 <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs")
#' grid <- nutrigon_build_grid(
#'   n_grid_points = 3, n_sim_per_point = 3, B_boot = 20,
#'   design = lee_drosophila, x_col = "protein_eaten", y_col = "carb_eaten",
#'   design_id_col = "treatment",
#'   snr = c(nutrigon_snr(fit1), nutrigon_snr(fit2)), seed = 99)
#' grid
#' }
#'
#' @seealso [nutrigon_check()], [nutrigon_snr()]
#' @export
nutrigon_build_grid <- function(n_grid_points    = 100,
                                n_sim_per_point  = 20,
                                B_boot           = 50,
                                design           = NULL,
                                x_col            = NULL,
                                y_col            = NULL,
                                design_id_col    = NULL,
                                snr              = NULL,
                                paired           = TRUE,
                                radius_quantiles = c(0.25, 0.95),
                                elongation_sampler = function() stats::runif(1, 1, 40),
                                rail_ratios      = c(0.1, 0.3, 1, 3, 10),
                                n_levels         = 3,
                                intake_sd        = 8,
                                response_sd      = 2,
                                n_cloud          = 100,
                                n_reps_cloud     = 3,
                                seed             = 1,
                                verbose          = TRUE,
                                n = NULL, n_sim = NULL,
                                theta_min = NULL, theta_max = NULL) {

  if (!is.null(n))     n_grid_points   <- n
  if (!is.null(n_sim)) n_sim_per_point <- n_sim
  if (!is.null(theta_min) || !is.null(theta_max))
    message("theta_min/theta_max are ignored: the angular support is ",
            "computed from the design.")

  use_design <- !is.null(design)
  if (!is.null(snr)) {
    if (!is.numeric(snr) || !(length(snr) %in% 1:2) ||
        any(!is.finite(snr)) || any(snr <= 0))
      stop("snr must be one or two positive numbers.")
    snr <- rep(snr, length.out = 2)
  }

  if (!use_design) {
    nominal   <- .build_rail_design(rail_ratios, n_levels, max_total = 200)
    theta_min <- atan(min(rail_ratios)) * 180 / pi
    theta_max <- atan(max(rail_ratios)) * 180 / pi
    width_scale <- 1
    gx <- nominal$nominal_x; gy <- nominal$nominal_y
  } else {
    if (is.null(x_col) || is.null(y_col) || is.null(design_id_col))
      stop("The design-matched mode needs x_col, y_col and design_id_col.")
    miss <- setdiff(c(x_col, y_col, design_id_col), names(design))
    if (length(miss))
      stop("Column(s) not found in `design`: ", paste(miss, collapse = ", "))
    d <- design[is.finite(design[[x_col]]) & is.finite(design[[y_col]]) &
                  !is.na(design[[design_id_col]]), , drop = FALSE]
    ids <- as.character(d[[design_id_col]])
    sd0 <- function(v) { s <- stats::sd(v); if (is.na(s)) 0 else s }
    sim_design <- data.frame(
      design_id = ids,
      mx = stats::ave(d[[x_col]], ids, FUN = mean),
      my = stats::ave(d[[y_col]], ids, FUN = mean),
      sx = stats::ave(d[[x_col]], ids, FUN = sd0),
      sy = stats::ave(d[[y_col]], ids, FUN = sd0),
      stringsAsFactors = FALSE)
    support   <- .angular_support_xy(d[[x_col]], d[[y_col]])
    theta_min <- support[["theta_min"]]; theta_max <- support[["theta_max"]]
    hull      <- .hull_polygon(d, x_col, y_col)
    radii     <- sqrt(d[[x_col]]^2 + d[[y_col]]^2)
    r_rng     <- stats::quantile(radii, radius_quantiles, names = FALSE)
    width_scale <- (max(radii) / 200)^2
    gx <- sim_design$mx; gy <- sim_design$my
  }

  draw_peak <- function() {
    if (!use_design) {
      th <- stats::runif(1, theta_min, theta_max)
      return(list(theta = th, r = NA_real_))
    }
    for (tries in seq_len(5000)) {
      th <- stats::runif(1, theta_min, theta_max)
      r  <- stats::runif(1, r_rng[1], r_rng[2])
      p  <- .point_at_angle(th, r)
      if (.point_in_hull(p[["x"]], p[["y"]], hull$x, hull$y))
        return(list(theta = th, r = r))
    }
    stop("Could not place a true optimum inside the design hull; ",
         "try widening radius_quantiles.")
  }

  # Draw every configuration first, so configurations depend only on seed
  cfgs <- withr::with_seed(seed, lapply(seq_len(n_grid_points), function(i) {
    list(pk1 = draw_peak(), pk2 = draw_peak(),
         elong1 = elongation_sampler(), elong2 = elongation_sampler(),
         orient1 = sample(c(TRUE, FALSE), 1), orient2 = sample(c(TRUE, FALSE), 1))
  }))

  rows <- vector("list", n_grid_points)
  for (i in seq_len(n_grid_points)) {
    cf <- cfgs[[i]]
    theta1 <- cf$pk1$theta; theta2 <- cf$pk2$theta
    r1 <- if (use_design) cf$pk1$r else 145
    r2 <- if (use_design) cf$pk2$r else 90
    peak1 <- .point_at_angle(theta1, r1); peak2 <- .point_at_angle(theta2, r2)
    base_width <- 50 * width_scale
    wX1 <- if (cf$orient1) base_width else base_width * cf$elong1
    wY1 <- if (cf$orient1) base_width * cf$elong1 else base_width
    wX2 <- if (cf$orient2) base_width else base_width * cf$elong2
    wY2 <- if (cf$orient2) base_width * cf$elong2 else base_width
    true_diff <- abs(theta1 - theta2)

    surf <- function(x, y, pk, wX, wY) -(x - pk[["x"]])^2 / wX - (y - pk[["y"]])^2 / wY
    if (is.null(snr)) {
      sd1 <- sd2 <- response_sd
    } else {
      sd1 <- max(stats::sd(surf(gx, gy, peak1, wX1, wY1)), 1e-8) / snr[1]
      sd2 <- max(stats::sd(surf(gx, gy, peak2, wX2, wY2)), 1e-8) / snr[2]
    }

    draw_intake <- function() {
      if (!use_design) {
        list(x = pmax(nominal$nominal_x + stats::rnorm(nrow(nominal), 0, intake_sd), 0.1),
             y = pmax(nominal$nominal_y + stats::rnorm(nrow(nominal), 0, intake_sd), 0.1),
             id = nominal$design_id)
      } else {
        n_obs <- nrow(sim_design)
        list(x = pmax(sim_design$mx + stats::rnorm(n_obs, 0, sim_design$sx), 0),
             y = pmax(sim_design$my + stats::rnorm(n_obs, 0, sim_design$sy), 0),
             id = sim_design$design_id)
      }
    }

    covered <- rep(NA, n_sim_per_point)
    for (s in seq_len(n_sim_per_point)) {
      covered[s] <- withr::with_seed(seed * 1e6 + i * 1000 + s, {
        it1 <- draw_intake()
        it2 <- if (paired) it1 else draw_intake()
        z1 <- surf(it1$x, it1$y, peak1, wX1, wY1) + stats::rnorm(length(it1$x), 0, sd1)
        z2 <- surf(it2$x, it2$y, peak2, wX2, wY2) + stats::rnorm(length(it2$x), 0, sd2)
        if (paired) {
          dat1 <- dat2 <- data.frame(x = it1$x, y = it1$y, z1 = z1, z2 = z2, id = it1$id)
        } else {
          dat1 <- data.frame(x = it1$x, y = it1$y, z1 = z1, id = it1$id)
          dat2 <- data.frame(x = it2$x, y = it2$y, z2 = z2, id = it2$id)
        }
        res <- tryCatch(
          suppressWarnings(nutrigon_compare(
            dat1, "x", "y", "z1", dat2, "x", "y", "z2",
            design_id_col1 = "id", design_id_col2 = "id",
            paired = paired, B = B_boot, n_cloud = n_cloud,
            n_reps_cloud = n_reps_cloud)),
          error = function(e) NULL)
        if (is.null(res) || (B_boot - res$n_failed) < max(8, 0.2 * B_boot) ||
            anyNA(res$ci)) NA else
          (res$ci[1] <= true_diff && true_diff <= res$ci[2])
      })
    }

    d1 <- min(theta1 - theta_min, theta_max - theta1)
    d2 <- min(theta2 - theta_min, theta_max - theta2)
    rows[[i]] <- data.frame(
      max_elongation    = max(cf$elong1, cf$elong2),
      min_edge_distance = min(d1, d2),
      angle_diff        = true_diff,
      coverage          = if (all(is.na(covered))) NA_real_ else mean(covered, na.rm = TRUE),
      n_valid_sims      = sum(!is.na(covered)),
      theta1 = theta1, theta2 = theta2, radius1 = r1, radius2 = r2)
    if (verbose)
      message(sprintf("  grid point %3d/%d: elong = %5.1f, edge = %5.1f, angle = %5.1f -> coverage = %.2f",
                      i, n_grid_points, rows[[i]]$max_elongation,
                      rows[[i]]$min_edge_distance, true_diff, rows[[i]]$coverage))
  }

  out <- do.call(rbind, rows)
  attr(out, "angular_support") <- c(theta_min = theta_min, theta_max = theta_max)
  attr(out, "grid_design") <- if (use_design) "user_design" else "default_rails"
  attr(out, "snr")    <- snr
  attr(out, "paired") <- paired
  out
}
