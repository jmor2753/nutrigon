# ============================================================================
# Internal helpers -- not exported, not documented for users.
# Ported from nutrigonometry_USER_FUNCTIONS.R (self-contained script).
# ============================================================================

# Point-in-polygon test (ray casting) - vectorised over test points
.point_in_hull <- function(px, py, hx, hy) {
  n <- length(hx)
  inside <- rep(FALSE, length(px))
  j <- n
  for (i in seq_len(n)) {
    xi <- hx[i]; yi <- hy[i]
    xj <- hx[j]; yj <- hy[j]
    cond <- ((yi > py) != (yj > py)) &
      (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
    inside <- xor(inside, cond)
    j <- i
  }
  inside
}

# Angle from origin (in degrees, always positive)
.theta_deg <- function(x, y) atan2(y, x) * 180 / pi

# Build the convex hull polygon of a data frame
.hull_polygon <- function(data, x_col, y_col) {
  x <- data[[x_col]]; y <- data[[y_col]]
  idx <- grDevices::chull(x, y)
  list(x = x[idx], y = y[idx],
       vertices = data.frame(x = x[idx], y = y[idx]))
}

# Sample uniformly inside the hull (rejection sampling)
.sample_in_hull <- function(hx, hy, n = 500) {
  xr <- range(hx); yr <- range(hy)
  out_x <- numeric(n); out_y <- numeric(n)
  filled <- 0L
  while (filled < n) {
    batch <- 5L * (n - filled)
    cx <- stats::runif(batch, xr[1], xr[2])
    cy <- stats::runif(batch, yr[1], yr[2])
    ok <- .point_in_hull(cx, cy, hx, hy)
    take <- min(sum(ok), n - filled)
    out_x[filled + seq_len(take)] <- cx[ok][seq_len(take)]
    out_y[filled + seq_len(take)] <- cy[ok][seq_len(take)]
    filled <- filled + take
  }
  data.frame(x = out_x, y = out_y)
}

# Fit a quadratic surface; return coefs, model object, and design hull
.fit_quadratic <- function(data, x_col, y_col, z_col) {
  d <- data.frame(.x = data[[x_col]], .y = data[[y_col]], .z = data[[z_col]])
  d <- d[stats::complete.cases(d), ]
  if (nrow(d) < 6L)
    stop("Need at least 6 complete observations to fit a quadratic surface.")
  mod <- stats::lm(.z ~ .x + .y + I(.x^2) + I(.y^2) + .x:.y, data = d)
  if (any(is.na(stats::coef(mod))))
    stop("Quadratic fit is rank-deficient (too few distinct design points?).")
  hull <- .hull_polygon(d, ".x", ".y")
  list(model = mod, data = d, hull = hull)
}

# Locate the optimum of a fitted quadratic (closed-form or MC fallback)
# Returns a list: x, y, used_fallback (logical)
.locate_optimum <- function(fit, peak = TRUE, n_cloud = 500, n_reps = 5,
                             region_quantile = NULL) {
  if (is.null(region_quantile))
    region_quantile <- if (peak) 0.95 else 0.05

  cf <- stats::coef(fit$model)
  a  <- unname(cf["I(.x^2)"])
  b  <- unname(cf["I(.y^2)"])
  e  <- unname(cf[".x:.y"])
  cc <- unname(cf[".x"])
  d  <- unname(cf[".y"])
  det <- 4 * a * b - e^2

  correct_sign <- if (peak) a < 0 else a > 0
  valid_det    <- det > 0

  if (valid_det && correct_sign) {
    xs <- (e * d - 2 * b * cc) / det
    ys <- (e * cc - 2 * a * d) / det
    if (.point_in_hull(xs, ys, fit$hull$x, fit$hull$y)) {
      return(list(x = xs, y = ys, used_fallback = FALSE))
    }
  }

  # Monte Carlo fallback: sample cloud inside hull, evaluate, take top region
  pool <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    cloud <- .sample_in_hull(fit$hull$x, fit$hull$y, n = n_cloud)
    names(cloud) <- c(".x", ".y")
    cloud$pred <- stats::predict(fit$model, newdata = cloud)
    cutoff <- stats::quantile(cloud$pred, region_quantile)
    region <- if (peak) cloud[cloud$pred >= cutoff, ] else
                         cloud[cloud$pred <= cutoff, ]
    pool[[i]] <- region
  }
  region <- do.call(rbind, pool)
  if (nrow(region) == 0L)
    stop("MC fallback returned empty region -- too few cloud points?")

  opt <- region[which.max(if (peak) region$pred else -region$pred), ]
  list(x = unname(opt$.x), y = unname(opt$.y), used_fallback = TRUE)
}

# Internal helper: build nominal design matrix for grid-building simulations
.build_nominal_design <- function(rail_ratios, n_levels, max_total) {
  total_intakes <- seq(max_total * 0.1, max_total, length.out = n_levels)
  rows <- list()
  for (r in rail_ratios) {
    angle <- atan(r)
    for (tot in total_intakes) {
      rows[[length(rows) + 1L]] <- data.frame(
        id = paste0("r", round(r, 3), "_t", round(tot)),
        nom_x = tot * cos(angle),
        nom_y = tot * sin(angle)
      )
    }
  }
  do.call(rbind, rows)
}

# Internal helper: simulate two correlated trait surfaces for grid building
.simulate_two_surfaces <- function(nominal, cfg, intake_sd, response_sd) {
  n <- nrow(nominal)
  x <- pmax(nominal$nom_x + stats::rnorm(n, 0, intake_sd), 0.1)
  y <- pmax(nominal$nom_y + stats::rnorm(n, 0, intake_sd), 0.1)
  angle1_rad <- cfg$theta1 * pi / 180
  angle2_rad <- cfg$theta2 * pi / 180
  cx1 <- 100 * cos(angle1_rad); cy1 <- 100 * sin(angle1_rad)
  cx2 <- 100 * cos(angle2_rad); cy2 <- 100 * sin(angle2_rad)
  z1 <- -(x - cx1)^2 / (50 * cfg$elong1) - (y - cy1)^2 * cfg$elong1 / 50 +
         stats::rnorm(n, 0, response_sd)
  z2 <- -(x - cx2)^2 / (50 * cfg$elong2) - (y - cy2)^2 * cfg$elong2 / 50 +
         stats::rnorm(n, 0, response_sd)
  data.frame(x = x, y = y, z1 = z1, z2 = z2, id = nominal$id)
}

# With LazyData: true (see DESCRIPTION), nutrigon_validation_grid is
# lazily bound in the package namespace and available directly -- no
# data() call needed. get0() with the package's namespace as envir is the
# standard way to reference it defensively from inside package code.
.load_default_validation_grid <- function() {
  get0("nutrigon_validation_grid", envir = asNamespace("nutrigon"))
}
