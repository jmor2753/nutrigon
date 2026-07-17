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

# NOTE: an earlier lm-only .fit_quadratic()/.locate_optimum() design
# lived here during development. They are superseded by .fit_surface()
# (internal-fit.R, which dispatches to lm/glm/lmer) and .locate_optimum()
# (now in internal-locate-optimum.R, reading coefficients/predictions
# through the .get_link_coefs()/.predict_link() dispatch layer instead of
# calling stats::coef()/stats::predict() on a bare lm object directly).


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
