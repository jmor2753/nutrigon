# ============================================================================
# Internal helpers -- not exported.
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

# Polar angle from the origin (degrees)
.theta_deg <- function(x, y) atan2(y, x) * 180 / pi

# Point at a given polar angle (degrees) and radius
.point_at_angle <- function(theta_deg, r) {
  th <- theta_deg * pi / 180
  c(x = r * cos(th), y = r * sin(th))
}

# Build the convex hull polygon of a data frame
.hull_polygon <- function(data, x_col, y_col) {
  x <- data[[x_col]]; y <- data[[y_col]]
  idx <- grDevices::chull(x, y)
  list(x = x[idx], y = y[idx],
       vertices = data.frame(x = x[idx], y = y[idx]))
}

# Angular support of a set of intake points: the range of polar angles
# (degrees) spanned by their convex hull, as seen from the origin. For a
# convex hull in the positive quadrant this range is attained at the hull
# vertices. Points at the origin (angle undefined) and non-finite values
# are dropped.
.angular_support_xy <- function(x, y) {
  keep <- is.finite(x) & is.finite(y) & !(x == 0 & y == 0)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 3L)
    stop("Need at least 3 non-origin points to compute the angular support.")
  idx <- grDevices::chull(x, y)
  a <- .theta_deg(x[idx], y[idx])
  c(theta_min = min(a), theta_max = max(a))
}

# Sample uniformly inside the hull (rejection sampling)
.sample_in_hull <- function(hx, hy, n = 500) {
  xr <- range(hx); yr <- range(hy)
  out_x <- numeric(n); out_y <- numeric(n)
  filled <- 0L
  attempts <- 0L
  while (filled < n) {
    attempts <- attempts + 1L
    if (attempts > 1000L)
      stop("Could not sample points inside the design hull; check that it ",
           "is a non-degenerate polygon.")
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

# Wilson score interval for a binomial proportion
.wilson_ci <- function(x, n, conf_level = 0.95) {
  if (!is.finite(n) || n <= 0) return(c(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(lower = max(0, centre - half), upper = min(1, centre + half))
}

# Generic rail design used by the default (non-design-matched) validation
# grid: fixed nutrient ratios (rails) x intake levels, one observation per
# design point.
.build_rail_design <- function(rail_ratios, n_levels, max_total,
                               min_total = 15) {
  out <- do.call(rbind, lapply(seq_along(rail_ratios), function(i) {
    r <- rail_ratios[i]
    t_vals <- if (n_levels == 1) max_total else
      seq(min_total, max_total, length.out = n_levels)
    data.frame(nominal_x = t_vals / (1 + r), nominal_y = t_vals * r / (1 + r),
               rail_id = i, level_id = seq_len(n_levels))
  }))
  out$design_id <- paste(out$rail_id, out$level_id, sep = "_")
  out
}

.RAIL_RATIOS_DEFAULT <- c(0.1, 0.3, 1, 3, 10)

# Shipped validation grid (generic 5-rail design). The angular support of
# that design is attached so that nutrigon_check() can warn when it is used
# for a result from a different design.
.load_default_validation_grid <- function() {
  g <- get0("nutrigon_validation_grid", envir = asNamespace("nutrigon"))
  if (!is.null(g) && is.null(attr(g, "angular_support")))
    attr(g, "angular_support") <- c(theta_min = 5.71, theta_max = 84.29)
  if (!is.null(g) && is.null(attr(g, "grid_design")))
    attr(g, "grid_design") <- "default_rails"
  g
}
