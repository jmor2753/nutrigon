#' Build a prediction grid inside the design hull for surface plotting
#'
#' Constructs an evenly-spaced grid over the range of the fitted data,
#' restricted to points inside the convex hull of the design, with model
#' predictions attached. Intended for use with `ggplot2::geom_raster()`.
#'
#' @param fit A `"nutrigon_fit"` object (from [nutrigon_fit()]).
#' @param n Grid resolution per axis (default 120).
#'
#' @return A data frame with columns named after the original `x_col` and
#'   `y_col`, plus `predicted`. Only points inside the convex hull are
#'   included.
#'
#' @examples
#' data(lee_drosophila)
#' fit_lifespan <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                               "lifespan", design_id_col = "treatment")
#' grid_lifespan <- nutrigon_grid(fit_lifespan, n = 40)
#' head(grid_lifespan)
#'
#' @seealso [nutrigon_hull()], [nutrigon_region_polygon()]
#' @export
nutrigon_grid <- function(fit, n = 120) {
  x_col <- fit$x_col; y_col <- fit$y_col
  d  <- fit$data
  xs <- seq(min(d$.x), max(d$.x), length.out = n)
  ys <- seq(min(d$.y), max(d$.y), length.out = n)
  gd <- expand.grid(.x = xs, .y = ys)
  inside <- .point_in_hull(gd$.x, gd$.y, fit$hull$x, fit$hull$y)
  gd <- gd[inside, ]
  gd$predicted <- stats::predict(fit$model, newdata = gd)
  names(gd)[1:2] <- c(x_col, y_col)
  gd
}

#' Design hull as a closed polygon for plotting
#'
#' Returns the convex hull of the design as a closed polygon data frame,
#' suitable for `ggplot2::geom_polygon()`.
#'
#' @param fit A `"nutrigon_fit"` object (from [nutrigon_fit()]).
#'
#' @return A data frame with columns named after the original `x_col` and
#'   `y_col`, with the first row repeated at the end to close the polygon.
#'
#' @examples
#' data(lee_drosophila)
#' fit_lifespan <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                               "lifespan", design_id_col = "treatment")
#' nutrigon_hull(fit_lifespan)
#'
#' @seealso [nutrigon_grid()], [nutrigon_region_polygon()]
#' @export
nutrigon_hull <- function(fit) {
  df <- data.frame(x = c(fit$hull$x, fit$hull$x[1]),
                   y = c(fit$hull$y, fit$hull$y[1]))
  names(df) <- c(fit$x_col, fit$y_col)
  df
}

#' Peak (or valley) region as a convex hull polygon for plotting
#'
#' Returns the convex hull of the top- (or bottom-) quantile region stored
#' in a `"nutrigon_fit"` object, as a closed polygon data frame suitable for
#' `ggplot2::geom_polygon()`.
#'
#' @param fit A `"nutrigon_fit"` object (from [nutrigon_fit()]).
#'
#' @return A data frame with the same column names as `fit$region`, closed
#'   into a polygon.
#'
#' @examples
#' data(lee_drosophila)
#' fit_lifespan <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                               "lifespan", design_id_col = "treatment")
#' nutrigon_region_polygon(fit_lifespan)
#'
#' @seealso [nutrigon_grid()], [nutrigon_hull()]
#' @export
nutrigon_region_polygon <- function(fit) {
  reg <- fit$region
  x_col <- names(reg)[1]; y_col <- names(reg)[2]
  idx <- grDevices::chull(reg[[x_col]], reg[[y_col]])
  reg_hull <- reg[c(idx, idx[1]), ]
  reg_hull
}
