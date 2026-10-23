# ============================================================================
# plot-helpers-generalized.R
#
# Fixes nutrigon_grid() (the only one of the three plot helpers that touched
# the fitted model directly, via stats::predict(fit$model, ...)) and adds
# GLM link-scale messaging to all three, since the values a user is about
# to plot are on the link scale, not the original trait's natural scale,
# whenever model_type = "glm" was used.
#
# nutrigon_hull() and nutrigon_region_polygon() only use fit$hull and
# fit$region, which are plain (x, y) data frames for every model_type. For
# glm fits they print a reminder that the region was defined on the link
# scale.
# ============================================================================

#' Build a prediction grid inside the design hull for surface plotting
#'
#' Constructs an evenly-spaced grid over the range of the fitted data,
#' restricted to points inside the convex hull of the design, with model
#' predictions attached. Intended for use with `ggplot2::geom_raster()`.
#'
#' @section GLM fits: when `fit` was produced with `model_type = "glm"`
#'   (see [nutrigon_fit()]), the `predicted` column is on the **link
#'   scale** (e.g. logit(p) for a binomial GLM, log(mu) for a Poisson/
#'   quasi-Poisson GLM), not the original response scale. This is
#'   deliberate and consistent with how the optimum itself was located
#'   (see `?nutrigon_fit`, Details). A message is printed reminding you of
#'   this every time `nutrigon_grid()` is called on a GLM fit; pass
#'   `quiet = TRUE` to suppress it. If you want a response-scale surface
#'   for plotting (e.g. predicted probabilities rather than logit-odds),
#'   back-transform the `predicted` column yourself with the appropriate
#'   inverse link, e.g. `plogis()` for a logit link or `exp()` for a log
#'   link -- see Examples.
#'
#' @param fit A `"nutrigon_fit"` object (from [nutrigon_fit()]).
#' @param n Grid resolution per axis (default 120).
#' @param quiet Suppress the link-scale reminder message for GLM fits
#'   (default `FALSE`).
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
#' \donttest{
#' # GLM example: predicted is on the logit scale; back-transform for plots
#' # of predicted probability rather than logit-odds
#' fit_glm <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
#'                          "dailyeggs", model_type = "glm",
#'                          family = quasipoisson())
#' grid_glm <- nutrigon_grid(fit_glm, n = 40)
#' grid_glm$predicted_response <- exp(grid_glm$predicted)  # log link -> mu
#' }
#'
#' @seealso [nutrigon_hull()], [nutrigon_region_polygon()], [nutrigon_fit()]
#' @export
nutrigon_grid <- function(fit, n = 120, quiet = FALSE) {
  if (!inherits(fit, "nutrigon_fit"))
    stop("fit must be a \"nutrigon_fit\" object (from nutrigon_fit()).")

  x_col <- fit$x_col; y_col <- fit$y_col
  d  <- fit$data
  xs <- seq(min(d$.x), max(d$.x), length.out = n)
  ys <- seq(min(d$.y), max(d$.y), length.out = n)
  gd <- expand.grid(.x = xs, .y = ys)
  inside <- .point_in_hull(gd$.x, gd$.y, fit$hull$x, fit$hull$y)
  gd <- gd[inside, ]

  # .predict_link() dispatches on fit$model's class (nutrigon_lm_fit /
  # nutrigon_glm_fit / nutrigon_lmer_fit) -- see internal-fit.R.
  gd$predicted <- .predict_link(fit$model, newdata = gd)

  if (fit$model_type == "glm" && !quiet) {
    message("nutrigon_grid(): `predicted` is on the LINK scale (family: ",
            .glm_family_name(fit$model), "), not the original response ",
            "scale. Back-transform with the appropriate inverse link if ",
            "you want response-scale predictions for plotting -- see ",
            "?nutrigon_grid. Pass quiet = TRUE to suppress this message.")
  }

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
  if (!inherits(fit, "nutrigon_fit"))
    stop("fit must be a \"nutrigon_fit\" object (from nutrigon_fit()).")
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
#' @section GLM fits: the region stored in `fit$region` (and therefore
#'   returned here) was extracted on the **link scale** for `model_type =
#'   "glm"` fits, consistent with [nutrigon_grid()] and the optimum-finding
#'   logic in [nutrigon_fit()]. The (x, y) coordinates returned are always
#'   in the original nutrient units regardless of model type -- it is only
#'   the *criterion* used to select which points fall inside the region
#'   (top/bottom quantile of the link-scale prediction) that differs. A
#'   reminder message is printed for GLM fits; pass `quiet = TRUE` to
#'   suppress it.
#'
#' @param fit A `"nutrigon_fit"` object (from [nutrigon_fit()]).
#' @param quiet Suppress the link-scale reminder message for GLM fits
#'   (default `FALSE`).
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
nutrigon_region_polygon <- function(fit, quiet = FALSE) {
  if (!inherits(fit, "nutrigon_fit"))
    stop("fit must be a \"nutrigon_fit\" object (from nutrigon_fit()).")

  reg <- fit$region
  x_col <- names(reg)[1]; y_col <- names(reg)[2]
  idx <- grDevices::chull(reg[[x_col]], reg[[y_col]])
  reg_hull <- reg[c(idx, idx[1]), ]

  if (fit$model_type == "glm" && !quiet) {
    message("nutrigon_region_polygon(): this region was selected using ",
            "the LINK-scale prediction (family: ",
            .glm_family_name(fit$model), "). The (x, y) coordinates ",
            "returned are in the original nutrient units. Pass quiet = ",
            "TRUE to suppress this message.")
  }

  reg_hull
}

# Internal helper: pull a readable family name out of a nutrigon_glm_fit
# wrapper for use in the messages above.
.glm_family_name <- function(model_obj) {
  if (!inherits(model_obj, "nutrigon_glm_fit")) return(NA_character_)
  fam <- stats::family(model_obj$fit)
  paste0(fam$family, "(link = \"", fam$link, "\")")
}
