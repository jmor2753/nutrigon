#' Compare two nutritional optima and bootstrap the angular separation
#'
#' Fits two quadratic response surfaces and computes the angular separation
#' theta between their nutritional optima, together with a bootstrap
#' percentile confidence interval. Supports paired resampling (both traits
#' measured on the same individuals) and unpaired resampling with optional
#' resampling by design point.
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
#' result <- nutrigon_compare(
#'   lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
#'   lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
#'   design_id_col1 = "treatment", design_id_col2 = "treatment",
#'   paired = TRUE, B = 50, seed = 1
#' )
#' print(result)
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
                              B              = 200,
                              n_cloud        = 500,
                              n_reps_cloud   = 5,
                              seed           = NULL) {

  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit({
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }, add = TRUE)
    } else {
      on.exit({
        if (exists(".Random.seed", envir = .GlobalEnv)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }, add = TRUE)
    }
    set.seed(seed)
  }

  # Point estimates
  fit1 <- .fit_quadratic(data1, x1_col, y1_col, z1_col)
  fit2 <- .fit_quadratic(data2, x2_col, y2_col, z2_col)
  opt1 <- .locate_optimum(fit1, peak1, n_cloud, n_reps_cloud)
  opt2 <- .locate_optimum(fit2, peak2, n_cloud, n_reps_cloud)

  theta1   <- .theta_deg(opt1$x, opt1$y)
  theta2   <- .theta_deg(opt2$x, opt2$y)
  theta_hat <- abs(theta1 - theta2)

  # Bootstrap
  resample_one <- function(data, x_col, y_col, z_col, id_col) {
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
        b1 <- resample_one(data1, x1_col, y1_col, z1_col, design_id_col1)
        b2 <- resample_one(data2, x2_col, y2_col, z2_col, design_id_col2)
      }
      f1 <- .fit_quadratic(b1, x1_col, y1_col, z1_col)
      f2 <- .fit_quadratic(b2, x2_col, y2_col, z2_col)
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
         z2_col     = z2_col),
    class = "nutrigon_comparison"
  )
}

#' @export
print.nutrigon_comparison <- function(x, ...) {
  cat("Nutrigonometry angular separation\n")
  cat(sprintf("  Trait 1 (%s): theta = %.2f degrees\n", x$z1_col, x$theta1))
  cat(sprintf("  Trait 2 (%s): theta = %.2f degrees\n", x$z2_col, x$theta2))
  cat(sprintf("  Angular separation: %.2f degrees\n",   x$theta_hat))
  cat(sprintf("  95%% CI: [%.2f, %.2f] degrees\n",
              x$ci["lower_95"], x$ci["upper_95"]))
  cat(sprintf("  Bootstrap: B=%d, %d failed resamples\n",
              x$B, x$n_fallback))
  invisible(x)
}
