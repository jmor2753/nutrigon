#' Equivalence test for alignment of two nutritional optima
#'
#' Because the angular separation is non-negative, its percentile interval
#' lies above zero even when the true separation is zero, so it cannot by
#' itself show that two optima are aligned. This function implements the
#' equivalence approach: given a smallest angular separation of biological
#' interest `delta`, chosen before looking at the data, it tests
#' \eqn{H_0: \theta \ge \delta} against \eqn{H_1: \theta < \delta} and
#' rejects \eqn{H_0} (concluding equivalence, i.e. alignment within
#' `delta`) when the upper `1 - alpha` percentile of the bootstrap
#' distribution of theta lies below `delta`.
#'
#' The estimator is biased upwards near zero, which inflates the upper
#' bound, so the procedure errs towards not declaring alignment. Choose
#' `delta` above the near-zero zone in which coverage is poor (about 5
#' degrees at realistic signal-to-noise ratios). The calibration of this
#' one-sided bound has not been validated separately by simulation.
#'
#' @param result A `"nutrigon_comparison"` object from [nutrigon_compare()].
#' @param delta Smallest angular separation of biological interest
#'   (degrees), specified in advance.
#' @param alpha Significance level (default 0.05; the bound is the upper
#'   limit of a `1 - 2 * alpha` two-sided percentile interval).
#'
#' @return A list with `upper_bound`, `delta`, `alpha`, `equivalent`
#'   (logical) and a plain-language `conclusion`.
#'
#' @references Schuirmann, D.J. (1987) A comparison of the two one-sided
#'   tests procedure and the power approach for assessing the equivalence
#'   of average bioavailability. \emph{Journal of Pharmacokinetics and
#'   Biopharmaceutics} 15: 657-680.
#'
#' @examples
#' data(lee_drosophila)
#' result <- nutrigon_compare(
#'   lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
#'   lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
#'   paired = TRUE, B = 50, seed = 1)
#' nutrigon_equivalence(result, delta = 10)
#'
#' @seealso [nutrigon_compare()], [nutrigon_check()]
#' @export
nutrigon_equivalence <- function(result, delta, alpha = 0.05) {
  if (!inherits(result, "nutrigon_comparison"))
    stop("result must be a \"nutrigon_comparison\" object from nutrigon_compare().")
  if (!is.numeric(delta) || length(delta) != 1L || delta <= 0)
    stop("delta must be a single positive number (degrees).")
  valid <- result$boot_diffs[!is.na(result$boot_diffs)]
  if (length(valid) < 10L) stop("Too few valid bootstrap replicates.")
  ub <- stats::quantile(valid, 1 - alpha, names = FALSE)
  eq <- ub < delta
  structure(list(
    upper_bound = ub, delta = delta, alpha = alpha, equivalent = eq,
    conclusion = if (eq)
      sprintf("Aligned within %.1f deg: upper %.0f%% bound = %.2f deg < delta.",
              delta, 100 * (1 - alpha), ub) else
      sprintf("Alignment within %.1f deg not shown: upper %.0f%% bound = %.2f deg >= delta (the data are not evidence of alignment).",
              delta, 100 * (1 - alpha), ub)),
    class = "nutrigon_equivalence")
}

#' @export
print.nutrigon_equivalence <- function(x, ...) {
  cat("Equivalence test for alignment of optima\n  ", x$conclusion, "\n", sep = "")
  invisible(x)
}
