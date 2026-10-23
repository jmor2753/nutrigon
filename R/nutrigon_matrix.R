#' Angular trade-off matrix for N surfaces, with joint bootstrap uncertainty
#'
#' Fits N response surfaces, locates each optimum, and returns the N x N
#' matrix of pairwise angular separations with percentile bootstrap
#' intervals. Each surface is bootstrapped independently (the surfaces are
#' assumed to come from independent samples, e.g. different traits measured
#' on different units or different environments), and in replicate `b` the
#' `b`-th bootstrap optimum of every surface is used in all pairs in which
#' that surface appears. Each replicate therefore yields a complete,
#' internally consistent matrix, so the dependence among entries (each
#' surface contributes to N - 1 of them) is preserved, and any summary of
#' the matrix can be given bootstrap uncertainty. Optionally, the
#' stability of an average-linkage hierarchical clustering is assessed in
#' this way.
#'
#' With two nutrients, every optimum is described by one polar angle and
#' the matrix is one-dimensional (`theta_ij = |theta_i - theta_j|`): the
#' ordering of the optimum angles (`optima`) carries all of its
#' information, and ordination methods such as multidimensional scaling
#' add no further dimension.
#'
#' @param surfaces A named list. Each element is a list with `data`,
#'   `x_col`, `y_col`, `z_col` and, optionally, `model_type` (default
#'   `"lm"`), `family`, `design_id_col`, `peak` (default `TRUE`),
#'   `re_formula`, `re_group_cols` and `extra_args` (as in
#'   [nutrigon_fit()]).
#' @param B Number of bootstrap replicates (default 200).
#' @param resample `"case"` (default) or `"stratified"` (within design
#'   points; needs `design_id_col` in every surface).
#' @param conf_level Confidence level of the percentile intervals.
#' @param cluster_k If not `NULL` (default 2), assess clustering stability
#'   with the dendrogram cut into `cluster_k` groups.
#' @param n_cloud,n_reps_cloud Monte Carlo fallback settings.
#' @param seed Optional random seed (set locally).
#'
#' @return An object of class `"nutrigon_matrix"`, a list with:
#'   \describe{
#'     \item{Theta}{N x N matrix of point estimates (degrees).}
#'     \item{lower, upper}{N x N matrices of interval limits.}
#'     \item{optima}{Data frame with each surface's optimum location, polar
#'       angle, bootstrap interval for the angle, and whether the Monte
#'       Carlo fallback was used.}
#'     \item{boot_theta}{B x N matrix of bootstrap polar angles.}
#'     \item{pairwise}{Long-format data frame of the N(N-1)/2 entries.}
#'     \item{clustering}{If requested: the point-estimate `hclust` object,
#'       the bootstrap cophenetic correlations, a co-clustering frequency
#'       matrix, and per-surface stability (proportion of replicates in
#'       which a surface clusters with the members of its point-estimate
#'       cluster).}
#'   }
#'
#' @examples
#' \donttest{
#' data(lee_drosophila)
#' surfaces <- list(
#'   lifespan = list(data = lee_drosophila, x_col = "protein_eaten",
#'                   y_col = "carb_eaten", z_col = "lifespan"),
#'   eggs     = list(data = lee_drosophila, x_col = "protein_eaten",
#'                   y_col = "carb_eaten", z_col = "dailyeggs"),
#'   lifetime = list(data = lee_drosophila, x_col = "protein_eaten",
#'                   y_col = "carb_eaten", z_col = "lifetimeegg"))
#' m <- nutrigon_matrix(surfaces, B = 50, seed = 1, cluster_k = NULL)
#' m
#' }
#'
#' @seealso [nutrigon_compare()]
#' @export
nutrigon_matrix <- function(surfaces, B = 200,
                            resample = c("case", "stratified"),
                            conf_level = 0.95, cluster_k = 2,
                            n_cloud = 500, n_reps_cloud = 5, seed = NULL) {
  resample <- match.arg(resample)
  if (!is.list(surfaces) || length(surfaces) < 2L || is.null(names(surfaces)) ||
      any(names(surfaces) == ""))
    stop("surfaces must be a named list of at least two surface specifications.")
  if (!is.null(seed)) withr::local_seed(seed)

  get <- function(s, nm, default = NULL) if (is.null(s[[nm]])) default else s[[nm]]
  fit_spec <- function(s, data) {
    extra <- get(s, "extra_args", list())
    if (!is.null(s$family)) extra$family <- s$family
    f <- do.call(.fit_surface,
                 c(list(data = data, x_col = s$x_col, y_col = s$y_col,
                        z_col = s$z_col, model_type = get(s, "model_type", "lm"),
                        re_formula = s$re_formula, re_group_cols = s$re_group_cols),
                   extra))
    .locate_optimum(f, get(s, "peak", TRUE), n_cloud, n_reps_cloud)
  }
  draw_rows <- function(s) {
    d <- s$data
    if (resample == "stratified") {
      if (is.null(s$design_id_col))
        stop("resample = \"stratified\" needs design_id_col in every surface.")
      idx <- split(seq_len(nrow(d)), d[[s$design_id_col]])
      unlist(lapply(idx, function(ii) ii[sample.int(length(ii), length(ii), replace = TRUE)]),
             use.names = FALSE)
    } else sample.int(nrow(d), nrow(d), replace = TRUE)
  }

  nm <- names(surfaces); N <- length(surfaces)
  pts <- lapply(surfaces, function(s) fit_spec(s, s$data))
  theta <- vapply(pts, function(o) .theta_deg(o$x, o$y), numeric(1))

  boot_theta <- matrix(NA_real_, B, N, dimnames = list(NULL, nm))
  for (j in seq_len(N)) {
    s <- surfaces[[j]]
    for (b in seq_len(B)) {
      boot_theta[b, j] <- tryCatch({
        o <- fit_spec(s, s$data[draw_rows(s), , drop = FALSE])
        .theta_deg(o$x, o$y)
      }, error = function(e) NA_real_)
    }
  }

  alpha <- 1 - conf_level
  Theta <- abs(outer(theta, theta, "-")); dimnames(Theta) <- list(nm, nm)
  lower <- upper <- matrix(0, N, N, dimnames = list(nm, nm))
  pair_rows <- list()
  for (i in seq_len(N - 1L)) for (j in (i + 1L):N) {
    d_b <- abs(boot_theta[, i] - boot_theta[, j]); d_b <- d_b[!is.na(d_b)]
    q <- if (length(d_b) >= 10L) stats::quantile(d_b, c(alpha / 2, 1 - alpha / 2), names = FALSE) else c(NA, NA)
    lower[i, j] <- lower[j, i] <- q[1]; upper[i, j] <- upper[j, i] <- q[2]
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      surface_1 = nm[i], surface_2 = nm[j], theta_deg = Theta[i, j],
      ci_lower = q[1], ci_upper = q[2], n_valid_boot = length(d_b))
  }

  optima <- data.frame(
    surface = nm,
    x = vapply(pts, `[[`, numeric(1), "x"),
    y = vapply(pts, `[[`, numeric(1), "y"),
    theta = theta,
    ci_lower = apply(boot_theta, 2, stats::quantile, alpha / 2, na.rm = TRUE),
    ci_upper = apply(boot_theta, 2, stats::quantile, 1 - alpha / 2, na.rm = TRUE),
    used_fallback = vapply(pts, `[[`, logical(1), "used_fallback"),
    row.names = NULL)

  clustering <- NULL
  if (!is.null(cluster_k) && N > cluster_k) {
    hc <- stats::hclust(stats::as.dist(Theta), method = "average")
    cl0 <- stats::cutree(hc, k = cluster_k)
    coph0 <- stats::cophenetic(hc)
    co <- matrix(0, N, N, dimnames = list(nm, nm)); coph_cor <- numeric(0)
    for (b in seq_len(B)) {
      th <- boot_theta[b, ]
      if (anyNA(th)) next
      D <- abs(outer(th, th, "-")); dimnames(D) <- list(nm, nm)
      hb <- stats::hclust(stats::as.dist(D), method = "average")
      clb <- stats::cutree(hb, k = cluster_k)
      co <- co + outer(clb, clb, "==")
      coph_cor <- c(coph_cor, stats::cor(as.vector(stats::cophenetic(hb)), as.vector(coph0)))
    }
    co <- co / max(1, length(coph_cor))
    stab <- vapply(nm, function(s) {
      mates <- setdiff(names(cl0)[cl0 == cl0[s]], s)
      if (!length(mates)) NA_real_ else mean(co[s, mates])
    }, numeric(1))
    clustering <- list(hclust = hc, k = cluster_k, clusters = cl0,
                       cophenetic_cor = coph_cor, co_clustering = co,
                       stability = stab)
  }

  structure(list(Theta = Theta, lower = lower, upper = upper, optima = optima,
                 boot_theta = boot_theta, pairwise = do.call(rbind, pair_rows),
                 clustering = clustering, B = B, resample = resample,
                 conf_level = conf_level),
            class = "nutrigon_matrix")
}

#' @export
print.nutrigon_matrix <- function(x, ...) {
  cat(sprintf("Angular trade-off matrix: %d surfaces, %d pairs, B = %d (%s resampling)\n",
              nrow(x$Theta), nrow(x$pairwise), x$B, x$resample))
  o <- x$optima[order(x$optima$theta), c("surface", "theta", "ci_lower", "ci_upper")]
  cat("Optimum angles (degrees), ordered:\n")
  print(o, row.names = FALSE, digits = 3)
  if (!is.null(x$clustering)) {
    cc <- x$clustering$cophenetic_cor
    cat(sprintf("Clustering (k = %d): cophenetic correlation with point estimate, median %.2f [%.2f, %.2f]\n",
                x$clustering$k, stats::median(cc),
                stats::quantile(cc, 0.025), stats::quantile(cc, 0.975)))
  }
  invisible(x)
}
