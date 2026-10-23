lee_result <- function(...) {
  data(lee_drosophila, envir = environment())
  nutrigon_compare(
    lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
    lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
    paired = TRUE, B = 30, seed = 1, ...)
}

test_that("nutrigon_compare returns optima, supports and bootstrap details", {
  res <- lee_result()
  expect_named(res$point1, c("x", "y"))
  expect_named(res$angular_support, c("surface1", "surface2"))
  s1 <- res$angular_support$surface1
  expect_true(s1[["theta_min"]] < s1[["theta_max"]])
  # Optima lie inside (or on) their surface's angular support
  expect_true(res$theta1 >= s1[["theta_min"]] - 1e-8 &&
                res$theta1 <= s1[["theta_max"]] + 1e-8)
  expect_length(res$boot_signed, 30)
  expect_equal(res$boot_diffs, abs(res$boot_signed))
  expect_true(res$fallback_rate >= 0 && res$fallback_rate <= 1)
  expect_identical(res$resample, "case")
  expect_equal(res$n_fallback, res$n_failed)
})

test_that("seeded comparisons are reproducible and leave the global RNG alone", {
  set.seed(123); before <- .Random.seed
  r1 <- lee_result(); r2 <- lee_result()
  expect_identical(.Random.seed, before)
  expect_equal(r1$ci, r2$ci)
})

test_that("conf_level names the interval", {
  res <- lee_result(conf_level = 0.9)
  expect_identical(names(res$ci), c("lower_90", "upper_90"))
})

test_that("paired comparisons require matching rows; stratified needs ids", {
  data(lee_drosophila)
  expect_error(
    nutrigon_compare(lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
                     lee_drosophila[-1, ], "protein_eaten", "carb_eaten", "dailyeggs",
                     paired = TRUE, B = 5),
    "same rows")
  expect_error(
    nutrigon_compare(lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
                     lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
                     resample = "stratified", B = 5),
    "design_id_col")
  res <- nutrigon_compare(
    lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
    lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
    design_id_col1 = "treatment", design_id_col2 = "treatment",
    paired = TRUE, resample = "stratified", B = 20, seed = 2)
  expect_identical(res$resample, "stratified")
})

test_that("nutrigon_check uses the result's own design support and reports a Wilson CI", {
  data(lee_drosophila)
  res <- lee_result()
  # Synthetic grid on the Lee support, with n_valid_sims
  set.seed(5)
  g <- data.frame(max_elongation = runif(40, 1, 40),
                  min_edge_distance = runif(40, 0, 20),
                  angle_diff = runif(40, 0, 60),
                  coverage = runif(40, 0.8, 1), n_valid_sims = 20L)
  attr(g, "angular_support") <- c(theta_min = min(res$angular_support$surface1[1],
                                                  res$angular_support$surface2[1]),
                                  theta_max = max(res$angular_support$surface1[2],
                                                  res$angular_support$surface2[2]))
  expect_no_warning(nutrigon_check(res, 4, 3, validation_grid = g, k = 5))
  chk <- nutrigon_check(res, 4, 3, validation_grid = g, k = 5)
  expect_s3_class(chk, "nutrigon_check")
  expect_true(all(chk$neighbours$angle_diff >= 5))
  expect_true(chk$coverage_ci[["lower"]] <= chk$estimated_coverage &&
                chk$estimated_coverage <= chk$coverage_ci[["upper"]])
  expect_equal(chk$n_sims_pooled, 100)
  # Edge distance is measured from the result's own support
  expected_edge <- with(res, min(theta1 - angular_support$surface1[["theta_min"]],
                                 angular_support$surface1[["theta_max"]] - theta1,
                                 theta2 - angular_support$surface2[["theta_min"]],
                                 angular_support$surface2[["theta_max"]] - theta2))
  expect_equal(chk$risk_factors$min_edge_distance, expected_edge)
})

test_that("nutrigon_check does not assess near-zero results", {
  res <- lee_result()
  res$theta1 <- res$theta2 + 1; res$theta_hat <- 1
  g <- data.frame(max_elongation = 1:20, min_edge_distance = 1:20,
                  angle_diff = seq(1, 40, length.out = 20), coverage = 0.9)
  chk <- nutrigon_check(res, 2, 2, validation_grid = g)
  expect_true(is.na(chk$estimated_coverage))
  expect_match(chk$flag, "NEAR-ZERO")
})

test_that("nutrigon_check warns for a grid built on a different design", {
  res <- lee_result()
  g <- data.frame(max_elongation = 1:20, min_edge_distance = 1:20,
                  angle_diff = seq(6, 60, length.out = 20), coverage = 0.9)
  attr(g, "angular_support") <- c(theta_min = 5.71, theta_max = 84.29)
  expect_warning(nutrigon_check(res, 2, 2, validation_grid = g), "angular support")
})

test_that("nutrigon_snr is positive for lm fits and errors for glm fits", {
  data(lee_drosophila)
  fit <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten", "lifespan")
  expect_gt(nutrigon_snr(fit), 0)
  fit_glm <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                          "dailyeggs", model_type = "glm",
                          family = quasipoisson())
  expect_error(nutrigon_snr(fit_glm), "Gaussian")
})

test_that("nutrigon_equivalence uses the upper bootstrap percentile", {
  res <- lee_result()
  eq <- nutrigon_equivalence(res, delta = 90)
  expect_true(eq$equivalent)
  expect_equal(eq$upper_bound,
               unname(quantile(res$boot_diffs, 0.95, na.rm = TRUE)))
  expect_false(nutrigon_equivalence(res, delta = 0.01)$equivalent)
  expect_error(nutrigon_equivalence(res, delta = -1))
})

test_that("nutrigon_build_grid design-matched mode returns a tagged grid", {
  skip_on_cran()
  data(lee_drosophila)
  g <- nutrigon_build_grid(n_grid_points = 2, n_sim_per_point = 2, B_boot = 12,
                           design = lee_drosophila, x_col = "protein_eaten",
                           y_col = "carb_eaten", design_id_col = "treatment",
                           snr = c(1, 0.7), seed = 3, verbose = FALSE)
  expect_equal(nrow(g), 2)
  expect_true(all(c("coverage", "n_valid_sims", "radius1", "radius2") %in% names(g)))
  expect_identical(attr(g, "grid_design"), "user_design")
  # Same seed, different SNR: identical configurations
  g2 <- nutrigon_build_grid(n_grid_points = 2, n_sim_per_point = 2, B_boot = 12,
                            design = lee_drosophila, x_col = "protein_eaten",
                            y_col = "carb_eaten", design_id_col = "treatment",
                            seed = 3, verbose = FALSE)
  expect_equal(g$theta1, g2$theta1)
  expect_equal(g$radius2, g2$radius2)
})

test_that("nutrigon_matrix preserves joint replicates", {
  skip_on_cran()
  data(lee_drosophila)
  sp <- list(
    a = list(data = lee_drosophila, x_col = "protein_eaten", y_col = "carb_eaten", z_col = "lifespan"),
    b = list(data = lee_drosophila, x_col = "protein_eaten", y_col = "carb_eaten", z_col = "dailyeggs"),
    c = list(data = lee_drosophila, x_col = "protein_eaten", y_col = "carb_eaten", z_col = "lifetimeegg"))
  m <- nutrigon_matrix(sp, B = 15, seed = 4, cluster_k = NULL)
  expect_equal(dim(m$Theta), c(3, 3))
  expect_equal(unname(diag(m$Theta)), c(0, 0, 0))
  expect_equal(m$Theta["a", "b"], abs(m$optima$theta[1] - m$optima$theta[2]))
  expect_equal(dim(m$boot_theta), c(15, 3))
  expect_equal(nrow(m$pairwise), 3)
})
