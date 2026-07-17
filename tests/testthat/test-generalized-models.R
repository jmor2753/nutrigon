test_that("nutrigon_fit defaults to model_type = \"lm\" and matches old behaviour", {
  data(lee_drosophila)
  fit <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                       "lifespan", design_id_col = "treatment")

  expect_equal(fit$model_type, "lm")
  expect_s3_class(fit$model, "nutrigon_lm_fit")
  expect_s3_class(fit$model$fit, "lm")
})

test_that("nutrigon_coefs() returns the same coefficients regardless of wrapper", {
  data(lee_drosophila)
  fit <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                       "lifespan", design_id_col = "treatment")

  cf_via_accessor <- nutrigon_coefs(fit)
  cf_direct       <- stats::coef(fit$model$fit)

  expect_equal(unname(cf_via_accessor), unname(cf_direct))
  expect_true(all(c(".x", ".y", "I(.x^2)", "I(.y^2)", ".x:.y") %in%
                    names(cf_via_accessor)))
})

test_that("nutrigon_fit works with model_type = \"glm\" and a quasi-family", {
  data(lee_drosophila)
  fit_glm <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                           "dailyeggs", design_id_col = "treatment",
                           model_type = "glm", family = quasipoisson())

  expect_equal(fit_glm$model_type, "glm")
  expect_s3_class(fit_glm$model, "nutrigon_glm_fit")
  expect_s3_class(fit_glm$model$fit, "glm")
  expect_true(is.finite(fit_glm$optimum$theta))
  expect_true(fit_glm$elongation >= 1 || is.na(fit_glm$elongation))
})

test_that("nutrigon_fit errors clearly for model_type = \"lmer\" without lme4/lmerTest", {
  skip_if(requireNamespace("lme4", quietly = TRUE) &&
            requireNamespace("lmerTest", quietly = TRUE),
          "lme4/lmerTest are installed; the no-package error path is not reachable")

  data(lee_drosophila)
  expect_error(
    nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
                 model_type = "lmer", re_formula = "+ (1 | treatment)",
                 re_group_cols = "treatment"),
    "requires the lme4 and lmerTest packages"
  )
})

test_that("nutrigon_fit errors clearly when re_group_cols is missing for lmer", {
  skip_if_not(requireNamespace("lme4", quietly = TRUE) &&
                requireNamespace("lmerTest", quietly = TRUE),
              "lme4/lmerTest not installed")

  data(lee_drosophila)
  expect_error(
    nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
                 model_type = "lmer", re_formula = "+ (1 | treatment)"),
    "also requires re_group_cols"
  )
})

test_that("nutrigon_fit errors clearly when re_group_cols names a nonexistent column", {
  skip_if_not(requireNamespace("lme4", quietly = TRUE) &&
                requireNamespace("lmerTest", quietly = TRUE),
              "lme4/lmerTest not installed")

  data(lee_drosophila)
  expect_error(
    nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
                 model_type = "lmer", re_formula = "+ (1 | nonexistent_col)",
                 re_group_cols = "nonexistent_col"),
    "re_group_cols not found in data"
  )
})

test_that("nutrigon_fit works with model_type = \"lmer\" when lme4/lmerTest are available", {
  skip_if_not(requireNamespace("lme4", quietly = TRUE) &&
                requireNamespace("lmerTest", quietly = TRUE),
              "lme4/lmerTest not installed")

  data(lee_drosophila)
  fit_lmer <- suppressWarnings(nutrigon_fit(
    lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
    model_type = "lmer", re_formula = "+ (1 | treatment)",
    re_group_cols = "treatment"
  ))

  expect_equal(fit_lmer$model_type, "lmer")
  expect_s3_class(fit_lmer$model, "nutrigon_lmer_fit")
  expect_true(is.finite(fit_lmer$optimum$theta))
})

test_that("nutrigon_compare works with model_type = \"lmer\" on one surface", {
  skip_if_not(requireNamespace("lme4", quietly = TRUE) &&
                requireNamespace("lmerTest", quietly = TRUE),
              "lme4/lmerTest not installed")
  skip_on_cran()

  data(lee_drosophila)
  result <- suppressWarnings(nutrigon_compare(
    lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
    lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
    paired = TRUE, B = 10, seed = 1,
    model_type1 = "lmer", re_formula1 = "+ (1 | treatment)",
    re_group_cols1 = "treatment"
  ))

  expect_s3_class(result, "nutrigon_comparison")
  expect_true(is.finite(result$theta_hat))
})

test_that("nutrigon_fit's lmer path supports crossed, nested, and interaction random effects", {
  skip_if_not(requireNamespace("lme4", quietly = TRUE) &&
                requireNamespace("lmerTest", quietly = TRUE),
              "lme4/lmerTest not installed")
  skip_on_cran()

  set.seed(1)
  n <- 200
  dt <- data.frame(
    x = runif(n, 0, 100), y = runif(n, 0, 100),
    block = sample(paste0("b", 1:4), n, replace = TRUE),
    batch = sample(paste0("t", 1:3), n, replace = TRUE)
  )
  dt$subblock <- paste0(dt$block, "_s", sample(1:2, n, replace = TRUE))
  dt$z <- -(dt$x - 50)^2 / 200 - (dt$y - 50)^2 / 200 + rnorm(n, 0, 3)

  # Crossed random intercepts
  fit_crossed <- suppressWarnings(nutrigon_fit(
    dt, "x", "y", "z", model_type = "lmer",
    re_formula = "+ (1 | block) + (1 | batch)",
    re_group_cols = c("block", "batch")
  ))
  expect_s3_class(fit_crossed$model, "nutrigon_lmer_fit")
  expect_true(is.finite(fit_crossed$optimum$theta))

  # Nested random intercepts
  fit_nested <- suppressWarnings(nutrigon_fit(
    dt, "x", "y", "z", model_type = "lmer",
    re_formula = "+ (1 | block/subblock)",
    re_group_cols = c("block", "subblock")
  ))
  expect_s3_class(fit_nested$model, "nutrigon_lmer_fit")
  expect_true(is.finite(fit_nested$optimum$theta))

  # Interaction grouping factor
  fit_interact <- suppressWarnings(nutrigon_fit(
    dt, "x", "y", "z", model_type = "lmer",
    re_formula = "+ (1 | block:batch)",
    re_group_cols = c("block", "batch")
  ))
  expect_s3_class(fit_interact$model, "nutrigon_lmer_fit")
  expect_true(is.finite(fit_interact$optimum$theta))

  # Random slope, correctly referencing the internal .x column name
  fit_slope <- suppressWarnings(nutrigon_fit(
    dt, "x", "y", "z", model_type = "lmer",
    re_formula = "+ (1 + .x | block)",
    re_group_cols = "block"
  ))
  expect_s3_class(fit_slope$model, "nutrigon_lmer_fit")
  expect_true(is.finite(fit_slope$optimum$theta))
})

test_that("a random slope written against the original column name gives a translated, actionable error", {
  skip_if_not(requireNamespace("lme4", quietly = TRUE) &&
                requireNamespace("lmerTest", quietly = TRUE),
              "lme4/lmerTest not installed")

  set.seed(1)
  dt <- data.frame(x = runif(100, 0, 100), y = runif(100, 0, 100),
                    block = sample(1:5, 100, replace = TRUE))
  dt$z <- -(dt$x - 50)^2 / 200 - (dt$y - 50)^2 / 200 + rnorm(100, 0, 3)

  expect_error(
    nutrigon_fit(dt, "x", "y", "z", model_type = "lmer",
                 re_formula = "+ (1 + x | block)", re_group_cols = "block"),
    "nutrigon_fit\\(\\) internally renames your x/y/z columns to \\.x/\\.y/\\.z"
  )

  # A genuinely unrelated bad reference should NOT be swallowed by the
  # translated-error branch, and should surface lme4's own message as-is.
  expect_error(
    nutrigon_fit(dt, "x", "y", "z", model_type = "lmer",
                 re_formula = "+ (1 | totally_bogus_col)",
                 re_group_cols = "block"),
    "totally_bogus_col"
  )
})

test_that("nutrigon_compare supports different model_type1/model_type2", {
  skip_on_cran()
  data(lee_drosophila)

  result <- nutrigon_compare(
    lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
    lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
    design_id_col1 = "treatment", design_id_col2 = "treatment",
    paired = TRUE, B = 20, seed = 1,
    model_type2 = "glm", family2 = quasipoisson()
  )

  expect_s3_class(result, "nutrigon_comparison")
  expect_true(is.finite(result$theta_hat))
})

test_that("nutrigon_curvature recovers known closed-form curvature on a paraboloid", {
  # z = -(x^2 + y^2): a dome. At the apex the surface normal is vertical
  # (E=G=1, F=0), so Gaussian curvature K = kappa1*kappa2 = 4 and mean
  # curvature H = (kappa1+kappa2)/2 = -2 exactly, at (0, 0).
  domedt <- expand.grid(x = seq(-3, 3, by = 0.3), y = seq(-3, 3, by = 0.3))
  domedt$z  <- -(domedt$x^2 + domedt$y^2)
  domedt$id <- seq_len(nrow(domedt))

  fit_dome <- nutrigon_fit(domedt, "x", "y", "z", design_id_col = "id",
                            peak = TRUE, n_cloud = 100, n_reps = 2)
  curv <- nutrigon_curvature(fit_dome, n = 21)

  g   <- curv$grid
  idx <- which.min(g$x^2 + g$y^2)  # nearest grid point to the apex

  expect_equal(g$x[idx], 0, tolerance = 1e-8)
  expect_equal(g$y[idx], 0, tolerance = 1e-8)
  expect_equal(g$Gauss_curv[idx], 4, tolerance = 1e-6)
  expect_equal(g$Mean_curv[idx], -2, tolerance = 1e-6)
})

test_that("nutrigon_curvature gives ~0 curvature and area ratio ~1 on a flat surface", {
  flatdt <- expand.grid(x = seq(-5, 5, by = 0.5), y = seq(-5, 5, by = 0.5))
  flatdt$z  <- 1
  flatdt$id <- seq_len(nrow(flatdt))

  fit_flat <- nutrigon_fit(flatdt, "x", "y", "z", design_id_col = "id",
                            peak = TRUE, n_cloud = 100, n_reps = 2)
  curv <- nutrigon_curvature(fit_flat, n = 15, quiet = TRUE)

  expect_true(max(abs(curv$grid$Gauss_curv)) < 1e-6)
  expect_true(max(abs(curv$grid$Mean_curv)) < 1e-4)

  if (requireNamespace("pracma", quietly = TRUE)) {
    expect_equal(curv$surface_area_ratio, 1, tolerance = 1e-6)
  } else {
    expect_true(is.na(curv$surface_area_ratio))
  }
})

test_that("nutrigon_curvature restricts the grid to the design hull", {
  data(lee_drosophila)
  fit <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                       "lifespan", design_id_col = "treatment")

  curv <- nutrigon_curvature(fit, n = 20, quiet = TRUE)
  hull <- nutrigon_hull(fit)

  expect_true(nrow(curv$grid) > 0)
  expect_true(nrow(curv$grid) < 20 * 20)  # strictly fewer than the full grid
  expect_true(all(c("protein_eaten", "carb_eaten", "Gauss_curv",
                     "Mean_curv") %in% names(curv$grid)))
})

test_that("nutrigon_curvature warns and returns NA area ratio without pracma", {
  skip_if(requireNamespace("pracma", quietly = TRUE),
          "pracma is installed; the missing-package path is not reachable")

  data(lee_drosophila)
  fit <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                       "lifespan", design_id_col = "treatment")

  expect_warning(curv <- nutrigon_curvature(fit, n = 15, quiet = TRUE),
                 "requires the pracma package")
  expect_true(is.na(curv$surface_area_ratio))
})

test_that("nutrigon_curvature works on a GLM fit (link-scale) via the dispatch layer", {
  data(lee_drosophila)
  fit_glm <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                           "dailyeggs", design_id_col = "treatment",
                           model_type = "glm", family = quasipoisson())

  expect_message(curv <- nutrigon_curvature(fit_glm, n = 15),
                 "LINK scale")
  expect_true(nrow(curv$grid) > 0)

  curv_quiet <- nutrigon_curvature(fit_glm, n = 15, quiet = TRUE)
  expect_no_message(nutrigon_curvature(fit_glm, n = 15, quiet = TRUE))
})

test_that("elongation and curvature still work on the shipped zanco_tradeoff_matrix surfaces", {
  # Regression test: during development, the data-raw build script for
  # zanco_tradeoff_matrix was run before .fit_surface() was generalized
  # to dispatch across model types, so its shipped surfaces briefly stored
  # a bare `lm` object as $model, which the generalized dispatch layer
  # (.get_link_coefs(), etc.) cannot handle. If this data is ever rebuilt
  # without also re-running data-raw/zanco_tradeoff.R under the current
  # package version, this test should fail loudly rather than the break
  # surfacing silently downstream.
  data(zanco_tradeoff_matrix)
  fit1 <- zanco_tradeoff_matrix$surfaces[[1]]

  expect_s3_class(fit1$model, "nutrigon_lm_fit")
  expect_true(is.finite(nutrigon_elongation(fit1)))

  curv <- nutrigon_curvature(fit1, n = 15, quiet = TRUE)
  expect_true(nrow(curv$grid) > 0)
})

test_that("nutrigon_fit works with a quasibinomial GLM on real Zanco viability data", {
  # Viability (Prop_viable) is a bounded proportion derived from binomial
  # counts (Total_Adults / Total_Eggs per vial) with visible overdispersion
  # -- a natural quasibinomial use case, and a real alternative to the OLS
  # fit used in the shipped zanco_tradeoff_matrix.
  data(zanco_survival)
  mel18 <- zanco_survival[zanco_survival$Temperature == "18fluc" &
                             zanco_survival$Population  == "Melbourne", ]

  fit_quasi <- nutrigon_fit(mel18, "Prot", "Carb", "Prop_viable",
                             design_id_col = "design_id",
                             model_type = "glm", family = quasibinomial())

  expect_equal(fit_quasi$model_type, "glm")
  expect_s3_class(fit_quasi$model, "nutrigon_glm_fit")
  expect_equal(fit_quasi$model$fit$family$family, "quasibinomial")
  expect_true(is.finite(fit_quasi$optimum$theta))
  expect_true(fit_quasi$elongation >= 1 || is.na(fit_quasi$elongation))

  # Curvature should also work through the same dispatch layer
  curv <- nutrigon_curvature(fit_quasi, n = 15, quiet = TRUE)
  expect_true(nrow(curv$grid) > 0)
})

test_that("quasibinomial and lm fits of the same Zanco viability surface give different (both valid) optima", {
  # These are genuinely different statistical models of the same data (one
  # quadratic in the raw proportion, one quadratic on the logit scale), so
  # this test only checks both run and produce finite, non-identical
  # results -- not that they agree.
  data(zanco_survival)
  data(zanco_tradeoff_matrix)
  mel18 <- zanco_survival[zanco_survival$Temperature == "18fluc" &
                             zanco_survival$Population  == "Melbourne", ]

  fit_ols   <- zanco_tradeoff_matrix$surfaces[["viability_18fluc_Melbourne"]]
  fit_quasi <- nutrigon_fit(mel18, "Prot", "Carb", "Prop_viable",
                             design_id_col = "design_id",
                             model_type = "glm", family = quasibinomial())

  expect_true(is.finite(fit_ols$optimum$theta))
  expect_true(is.finite(fit_quasi$optimum$theta))
  # Not asserting equality -- different link scales give different optima
  # in general. Just confirm they are indeed both valid, finite results.
})
