test_that("nutrigon_fit finds a sensible optimum on Lee data", {
  data(lee_drosophila)
  fit <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                       "lifespan", design_id_col = "treatment")

  expect_s3_class(fit, "nutrigon_fit")
  expect_true(all(c("optimum", "model", "data", "hull", "region",
                     "elongation") %in% names(fit)))
  expect_true(is.finite(fit$optimum$theta))
  expect_true(fit$elongation >= 1 || is.na(fit$elongation))
})

test_that("nutrigon_compare returns a valid CI structure", {
  data(lee_drosophila)
  set.seed(1)
  result <- nutrigon_compare(
    lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
    lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
    design_id_col1 = "treatment", design_id_col2 = "treatment",
    paired = TRUE, B = 30, seed = 1
  )

  expect_s3_class(result, "nutrigon_comparison")
  expect_true(is.finite(result$theta_hat))
  expect_length(result$ci, 2)
  expect_identical(names(result$ci), c("lower_95", "upper_95"))
  expect_true(result$ci["lower_95"] <= result$ci["upper_95"] ||
                is.na(result$ci["lower_95"]))
})

test_that("nutrigon_elongation matches fit$elongation", {
  data(lee_drosophila)
  fit <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                       "lifespan", design_id_col = "treatment")
  expect_equal(nutrigon_elongation(fit), fit$elongation)
})

test_that("nutrigon_check works with an explicitly supplied validation grid", {
  data(lee_drosophila)
  data(nutrigon_validation_grid)

  fit1 <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                        "lifespan", design_id_col = "treatment")
  fit2 <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                        "dailyeggs", design_id_col = "treatment")
  result <- nutrigon_compare(
    lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
    lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
    design_id_col1 = "treatment", design_id_col2 = "treatment",
    paired = TRUE, B = 30, seed = 1
  )

  # The shipped grid uses a generic design, so a design-mismatch warning
  # is expected for the Lee et al. design.
  check <- suppressWarnings(
    nutrigon_check(result,
                   elong1 = nutrigon_elongation(fit1),
                   elong2 = nutrigon_elongation(fit2),
                   k = 10,
                   validation_grid = nutrigon_validation_grid))

  expect_true(check$estimated_coverage >= 0 && check$estimated_coverage <= 1)
  expect_true(is.character(check$flag))
  expect_equal(check$n_neighbours, 10)
})

test_that("nutrigon_check defaults to the shipped grid without rebuilding", {
  data(lee_drosophila)
  data(nutrigon_validation_grid)

  fit1 <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                        "lifespan", design_id_col = "treatment")
  fit2 <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                        "dailyeggs", design_id_col = "treatment")
  result <- nutrigon_compare(
    lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
    lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
    design_id_col1 = "treatment", design_id_col2 = "treatment",
    paired = TRUE, B = 30, seed = 1
  )

  # No validation_grid supplied: uses the shipped grid (no rebuild message)
  # and warns that it was simulated on a different design.
  expect_warning(
    expect_no_message(
      check_default <- nutrigon_check(result,
                                      elong1 = nutrigon_elongation(fit1),
                                      elong2 = nutrigon_elongation(fit2))
    ),
    "different|angular support"
  )

  check_explicit <- suppressWarnings(
    nutrigon_check(result,
                   elong1 = nutrigon_elongation(fit1),
                   elong2 = nutrigon_elongation(fit2),
                   validation_grid = nutrigon_validation_grid))

  # Defaulting should give identical results to passing the shipped grid
  # explicitly (same object, same neighbours).
  expect_equal(check_default$estimated_coverage, check_explicit$estimated_coverage)
  expect_equal(check_default$n_neighbours, check_explicit$n_neighbours)
})

test_that("nutrigon_check with validation_grid = NULL rebuilds from scratch and messages", {
  skip_on_cran()
  data(lee_drosophila)

  fit1 <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                        "lifespan", design_id_col = "treatment")
  fit2 <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                        "dailyeggs", design_id_col = "treatment")
  result <- nutrigon_compare(
    lee_drosophila, "protein_eaten", "carb_eaten", "lifespan",
    lee_drosophila, "protein_eaten", "carb_eaten", "dailyeggs",
    design_id_col1 = "treatment", design_id_col2 = "treatment",
    paired = TRUE, B = 30, seed = 1
  )

  # Explicit NULL rebuilds a small generic grid; it warns because the
  # generic design differs from the Lee et al. design.
  expect_message(
    check_fresh <- suppressWarnings(
      nutrigon_check(result,
                     elong1 = nutrigon_elongation(fit1),
                     elong2 = nutrigon_elongation(fit2),
                     validation_grid = NULL,
                     k = 3)),
    "Building a fresh validation grid"
  )
  expect_true(check_fresh$estimated_coverage >= 0 &&
                check_fresh$estimated_coverage <= 1)
})

test_that("plot helpers return correctly named data frames", {
  data(lee_drosophila)
  fit <- nutrigon_fit(lee_drosophila, "protein_eaten", "carb_eaten",
                       "lifespan", design_id_col = "treatment")

  g <- nutrigon_grid(fit, n = 20)
  expect_true(all(c("protein_eaten", "carb_eaten", "predicted") %in% names(g)))

  h <- nutrigon_hull(fit)
  expect_true(all(c("protein_eaten", "carb_eaten") %in% names(h)))
  first_row <- as.numeric(h[1, ])
  last_row  <- as.numeric(h[nrow(h), ])
  expect_equal(first_row, last_row)

  r <- nutrigon_region_polygon(fit)
  expect_equal(ncol(r), ncol(fit$region))
})

test_that("nutrigon_fit errors informatively on too little data", {
  tiny <- data.frame(x = 1:5, y = 1:5, z = 1:5)
  expect_error(nutrigon_fit(tiny, "x", "y", "z"),
               "at least 6 complete observations")
})
