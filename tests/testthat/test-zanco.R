test_that("zanco_survival and zanco_wingsize load with expected structure", {
  data(zanco_survival)
  data(zanco_wingsize)

  expect_s3_class(zanco_survival, "data.frame")
  expect_s3_class(zanco_wingsize, "data.frame")

  expect_true(all(c("Temperature", "Population", "Prot", "Carb",
                     "Prop_viable", "design_id") %in% names(zanco_survival)))
  expect_true(all(c("Temperature", "Population", "Prot", "Carb",
                     "WingSize", "design_id") %in% names(zanco_wingsize)))

  # constant-25 design point should already be excluded
  expect_false("25" %in% zanco_survival$Temperature)
  expect_false("25" %in% zanco_wingsize$Temperature)

  expect_setequal(unique(zanco_survival$Temperature),
                   c("18fluc", "25fluc", "28fluc"))
  expect_setequal(unique(zanco_survival$Population),
                   c("Melbourne", "Townsville"))
})

test_that("zanco_tradeoff_matrix has the expected structure and is symmetric", {
  data(zanco_tradeoff_matrix)

  expect_true(all(c("Theta", "ci_lower", "ci_upper", "pairwise",
                     "surfaces", "build_settings") %in%
                    names(zanco_tradeoff_matrix)))

  Theta <- zanco_tradeoff_matrix$Theta
  expect_equal(dim(Theta), c(12L, 12L))
  expect_equal(Theta, t(Theta))
  expect_true(all(diag(Theta) == 0))

  expect_equal(nrow(zanco_tradeoff_matrix$pairwise), choose(12, 2))
  expect_length(zanco_tradeoff_matrix$surfaces, 12)
  expect_true(all(vapply(zanco_tradeoff_matrix$surfaces,
                          inherits, logical(1), "nutrigon_fit")))
})

test_that("zanco_tradeoff_matrix pairwise table is consistent with the Theta matrix", {
  data(zanco_tradeoff_matrix)
  Theta <- zanco_tradeoff_matrix$Theta
  pw    <- zanco_tradeoff_matrix$pairwise

  # spot-check a handful of rows against the matrix
  set.seed(1)
  idx <- sample(nrow(pw), 5)
  for (k in idx) {
    r <- pw[k, ]
    expect_equal(Theta[r$surface_1, r$surface_2], r$theta_deg,
                 tolerance = 1e-6)
  }
})

test_that("hierarchical clustering and MDS run on the shipped Theta matrix", {
  data(zanco_tradeoff_matrix)
  Theta <- zanco_tradeoff_matrix$Theta

  dist_obj <- as.dist(Theta)
  hc <- hclust(dist_obj, method = "average")
  expect_s3_class(hc, "hclust")
  expect_length(hc$labels, 12)

  mds <- cmdscale(dist_obj, k = 2, eig = TRUE)
  expect_equal(dim(mds$points), c(12L, 2L))
})

test_that("plot helpers work on surfaces stored in zanco_tradeoff_matrix", {
  data(zanco_tradeoff_matrix)
  fit <- zanco_tradeoff_matrix$surfaces[[1]]

  g <- nutrigon_grid(fit, n = 20)
  expect_true(nrow(g) > 0)

  h <- nutrigon_hull(fit)
  first_row <- as.numeric(h[1, ])
  last_row  <- as.numeric(h[nrow(h), ])
  expect_equal(first_row, last_row)

  r <- nutrigon_region_polygon(fit)
  expect_true(nrow(r) > 0)
})

test_that("an unpaired nutrigon_compare() call on Zanco data reproduces a shipped pair", {
  skip_on_cran()  # bootstrap-based, kept light with small B for CI speed

  data(zanco_survival)
  data(zanco_wingsize)
  data(zanco_tradeoff_matrix)

  sub_s <- zanco_survival[zanco_survival$Temperature == "18fluc" &
                           zanco_survival$Population  == "Melbourne", ]
  sub_w <- zanco_wingsize[zanco_wingsize$Temperature == "18fluc" &
                           zanco_wingsize$Population  == "Melbourne", ]

  res <- nutrigon_compare(
    sub_s, "Prot", "Carb", "Prop_viable",
    sub_w, "Prot", "Carb", "WingSize",
    design_id_col1 = "design_id", design_id_col2 = "design_id",
    paired = FALSE, B = 30, seed = 0
  )

  shipped <- zanco_tradeoff_matrix$pairwise
  ref <- shipped[shipped$surface_1 == "viability_18fluc_Melbourne" &
                   shipped$surface_2 == "wingsize_18fluc_Melbourne", ]

  # Point estimate is deterministic given the data (bootstrap CI width may
  # differ slightly with smaller B here, but theta_hat should match closely)
  expect_equal(res$theta_hat, ref$theta_deg, tolerance = 0.5)
})
