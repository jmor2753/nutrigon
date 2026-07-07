# data-raw/zanco_tradeoff.R
#
# Builds three package data objects from the raw Zanco et al. (2025) CSVs
# and a pre-computed pairwise comparison table:
#   - zanco_survival          cleaned survival/viability data
#   - zanco_wingsize          cleaned wing size data
#   - zanco_tradeoff_matrix   12x12 pairwise angular trade-off matrix
#                             (Theta) with 95% CI bounds, the flat pairwise
#                             table, and the 12 fitted nutrigon_fit() surfaces
#
# The pairwise comparisons (66 unpaired bootstrap comparisons at B = 200)
# are expensive -- several minutes -- so they are computed once and the
# resulting flat table is checked in at data-raw/pairwise_matrix_zanco.rda
# (object `pairwise_df`) rather than recomputed here. This script re-fits
# the 12 surfaces (fast, seconds) and assembles the full matrix from that
# pre-computed table.
#
# To regenerate pairwise_matrix_zanco.rda from scratch (e.g. after changing
# nutrigon_compare()), see the pairwise loop commented out at the bottom of
# this script.
#
# Run once with: source("data-raw/zanco_tradeoff.R") from the package root,
# with the package loaded (devtools::load_all() or library(nutrigon)).

surv <- read.csv("data-raw/dataset_zanco_survival.csv", stringsAsFactors = FALSE)
wing <- read.csv("data-raw/dataset_zanco_wingsize.csv", stringsAsFactors = FALSE)

# Exclude constant-25 (single design point; cannot fit a quadratic surface)
TEMPS <- c("18fluc", "25fluc", "28fluc")
POPS  <- c("Melbourne", "Townsville")
surv  <- surv[surv$Temperature %in% TEMPS, ]
wing  <- wing[wing$Temperature %in% TEMPS, ]

# Add a design ID column (Ratio x Food combination), used for
# within-design-point bootstrap resampling
surv$design_id <- paste(surv$Ratio, surv$Food, sep = "_")
wing$design_id <- paste(wing$Ratio, wing$Food, sep = "_")

zanco_survival <- surv
zanco_wingsize <- wing

# ---- Fit all 12 surfaces (fast) ---------------------------------------------
surfaces <- list()
for (temp in TEMPS) {
  for (pop in POPS) {
    sub_s <- surv[surv$Temperature == temp & surv$Population == pop, ]
    key_s <- paste("viability", temp, pop, sep = "_")
    surfaces[[key_s]] <- nutrigon_fit(
      sub_s, "Prot", "Carb", "Prop_viable",
      design_id_col = "design_id", peak = TRUE
    )

    sub_w <- wing[wing$Temperature == temp & wing$Population == pop, ]
    key_w <- paste("wingsize", temp, pop, sep = "_")
    surfaces[[key_w]] <- nutrigon_fit(
      sub_w, "Prot", "Carb", "WingSize",
      design_id_col = "design_id", peak = TRUE
    )
  }
}

# ---- Load pre-computed pairwise comparisons and assemble matrices ----------
load("data-raw/pairwise_matrix_zanco.rda")  # provides pairwise_df

surface_names <- names(surfaces)
N <- length(surface_names)
Theta_matrix <- matrix(0, N, N, dimnames = list(surface_names, surface_names))
ci_lower_mat <- matrix(NA, N, N, dimnames = list(surface_names, surface_names))
ci_upper_mat <- matrix(NA, N, N, dimnames = list(surface_names, surface_names))

for (k in seq_len(nrow(pairwise_df))) {
  r <- pairwise_df[k, ]
  i <- r$surface_1; j <- r$surface_2
  Theta_matrix[i, j] <- Theta_matrix[j, i] <- r$theta_deg
  ci_lower_mat[i, j] <- ci_lower_mat[j, i] <- r$ci_lower
  ci_upper_mat[i, j] <- ci_upper_mat[j, i] <- r$ci_upper
}

zanco_tradeoff_matrix <- list(
  Theta          = Theta_matrix,
  ci_lower       = ci_lower_mat,
  ci_upper       = ci_upper_mat,
  pairwise       = pairwise_df,
  surfaces       = surfaces,
  build_settings = list(B = 200, temps = TEMPS, pops = POPS)
)

usethis_available <- requireNamespace("usethis", quietly = TRUE)
if (usethis_available) {
  usethis::use_data(zanco_survival, overwrite = TRUE)
  usethis::use_data(zanco_wingsize, overwrite = TRUE)
  usethis::use_data(zanco_tradeoff_matrix, overwrite = TRUE)
} else {
  save(zanco_survival, file = "data/zanco_survival.rda", compress = "bzip2")
  save(zanco_wingsize, file = "data/zanco_wingsize.rda", compress = "bzip2")
  save(zanco_tradeoff_matrix, file = "data/zanco_tradeoff_matrix.rda", compress = "bzip2")
}

# ---- To regenerate pairwise_matrix_zanco.rda from scratch -------------------
# (expensive: ~7 sec/pair x 66 pairs at B = 200; several minutes total)
#
# get_raw <- function(key) {
#   parts <- strsplit(key, "_")[[1]]
#   trait <- parts[1]; temp <- parts[2]; pop <- parts[3]
#   if (trait == "viability") {
#     surv[surv$Temperature == temp & surv$Population == pop, ]
#   } else {
#     wing[wing$Temperature == temp & wing$Population == pop, ]
#   }
# }
# get_zcol <- function(key) if (startsWith(key, "viability")) "Prop_viable" else "WingSize"
# parse_key <- function(k) {
#   p <- strsplit(k, "_")[[1]]
#   list(trait = p[1], temp = p[2], pop = p[3])
# }
#
# pair_results <- list(); pair_id <- 0L
# for (i in seq_len(N - 1L)) {
#   for (j in seq(i + 1L, N)) {
#     si_name <- surface_names[i]; sj_name <- surface_names[j]
#     si <- surfaces[[si_name]]; sj <- surfaces[[sj_name]]
#     raw_i <- get_raw(si_name); raw_j <- get_raw(sj_name)
#     z_i <- get_zcol(si_name); z_j <- get_zcol(sj_name)
#     res <- tryCatch(
#       nutrigon_compare(raw_i, "Prot", "Carb", z_i, raw_j, "Prot", "Carb", z_j,
#                        design_id_col1 = "design_id", design_id_col2 = "design_id",
#                        paired = FALSE, B = 200, seed = pair_id),
#       error = function(e) NULL
#     )
#     theta_pt <- if (is.null(res)) abs(si$optimum$theta - sj$optimum$theta) else res$theta_hat
#     ci_l <- if (is.null(res)) NA_real_ else res$ci["lower_95"]
#     ci_u <- if (is.null(res)) NA_real_ else res$ci["upper_95"]
#     pair_id <- pair_id + 1L
#     mi <- parse_key(si_name); mj <- parse_key(sj_name)
#     pair_results[[pair_id]] <- data.frame(
#       surface_1 = si_name, surface_2 = sj_name,
#       trait_1 = mi$trait, trait_2 = mj$trait,
#       temperature_1 = mi$temp, temperature_2 = mj$temp,
#       population_1 = mi$pop, population_2 = mj$pop,
#       theta_deg = round(theta_pt, 2), ci_lower = round(ci_l, 2), ci_upper = round(ci_u, 2),
#       stringsAsFactors = FALSE
#     )
#   }
# }
# pairwise_df <- do.call(rbind, pair_results)
# save(pairwise_df, file = "data-raw/pairwise_matrix_zanco.rda")
