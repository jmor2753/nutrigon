# data-raw/lee_drosophila.R
# Builds data/lee_drosophila.rda from the raw Lee et al. (2008) CSV.
# Run once with: source("data-raw/lee_drosophila.R") from the package root.

lee_raw <- read.csv("data-raw/dataset_lee.csv", na.strings = ".")

lee_drosophila <- lee_raw[complete.cases(
  lee_raw[, c("protein_eaten", "carb_eaten", "lifespan", "dailyeggs", "treatment")]
), ]
rownames(lee_drosophila) <- NULL

usethis_available <- requireNamespace("usethis", quietly = TRUE)
if (usethis_available) {
  usethis::use_data(lee_drosophila, overwrite = TRUE)
} else {
  save(lee_drosophila, file = "data/lee_drosophila.rda", compress = "bzip2")
}
