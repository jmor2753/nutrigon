#' Lee et al. (2008) Drosophila lifespan and fecundity data
#'
#' Nutrient intake, lifespan, and daily egg production for female
#' \emph{Drosophila melanogaster} reared across a geometric array of
#' protein:carbohydrate diets. Used throughout the package examples and
#' vignette as a worked, paired comparison (both traits measured on the
#' same individual flies).
#'
#' @format A data frame with 973 rows (complete cases) and 11 columns:
#' \describe{
#'   \item{carb_eaten}{Carbohydrate eaten (mg).}
#'   \item{protein_eaten}{Protein eaten (mg).}
#'   \item{lifespan}{Lifespan (days).}
#'   \item{lifetimeegg}{Lifetime egg count.}
#'   \item{dailyeggs}{Daily egg production rate.}
#'   \item{treatment}{Design-point identifier (diet treatment code); used
#'     for within-design-point bootstrap resampling.}
#'   \item{replicate}{Replicate identifier within treatment.}
#'   \item{Ratio}{Nominal protein:carbohydrate ratio, as a string, e.g.
#'     `"(0:1)"`.}
#'   \item{Food}{Total food concentration.}
#'   \item{P_fixed}{Fixed protein concentration of the diet.}
#'   \item{C_fixed}{Fixed carbohydrate concentration of the diet.}
#' }
#'
#' @source Lee, K.P., Simpson, S.J., Clissold, F.J., Brooks, R., Ballard,
#'   J.W.O., Taylor, P.W., Soran, N. & Raubenheimer, D. (2008). Lifespan and
#'   reproduction in Drosophila: new insights from nutritional geometry.
#'   \emph{PNAS}, 105(7), 2498-2503.
#'   \doi{10.1073/pnas.0710787105}
#'
#'   Included with the original authors' permission. Please cite the
#'   original paper when using these data.
"lee_drosophila"


#' Pre-built validation grid (seed = 2025)
#'
#' A validation grid built by [nutrigon_build_grid()]-equivalent simulation
#' with `seed = 2025`, provided so that [nutrigon_check()] can be used
#' immediately without waiting for a fresh simulation. See
#' [nutrigon_build_grid()] for the meaning of each column and the
#' simulation design used to construct it.
#'
#' @format A data frame with 1000 rows and columns:
#' \describe{
#'   \item{max_elongation}{Maximum curvature elongation ratio across the
#'     two simulated surfaces in that configuration.}
#'   \item{min_edge_distance}{Minimum distance (degrees) from either
#'     surface's optimum to the angular boundary of the design.}
#'   \item{angle_diff}{True angular separation (degrees) between the two
#'     simulated optima.}
#'   \item{coverage}{Ground-truth empirical bootstrap CI coverage for that
#'     configuration, from repeated simulation.}
#'   \item{n_valid_sims}{Number of simulated datasets that produced a valid
#'     (non-error) bootstrap result for that configuration, out of the
#'     `n_sim_per_point` attempted.}
#' }
#'
#' Build settings are stored in the `"build_settings"` attribute:
#' `n_grid_points = 1000`, `n_sim_per_point = 20`, `B_boot = 100`,
#' `n_rails = 5`, `seed = 2025`.
#'
#' @seealso [nutrigon_check()], [nutrigon_build_grid()]
"nutrigon_validation_grid"


#' Zanco et al. (2025) Drosophila viability data
#'
#' Egg-to-adult viability (proportion of eggs eclosing to adulthood) for
#' \emph{Drosophila melanogaster}, measured per vial, across three
#' temperature regimes and two locally-adapted populations, on a geometric
#' array of protein:carbohydrate diets. Used in the package vignette
#' alongside [zanco_wingsize] as a worked, unpaired comparison (viability
#' is measured per vial; wing size is measured per individual fly, so the
#' two traits do not share a sampling unit).
#'
#' The single constant-25 degree design point (a single design point with
#' no angular design, present in the raw data) has already been excluded,
#' since a quadratic surface cannot be fit to it.
#'
#' @format A data frame with 1495 rows and 15 columns:
#' \describe{
#'   \item{Code}{Integer vial identifier.}
#'   \item{Block}{Integer experimental block.}
#'   \item{Temperature}{Temperature regime: `"18fluc"`, `"25fluc"`, or
#'     `"28fluc"` (fluctuating around 18, 25, or 28 degrees C).}
#'   \item{Population}{`"Melbourne"` or `"Townsville"`.}
#'   \item{Food}{Total food concentration.}
#'   \item{Ratio}{Nominal protein:carbohydrate ratio, as a string, e.g.
#'     `"(1:8)"`.}
#'   \item{Replicate}{Replicate identifier within design point (character).}
#'   \item{Ratio_Num}{Nominal protein:carbohydrate ratio, as a number
#'     (protein fraction).}
#'   \item{Prot, Carb}{Protein and carbohydrate concentration of the diet.}
#'   \item{Total_Eggs}{Number of eggs placed in the vial (out of which
#'     adults may eclose).}
#'   \item{Total_Adults}{Number of adults that eclosed.}
#'   \item{dead_adults}{Number of eggs that did not eclose as adults.}
#'   \item{Prop_viable}{Proportion of eggs surviving to adulthood in that
#'     vial (the trait of interest).}
#'   \item{design_id}{Design-point identifier (Ratio x Food combination);
#'     used for within-design-point bootstrap resampling.}
#' }
#'
#' @source Zanco, B., Morimoto, J., Cockerell, F., Mirth, C. & Sgro, C.M.
#'   (2025). Nutritional optima for life-history traits vary with
#'   temperature and across locally-adapted populations. \emph{Journal of
#'   Insect Physiology}, 163, 104815. \doi{10.1016/j.jinsphys.2025.104815}
#'
#' @seealso [zanco_wingsize], [zanco_tradeoff_matrix]
"zanco_survival"


#' Zanco et al. (2025) Drosophila wing size data
#'
#' Wing size for individual female \emph{Drosophila melanogaster}, measured
#' across the same three temperature regimes and two populations as
#' [zanco_survival], on the same geometric array of protein:carbohydrate
#' diets. Wing size is measured per individual fly (unlike viability, which
#' is measured per vial), so comparisons between the two traits are
#' unpaired.
#'
#' The single constant-25 degree design point has already been excluded,
#' since a quadratic surface cannot be fit to it.
#'
#' @format A data frame with 2073 rows and 9 columns:
#' \describe{
#'   \item{Food}{Total food concentration.}
#'   \item{Ratio}{Nominal protein:carbohydrate ratio, as a string, e.g.
#'     `"(1:8)"`.}
#'   \item{Replicate}{Integer replicate identifier within design point.}
#'   \item{Prot, Carb}{Protein and carbohydrate concentration of the diet.}
#'   \item{Population}{`"Melbourne"` or `"Townsville"`.}
#'   \item{Temperature}{Temperature regime: `"18fluc"`, `"25fluc"`, or
#'     `"28fluc"`.}
#'   \item{WingSize}{Wing size of the individual fly (the trait of
#'     interest).}
#'   \item{design_id}{Design-point identifier (Ratio x Food combination);
#'     used for within-design-point bootstrap resampling.}
#' }
#'
#' @source Zanco, B., Morimoto, J., Cockerell, F., Mirth, C. & Sgro, C.M.
#'   (2025). Nutritional optima for life-history traits vary with
#'   temperature and across locally-adapted populations. \emph{Journal of
#'   Insect Physiology}, 163, 104815. \doi{10.1016/j.jinsphys.2025.104815}
#'
#' @seealso [zanco_survival], [zanco_tradeoff_matrix]
"zanco_wingsize"


#' Pre-computed pairwise angular trade-off matrix for the Zanco example
#'
#' The full N x N pairwise angular trade-off matrix across all 12 surfaces
#' fitted to [zanco_survival] and [zanco_wingsize] (viability and wing
#' size, each at 3 temperatures x 2 populations), together with the fitted
#' surfaces themselves. The 66 pairwise comparisons are unpaired bootstrap
#' comparisons (`nutrigon_compare(..., paired = FALSE, B = 200)`) and are
#' expensive to compute (several minutes), so this object is shipped
#' pre-built rather than recomputed at package load or in the vignette.
#'
#' @format A list with components:
#' \describe{
#'   \item{Theta}{A 12x12 numeric matrix of angular separations (degrees),
#'     named by surface (`"<trait>_<temperature>_<population>"`, e.g.
#'     `"viability_18fluc_Melbourne"`).}
#'   \item{ci_lower, ci_upper}{Matching 12x12 matrices of 95\% bootstrap CI
#'     bounds.}
#'   \item{pairwise}{The same information as a flat data frame with one row
#'     per pair (66 rows), including `surface_1`, `surface_2`, `trait_1`,
#'     `trait_2`, `temperature_1`, `temperature_2`, `population_1`,
#'     `population_2`, `theta_deg`, `ci_lower`, `ci_upper`.}
#'   \item{surfaces}{A named list of the 12 `"nutrigon_fit"` objects
#'     underlying the matrix, for use with the plot helpers
#'     ([nutrigon_grid()], [nutrigon_hull()], [nutrigon_region_polygon()]).}
#'   \item{build_settings}{List recording `B` (bootstrap resamples),
#'     `temps`, and `pops` used to build the matrix.}
#' }
#'
#' @seealso [zanco_survival], [zanco_wingsize], [nutrigon_compare()]
"zanco_tradeoff_matrix"
