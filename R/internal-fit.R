# ============================================================================
# internal-fit.R
#
# Model-fitting layer for nutrigon_fit(): a dispatcher supporting lm(),
# glm() (any family, including quasi-families) and lmerTest::lmer().
#
# DESIGN NOTE (read before touching this file):
# -----------------------------------------------------------------------
# The closed-form stationary-point algebra in .locate_optimum() operates on
# a plain quadratic in two variables:
#
#   linear_predictor(.x, .y) = beta0 + b_x .x + b_y .y
#                               + a_xx .x^2 + a_yy .y^2 + a_xy .x:.y
#
# This is exact for lm() (identity link, response scale = linear predictor
# scale). For glm() it is only exact on the LINK scale -- e.g. for a
# binomial GLM with a logit link, the surface that is genuinely quadratic
# in (.x, .y) is logit(p), not p itself. Locating the stationary point of
# the quadratic therefore means differentiating the LINEAR PREDICTOR, and
# admissibility (interior + correct curvature sign) is also evaluated on
# the link scale. Because the link function is monotonic, the location of
# the constrained optimum is the same whether you maximise the link-scale
# quadratic or the response-scale prediction -- so this loses nothing, and
# choosing the link scale is what makes the SAME closed-form algebra reused
# without modification.
#
# For lmer() (via lme4), "coefficients" means FIXED EFFECTS ONLY
# (lme4::fixef()), and predictions used to locate the population-level
# optimum use re.form = NA (average over random effects, not conditional
# on any particular group/subject). A random-effects structure is supplied
# by the user via the `re_formula` argument.
#
# Every model type below must implement two S3 methods:
#   .get_link_coefs(fit_obj)          -> named numeric vector with the five
#                                         canonical names: ".x", ".y",
#                                         "I(.x^2)", "I(.y^2)", ".x:.y"
#   .predict_link(fit_obj, newdata)   -> numeric vector of link-scale (or
#                                         identity-scale for lm) predictions
#
# Both are dispatched on the *nutrigon internal wrapper class*
# (nutrigon_lm_fit / nutrigon_glm_fit / nutrigon_lmer_fit), not on the raw
# model object's own class, so that the wrapper can carry model-specific
# extras (e.g. the family, for glm) alongside the fitted object without
# polluting the fitted object itself.
# ============================================================================


# ============================================================================
# .fit_surface(): dispatcher.
# ============================================================================

#' Fit the quadratic surface (internal dispatcher)
#'
#' @param data       Data frame.
#' @param x_col,y_col,z_col  Column names .
#' @param model_type One of "lm", "glm", "lmer".
#' @param re_formula For model_type = "lmer" only: the random-effects part
#'                    of the formula, as a raw string starting with "+",
#'                    e.g. "+ (1 | block)" or "+ (1 + .x | block)". Appended
#'                    to the fixed-effects quadratic formula.
#' @param re_group_cols For model_type = "lmer" only: character vector of
#'                    original column name(s) in `data` referenced by
#'                    `re_formula` as grouping factors (e.g. "block" for
#'                    re_formula = "+ (1 | block)"). Required whenever
#'                    re_formula is supplied, because internally the data
#'                    is reduced to just .x/.y/.z before fitting -- any
#'                    column re_formula needs must be explicitly carried
#'                    through via this argument, or lme4::lmer() will fail
#'                    with "object '<name>' not found".
#' @param ...        Forwarded to the underlying fitter: `family = ` and
#'                    `weights = ` for glm(); `weights = `, `control = `,
#'                    etc. for lmer(). Ignored for lm() (accepted for
#'                    interface uniformity, but lm() takes no such args
#'                    here beyond what's already fixed).
#'
#' @return list(model_obj = <wrapper>, data = <cleaned data>, hull = <hull>)
#'   where model_obj is a nutrigon_*_fit wrapper object (S3-dispatchable),
#'   NOT the raw lm/glm/lmer object directly -- callers should never call
#'   stats::coef()/stats::predict() on it themselves; always go through
#'   .get_link_coefs() / .predict_link().
#' @noRd
.fit_surface <- function(data, x_col, y_col, z_col,
                          model_type = c("lm", "glm", "lmer"),
                          re_formula = NULL,
                          re_group_cols = NULL,
                          ...) {
  model_type <- match.arg(model_type)

  d <- data.frame(.x = data[[x_col]], .y = data[[y_col]], .z = data[[z_col]])

  # lmer's re_formula references grouping columns by their ORIGINAL names
  # (e.g. "block" in "+ (1 | block)"), but the .x/.y/.z data frame above
  # drops every column except the three renamed ones. Without this, lmer()
  # fails downstream with "object '<name>' not found" -- carry the
  # requested grouping column(s) through by their original names so the
  # random-effects formula can find them.
  if (model_type == "lmer" && !is.null(re_formula)) {
    if (is.null(re_group_cols) || !length(re_group_cols))
      stop("model_type = \"lmer\" with re_formula also requires ",
           "re_group_cols: the name(s) of the column(s) in `data` that ",
           "re_formula's grouping factor(s) refer to, e.g. ",
           "re_group_cols = \"block\" for re_formula = \"+ (1 | block)\".",
           call. = FALSE)
    missing_cols <- setdiff(re_group_cols, names(data))
    if (length(missing_cols))
      stop("re_group_cols not found in data: ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    for (col in re_group_cols) d[[col]] <- data[[col]]
  }

  d <- d[stats::complete.cases(d), ]
  if (nrow(d) < 6L)
    stop("Need at least 6 complete observations to fit a quadratic surface.")

  quad_rhs <- ".x + .y + I(.x^2) + I(.y^2) + .x:.y"

  model_obj <- switch(
    model_type,
    lm = {
      mod <- stats::lm(stats::as.formula(paste(".z ~", quad_rhs)), data = d)
      if (any(is.na(stats::coef(mod))))
        stop("Quadratic fit is rank-deficient (too few distinct design points?).")
      structure(list(fit = mod), class = "nutrigon_lm_fit")
    },
    glm = {
      mod <- stats::glm(stats::as.formula(paste(".z ~", quad_rhs)), data = d, ...)
      if (any(is.na(stats::coef(mod))))
        stop("Quadratic fit is rank-deficient (too few distinct design points?).")
      structure(list(fit = mod), class = "nutrigon_glm_fit")
    },
    lmer = {
      # lme4 and lmerTest are Suggests, not Imports: most users
      # never fit mixed models with nutrigon, so they shouldn't have to
      # install lme4's compiled-code dependency chain just to use
      # nutrigon_fit()/nutrigon_compare() with model_type = "lm"/"glm".
      # lmerTest::lmer() is used rather than lme4::lmer() so that a raw
      # fitted model (fit$model$fit) already carries Satterthwaite/KR
      # degrees of freedom for summary()/anova() if a user inspects it
      # directly -- fixef() and predict() behave identically either way,
      # since lmerTest's lmer() returns an object that IS-A merMod with
      # extra attributes, not a different class hierarchy.
      if (!requireNamespace("lme4", quietly = TRUE) ||
          !requireNamespace("lmerTest", quietly = TRUE)) {
        stop("model_type = \"lmer\" requires the lme4 and lmerTest ",
             "packages. Install them with:\n",
             "  install.packages(c(\"lme4\", \"lmerTest\"))",
             call. = FALSE)
      }
      if (is.null(re_formula))
        stop("model_type = \"lmer\" requires re_formula, e.g. ",
             "re_formula = \"+ (1 | block)\".")
      full_formula <- stats::as.formula(paste(".z ~", quad_rhs, re_formula))
      mod <- tryCatch(
        lmerTest::lmer(full_formula, data = d, ...),
        error = function(e) {
          msg <- conditionMessage(e)
          # A random-slope term written against the ORIGINAL column name
          # (e.g. "+ (1 + x | block)" when x_col = "x") fails inside
          # lme4's model.frame() with "object '<name>' not found", because
          # the fitting data frame only ever has .x/.y/.z (plus
          # re_group_cols), never the original x_col/y_col names. Detect
          # that specific failure and point the user at the fix, rather
          # than surfacing lme4's internal object-not-found error as-is.
          obj_not_found <- regmatches(
            msg, regexpr("(?<=object ')[^']+(?=' not found)", msg, perl = TRUE)
          )
          if (length(obj_not_found) == 1 &&
                obj_not_found %in% c(x_col, y_col, z_col)) {
            stop("re_formula referenced \"", obj_not_found, "\", but ",
                 "nutrigon_fit() internally renames your x/y/z columns ",
                 "to .x/.y/.z before fitting. For random SLOPES, write ",
                 "re_formula using \".x\"/\".y\" instead of your original ",
                 "column names, e.g. \"+ (1 + .x | block)\" rather than ",
                 "\"+ (1 + ", obj_not_found, " | block)\". Random ",
                 "intercepts alone (e.g. \"+ (1 | block)\") are ",
                 "unaffected by this.", call. = FALSE)
          }
          stop(e)
        }
      )
      if (any(is.na(lme4::fixef(mod))))
        stop("Quadratic fit is rank-deficient in its fixed effects ",
             "(too few distinct design points, or overparameterised ",
             "random-effects structure?).")
      structure(list(fit = mod), class = "nutrigon_lmer_fit")
    }
  )

  hull <- .hull_polygon(d, ".x", ".y")
  list(model_obj = model_obj, data = d, hull = hull, model_type = model_type)
}


# ============================================================================
# S3 dispatch: .get_link_coefs()
# ============================================================================

.get_link_coefs <- function(model_obj) UseMethod(".get_link_coefs")

.get_link_coefs.nutrigon_lm_fit <- function(model_obj) {
  stats::coef(model_obj$fit)
}

.get_link_coefs.nutrigon_glm_fit <- function(model_obj) {
  # coef() on a glm object already returns coefficients on the LINK scale
  # (this is the linear-predictor parameterisation glm() estimates
  # directly, regardless of family) -- no extra transformation needed.
  stats::coef(model_obj$fit)
}

.get_link_coefs.nutrigon_lmer_fit <- function(model_obj) {
  # lme4 availability is checked in .fit_surface() before an lmer fit.
  lme4::fixef(model_obj$fit)
}


# ============================================================================
# S3 dispatch: .predict_link()
# ============================================================================

.predict_link <- function(model_obj, newdata) UseMethod(".predict_link")

.predict_link.nutrigon_lm_fit <- function(model_obj, newdata) {
  stats::predict(model_obj$fit, newdata = newdata)
}

.predict_link.nutrigon_glm_fit <- function(model_obj, newdata) {
  # type = "link": predictions on the SAME scale the quadratic coefficients
  # were estimated on. This is deliberate -- see the design note at the top
  # of this file. Do not change to type = "response" without also changing
  # the closed-form algebra in .locate_optimum(), which assumes a plain
  # quadratic in (.x, .y).
  stats::predict(model_obj$fit, newdata = newdata, type = "link")
}

.predict_link.nutrigon_lmer_fit <- function(model_obj, newdata) {
  # re.form = NA: population-level prediction, ignoring random effects --
  # this is what "the nutritional optimum" should mean when random effects
  # represent blocks/individuals/replicates rather than a factor of
  # substantive biological interest. If you need optimum for a specific
  # group, fit a separate model per group instead of using re.form = NULL
  # here, since a group-conditional optimum is a different estimand from
  # what nutrigon_compare()'s bootstrap machinery assumes.
  stats::predict(model_obj$fit, newdata = newdata, re.form = NA)
}


# ============================================================================
# Convenience: does this model_obj's coefficient vector contain everything
# .locate_optimum() needs? (guards against e.g. a user-supplied re_formula
# that somehow drops a fixed-effect term via singular fit / accidental
# formula typo -- fails loudly here rather than producing NA arithmetic
# silently inside the stationary-point algebra)
# ============================================================================

.check_quadratic_coefs <- function(cf) {
  required <- c(".x", ".y", "I(.x^2)", "I(.y^2)", ".x:.y")
  missing <- setdiff(required, names(cf))
  if (length(missing) > 0)
    stop("Fitted model is missing expected quadratic term(s): ",
         paste(missing, collapse = ", "),
         ". Check that x_col/y_col are not aliased with another predictor ",
         "(e.g. via a random-effects term using the same variable).")
  invisible(TRUE)
}
