# ================================================================
# Necessary Condition Analysis (NCA) with psychological data
# Simulation blocks 1A-1E
# ---------------------------------------------------------------
# Purpose:
#   Code for the one-factor-at-a-time simulation blocks proposed in
#   the methods paper on NCA and psychological measurement.
#
# Blocks:
#   1A = score construction
#   1B = reliability / loading strength
#   1C = ordinal scaling (number of categories + threshold spacing)
#   1D = skewness of latent variables
#   1E = range restriction / sample composition
#
# Notes:
#   - This script is written to be modular and easy to extend.
#   - It is conservative in its package use: tidyverse, lavaan, NCA.
#   - The extraction of some NCA test fields is deliberately defensive,
#     because the exact internal structure of model$tests can vary.
#   - The full NCA model object is preserved in list-columns whenever you
#     want to inspect bottleneck tables or raw test outputs later.
# ================================================================

# -----------------------------
# 0. Packages
# -----------------------------

required_packages <- c("dplyr", "purrr", "tibble", "tidyr", "ggplot2", "lavaan", "NCA")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(ggplot2)

# -----------------------------
# 1. Small utility helpers
# -----------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

safe_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  if (length(out) == 0) return(NA_real_)
  out[[1]]
}

standardize <- function(x) {
  as.numeric(scale(x))
}

sample_skewness <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA_real_)
  m <- mean(x)
  s <- stats::sd(x)
  if (is.na(s) || s == 0) return(NA_real_)
  mean(((x - m) / s)^3)
}

# -----------------------------
# 2. Latent-variable generation
# -----------------------------

# Rank-based monotone remapping.
# This preserves ordering while changing marginal shape.
induce_skew <- function(x,
                        kind = c("none", "positive", "negative"),
                        strength = 6) {
  kind <- match.arg(kind)

  if (kind == "none") {
    return(standardize(x))
  }

  u <- rank(x, ties.method = "average") / (length(x) + 1)
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)

  # Beta distributions provide convenient bounded skewed shapes.
  # positive  = right-skewed  (mass near lower values)
  # negative  = left-skewed   (mass near higher values)
  z <- if (kind == "positive") {
    stats::qbeta(u, shape1 = 2, shape2 = strength)
  } else {
    stats::qbeta(u, shape1 = strength, shape2 = 2)
  }

  standardize(z)
}

simulate_latent_pair <- function(n = 1000,
                                 beta = 0.40,
                                 x_skew = c("none", "positive", "negative"),
                                 y_skew = c("none", "positive", "negative"),
                                 residual_sd = NULL,
                                 seed = NULL) {
  x_skew <- match.arg(x_skew)
  y_skew <- match.arg(y_skew)

  if (!is.null(seed)) set.seed(seed)

  if (is.null(residual_sd)) {
    residual_sd <- sqrt(1 - beta^2)
  }

  eta_x_raw <- stats::rnorm(n)
  eta_x <- induce_skew(eta_x_raw, kind = x_skew)

  eps_raw <- stats::rnorm(n, mean = 0, sd = residual_sd)
  eps <- standardize(eps_raw) * residual_sd

  eta_y_linear <- beta * eta_x + eps
  eta_y <- induce_skew(eta_y_linear, kind = y_skew)

  tibble(
    id = seq_len(n),
    eta_x = eta_x,
    eta_y = eta_y
  )
}

# -----------------------------
# 3. Ordinal item generation
# -----------------------------

# Threshold spacing for the optional ordinal-scaling manipulation.
# symmetric = approximately equally spaced category probabilities
# skewed    = asymmetric thresholds, yielding uneven category widths
make_thresholds <- function(n_cat = 5,
                            spacing = c("symmetric", "skewed"),
                            skew_direction = c("lower", "upper")) {
  spacing <- match.arg(spacing)
  skew_direction <- match.arg(skew_direction)

  if (n_cat < 2) stop("n_cat must be at least 2.")

  probs_sym <- seq(1 / n_cat, (n_cat - 1) / n_cat, by = 1 / n_cat)

  if (spacing == "symmetric") {
    probs <- probs_sym
  } else {
    # Asymmetric cumulative probabilities using beta quantiles.
    # lower = more observations in lower response categories
    # upper = more observations in upper response categories
    u <- probs_sym
    probs <- if (skew_direction == "lower") {
      stats::qbeta(u, shape1 = 2, shape2 = 5)
    } else {
      stats::qbeta(u, shape1 = 5, shape2 = 2)
    }
  }

  stats::qnorm(probs)
}

latent_to_ordinal <- function(y_star, thresholds) {
  as.integer(cut(y_star,
                 breaks = c(-Inf, thresholds, Inf),
                 labels = FALSE,
                 right = TRUE))
}

simulate_scale_items <- function(eta,
                                 prefix,
                                 n_items = 6,
                                 loading = 0.70,
                                 n_cat = 5,
                                 threshold_spacing = c("symmetric", "skewed"),
                                 threshold_skew_direction = c("lower", "upper"),
                                 heterogeneous_loadings = FALSE,
                                 loading_sd = 0.05,
                                 seed = NULL) {
  threshold_spacing <- match.arg(threshold_spacing)
  threshold_skew_direction <- match.arg(threshold_skew_direction)

  if (!is.null(seed)) set.seed(seed)

  if (length(loading) == 1) {
    lambda <- rep(loading, n_items)
  } else {
    lambda <- loading
  }

  if (heterogeneous_loadings) {
    lambda <- pmin(pmax(stats::rnorm(n_items, mean = mean(lambda), sd = loading_sd), 0.20), 0.95)
  }

  thresholds <- make_thresholds(
    n_cat = n_cat,
    spacing = threshold_spacing,
    skew_direction = threshold_skew_direction
  )

  out <- vector("list", length = n_items)

  for (j in seq_len(n_items)) {
    err_sd <- sqrt(max(1 - lambda[j]^2, 1e-6))
    y_star <- lambda[j] * eta + stats::rnorm(length(eta), mean = 0, sd = err_sd)
    out[[j]] <- latent_to_ordinal(y_star, thresholds)
  }

  names(out) <- paste0(prefix, seq_len(n_items))
  out <- as_tibble(out)
  out
}

# -----------------------------
# 4. Score construction
# -----------------------------

make_lavaan_factor_scores <- function(dat, items, factor_name = "F") {
  model_syntax <- paste0(factor_name, " =~ ", paste(items, collapse = " + "))

  dat_ord <- dat
  for (nm in items) {
    dat_ord[[nm]] <- ordered(dat_ord[[nm]], levels = sort(unique(dat_ord[[nm]])))
  }

  fit <- tryCatch(
    lavaan::cfa(
      model = model_syntax,
      data = dat_ord,
      ordered = items,
      std.lv = TRUE,
      estimator = "WLSMV",
      parameterization = "theta"
    ),
    error = function(e) NULL,
    warning = function(w) invokeRestart("muffleWarning")
  )

  if (is.null(fit)) {
    return(list(scores = rep(NA_real_, nrow(dat)), fit = NULL))
  }

  scores <- tryCatch(
    as.numeric(lavaan::lavPredict(fit)[, 1]),
    error = function(e) rep(NA_real_, nrow(dat))
  )

  list(scores = scores, fit = fit)
}

construct_scores <- function(dat, x_items, y_items) {
  dat <- dat %>%
    mutate(
      x_sum  = rowSums(across(all_of(x_items))),
      y_sum  = rowSums(across(all_of(y_items))),
      x_mean = rowMeans(across(all_of(x_items))),
      y_mean = rowMeans(across(all_of(y_items)))
    )

  x_fs <- make_lavaan_factor_scores(dat, x_items, factor_name = "FX")
  y_fs <- make_lavaan_factor_scores(dat, y_items, factor_name = "FY")

  dat <- dat %>%
    mutate(
      x_factor = x_fs$scores,
      y_factor = y_fs$scores
    )

  list(
    data = dat,
    fit_x = x_fs$fit,
    fit_y = y_fs$fit
  )
}

# -----------------------------
# 5. Sample / range restriction
# -----------------------------

apply_range_restriction <- function(dat,
                                    restriction = c(
                                      "full",
                                      "x_upper",
                                      "x_middle",
                                      "y_upper",
                                      "y_middle"
                                    ),
                                    lower_q = 0.30,
                                    upper_q = 0.70) {
  restriction <- match.arg(restriction)

  if (restriction == "full") return(dat)

  if (restriction == "x_upper") {
    cut <- stats::quantile(dat$eta_x, probs = lower_q, na.rm = TRUE)
    return(filter(dat, eta_x > cut))
  }

  if (restriction == "x_middle") {
    cuts <- stats::quantile(dat$eta_x, probs = c(lower_q, upper_q), na.rm = TRUE)
    return(filter(dat, eta_x > cuts[[1]], eta_x < cuts[[2]]))
  }

  if (restriction == "y_upper") {
    cut <- stats::quantile(dat$eta_y, probs = lower_q, na.rm = TRUE)
    return(filter(dat, eta_y > cut))
  }

  if (restriction == "y_middle") {
    cuts <- stats::quantile(dat$eta_y, probs = c(lower_q, upper_q), na.rm = TRUE)
    return(filter(dat, eta_y > cuts[[1]], eta_y < cuts[[2]]))
  }
}

# -----------------------------
# 6. NCA helpers
# -----------------------------

extract_named_numeric <- function(x, candidates) {
  if (is.null(x)) return(NA_real_)

  nms <- names(x)
  if (is.null(nms)) return(NA_real_)

  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_real_)

  safe_numeric(x[[hit[[1]]]])
}

extract_nca_test_fields <- function(test_object) {
  # Internal structure may differ across versions; this tries common names.
  if (is.null(test_object)) {
    return(tibble(
      p_value = NA_real_,
      p_accuracy = NA_real_,
      significant = NA
    ))
  }

  p_value <- extract_named_numeric(
    test_object,
    c("p", "p.value", "p_value", "P-value", "p-value", "P")
  )

  p_accuracy <- extract_named_numeric(
    test_object,
    c("p.accuracy", "p_accuracy", "accuracy", "P accuracy", "p-accuracy")
  )

  significant_raw <- NULL
  if (!is.null(names(test_object))) {
    cand <- c("significant", "sig", "Significant")
    hit <- cand[cand %in% names(test_object)]
    if (length(hit) > 0) significant_raw <- test_object[[hit[[1]]]]
  }

  significant <- if (is.null(significant_raw)) {
    NA
  } else {
    as.logical(significant_raw[[1]])
  }

  tibble(
    p_value = p_value,
    p_accuracy = p_accuracy,
    significant = significant
  )
}

extract_nca_metrics <- function(model,
                                condition_name,
                                ceiling = "ce_fdh") {
  effect_size <- tryCatch(
    safe_numeric(NCA::nca_extract(model, x = condition_name, ceiling = ceiling, param = "Effect size")),
    error = function(e) NA_real_
  )

  ceiling_zone <- tryCatch(
    safe_numeric(NCA::nca_extract(model, x = condition_name, ceiling = ceiling, param = "Ceiling zone")),
    error = function(e) NA_real_
  )

  scope <- tryCatch(
    safe_numeric(NCA::nca_extract(model, x = condition_name, ceiling = ceiling, param = "Scope")),
    error = function(e) NA_real_
  )

  ceiling_accuracy <- tryCatch(
    safe_numeric(NCA::nca_extract(model, x = condition_name, ceiling = ceiling, param = "Ceiling accuracy")),
    error = function(e) NA_real_
  )

  fit <- tryCatch(
    safe_numeric(NCA::nca_extract(model, x = condition_name, ceiling = ceiling, param = "Fit")),
    error = function(e) NA_real_
  )

  slope <- tryCatch(
    safe_numeric(NCA::nca_extract(model, x = condition_name, ceiling = ceiling, param = "Slope")),
    error = function(e) NA_real_
  )

  intercept <- tryCatch(
    safe_numeric(NCA::nca_extract(model, x = condition_name, ceiling = ceiling, param = "Intercept")),
    error = function(e) NA_real_
  )

  test_object <- tryCatch(model$tests[[condition_name]][[ceiling]], error = function(e) NULL)
  test_tbl <- extract_nca_test_fields(test_object)

  bottleneck_table <- tryCatch(model$bottlenecks[[ceiling]], error = function(e) NULL)

  tibble(
    effect_size = effect_size,
    ceiling_zone = ceiling_zone,
    scope = scope,
    ceiling_accuracy = ceiling_accuracy,
    fit = fit,
    slope = slope,
    intercept = intercept,
    p_value = test_tbl$p_value,
    p_accuracy = test_tbl$p_accuracy,
    significant = test_tbl$significant,
    bottleneck_table = list(bottleneck_table),
    nca_model = list(model)
  )
}

run_one_nca <- function(dat,
                        x_var,
                        y_var,
                        ceiling = "ce_fdh",
                        test_rep = 1000,
                        steps = 10,
                        cutoff = 0) {
  tmp <- dat %>% select(all_of(c(x_var, y_var))) %>% tidyr::drop_na()
  names(tmp) <- c("X", "Y")

  model <- NCA::nca_analysis(
    data = tmp,
    x = "X",
    y = "Y",
    ceilings = ceiling,
    steps = steps,
    cutoff = cutoff,
    test.rep = test_rep
  )

  metrics <- extract_nca_metrics(model, condition_name = "X", ceiling = ceiling)

  metrics %>%
    mutate(
      n_used = nrow(tmp),
      correlation = stats::cor(tmp$X, tmp$Y, use = "pairwise.complete.obs"),
      x_mean_obs = mean(tmp$X),
      y_mean_obs = mean(tmp$Y),
      x_sd_obs = stats::sd(tmp$X),
      y_sd_obs = stats::sd(tmp$Y),
      x_skew_obs = sample_skewness(tmp$X),
      y_skew_obs = sample_skewness(tmp$Y)
    )
}

# -----------------------------
# 7. General replication engine
# -----------------------------

run_single_replication <- function(rep_id,
                                   n = 1000,
                                   beta = 0.40,
                                   n_items = 6,
                                   loading = 0.70,
                                   n_cat = 5,
                                   threshold_spacing = c("symmetric", "skewed"),
                                   threshold_skew_direction = c("lower", "upper"),
                                   x_skew = c("none", "positive", "negative"),
                                   y_skew = c("none", "positive", "negative"),
                                   restriction = c("full", "x_upper", "x_middle", "y_upper", "y_middle"),
                                   ceiling = "ce_fdh",
                                   test_rep = 1000,
                                   steps = 10,
                                   cutoff = 0,
                                   score_types = c("sum", "mean", "factor"),
                                   seed = NULL) {
  threshold_spacing <- match.arg(threshold_spacing)
  threshold_skew_direction <- match.arg(threshold_skew_direction)
  x_skew <- match.arg(x_skew)
  y_skew <- match.arg(y_skew)
  restriction <- match.arg(restriction)

  # Reproducibility by replication id if desired.
  if (!is.null(seed)) set.seed(seed + rep_id)

  latent <- simulate_latent_pair(
    n = n,
    beta = beta,
    x_skew = x_skew,
    y_skew = y_skew
  )

  x_items <- simulate_scale_items(
    eta = latent$eta_x,
    prefix = "x",
    n_items = n_items,
    loading = loading,
    n_cat = n_cat,
    threshold_spacing = threshold_spacing,
    threshold_skew_direction = threshold_skew_direction
  )

  y_items <- simulate_scale_items(
    eta = latent$eta_y,
    prefix = "y",
    n_items = n_items,
    loading = loading,
    n_cat = n_cat,
    threshold_spacing = threshold_spacing,
    threshold_skew_direction = threshold_skew_direction
  )

  dat <- bind_cols(latent, x_items, y_items)
  dat <- apply_range_restriction(dat, restriction = restriction)

  scored <- construct_scores(
    dat = dat,
    x_items = names(x_items),
    y_items = names(y_items)
  )

  dat_scored <- scored$data

  score_map <- list(
    sum = c("x_sum", "y_sum"),
    mean = c("x_mean", "y_mean"),
    factor = c("x_factor", "y_factor")
  )

  results <- purrr::map_dfr(score_types, function(st) {
    vars <- score_map[[st]]

    run_one_nca(
      dat = dat_scored,
      x_var = vars[[1]],
      y_var = vars[[2]],
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff
    ) %>%
      mutate(score_type = st)
  })

  results %>%
    mutate(
      rep_id = rep_id,
      n_input = n,
      beta = beta,
      n_items = n_items,
      loading = loading,
      n_cat = n_cat,
      threshold_spacing = threshold_spacing,
      threshold_skew_direction = threshold_skew_direction,
      x_skew = x_skew,
      y_skew = y_skew,
      restriction = restriction,
      ceiling = ceiling,
      test_rep = test_rep
    )
}

run_replicated_condition <- function(reps = 100,
                                     progress = TRUE,
                                     ...) {
  out <- vector("list", length = reps)

  for (r in seq_len(reps)) {
    if (progress && (r %% max(1, floor(reps / 10)) == 0 || r == 1 || r == reps)) {
      message("Running replication ", r, " / ", reps)
    }

    out[[r]] <- run_single_replication(rep_id = r, ...)
  }

  bind_rows(out)
}

# -----------------------------
# 8. Block-specific wrappers
# -----------------------------

# ------------------------------------------------
# 1A. Score construction
# ------------------------------------------------
# Same latent relation, same measurement model, same sample.
# The only comparison is across score types in the output.

simulate_1A_score_construction <- function(reps = 500,
                                           n = 1000,
                                           beta = 0.40,
                                           n_items = 6,
                                           loading = 0.70,
                                           n_cat = 5,
                                           threshold_spacing = "symmetric",
                                           threshold_skew_direction = "lower",
                                           x_skew = "none",
                                           y_skew = "none",
                                           restriction = "full",
                                           ceiling = "ce_fdh",
                                           test_rep = 1000,
                                           steps = 10,
                                           cutoff = 0,
                                           seed = 1234,
                                           progress = TRUE) {
  run_replicated_condition(
    reps = reps,
    progress = progress,
    n = n,
    beta = beta,
    n_items = n_items,
    loading = loading,
    n_cat = n_cat,
    threshold_spacing = threshold_spacing,
    threshold_skew_direction = threshold_skew_direction,
    x_skew = x_skew,
    y_skew = y_skew,
    restriction = restriction,
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff,
    seed = seed
  ) %>%
    mutate(block = "1A_score_construction")
}

# ------------------------------------------------
# 1B. Reliability / loading strength
# ------------------------------------------------

simulate_1B_reliability <- function(reps = 500,
                                    loading_levels = c(0.50, 0.70, 0.85),
                                    n = 1000,
                                    beta = 0.40,
                                    n_items = 6,
                                    n_cat = 5,
                                    threshold_spacing = "symmetric",
                                    threshold_skew_direction = "lower",
                                    x_skew = "none",
                                    y_skew = "none",
                                    restriction = "full",
                                    ceiling = "ce_fdh",
                                    test_rep = 1000,
                                    steps = 10,
                                    cutoff = 0,
                                    seed = 2234,
                                    progress = TRUE) {
  purrr::map_dfr(loading_levels, function(ld) {
    run_replicated_condition(
      reps = reps,
      progress = progress,
      n = n,
      beta = beta,
      n_items = n_items,
      loading = ld,
      n_cat = n_cat,
      threshold_spacing = threshold_spacing,
      threshold_skew_direction = threshold_skew_direction,
      x_skew = x_skew,
      y_skew = y_skew,
      restriction = restriction,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      seed = seed + round(ld * 100)
    )
  }) %>%
    mutate(
      block = "1B_reliability",
      reliability_condition = case_when(
        abs(loading - 0.50) < 1e-8 ~ "low",
        abs(loading - 0.70) < 1e-8 ~ "medium",
        abs(loading - 0.85) < 1e-8 ~ "high",
        TRUE ~ paste0("loading_", loading)
      )
    )
}

# ------------------------------------------------
# 1C. Ordinal scaling
# ------------------------------------------------
# Includes the optional threshold-spacing manipulation:
#   - symmetric thresholds
#   - skewed thresholds

simulate_1C_ordinal_scaling <- function(reps = 500,
                                        n_cat_levels = c(5, 7),
                                        threshold_spacings = c("symmetric", "skewed"),
                                        threshold_skew_direction = "lower",
                                        n = 1000,
                                        beta = 0.40,
                                        n_items = 6,
                                        loading = 0.70,
                                        x_skew = "none",
                                        y_skew = "none",
                                        restriction = "full",
                                        ceiling = "ce_fdh",
                                        test_rep = 1000,
                                        steps = 10,
                                        cutoff = 0,
                                        seed = 3234,
                                        progress = TRUE) {
  design <- tidyr::crossing(
    n_cat = n_cat_levels,
    threshold_spacing = threshold_spacings
  )

  purrr::pmap_dfr(design, function(n_cat, threshold_spacing) {
    run_replicated_condition(
      reps = reps,
      progress = progress,
      n = n,
      beta = beta,
      n_items = n_items,
      loading = loading,
      n_cat = n_cat,
      threshold_spacing = threshold_spacing,
      threshold_skew_direction = threshold_skew_direction,
      x_skew = x_skew,
      y_skew = y_skew,
      restriction = restriction,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      seed = seed + n_cat
    )
  }) %>%
    mutate(block = "1C_ordinal_scaling")
}

# ------------------------------------------------
# 1D. Skewness
# ------------------------------------------------

simulate_1D_skewness <- function(reps = 500,
                                 skew_conditions = NULL,
                                 n = 1000,
                                 beta = 0.40,
                                 n_items = 6,
                                 loading = 0.70,
                                 n_cat = 5,
                                 threshold_spacing = "symmetric",
                                 threshold_skew_direction = "lower",
                                 restriction = "full",
                                 ceiling = "ce_fdh",
                                 test_rep = 1000,
                                 steps = 10,
                                 cutoff = 0,
                                 seed = 4234,
                                 progress = TRUE) {
  if (is.null(skew_conditions)) {
    skew_conditions <- tibble(
      skew_label = c(
        "symmetric_X_symmetric_Y",
        "negative_X_symmetric_Y",
        "symmetric_X_positive_Y",
        "negative_X_positive_Y"
      ),
      x_skew = c("none", "negative", "none", "negative"),
      y_skew = c("none", "none", "positive", "positive")
    )
  }

  purrr::pmap_dfr(skew_conditions, function(skew_label, x_skew, y_skew) {
    run_replicated_condition(
      reps = reps,
      progress = progress,
      n = n,
      beta = beta,
      n_items = n_items,
      loading = loading,
      n_cat = n_cat,
      threshold_spacing = threshold_spacing,
      threshold_skew_direction = threshold_skew_direction,
      x_skew = x_skew,
      y_skew = y_skew,
      restriction = restriction,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      seed = seed + which(skew_conditions$skew_label == skew_label)
    ) %>%
      mutate(skew_label = skew_label)
  }) %>%
    mutate(block = "1D_skewness")
}

# ------------------------------------------------
# 1E. Range restriction / sample composition
# ------------------------------------------------

simulate_1E_range_restriction <- function(reps = 500,
                                          restrictions = c("full", "x_upper", "x_middle", "y_upper"),
                                          n = 1000,
                                          beta = 0.40,
                                          n_items = 6,
                                          loading = 0.70,
                                          n_cat = 5,
                                          threshold_spacing = "symmetric",
                                          threshold_skew_direction = "lower",
                                          x_skew = "none",
                                          y_skew = "none",
                                          ceiling = "ce_fdh",
                                          test_rep = 1000,
                                          steps = 10,
                                          cutoff = 0,
                                          seed = 5234,
                                          progress = TRUE) {
  purrr::map_dfr(restrictions, function(rr) {
    run_replicated_condition(
      reps = reps,
      progress = progress,
      n = n,
      beta = beta,
      n_items = n_items,
      loading = loading,
      n_cat = n_cat,
      threshold_spacing = threshold_spacing,
      threshold_skew_direction = threshold_skew_direction,
      x_skew = x_skew,
      y_skew = y_skew,
      restriction = rr,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      seed = seed + match(rr, restrictions)
    )
  }) %>%
    mutate(block = "1E_range_restriction")
}

# -----------------------------
# 9. Summaries and plotting
# -----------------------------

summarise_simulation_results <- function(results) {
  results %>%
    group_by(block, score_type, loading, n_cat, threshold_spacing, x_skew, y_skew, restriction) %>%
    summarise(
      n_replications = dplyr::n(),
      mean_effect_size = mean(effect_size, na.rm = TRUE),
      sd_effect_size = stats::sd(effect_size, na.rm = TRUE),
      mean_correlation = mean(correlation, na.rm = TRUE),
      mean_p_value = mean(p_value, na.rm = TRUE),
      prop_significant = mean(significant %in% TRUE, na.rm = TRUE),
      mean_ceiling_accuracy = mean(ceiling_accuracy, na.rm = TRUE),
      mean_x_skew_obs = mean(x_skew_obs, na.rm = TRUE),
      mean_y_skew_obs = mean(y_skew_obs, na.rm = TRUE),
      mean_n_used = mean(n_used, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_effect_size_density <- function(results, facet_var = NULL) {
  p <- ggplot(results, aes(x = effect_size, colour = score_type)) +
    geom_density(na.rm = TRUE) +
    labs(
      x = "NCA effect size",
      y = "Density",
      colour = "Score type"
    ) +
    theme_minimal(base_size = 12)

  if (!is.null(facet_var)) {
    p <- p + facet_wrap(stats::as.formula(paste("~", facet_var)), scales = "free_y")
  }

  p
}

plot_mean_effects <- function(summary_tbl, x_var) {
  ggplot(summary_tbl, aes_string(x = x_var, y = "mean_effect_size", colour = "score_type", group = "score_type")) +
    geom_line() +
    geom_point() +
    labs(
      x = x_var,
      y = "Mean NCA effect size",
      colour = "Score type"
    ) +
    theme_minimal(base_size = 12)
}

# -----------------------------
# 10. Example runs (commented)
# -----------------------------

# ---- 1A. Score construction ----
res_1A <- simulate_1A_score_construction(reps = 10, n = 1000, test_rep = 1000)
sum_1A <- summarise_simulation_results(res_1A)
print(sum_1A)
plot_effect_size_density(res_1A)

# ---- 1B. Reliability ----
# res_1B <- simulate_1B_reliability(reps = 100, loading_levels = c(0.50, 0.70, 0.85), test_rep = 500)
# sum_1B <- summarise_simulation_results(res_1B)
# print(sum_1B)
# plot_mean_effects(sum_1B, "loading")

# ---- 1C. Ordinal scaling ----
# res_1C <- simulate_1C_ordinal_scaling(
#   reps = 100,
#   n_cat_levels = c(5, 7),
#   threshold_spacings = c("symmetric", "skewed"),
#   test_rep = 500
# )
# sum_1C <- summarise_simulation_results(res_1C)
# print(sum_1C)
# plot_effect_size_density(res_1C, facet_var = "threshold_spacing")

# ---- 1D. Skewness ----
# res_1D <- simulate_1D_skewness(reps = 100, test_rep = 500)
# sum_1D <- summarise_simulation_results(res_1D)
# print(sum_1D)
# plot_effect_size_density(res_1D, facet_var = "skew_label")

# ---- 1E. Range restriction ----
# res_1E <- simulate_1E_range_restriction(
#   reps = 100,
#   restrictions = c("full", "x_upper", "x_middle", "y_upper"),
#   test_rep = 500
# )
# sum_1E <- summarise_simulation_results(res_1E)
# print(sum_1E)
# plot_effect_size_density(res_1E, facet_var = "restriction")

# -----------------------------
# 11. Optional saving helpers
# -----------------------------

save_simulation_bundle <- function(results, file_stub) {
  saveRDS(results, paste0(file_stub, ".rds"))
  write.csv(results %>% select(-bottleneck_table, -nca_model), paste0(file_stub, ".csv"), row.names = FALSE)
}

