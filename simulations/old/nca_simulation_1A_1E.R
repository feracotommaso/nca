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

# Threshold systems for the ordinal-scaling manipulation.
# central  = middle response categories are more likely than extremes
# shifted  = same bell-shaped profile, but all thresholds are moved left/right
#            to induce gradual skew in observed responses.
make_thresholds <- function(n_cat = 5,
                            spacing = c("central", "shifted"),
                            shift_direction = c("none", "lower", "upper"),
                            shift_magnitude = 0.50,
                            central_sd = 1.2) {
  spacing <- match.arg(spacing)
  shift_direction <- match.arg(shift_direction)
  
  if (n_cat < 2) stop("n_cat must be at least 2.")
  if (central_sd <= 0) stop("central_sd must be > 0.")
  
  # Build a bell-shaped target distribution over response categories so that
  # middle categories are endorsed more often than extremes.
  centers <- seq(
    from = -(n_cat - 1) / 2,
    to = (n_cat - 1) / 2,
    length.out = n_cat
  )
  
  probs <- stats::dnorm(centers, mean = 0, sd = central_sd)
  probs <- probs / sum(probs)
  
  cum_probs <- cumsum(probs)[seq_len(n_cat - 1)]
  thresholds <- stats::qnorm(cum_probs)
  
  # Shift thresholds to change endorsement difficulty while preserving
  # the overall central-concentration shape.
  if (spacing == "shifted") {
    if (shift_direction == "lower") {
      thresholds <- thresholds - shift_magnitude
    } else if (shift_direction == "upper") {
      thresholds <- thresholds + shift_magnitude
    }
  }
  
  thresholds
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
                                 spacing = c("central", "shifted"),
                                 shift_direction = c("none", "lower", "upper"),
                                 shift_magnitude = 0.50,
                                 central_sd = 1.2,
                                 heterogeneous_loadings = FALSE,
                                 loading_sd = 0.05,
                                 seed = NULL) {
  spacing <- match.arg(spacing)
  shift_direction <- match.arg(shift_direction)
  
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
    spacing = spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd
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

construct_scores <- function(dat,
                             x_items,
                             y_items,
                             score_types = c("sum", "mean", "factor")) {
  score_types <- unique(score_types)

  need_sum_or_mean <- any(score_types %in% c("sum", "mean"))
  need_factor <- any(score_types %in% "factor")

  if (need_sum_or_mean) {
    dat <- dat %>%
      mutate(
        x_sum  = rowSums(across(all_of(x_items))),
        y_sum  = rowSums(across(all_of(y_items))),
        x_mean = rowMeans(across(all_of(x_items))),
        y_mean = rowMeans(across(all_of(y_items)))
      )
  }

  x_fit <- NULL
  y_fit <- NULL

  if (need_factor) {
    x_fs <- make_lavaan_factor_scores(dat, x_items, factor_name = "FX")
    y_fs <- make_lavaan_factor_scores(dat, y_items, factor_name = "FY")

    dat <- dat %>%
      mutate(
        x_factor = x_fs$scores,
        y_factor = y_fs$scores
      )

    x_fit <- x_fs$fit
    y_fit <- y_fs$fit
  }

  list(
    data = dat,
    fit_x = x_fit,
    fit_y = y_fit
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
  if (is.null(test_object)) {
    return(tibble(
      p_value = NA_real_,
      p_accuracy = NA_real_,
      threshold_value = NA_real_,
      observed_test_stat = NA_real_,
      significant = NA
    ))
  }

  p_value <- safe_numeric(test_object$p_value)
  p_accuracy <- safe_numeric(test_object$test.params$p_accuracy)
  threshold_value <- safe_numeric(test_object$threshold.value)
  observed_test_stat <- safe_numeric(test_object$observed)

  significant <- if (is.na(p_value)) {
    NA
  } else {
    p_value < 0.05
  }

  tibble(
    p_value = p_value,
    p_accuracy = p_accuracy,
    threshold_value = threshold_value,
    observed_test_stat = observed_test_stat,
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
    threshold_value = test_tbl$threshold_value,
    observed_test_stat = test_tbl$observed_test_stat,
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

  n_used <- nrow(tmp)
  x_unique <- dplyr::n_distinct(tmp$X)
  y_unique <- dplyr::n_distinct(tmp$Y)
  corr_xy <- if (n_used >= 2 && stats::sd(tmp$X) > 0 && stats::sd(tmp$Y) > 0) {
    stats::cor(tmp$X, tmp$Y, use = "pairwise.complete.obs")
  } else {
    NA_real_
  }

  # Guard against degenerate datasets before calling NCA.
  if (n_used < 10 || x_unique < 2 || y_unique < 2 || is.na(corr_xy)) {
    return(
      tibble(
        effect_size = NA_real_,
        ceiling_zone = NA_real_,
        scope = NA_real_,
        ceiling_accuracy = NA_real_,
        fit = NA_real_,
        slope = NA_real_,
        intercept = NA_real_,
        p_value = NA_real_,
        p_accuracy = NA_real_,
        significant = NA,
        bottleneck_table = list(NULL),
        nca_model = list(NULL),
        test_failed = TRUE,
        test_error = "Degenerate dataset for NCA (too few cases or unique values).",
        n_used = n_used,
        x_unique = x_unique,
        y_unique = y_unique,
        correlation = corr_xy,
        x_mean_obs = mean(tmp$X),
        y_mean_obs = mean(tmp$Y),
        x_sd_obs = stats::sd(tmp$X),
        y_sd_obs = stats::sd(tmp$Y),
        x_skew_obs = sample_skewness(tmp$X),
        y_skew_obs = sample_skewness(tmp$Y)
      )
    )
  }

  analysis_args <- list(
    data = tmp,
    x = "X",
    y = "Y",
    ceilings = ceiling,
    steps = steps,
    cutoff = cutoff
  )

  model <- tryCatch(
    do.call(NCA::nca_analysis, c(analysis_args, list(test.rep = test_rep))),
    error = function(e) structure(list(error = e), class = "nca_error")
  )

  test_failed <- inherits(model, "nca_error")
  test_error <- if (test_failed) conditionMessage(model$error) else NA_character_

  # Fallback: rerun without the permutation test so the simulation can continue
  # and still return effect-size information.
  if (test_failed) {
    model <- tryCatch(
      do.call(NCA::nca_analysis, c(analysis_args, list(test.rep = 0))),
      error = function(e) structure(list(error = e), class = "nca_error")
    )
  }

  if (inherits(model, "nca_error")) {
    return(
      tibble(
        effect_size = NA_real_,
        ceiling_zone = NA_real_,
        scope = NA_real_,
        ceiling_accuracy = NA_real_,
        fit = NA_real_,
        slope = NA_real_,
        intercept = NA_real_,
        p_value = NA_real_,
        p_accuracy = NA_real_,
        significant = NA,
        bottleneck_table = list(NULL),
        nca_model = list(NULL),
        test_failed = TRUE,
        test_error = conditionMessage(model$error),
        n_used = n_used,
        x_unique = x_unique,
        y_unique = y_unique,
        correlation = corr_xy,
        x_mean_obs = mean(tmp$X),
        y_mean_obs = mean(tmp$Y),
        x_sd_obs = stats::sd(tmp$X),
        y_sd_obs = stats::sd(tmp$Y),
        x_skew_obs = sample_skewness(tmp$X),
        y_skew_obs = sample_skewness(tmp$Y)
      )
    )
  }

  metrics <- extract_nca_metrics(model, condition_name = "X", ceiling = ceiling)

  metrics %>%
    mutate(
      test_failed = test_failed,
      test_error = test_error,
      n_used = n_used,
      x_unique = x_unique,
      y_unique = y_unique,
      correlation = corr_xy,
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
                                   spacing = c("central", "shifted"),
                                   shift_direction = c("none", "lower", "upper"),
                                   shift_magnitude = 0.50,
                                   central_sd = 1.2,
                                   x_skew = c("none", "positive", "negative"),
                                   y_skew = c("none", "positive", "negative"),
                                   restriction = c("full", "x_upper", "x_middle", "y_upper", "y_middle"),
                                   ceiling = "ce_fdh",
                                   test_rep = 1000,
                                   steps = 10,
                                   cutoff = 0,
                                   score_types = "sum",
                                   seed = NULL) {
  spacing <- match.arg(spacing)
  shift_direction <- match.arg(shift_direction)
  x_skew <- match.arg(x_skew)
  y_skew <- match.arg(y_skew)
  restriction <- match.arg(restriction)

  score_types <- unique(score_types)
  score_map <- list(
    sum = c("x_sum", "y_sum"),
    mean = c("x_mean", "y_mean"),
    factor = c("x_factor", "y_factor")
  )

  invalid_score_types <- setdiff(score_types, names(score_map))
  if (length(invalid_score_types) > 0) {
    stop(
      "Unknown score_types: ",
      paste(invalid_score_types, collapse = ", "),
      ". Valid options are: ",
      paste(names(score_map), collapse = ", ")
    )
  }

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
    spacing = spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd
  )

  y_items <- simulate_scale_items(
    eta = latent$eta_y,
    prefix = "y",
    n_items = n_items,
    loading = loading,
    n_cat = n_cat,
    spacing = spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd
  )

  dat <- bind_cols(latent, x_items, y_items)
  dat <- apply_range_restriction(dat, restriction = restriction)

  scored <- construct_scores(
    dat = dat,
    x_items = names(x_items),
    y_items = names(y_items),
    score_types = score_types
  )

  dat_scored <- scored$data

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
      spacing = spacing,
      shift_direction = shift_direction,
      shift_magnitude = shift_magnitude,
      central_sd = central_sd,
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
                                           spacing = "central",
                                           shift_direction = c("none"),
                                           shift_magnitude = 0.50,
                                           central_sd = 1.2,
                                           x_skew = "none",
                                           y_skew = "none",
                                           restriction = "full",
                                           ceiling = "ce_fdh",
                                           test_rep = 1000,
                                           steps = 10,
                                           cutoff = 0,
                                           score_types = "sum",
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
    spacing = spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd,
    x_skew = x_skew,
    y_skew = y_skew,
    restriction = restriction,
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff,
    score_types = score_types,
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
                                    spacing = "central",
                                    shift_direction = c("none"),
                                    shift_magnitude = 0.50,
                                    central_sd = 1.2,
                                    x_skew = "none",
                                    y_skew = "none",
                                    restriction = "full",
                                    ceiling = "ce_fdh",
                                    test_rep = 1000,
                                    steps = 10,
                                    cutoff = 0,
                                    score_types = "sum",
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
      spacing = spacing,
      shift_direction = shift_direction,
      shift_magnitude = shift_magnitude,
      central_sd = central_sd,
      x_skew = x_skew,
      y_skew = y_skew,
      restriction = restriction,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      score_types = score_types,
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
                                        spacing = c("central", "shifted"),
                                        shift_direction = c("lower"),
                                        shift_magnitude = 0.50,
                                        central_sd = 1.2,
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
                                        score_types = "sum",
                                        seed = 3234,
                                        progress = TRUE) {
  design <- tidyr::crossing(
    n_cat = n_cat_levels,
    spacing = spacing
  )

  purrr::pmap_dfr(design, function(n_cat, spacing) {
    run_replicated_condition(
      reps = reps,
      progress = progress,
      n = n,
      beta = beta,
      n_items = n_items,
      loading = loading,
      n_cat = n_cat,
      spacing = spacing,
      shift_direction = shift_direction,
      shift_magnitude = shift_magnitude,
      central_sd = central_sd,
      x_skew = x_skew,
      y_skew = y_skew,
      restriction = restriction,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      score_types = score_types,
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
                                 spacing = "central",
                                 shift_direction = c("none"),
                                 shift_magnitude = 0.50,
                                 central_sd = 1.2,
                                 restriction = "full",
                                 ceiling = "ce_fdh",
                                 test_rep = 1000,
                                 steps = 10,
                                 cutoff = 0,
                                 score_types = "sum",
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
      spacing = spacing,
      shift_direction = shift_direction,
      shift_magnitude = shift_magnitude,
      central_sd = central_sd,
      x_skew = x_skew,
      y_skew = y_skew,
      restriction = restriction,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      score_types = score_types,
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
                                          spacing = "central",
                                          shift_direction = c("none"),
                                          shift_magnitude = 0.50,
                                          central_sd = 1.2,
                                          x_skew = "none",
                                          y_skew = "none",
                                          ceiling = "ce_fdh",
                                          test_rep = 1000,
                                          steps = 10,
                                          cutoff = 0,
                                          score_types = "sum",
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
      spacing = spacing,
      shift_direction = shift_direction,
      shift_magnitude = shift_magnitude,
      central_sd = central_sd,
      x_skew = x_skew,
      y_skew = y_skew,
      restriction = rr,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      score_types = score_types,
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
    group_by(block, score_type, loading, n_cat, spacing, x_skew, y_skew, restriction) %>%
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

# N = 1000
# # ---- 1A. Score construction ----
# res_1A <- simulate_1A_score_construction(
#   reps = 500, 
#   n = N, 
#   test_rep = 0, 
#   score_types = c("sum","factor")
# )
# sum_1A <- summarise_simulation_results(res_1A)
# print(sum_1A)
# plot_effect_size_density(res_1A)
# 
# # ---- 1B. Reliability ----
# res_1B <- simulate_1B_reliability( # Lower rel > effect 12-9-6
#   reps = 500, 
#   n = N, 
#   test_rep = 0, 
#   loading_levels = c(0.50, 0.70, 0.85)
# )
# sum_1B <- summarise_simulation_results(res_1B)
# print(sum_1B)
# plot_mean_effects(sum_1B, "loading")
# 
# # ---- 1C. Ordinal scaling ----
# res_1C <- simulate_1C_ordinal_scaling( # 5items .099; .094(s) - 7items .152; .146(s) 
#   reps = 500, 
#   n = N, 
#   test_rep = 0, 
#   n_cat_levels = c(5, 7),
#   spacing = c("central", "shifted")
# )
# sum_1C <- summarise_simulation_results(res_1C)
# print(sum_1C)
# plot_effect_size_density(res_1C, facet_var = "spacing")
# 
# # ---- 1D. Skewness ----
# res_1D <- simulate_1D_skewness(
#   reps = 500, 
#   n = N, 
#   test_rep = 0)
# sum_1D <- summarise_simulation_results(res_1D)
# print(sum_1D)
# plot_effect_size_density(res_1D, facet_var = "skew_label")
# 
# # ---- 1E. Range restriction ----
# res_1E <- simulate_1E_range_restriction(
#   reps = 500, 
#   n = N, 
#   test_rep = 0, 
#   restrictions = c("full", "x_upper", "x_middle", "y_upper")
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

