# ================================================================ #
# Necessary Condition Analysis (NCA) with psychological data -----
# Unified paired/common-source simulation script
# --------------------------------------------------------------- #
# Purpose:
#   Self-contained simulation framework for the methods paper on
#   NCA and psychological measurement.
#
# Core design:
#   - In each replication, generate ONE latent paired dataset from a
#     common random source.
#   - Estimate latent-reference NCA on the latent variables.
#   - Derive multiple observed-score conditions from that same latent
#     source, varying one feature at a time.
#   - Compare observed NCA estimates against latent-reference NCA.
#
# Outer factors:
#   - sample size (n)
#   - latent association strength (beta)
#
# Inner blocks:
#   1A = score construction      (sum vs factor)
#   1B = reliability/loading     (0.50, 0.70, 0.85)
#   1C = response categories     (4, 5, 7)
#   1D = threshold shift         (lower/upper x 0.50 0.70 1.00)
#   1E = range restriction       (x_lower, x_upper, y_lower, y_upper)
#
# Notes:
#   - Block 1A is the only block that varies score type.
#   - All other blocks use sum scores by default.
#   - The script stores distribution descriptors of the final analyzed
#     X and Y scores at each replication: mean, SD, skewness, kurtosis.
#   - Example runs are commented out so the file is source-safe.
# ================================================================ #

# ----------------------------- #
# 0. Packages -------------------
# ----------------------------- #

required_packages <- c(
  "dplyr", "purrr", "tibble", "tidyr", "ggplot2",
  "lavaan", "NCA", "furrr", "future"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

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
library(furrr)
library(future)

# ----------------------------- #
# 1. Small utility helpers ------
# ----------------------------- #

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

sample_kurtosis <- function(x, na.rm = TRUE, excess = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n < 4) return(NA_real_)
  m <- mean(x)
  s <- stats::sd(x)
  if (is.na(s) || s == 0) return(NA_real_)
  k <- mean(((x - m) / s)^4)
  if (excess) k - 3 else k
}

# ----------------------------- #
# 2. Design validators ----------
# ----------------------------- #

.validate_score_types <- function(score_types) {
  score_types <- unique(score_types)
  allowed <- c("sum", "factor")
  bad <- setdiff(score_types, allowed)
  if (length(bad) > 0) {
    stop(
      "Unsupported score_types: ",
      paste(bad, collapse = ", "),
      ". Valid options are: ",
      paste(allowed, collapse = ", ")
    )
  }
  score_types
}

.validate_restrictions <- function(restrictions) {
  restrictions <- unique(restrictions)
  allowed <- c("full", "x_lower", "x_upper", "y_lower", "y_upper")
  bad <- setdiff(restrictions, allowed)
  if (length(bad) > 0) {
    stop(
      "Unsupported restrictions: ",
      paste(bad, collapse = ", "),
      ". Valid options are: ",
      paste(allowed, collapse = ", ")
    )
  }
  restrictions
}

# -------------------------------- #
# 3. Latent-variable generation ----
# -------------------------------- #

simulate_latent_from_common_source <- function(n = 1000,
                                               beta = 0.40,
                                               residual_sd = NULL,
                                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (is.null(residual_sd)) {
    residual_sd <- sqrt(1 - beta^2)
  }

  z_x <- stats::rnorm(n)
  z_eps <- stats::rnorm(n)

  eta_x <- standardize(z_x)
  eps <- standardize(z_eps) * residual_sd
  eta_y <- beta * eta_x + eps

  tibble(
    id = seq_len(n),
    eta_x = eta_x,
    eta_y = eta_y
  )
}

# ----------------------------- #
# 4. Ordinal item generation ----
# ----------------------------- #

make_thresholds <- function(n_cat = 5,
                            spacing = c("central", "shifted"),
                            shift_direction = c("none", "lower", "upper"),
                            shift_magnitude = 0.50,
                            central_sd = 1.2) {
  spacing <- match.arg(spacing)
  shift_direction <- match.arg(shift_direction)

  if (n_cat < 2) stop("n_cat must be at least 2.")
  if (central_sd <= 0) stop("central_sd must be > 0.")
  if (shift_magnitude < 0) stop("shift_magnitude must be >= 0.")

  centers <- seq(
    from = -(n_cat - 1) / 2,
    to = (n_cat - 1) / 2,
    length.out = n_cat
  )

  probs <- stats::dnorm(centers, mean = 0, sd = central_sd)
  probs <- probs / sum(probs)

  cum_probs <- cumsum(probs)[seq_len(n_cat - 1)]
  thresholds <- stats::qnorm(cum_probs)

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
  as.integer(cut(
    y_star,
    breaks = c(-Inf, thresholds, Inf),
    labels = FALSE,
    right = TRUE
  ))
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
    lambda <- pmin(
      pmax(stats::rnorm(n_items, mean = mean(lambda), sd = loading_sd), 0.20),
      0.95
    )
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
  as_tibble(out)
}

# ----------------------------- #
# 5. Score construction ---------
# ----------------------------- #

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
                             score_types = c("sum", "factor")) {
  score_types <- .validate_score_types(score_types)

  dat <- dat %>%
    mutate(
      x_sum = rowSums(across(all_of(x_items))),
      y_sum = rowSums(across(all_of(y_items)))
    )

  fit_x <- NULL
  fit_y <- NULL

  if ("factor" %in% score_types) {
    x_fs <- make_lavaan_factor_scores(dat, x_items, factor_name = "FX")
    y_fs <- make_lavaan_factor_scores(dat, y_items, factor_name = "FY")

    dat <- dat %>%
      mutate(
        x_factor = x_fs$scores,
        y_factor = y_fs$scores
      )

    fit_x <- x_fs$fit
    fit_y <- y_fs$fit
  }

  list(
    data = dat,
    fit_x = fit_x,
    fit_y = fit_y
  )
}

# -------------------------------------------- #
# 6. Range restriction / sample composition ----
# -------------------------------------------- #

apply_range_restriction <- function(dat,
                                    restriction = c("full", "x_lower", "x_upper", "y_lower", "y_upper"),
                                    lower_q = 0.30,
                                    upper_q = 0.70) {
  restriction <- match.arg(restriction)

  if (restriction == "full") return(dat)

  if (restriction == "x_upper") {
    cut <- stats::quantile(dat$eta_x, probs = lower_q, na.rm = TRUE)
    return(filter(dat, eta_x > cut))
  }

  if (restriction == "x_lower") {
    cut <- stats::quantile(dat$eta_x, probs = upper_q, na.rm = TRUE)
    return(filter(dat, eta_x < cut))
  }

  if (restriction == "y_upper") {
    cut <- stats::quantile(dat$eta_y, probs = lower_q, na.rm = TRUE)
    return(filter(dat, eta_y > cut))
  }

  if (restriction == "y_lower") {
    cut <- stats::quantile(dat$eta_y, probs = upper_q, na.rm = TRUE)
    return(filter(dat, eta_y < cut))
  }
}

# ----------------------------- #
# 7. NCA helpers ----------------
# ----------------------------- #

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

  significant <- if (is.na(p_value)) NA else p_value < 0.05

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

  x_mean_obs <- if (n_used > 0) mean(tmp$X) else NA_real_
  y_mean_obs <- if (n_used > 0) mean(tmp$Y) else NA_real_
  x_sd_obs <- if (n_used > 1) stats::sd(tmp$X) else NA_real_
  y_sd_obs <- if (n_used > 1) stats::sd(tmp$Y) else NA_real_
  x_skew_obs <- sample_skewness(tmp$X)
  y_skew_obs <- sample_skewness(tmp$Y)
  x_kurtosis_obs <- sample_kurtosis(tmp$X)
  y_kurtosis_obs <- sample_kurtosis(tmp$Y)

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
        threshold_value = NA_real_,
        observed_test_stat = NA_real_,
        significant = NA,
        bottleneck_table = list(NULL),
        nca_model = list(NULL),
        test_failed = TRUE,
        test_error = "Degenerate dataset for NCA (too few cases or unique values).",
        n_used = n_used,
        x_unique = x_unique,
        y_unique = y_unique,
        correlation = corr_xy,
        x_mean_obs = x_mean_obs,
        y_mean_obs = y_mean_obs,
        x_sd_obs = x_sd_obs,
        y_sd_obs = y_sd_obs,
        x_skew_obs = x_skew_obs,
        y_skew_obs = y_skew_obs,
        x_kurtosis_obs = x_kurtosis_obs,
        y_kurtosis_obs = y_kurtosis_obs
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
        threshold_value = NA_real_,
        observed_test_stat = NA_real_,
        significant = NA,
        bottleneck_table = list(NULL),
        nca_model = list(NULL),
        test_failed = TRUE,
        test_error = conditionMessage(model$error),
        n_used = n_used,
        x_unique = x_unique,
        y_unique = y_unique,
        correlation = corr_xy,
        x_mean_obs = x_mean_obs,
        y_mean_obs = y_mean_obs,
        x_sd_obs = x_sd_obs,
        y_sd_obs = y_sd_obs,
        x_skew_obs = x_skew_obs,
        y_skew_obs = y_skew_obs,
        x_kurtosis_obs = x_kurtosis_obs,
        y_kurtosis_obs = y_kurtosis_obs
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
      x_mean_obs = x_mean_obs,
      y_mean_obs = y_mean_obs,
      x_sd_obs = x_sd_obs,
      y_sd_obs = y_sd_obs,
      x_skew_obs = x_skew_obs,
      y_skew_obs = y_skew_obs,
      x_kurtosis_obs = x_kurtosis_obs,
      y_kurtosis_obs = y_kurtosis_obs
    )
}

get_reference_metrics <- function(dat,
                                  x_var = "eta_x",
                                  y_var = "eta_y",
                                  ceiling = "ce_fdh",
                                  test_rep = 0,
                                  steps = 10,
                                  cutoff = 0) {
  ref <- run_one_nca(
    dat = dat,
    x_var = x_var,
    y_var = y_var,
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff
  )

  list(
    effect_size = ref$effect_size[[1]],
    correlation = ref$correlation[[1]],
    n_used = ref$n_used[[1]],
    x_mean_obs = ref$x_mean_obs[[1]],
    y_mean_obs = ref$y_mean_obs[[1]],
    x_sd_obs = ref$x_sd_obs[[1]],
    y_sd_obs = ref$y_sd_obs[[1]],
    x_skew_obs = ref$x_skew_obs[[1]],
    y_skew_obs = ref$y_skew_obs[[1]],
    x_kurtosis_obs = ref$x_kurtosis_obs[[1]],
    y_kurtosis_obs = ref$y_kurtosis_obs[[1]],
    nca_model = ref$nca_model[[1]],
    bottleneck_table = ref$bottleneck_table[[1]]
  )
}

# ------------------------------- #
# 8. Observed-condition engine ----
# ------------------------------- #

run_observed_condition <- function(latent_full,
                                   latent_selected,
                                   block,
                                   condition_label,
                                   condition_value = NA_character_,
                                   score_types = "sum",
                                   n_items = 6,
                                   loading = 0.70,
                                   n_cat = 5,
                                   spacing = "central",
                                   shift_direction = "none",
                                   shift_magnitude = 0,
                                   central_sd = 1.2,
                                   ceiling = "ce_fdh",
                                   test_rep = 0,
                                   steps = 10,
                                   cutoff = 0,
                                   base_reference,
                                   condition_full_reference,
                                   condition_selected_reference) {
  score_types <- .validate_score_types(score_types)

  x_items <- simulate_scale_items(
    eta = latent_selected$eta_x,
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
    eta = latent_selected$eta_y,
    prefix = "y",
    n_items = n_items,
    loading = loading,
    n_cat = n_cat,
    spacing = spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd
  )

  dat_obs <- dplyr::bind_cols(latent_selected, x_items, y_items)
  scored <- construct_scores(
    dat = dat_obs,
    x_items = names(x_items),
    y_items = names(y_items),
    score_types = score_types
  )
  dat_scored <- scored$data

  score_map <- list(
    sum = c("x_sum", "y_sum"),
    factor = c("x_factor", "y_factor")
  )

  purrr::map_dfr(score_types, function(st) {
    vars <- score_map[[st]]

    out <- run_one_nca(
      dat = dat_scored,
      x_var = vars[[1]],
      y_var = vars[[2]],
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff
    )

    out %>%
      dplyr::mutate(
        block = block,
        condition_label = condition_label,
        condition_value = condition_value,
        score_type = st,
        latent_effect_base_full = base_reference$effect_size,
        latent_correlation_base_full = base_reference$correlation,
        latent_effect_condition_full = condition_full_reference$effect_size,
        latent_correlation_condition_full = condition_full_reference$correlation,
        latent_effect_condition_selected = condition_selected_reference$effect_size,
        latent_correlation_condition_selected = condition_selected_reference$correlation,
        delta_effect_from_base_full = effect_size - latent_effect_base_full,
        delta_effect_from_condition_full = effect_size - latent_effect_condition_full,
        delta_effect_from_condition_selected = effect_size - latent_effect_condition_selected,
        selection_effect_on_latent_d = latent_effect_condition_selected - latent_effect_condition_full,
        latent_n_full = nrow(latent_full),
        latent_n_selected = nrow(latent_selected),
        n_items = n_items,
        loading = loading,
        n_cat = n_cat,
        spacing = spacing,
        shift_direction = shift_direction,
        shift_magnitude = shift_magnitude,
        central_sd = central_sd
      )
  })
}

# ------------------------------------ #
# 9. One unified paired replication ----
# ------------------------------------ #

run_unified_paired_replication <- function(rep_id,
                                           n = 1000,
                                           beta = 0.40,
                                           n_items = 6,
                                           baseline_loading = 0.70,
                                           baseline_n_cat = 5,
                                           baseline_spacing = "central",
                                           baseline_shift_direction = "none",
                                           baseline_shift_magnitude = 0,
                                           baseline_central_sd = 1.2,
                                           baseline_restriction = "full",
                                           score_types_1A = c("sum", "factor"),
                                           loading_levels_1B = c(0.50, 0.70, 0.85),
                                           n_cat_levels_1C = c(4, 5, 7),
                                           shift_directions_1D = c("lower", "upper"),
                                           shift_magnitudes_1D = c(0.50, 0.70, 1.00),
                                           restrictions_1E = c("x_lower", "x_upper", "y_lower", "y_upper"),
                                           restriction_lower_q = 0.30,
                                           restriction_upper_q = 0.70,
                                           ceiling = "ce_fdh",
                                           test_rep = 0,
                                           steps = 10,
                                           cutoff = 0,
                                           seed = NULL) {
  score_types_1A <- .validate_score_types(score_types_1A)
  restrictions_1E <- .validate_restrictions(restrictions_1E)

  if (!is.null(seed)) set.seed(seed + rep_id)

  latent_base_full <- simulate_latent_from_common_source(
    n = n,
    beta = beta
  )

  latent_base_selected <- apply_range_restriction(
    latent_base_full,
    restriction = baseline_restriction,
    lower_q = restriction_lower_q,
    upper_q = restriction_upper_q
  )

  ref_base_full <- get_reference_metrics(
    dat = latent_base_full,
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff
  )

  ref_base_selected <- get_reference_metrics(
    dat = latent_base_selected,
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff
  )

  # 1A. Score construction: same observed dataset, different scores.
  res_1A <- run_observed_condition(
    latent_full = latent_base_full,
    latent_selected = latent_base_selected,
    block = "1A_score_construction",
    condition_label = "baseline_scores",
    condition_value = "baseline",
    score_types = score_types_1A,
    n_items = n_items,
    loading = baseline_loading,
    n_cat = baseline_n_cat,
    spacing = baseline_spacing,
    shift_direction = baseline_shift_direction,
    shift_magnitude = baseline_shift_magnitude,
    central_sd = baseline_central_sd,
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff,
    base_reference = ref_base_full,
    condition_full_reference = ref_base_full,
    condition_selected_reference = ref_base_selected
  )

  # 1B. Reliability/loading: same latent base, sum score only.
  res_1B <- purrr::map_dfr(loading_levels_1B, function(ld) {
    run_observed_condition(
      latent_full = latent_base_full,
      latent_selected = latent_base_selected,
      block = "1B_reliability",
      condition_label = paste0("loading_", formatC(ld, format = "f", digits = 2)),
      condition_value = as.character(ld),
      score_types = "sum",
      n_items = n_items,
      loading = ld,
      n_cat = baseline_n_cat,
      spacing = baseline_spacing,
      shift_direction = baseline_shift_direction,
      shift_magnitude = baseline_shift_magnitude,
      central_sd = baseline_central_sd,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      base_reference = ref_base_full,
      condition_full_reference = ref_base_full,
      condition_selected_reference = ref_base_selected
    )
  })

  # 1C. Number of response categories: same latent base, sum score only.
  res_1C <- purrr::map_dfr(n_cat_levels_1C, function(n_cat_current) {
    run_observed_condition(
      latent_full = latent_base_full,
      latent_selected = latent_base_selected,
      block = "1C_response_categories",
      condition_label = paste0("n_cat_", n_cat_current),
      condition_value = as.character(n_cat_current),
      score_types = "sum",
      n_items = n_items,
      loading = baseline_loading,
      n_cat = n_cat_current,
      spacing = "central",
      shift_direction = "none",
      shift_magnitude = 0,
      central_sd = baseline_central_sd,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      base_reference = ref_base_full,
      condition_full_reference = ref_base_full,
      condition_selected_reference = ref_base_selected
    )
  })

  # 1D. Threshold shift: fixed category count, sum score only.
  design_1D <- tidyr::crossing(
    shift_direction = shift_directions_1D,
    shift_magnitude = shift_magnitudes_1D
  )

  res_1D <- purrr::pmap_dfr(design_1D, function(shift_direction, shift_magnitude) {
    label <- paste0(
      shift_direction,
      "_",
      formatC(shift_magnitude, format = "f", digits = 2)
    )

    run_observed_condition(
      latent_full = latent_base_full,
      latent_selected = latent_base_selected,
      block = "1D_threshold_shift",
      condition_label = label,
      condition_value = label,
      score_types = "sum",
      n_items = n_items,
      loading = baseline_loading,
      n_cat = baseline_n_cat,
      spacing = "shifted",
      shift_direction = shift_direction,
      shift_magnitude = shift_magnitude,
      central_sd = baseline_central_sd,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      base_reference = ref_base_full,
      condition_full_reference = ref_base_full,
      condition_selected_reference = ref_base_selected
    )
  })

  # 1E. Range restriction: same latent base; compare to selected latent reference.
  res_1E <- purrr::map_dfr(restrictions_1E, function(rr) {
    latent_rr_selected <- apply_range_restriction(
      latent_base_full,
      restriction = rr,
      lower_q = restriction_lower_q,
      upper_q = restriction_upper_q
    )

    ref_rr_selected <- get_reference_metrics(
      dat = latent_rr_selected,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff
    )

    run_observed_condition(
      latent_full = latent_base_full,
      latent_selected = latent_rr_selected,
      block = "1E_range_restriction",
      condition_label = rr,
      condition_value = rr,
      score_types = "sum",
      n_items = n_items,
      loading = baseline_loading,
      n_cat = baseline_n_cat,
      spacing = baseline_spacing,
      shift_direction = baseline_shift_direction,
      shift_magnitude = baseline_shift_magnitude,
      central_sd = baseline_central_sd,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      base_reference = ref_base_full,
      condition_full_reference = ref_base_full,
      condition_selected_reference = ref_rr_selected
    )
  })

  bind_rows(res_1A, res_1B, res_1C, res_1D, res_1E) %>%
    mutate(
      rep_id = rep_id,
      n_input = n,
      beta = beta,
      baseline_loading = baseline_loading,
      baseline_n_cat = baseline_n_cat,
      baseline_spacing = baseline_spacing,
      baseline_shift_direction = baseline_shift_direction,
      baseline_shift_magnitude = baseline_shift_magnitude,
      baseline_central_sd = baseline_central_sd,
      baseline_restriction = baseline_restriction,
      restriction_lower_q = restriction_lower_q,
      restriction_upper_q = restriction_upper_q,
      ceiling = ceiling,
      test_rep = test_rep
    )
}

# ----------------------------- #
# 10. Replicated runs -----------
# ----------------------------- #

simulate_unified_paired <- function(reps = 100,
                                    progress = TRUE,
                                    ...) {
  out <- vector("list", length = reps)

  for (r in seq_len(reps)) {
    if (progress && (r %% max(1, floor(reps / 10)) == 0 || r == 1 || r == reps)) {
      message("Running unified paired replication ", r, " / ", reps)
    }

    out[[r]] <- run_unified_paired_replication(rep_id = r, ...)
  }

  bind_rows(out)
}

simulate_unified_paired <- function(reps = 100,
                                    progress = TRUE,
                                    seed = NULL,
                                    parallel = FALSE,
                                    ...) {
  rep_ids <- seq_len(reps)
  
  if (!parallel) {
    out <- vector("list", length = reps)
    
    for (r in rep_ids) {
      if (progress && (r %% max(1, floor(reps / 10)) == 0 || r == 1 || r == reps)) {
        message("Running unified paired replication ", r, " / ", reps)
      }
      
      out[[r]] <- run_unified_paired_replication(
        rep_id = r,
        seed = seed,
        ...
      )
    }
    
    return(bind_rows(out))
  }
  
  furrr::future_map_dfr(
    rep_ids,
    function(r) {
      run_unified_paired_replication(
        rep_id = r,
        seed = seed,
        ...
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
}

simulate_unified_paired_across_design <- function(n_values = c(250, 500, 1000),
                                                  beta_values = c(0.20, 0.40, 0.60),
                                                  reps = 100,
                                                  progress = TRUE,
                                                  seed = NULL,
                                                  parallel_inner = FALSE,
                                                  ...) {
  design <- tidyr::crossing(
    sample_size = n_values,
    beta = beta_values
  ) %>%
    mutate(design_id = row_number())

  purrr::pmap_dfr(design, function(sample_size, beta, design_id) {
    if (progress) {
      message(
        "Running design cell ", design_id, " / ", nrow(design),
        ": n = ", sample_size,
        ", beta = ", beta
      )
    }

    cell_seed <- if (is.null(seed)) NULL else seed + design_id * 100000

    simulate_unified_paired(
      reps = reps,
      progress = progress,
      n = sample_size,
      beta = beta,
      seed = cell_seed,
      parallel = parallel_inner,
      ...
    ) %>%
      mutate(
        sample_size = sample_size,
        design_beta = beta,
        design_id = design_id
      )
  })
}

# ----------------------------- #
# 11. Summaries -----------------
# ----------------------------- #

summarise_unified_paired_results <- function(results) {
  results %>%
    group_by(block, condition_label, score_type, n_input, beta) %>%
    summarise(
      n_replications = dplyr::n(),
      mean_effect_size = mean(effect_size, na.rm = TRUE),
      sd_effect_size = stats::sd(effect_size, na.rm = TRUE),
      mean_correlation = mean(correlation, na.rm = TRUE),
      mean_latent_effect_base_full = mean(latent_effect_base_full, na.rm = TRUE),
      mean_latent_effect_condition_full = mean(latent_effect_condition_full, na.rm = TRUE),
      mean_latent_effect_condition_selected = mean(latent_effect_condition_selected, na.rm = TRUE),
      mean_delta_from_base_full = mean(delta_effect_from_base_full, na.rm = TRUE),
      mean_delta_from_condition_full = mean(delta_effect_from_condition_full, na.rm = TRUE),
      mean_delta_from_condition_selected = mean(delta_effect_from_condition_selected, na.rm = TRUE),
      mean_selection_effect_on_latent_d = mean(selection_effect_on_latent_d, na.rm = TRUE),
      mean_x_mean_obs = mean(x_mean_obs, na.rm = TRUE),
      mean_y_mean_obs = mean(y_mean_obs, na.rm = TRUE),
      mean_x_sd_obs = mean(x_sd_obs, na.rm = TRUE),
      mean_y_sd_obs = mean(y_sd_obs, na.rm = TRUE),
      mean_x_skew_obs = mean(x_skew_obs, na.rm = TRUE),
      mean_y_skew_obs = mean(y_skew_obs, na.rm = TRUE),
      mean_x_kurtosis_obs = mean(x_kurtosis_obs, na.rm = TRUE),
      mean_y_kurtosis_obs = mean(y_kurtosis_obs, na.rm = TRUE),
      mean_n_used = mean(n_used, na.rm = TRUE),
      prop_significant = mean(significant %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(n_input, beta, block, condition_label, score_type)
}

summarise_unified_paired_by_block <- function(results) {
  results %>%
    group_by(n_input, beta, block, condition_label, score_type) %>%
    summarise(
      n_replications = dplyr::n(),
      mean_effect_size = mean(effect_size, na.rm = TRUE),
      sd_effect_size = stats::sd(effect_size, na.rm = TRUE),
      mean_delta_from_base_full = mean(delta_effect_from_base_full, na.rm = TRUE),
      mean_delta_from_condition_selected = mean(delta_effect_from_condition_selected, na.rm = TRUE),
      mean_x_skew_obs = mean(x_skew_obs, na.rm = TRUE),
      mean_y_skew_obs = mean(y_skew_obs, na.rm = TRUE),
      mean_x_kurtosis_obs = mean(x_kurtosis_obs, na.rm = TRUE),
      mean_y_kurtosis_obs = mean(y_kurtosis_obs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(n_input, beta, block, condition_label, score_type)
}

# ----------------------------- #
# 12. Plot helpers --------------
# ----------------------------- #

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

# ----------------------------- #
# 13. Saving helpers ------------
# ----------------------------- #

save_simulation_bundle <- function(results, file_stub) {
  saveRDS(results, paste0(file_stub, ".rds"))
  write.csv(
    results %>% select(-bottleneck_table, -nca_model),
    paste0(file_stub, ".csv"),
    row.names = FALSE
  )
}

# ----------------------------------------- #
# 14. Latent-only sample-size simulation ----
# ----------------------------------------- #

run_latent_n_only_replication <- function(rep_id,
                                          n = 1000,
                                          beta = 0.40,
                                          ceiling = "ce_fdh",
                                          test_rep = 0,
                                          steps = 10,
                                          cutoff = 0,
                                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed + rep_id)
  
  dat_latent <- simulate_latent_from_common_source(
    n = n,
    beta = beta
  )
  
  out <- run_one_nca(
    dat = dat_latent,
    x_var = "eta_x",
    y_var = "eta_y",
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff
  )
  
  out %>%
    mutate(
      rep_id = rep_id,
      block = "latent_n_only",
      condition_label = paste0("n_", n),
      condition_value = as.character(n),
      score_type = "latent",
      n_input = n,
      beta = beta,
      latent_effect_base_full = effect_size,
      latent_correlation_base_full = correlation,
      latent_effect_condition_full = effect_size,
      latent_correlation_condition_full = correlation,
      latent_effect_condition_selected = effect_size,
      latent_correlation_condition_selected = correlation,
      delta_effect_from_base_full = 0,
      delta_effect_from_condition_full = 0,
      delta_effect_from_condition_selected = 0,
      selection_effect_on_latent_d = 0,
      latent_n_full = n_used,
      latent_n_selected = n_used
    )
}

simulate_latent_n_only <- function(reps = 100,
                                   n = 1000,
                                   beta = 0.40,
                                   progress = TRUE,
                                   seed = NULL,
                                   parallel = FALSE,
                                   ceiling = "ce_fdh",
                                   test_rep = 0,
                                   steps = 10,
                                   cutoff = 0) {
  rep_ids <- seq_len(reps)
  
  if (!parallel) {
    out <- vector("list", length = reps)
    
    for (r in rep_ids) {
      if (progress && (r %% max(1, floor(reps / 10)) == 0 || r == 1 || r == reps)) {
        message("Running latent-only replication ", r, " / ", reps, " for n = ", n)
      }
      
      out[[r]] <- run_latent_n_only_replication(
        rep_id = r,
        n = n,
        beta = beta,
        ceiling = ceiling,
        test_rep = test_rep,
        steps = steps,
        cutoff = cutoff,
        seed = seed
      )
    }
    
    return(bind_rows(out))
  }
  
  furrr::future_map_dfr(
    rep_ids,
    function(r) {
      run_latent_n_only_replication(
        rep_id = r,
        n = n,
        beta = beta,
        ceiling = ceiling,
        test_rep = test_rep,
        steps = steps,
        cutoff = cutoff,
        seed = seed
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
}

simulate_latent_n_only_across_n <- function(n_values = c(250, 500, 1000, 2000, 3000, 5000, 7000),
                                            reps = 100,
                                            beta = 0.40,
                                            progress = TRUE,
                                            seed = NULL,
                                            parallel_inner = FALSE,
                                            ceiling = "ce_fdh",
                                            test_rep = 0,
                                            steps = 10,
                                            cutoff = 0) {
  design <- tibble(
    sample_size = unique(n_values),
    design_id = seq_along(unique(n_values))
  )
  
  purrr::pmap_dfr(design, function(sample_size, design_id) {
    if (progress) {
      message(
        "Running latent-only design cell ", design_id, " / ", nrow(design),
        ": n = ", sample_size
      )
    }
    
    cell_seed <- if (is.null(seed)) NULL else seed + design_id * 100000
    
    simulate_latent_n_only(
      reps = reps,
      n = sample_size,
      beta = beta,
      progress = progress,
      seed = cell_seed,
      parallel = parallel_inner,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff
    ) %>%
      mutate(
        sample_size = sample_size,
        design_id = design_id
      )
  })
}

summarise_latent_n_only <- function(results) {
  results %>%
    group_by(n_input, beta) %>%
    summarise(
      n_replications = n(),
      mean_effect_size = mean(effect_size, na.rm = TRUE),
      sd_effect_size = sd(effect_size, na.rm = TRUE),
      se_effect_size = sd_effect_size / sqrt(n_replications),
      ci_low = mean_effect_size - 1.96 * se_effect_size,
      ci_high = mean_effect_size + 1.96 * se_effect_size,
      mean_correlation = mean(correlation, na.rm = TRUE),
      mean_x_skew = mean(x_skew_obs, na.rm = TRUE),
      mean_y_skew = mean(y_skew_obs, na.rm = TRUE),
      mean_x_kurtosis = mean(x_kurtosis_obs, na.rm = TRUE),
      mean_y_kurtosis = mean(y_kurtosis_obs, na.rm = TRUE),
      prop_significant = mean(significant %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(n_input)
}


# ----------------------------- #
# 15. Main Runs (commented) -----
# ----------------------------- #

# future::plan(future::multisession, workers = 8)#parallelly::availableCores() - 2)

# # One design cell
# startT <- Sys.time()
# res_one <- simulate_unified_paired(
#   reps = 4,
#   n = 500,
#   beta = 0.40,
#   test_rep = 0,
#   progress = TRUE,
#   parallel = FALSE
# )
# endT <- Sys.time()
# endT-startT
# tab_one <- summarise_unified_paired_results(res_one)
# print(tab_one)

# # Across outer design factors (2 min each)
# startT <- Sys.time()
# res_grid <- simulate_unified_paired_across_design(
#   n_values = c(250, 500, 1000),
#   beta_values = c(0.20, 0.40, 0.60),
#   reps = 1000,
#   test_rep = 0,
#   progress = TRUE,
#   parallel_inner = TRUE,
#   seed = 20260411
# )
# endT <- Sys.time()
# endT-startT
# 
# save_simulation_bundle(
#   results = res_grid,
#   file_stub = "simulations/simulation_results"
# )
# tab_sum <- summarise_unified_paired_results(res_grid)
# write.csv(
#   tab_sum,
#   "simulations/results_table.csv",
#   row.names = FALSE
# )

# ----------------------------- #
# 16. Sample size runs ----------
# ----------------------------- #

# res_latent_n <- simulate_latent_n_only_across_n(
#   n_values = c(250, 500, 1000, 2000, 4000, 
#                6000, 8000, 10000),
#   reps = 1000,
#   beta = 0.40,
#   test_rep = 0,
#   progress = TRUE,
#   parallel_inner = TRUE,
#   seed = 20260412
# )
# tab_latent_n <- summarise_latent_n_only(res_latent_n)
# 
# save_simulation_bundle(
#   results = res_latent_n,
#   file_stub = "simulations/simulation_results_n"
# )
# write.csv(
#   tab_latent_n,
#   "simulations/results_table_n.csv",
#   row.names = FALSE
# )
# 
# res_latent_n %>%
#   group_by(n_input) %>%
#   summarise(
#     mean_effect = mean(effect_size, na.rm = TRUE),
#     se = sd(effect_size, na.rm = TRUE) / sqrt(n()),
#     ci_low = mean_effect - 1.96 * se,
#     ci_high = mean_effect + 1.96 * se,
#     .groups = "drop"
#   ) %>%
#   ggplot(aes(x = n_input, y = mean_effect)) +
#   geom_point(size = 2.8) +
#   geom_line() +
#   geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 80) +
#   labs(
#     x = "Sample size",
#     y = "Mean latent NCA effect size"
#   ) +
#   theme_minimal(base_size = 12)
