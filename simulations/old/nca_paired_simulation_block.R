# ================================================================
# Paired simulation block for NCA with psychological data
# ---------------------------------------------------------------
# Purpose:
#   Add-on module to be sourced AFTER the main simulation script.
#   It implements a paired / common-source design:
#     - one latent dataset per replication
#     - one latent-reference NCA on the full latent data
#     - optional latent-reference NCA on the selected latent data
#     - multiple psychometric modifications derived from the SAME
#       latent source within replication
#
# Main use:
#   Quantify how much observed-score NCA estimates move away from the
#   latent-reference NCA estimate when psychometric features change.
# ================================================================

# -----------------------------
# 0. Checks
# -----------------------------
source("nca_simulation_1A_1E.R")
required_main_functions <- c(
  "simulate_latent_pair",
  "apply_range_restriction",
  "run_one_nca"
)

missing_main_functions <- required_main_functions[
  !vapply(required_main_functions, exists, logical(1), mode = "function")
]

if (length(missing_main_functions) > 0) {
  stop(
    "Please source the main simulation script first. Missing functions: ",
    paste(missing_main_functions, collapse = ", ")
  )
}

# -----------------------------
# 1. Small helpers
# -----------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

safe_match_arg <- function(x, choices) {
  if (length(x) == 1 && x %in% choices) return(x)
  match.arg(x, choices)
}

safe_run_one_nca <- function(dat,
                             x_var,
                             y_var,
                             ceiling = "ce_fdh",
                             test_rep = 0,
                             steps = 10,
                             cutoff = 0,
                             keep_models = FALSE) {
  out <- tryCatch(
    run_one_nca(
      dat = dat,
      x_var = x_var,
      y_var = y_var,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff
    ),
    error = function(e) {
      tibble::tibble(
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
        n_used = NA_integer_,
        correlation = NA_real_,
        x_mean_obs = NA_real_,
        y_mean_obs = NA_real_,
        x_sd_obs = NA_real_,
        y_sd_obs = NA_real_,
        x_skew_obs = NA_real_,
        y_skew_obs = NA_real_,
        nca_error = conditionMessage(e)
      )
    }
  )

  if (!"nca_error" %in% names(out)) {
    out$nca_error <- NA_character_
  }

  if (!keep_models) {
    out <- out %>%
      dplyr::select(-dplyr::any_of(c("bottleneck_table", "nca_model")))
  }

  out
}

extract_reference_core <- function(df, prefix) {
  keep <- c(
    "effect_size",
    "correlation",
    "ceiling_accuracy",
    "n_used",
    "x_skew_obs",
    "y_skew_obs",
    "nca_error"
  )

  df %>%
    dplyr::select(dplyr::any_of(keep)) %>%
    dplyr::rename_with(~ paste0(prefix, .x))
}

# -----------------------------
# 2. Threshold and item helpers
# -----------------------------

# Uses the current make_thresholds() if available; otherwise falls back
# to a local implementation consistent with the latest discussion.
make_thresholds_paired <- function(n_cat = 5,
                                   threshold_spacing = c("central", "shifted"),
                                   shift_direction = c("none", "lower", "upper"),
                                   shift_magnitude = 0,
                                   central_sd = 1.2) {
  threshold_spacing <- safe_match_arg(threshold_spacing, c("central", "shifted"))
  shift_direction <- safe_match_arg(shift_direction, c("none", "lower", "upper"))

  if (exists("make_thresholds", mode = "function")) {
    fn <- get("make_thresholds", mode = "function")
    fml <- names(formals(fn))

    args <- list(n_cat = n_cat)

    # argument name for spacing is stable in all versions we used
    if ("spacing" %in% fml) {
      args$spacing <- threshold_spacing
    } else if ("threshold_spacing" %in% fml) {
      args$threshold_spacing <- threshold_spacing
    }

    if ("shift_direction" %in% fml) {
      args$shift_direction <- shift_direction
    } else if ("skew_direction" %in% fml) {
      # best effort fallback for older code paths
      args$skew_direction <- if (shift_direction == "none") "lower" else shift_direction
    } else if ("threshold_skew_direction" %in% fml) {
      args$threshold_skew_direction <- if (shift_direction == "none") "lower" else shift_direction
    }

    if ("shift_magnitude" %in% fml) {
      args$shift_magnitude <- shift_magnitude
    } else if ("shift_size" %in% fml) {
      args$shift_size <- shift_magnitude
    }

    if ("central_sd" %in% fml) {
      args$central_sd <- central_sd
    } else if ("latent_sd_for_shape" %in% fml) {
      args$latent_sd_for_shape <- central_sd
    }

    out <- tryCatch(do.call(fn, args), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }

  # Fallback implementation
  centers <- seq(
    from = -(n_cat - 1) / 2,
    to   =  (n_cat - 1) / 2,
    length.out = n_cat
  )

  probs <- stats::dnorm(centers, mean = 0, sd = central_sd)
  probs <- probs / sum(probs)
  cum_probs <- cumsum(probs)[seq_len(n_cat - 1)]
  thresholds <- stats::qnorm(cum_probs)

  if (threshold_spacing == "shifted") {
    if (shift_direction == "lower") thresholds <- thresholds - shift_magnitude
    if (shift_direction == "upper") thresholds <- thresholds + shift_magnitude
  }

  thresholds
}

latent_to_ordinal_paired <- function(y_star, thresholds) {
  as.integer(cut(y_star,
                 breaks = c(-Inf, thresholds, Inf),
                 labels = FALSE,
                 right = TRUE))
}

build_items_from_latent_source <- function(eta,
                                           error_mat,
                                           prefix,
                                           loading = 0.70,
                                           n_cat = 5,
                                           threshold_spacing = "central",
                                           shift_direction = "none",
                                           shift_magnitude = 0,
                                           central_sd = 1.2) {
  n_items <- ncol(error_mat)

  lambda <- if (length(loading) == 1) {
    rep(loading, n_items)
  } else {
    loading
  }

  thresholds <- make_thresholds_paired(
    n_cat = n_cat,
    threshold_spacing = threshold_spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd
  )

  out <- vector("list", n_items)

  for (j in seq_len(n_items)) {
    err_sd <- sqrt(max(1 - lambda[j]^2, 1e-6))
    y_star <- lambda[j] * eta + err_sd * error_mat[, j]
    out[[j]] <- latent_to_ordinal_paired(y_star, thresholds)
  }

  names(out) <- paste0(prefix, seq_len(n_items))
  out <- tibble::as_tibble(out)
  out
}

# -----------------------------
# 3. Score construction helper
# -----------------------------

construct_selected_scores <- function(dat,
                                      x_items,
                                      y_items,
                                      score_types = "sum") {
  score_types <- unique(score_types)

  # Manual sums / means to avoid unnecessary factor-score computation
  if (any(score_types %in% c("sum", "mean"))) {
    dat <- dat %>%
      dplyr::mutate(
        x_sum  = rowSums(dplyr::across(dplyr::all_of(x_items))),
        y_sum  = rowSums(dplyr::across(dplyr::all_of(y_items))),
        x_mean = rowMeans(dplyr::across(dplyr::all_of(x_items))),
        y_mean = rowMeans(dplyr::across(dplyr::all_of(y_items)))
      )
  }

  fit_x <- NULL
  fit_y <- NULL

  if ("factor" %in% score_types) {
    if (!exists("make_lavaan_factor_scores", mode = "function")) {
      stop("The main script does not expose make_lavaan_factor_scores().")
    }

    x_fs <- make_lavaan_factor_scores(dat, x_items, factor_name = "FX")
    y_fs <- make_lavaan_factor_scores(dat, y_items, factor_name = "FY")

    dat <- dat %>%
      dplyr::mutate(
        x_factor = x_fs$scores,
        y_factor = y_fs$scores
      )

    fit_x <- x_fs$fit
    fit_y <- y_fs$fit
  }

  list(data = dat, fit_x = fit_x, fit_y = fit_y)
}

# -----------------------------
# 4. Condition grid builders
# -----------------------------

make_paired_condition_grid <- function(loading = 0.70,
                                       n_cat = 5,
                                       threshold_spacing = "central",
                                       shift_direction = "none",
                                       shift_magnitude = 0,
                                       central_sd = 1.2,
                                       restriction = "full") {
  tidyr::crossing(
    loading = loading,
    n_cat = n_cat,
    threshold_spacing = threshold_spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd,
    restriction = restriction
  ) %>%
    dplyr::mutate(condition_id = dplyr::row_number()) %>%
    dplyr::select(condition_id, dplyr::everything())
}

make_block_1B_grid <- function(loading_levels = c(0.50, 0.70, 0.85),
                               n_cat = 5,
                               threshold_spacing = "central",
                               shift_direction = "none",
                               shift_magnitude = 0,
                               central_sd = 1.2,
                               restriction = "full") {
  make_paired_condition_grid(
    loading = loading_levels,
    n_cat = n_cat,
    threshold_spacing = threshold_spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd,
    restriction = restriction
  ) %>%
    dplyr::mutate(
      reliability_condition = dplyr::case_when(
        abs(loading - 0.50) < 1e-8 ~ "low",
        abs(loading - 0.70) < 1e-8 ~ "medium",
        abs(loading - 0.85) < 1e-8 ~ "high",
        TRUE ~ paste0("loading_", loading)
      )
    )
}

make_block_1C_grid <- function(n_cat_levels = c(5, 7),
                               shift_magnitudes = c(0.25, 0.50),
                               central_sd = 1.2,
                               restriction = "full") {
  baseline <- tidyr::crossing(
    n_cat = n_cat_levels,
    threshold_spacing = "central",
    shift_direction = "none",
    shift_magnitude = 0,
    central_sd = central_sd,
    restriction = restriction
  )

  shifted <- tidyr::crossing(
    n_cat = n_cat_levels,
    threshold_spacing = "shifted",
    shift_direction = c("lower", "upper"),
    shift_magnitude = shift_magnitudes,
    central_sd = central_sd,
    restriction = restriction
  )

  dplyr::bind_rows(baseline, shifted) %>%
    dplyr::mutate(condition_id = dplyr::row_number()) %>%
    dplyr::select(condition_id, dplyr::everything())
}

make_block_1E_grid <- function(restrictions = c("full", "x_upper", "x_middle", "y_upper"),
                               loading = 0.70,
                               n_cat = 5,
                               threshold_spacing = "central",
                               shift_direction = "none",
                               shift_magnitude = 0,
                               central_sd = 1.2) {
  make_paired_condition_grid(
    loading = loading,
    n_cat = n_cat,
    threshold_spacing = threshold_spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd,
    restriction = restrictions
  )
}

# -----------------------------
# 5. Core paired engine
# -----------------------------

run_paired_replication <- function(rep_id,
                                   condition_grid,
                                   n = 1000,
                                   beta = 0.40,
                                   n_items = 6,
                                   score_types = "sum",
                                   x_skew = "none",
                                   y_skew = "none",
                                   ceiling = "ce_fdh",
                                   test_rep = 0,
                                   steps = 10,
                                   cutoff = 0,
                                   seed = NULL,
                                   keep_models = FALSE) {
  score_types <- unique(score_types)

  if (!is.null(seed)) set.seed(seed + rep_id)

  latent <- simulate_latent_pair(
    n = n,
    beta = beta,
    x_skew = x_skew,
    y_skew = y_skew
  )

  # Common random numbers within replication: same base item errors for
  # all psychometric conditions derived from this latent source.
  err_x <- matrix(stats::rnorm(n * n_items), nrow = n, ncol = n_items)
  err_y <- matrix(stats::rnorm(n * n_items), nrow = n, ncol = n_items)

  latent_full <- safe_run_one_nca(
    dat = latent,
    x_var = "eta_x",
    y_var = "eta_y",
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff,
    keep_models = keep_models
  ) %>%
    extract_reference_core(prefix = "latent_full_")

  purrr::pmap_dfr(condition_grid, function(...) {
    cond <- list(...)

    x_items <- build_items_from_latent_source(
      eta = latent$eta_x,
      error_mat = err_x,
      prefix = "x",
      loading = cond$loading,
      n_cat = cond$n_cat,
      threshold_spacing = cond$threshold_spacing,
      shift_direction = cond$shift_direction,
      shift_magnitude = cond$shift_magnitude,
      central_sd = cond$central_sd
    )

    y_items <- build_items_from_latent_source(
      eta = latent$eta_y,
      error_mat = err_y,
      prefix = "y",
      loading = cond$loading,
      n_cat = cond$n_cat,
      threshold_spacing = cond$threshold_spacing,
      shift_direction = cond$shift_direction,
      shift_magnitude = cond$shift_magnitude,
      central_sd = cond$central_sd
    )

    dat_full <- dplyr::bind_cols(latent, x_items, y_items)
    dat_selected <- apply_range_restriction(dat_full, restriction = cond$restriction)

    latent_selected <- safe_run_one_nca(
      dat = dat_selected,
      x_var = "eta_x",
      y_var = "eta_y",
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      keep_models = keep_models
    ) %>%
      extract_reference_core(prefix = "latent_selected_")

    scored <- construct_selected_scores(
      dat = dat_selected,
      x_items = names(x_items),
      y_items = names(y_items),
      score_types = score_types
    )

    dat_scored <- scored$data

    score_map <- list(
      sum = c("x_sum", "y_sum"),
      mean = c("x_mean", "y_mean"),
      factor = c("x_factor", "y_factor")
    )

    observed <- purrr::map_dfr(score_types, function(st) {
      vars <- score_map[[st]]
      safe_run_one_nca(
        dat = dat_scored,
        x_var = vars[[1]],
        y_var = vars[[2]],
        ceiling = ceiling,
        test_rep = test_rep,
        steps = steps,
        cutoff = cutoff,
        keep_models = keep_models
      ) %>%
        dplyr::mutate(score_type = st)
    })

    dplyr::bind_cols(observed, latent_full, latent_selected) %>%
      dplyr::mutate(
        rep_id = rep_id,
        condition_id = cond$condition_id,
        loading = cond$loading,
        n_cat = cond$n_cat,
        threshold_spacing = cond$threshold_spacing,
        shift_direction = cond$shift_direction,
        shift_magnitude = cond$shift_magnitude,
        central_sd = cond$central_sd,
        restriction = cond$restriction,
        n_input = n,
        beta = beta,
        n_items = n_items,
        x_skew = x_skew,
        y_skew = y_skew,
        ceiling = ceiling,
        test_rep = test_rep,
        delta_effect_from_latent_full = effect_size - latent_full_effect_size,
        delta_effect_from_latent_selected = effect_size - latent_selected_effect_size,
        selection_effect_on_latent_d = latent_selected_effect_size - latent_full_effect_size,
        delta_cor_from_latent_full = correlation - latent_full_correlation,
        delta_cor_from_latent_selected = correlation - latent_selected_correlation
      )
  })
}

run_paired_condition_grid <- function(reps = 100,
                                      condition_grid,
                                      progress = TRUE,
                                      ...) {
  out <- vector("list", length = reps)

  for (r in seq_len(reps)) {
    if (progress && (r %% max(1, floor(reps / 10)) == 0 || r == 1 || r == reps)) {
      message("Running paired replication ", r, " / ", reps)
    }

    out[[r]] <- run_paired_replication(rep_id = r, condition_grid = condition_grid, ...)
  }

  dplyr::bind_rows(out)
}

# -----------------------------
# 6. Block wrappers (paired design)
# -----------------------------

simulate_1A_paired <- function(reps = 500,
                               n = 1000,
                               beta = 0.40,
                               n_items = 6,
                               score_types = c("sum", "mean", "factor"),
                               loading = 0.70,
                               n_cat = 5,
                               threshold_spacing = "central",
                               shift_direction = "none",
                               shift_magnitude = 0,
                               central_sd = 1.2,
                               x_skew = "none",
                               y_skew = "none",
                               restriction = "full",
                               ceiling = "ce_fdh",
                               test_rep = 0,
                               steps = 10,
                               cutoff = 0,
                               seed = 1234,
                               progress = TRUE,
                               keep_models = FALSE) {
  grid <- make_paired_condition_grid(
    loading = loading,
    n_cat = n_cat,
    threshold_spacing = threshold_spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd,
    restriction = restriction
  )

  run_paired_condition_grid(
    reps = reps,
    condition_grid = grid,
    progress = progress,
    n = n,
    beta = beta,
    n_items = n_items,
    score_types = score_types,
    x_skew = x_skew,
    y_skew = y_skew,
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff,
    seed = seed,
    keep_models = keep_models
  ) %>%
    dplyr::mutate(block = "1A_paired_score_construction")
}

simulate_1B_paired <- function(reps = 500,
                               loading_levels = c(0.50, 0.70, 0.85),
                               n = 1000,
                               beta = 0.40,
                               n_items = 6,
                               score_types = "sum",
                               n_cat = 5,
                               threshold_spacing = "central",
                               shift_direction = "none",
                               shift_magnitude = 0,
                               central_sd = 1.2,
                               x_skew = "none",
                               y_skew = "none",
                               restriction = "full",
                               ceiling = "ce_fdh",
                               test_rep = 0,
                               steps = 10,
                               cutoff = 0,
                               seed = 2234,
                               progress = TRUE,
                               keep_models = FALSE) {
  grid <- make_block_1B_grid(
    loading_levels = loading_levels,
    n_cat = n_cat,
    threshold_spacing = threshold_spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd,
    restriction = restriction
  )

  run_paired_condition_grid(
    reps = reps,
    condition_grid = grid,
    progress = progress,
    n = n,
    beta = beta,
    n_items = n_items,
    score_types = score_types,
    x_skew = x_skew,
    y_skew = y_skew,
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff,
    seed = seed,
    keep_models = keep_models
  ) %>%
    dplyr::mutate(block = "1B_paired_reliability")
}

simulate_1C_paired <- function(reps = 500,
                               n_cat_levels = c(5, 7),
                               shift_magnitudes = c(0.25, 0.50),
                               n = 1000,
                               beta = 0.40,
                               n_items = 6,
                               score_types = "sum",
                               loading = 0.70,
                               central_sd = 1.2,
                               x_skew = "none",
                               y_skew = "none",
                               restriction = "full",
                               ceiling = "ce_fdh",
                               test_rep = 0,
                               steps = 10,
                               cutoff = 0,
                               seed = 3234,
                               progress = TRUE,
                               keep_models = FALSE) {
  grid <- make_block_1C_grid(
    n_cat_levels = n_cat_levels,
    shift_magnitudes = shift_magnitudes,
    central_sd = central_sd,
    restriction = restriction
  ) %>%
    dplyr::mutate(loading = loading)

  run_paired_condition_grid(
    reps = reps,
    condition_grid = grid,
    progress = progress,
    n = n,
    beta = beta,
    n_items = n_items,
    score_types = score_types,
    x_skew = x_skew,
    y_skew = y_skew,
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff,
    seed = seed,
    keep_models = keep_models
  ) %>%
    dplyr::mutate(block = "1C_paired_ordinal_scaling")
}

simulate_1D_paired <- function(reps = 500,
                               skew_conditions = NULL,
                               n = 1000,
                               beta = 0.40,
                               n_items = 6,
                               score_types = "sum",
                               loading = 0.70,
                               n_cat = 5,
                               threshold_spacing = "central",
                               shift_direction = "none",
                               shift_magnitude = 0,
                               central_sd = 1.2,
                               restriction = "full",
                               ceiling = "ce_fdh",
                               test_rep = 0,
                               steps = 10,
                               cutoff = 0,
                               seed = 4234,
                               progress = TRUE,
                               keep_models = FALSE) {
  if (is.null(skew_conditions)) {
    skew_conditions <- tibble::tibble(
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

  grid <- make_paired_condition_grid(
    loading = loading,
    n_cat = n_cat,
    threshold_spacing = threshold_spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd,
    restriction = restriction
  )

  purrr::pmap_dfr(skew_conditions, function(skew_label, x_skew, y_skew) {
    run_paired_condition_grid(
      reps = reps,
      condition_grid = grid,
      progress = progress,
      n = n,
      beta = beta,
      n_items = n_items,
      score_types = score_types,
      x_skew = x_skew,
      y_skew = y_skew,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      seed = seed + which(skew_conditions$skew_label == skew_label),
      keep_models = keep_models
    ) %>%
      dplyr::mutate(skew_label = skew_label)
  }) %>%
    dplyr::mutate(block = "1D_paired_skewness")
}

simulate_1E_paired <- function(reps = 500,
                               restrictions = c("full", "x_upper", "x_middle", "y_upper"),
                               n = 1000,
                               beta = 0.40,
                               n_items = 6,
                               score_types = "sum",
                               loading = 0.70,
                               n_cat = 5,
                               threshold_spacing = "central",
                               shift_direction = "none",
                               shift_magnitude = 0,
                               central_sd = 1.2,
                               x_skew = "none",
                               y_skew = "none",
                               ceiling = "ce_fdh",
                               test_rep = 0,
                               steps = 10,
                               cutoff = 0,
                               seed = 5234,
                               progress = TRUE,
                               keep_models = FALSE) {
  grid <- make_block_1E_grid(
    restrictions = restrictions,
    loading = loading,
    n_cat = n_cat,
    threshold_spacing = threshold_spacing,
    shift_direction = shift_direction,
    shift_magnitude = shift_magnitude,
    central_sd = central_sd
  )

  run_paired_condition_grid(
    reps = reps,
    condition_grid = grid,
    progress = progress,
    n = n,
    beta = beta,
    n_items = n_items,
    score_types = score_types,
    x_skew = x_skew,
    y_skew = y_skew,
    ceiling = ceiling,
    test_rep = test_rep,
    steps = steps,
    cutoff = cutoff,
    seed = seed,
    keep_models = keep_models
  ) %>%
    dplyr::mutate(block = "1E_paired_range_restriction")
}

# -----------------------------
# 7. Summaries for paired design
# -----------------------------

summarise_paired_results <- function(results) {
  results %>%
    dplyr::group_by(
      block,
      score_type,
      loading,
      n_cat,
      threshold_spacing,
      shift_direction,
      shift_magnitude,
      central_sd,
      x_skew,
      y_skew,
      restriction
    ) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      mean_effect_size = mean(effect_size, na.rm = TRUE),
      sd_effect_size = stats::sd(effect_size, na.rm = TRUE),
      mean_latent_full_d = mean(latent_full_effect_size, na.rm = TRUE),
      mean_latent_selected_d = mean(latent_selected_effect_size, na.rm = TRUE),
      mean_delta_full = mean(delta_effect_from_latent_full, na.rm = TRUE),
      mean_abs_delta_full = mean(abs(delta_effect_from_latent_full), na.rm = TRUE),
      mean_delta_selected = mean(delta_effect_from_latent_selected, na.rm = TRUE),
      mean_abs_delta_selected = mean(abs(delta_effect_from_latent_selected), na.rm = TRUE),
      mean_selection_shift = mean(selection_effect_on_latent_d, na.rm = TRUE),
      mean_correlation = mean(correlation, na.rm = TRUE),
      mean_latent_full_cor = mean(latent_full_correlation, na.rm = TRUE),
      mean_n_used = mean(n_used, na.rm = TRUE),
      .groups = "drop"
    )
}

# -----------------------------
# 8. Example runs (commented)
# -----------------------------

# source("nca_simulation_1A_1E.R")
# source("nca_paired_simulation_block.R")
#
# # 1A paired: same latent source, compare score constructions
# res_1A_p <- simulate_1A_paired(
#   reps = 20,
#   n = 1000,
#   score_types = c("sum", "factor"),
#   test_rep = 0
# )
# sum_1A_p <- summarise_paired_results(res_1A_p)
# print(sum_1A_p)
#
# # 1C paired: same latent source, compare ordinal scaling conditions
# res_1C_p <- simulate_1C_paired(
#   reps = 20,
#   n_cat_levels = c(5, 7),
#   shift_magnitudes = c(0.25, 0.50),
#   score_types = "sum",
#   test_rep = 0
# )
# sum_1C_p <- summarise_paired_results(res_1C_p)
# print(sum_1C_p)
#
# # 1E paired: decompose selection vs measurement effects
# res_1E_p <- simulate_1E_paired(
#   reps = 20,
#   restrictions = c("full", "x_upper", "x_middle", "y_upper"),
#   score_types = "sum",
#   test_rep = 0
# )
# sum_1E_p <- summarise_paired_results(res_1E_p)
# print(sum_1E_p)
