# ================================================================
# Unified paired/common-source simulation block for NCA
# ---------------------------------------------------------------
# Usage:
#   source("nca_simulation_1A_1E.R")
#   source("nca_all_blocks_unified_paired.R")
#
# Purpose:
#   In each replication, generate ONE common latent source dataset,
#   estimate the latent-reference NCA, and then derive all 1A-1E
#   psychometric conditions from that same source.
#
# Design logic:
#   - same latent source per replication
#   - all blocks branch from that same source
#   - block 1D uses the same underlying random source but applies
#     alternative monotone remappings to create skewed latent margins
#   - block 1E also stores the latent NCA after selection/restriction
#
# Main outputs:
#   - observed NCA estimates for every condition
#   - latent references attached to every row:
#       * latent_effect_base_full
#       * latent_effect_condition_full
#       * latent_effect_condition_selected
#   - deltas relative to those references
# ================================================================

# -----------------------------
# 0. Defensive checks
# -----------------------------

required_upstream <- c(
  "run_one_nca", "apply_range_restriction", "induce_skew", "standardize",
  "sample_skewness"
)

missing_upstream <- required_upstream[!vapply(required_upstream, exists, logical(1), inherits = TRUE)]
if (length(missing_upstream) > 0) {
  stop(
    "Please source your main simulation script before this add-on. Missing objects: ",
    paste(missing_upstream, collapse = ", ")
  )
}

# -----------------------------
# 1. Local helpers
# -----------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

.validate_score_types <- function(score_types) {
  score_types <- unique(score_types)
  allowed <- c("sum", "mean", "factor")
  bad <- setdiff(score_types, allowed)
  if (length(bad) > 0) {
    stop("Unsupported score_types: ", paste(bad, collapse = ", "))
  }
  score_types
}

make_thresholds_unified <- function(n_cat = 5,
                                    spacing = c("central", "shifted"),
                                    shift_direction = c("none", "lower", "upper"),
                                    shift_magnitude = 0.35,
                                    central_sd = 1.2) {
  spacing <- match.arg(spacing)
  shift_direction <- match.arg(shift_direction)

  if (n_cat < 2) stop("n_cat must be at least 2.")
  if (central_sd <= 0) stop("central_sd must be > 0.")
  if (shift_magnitude < 0) stop("shift_magnitude must be >= 0.")

  centers <- seq(
    from = -(n_cat - 1) / 2,
    to   =  (n_cat - 1) / 2,
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

latent_to_ordinal_unified <- function(y_star, thresholds) {
  as.integer(cut(
    y_star,
    breaks = c(-Inf, thresholds, Inf),
    labels = FALSE,
    right = TRUE
  ))
}

simulate_scale_items_unified <- function(eta,
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

  thresholds <- make_thresholds_unified(
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
    out[[j]] <- latent_to_ordinal_unified(y_star, thresholds)
  }

  names(out) <- paste0(prefix, seq_len(n_items))
  out <- as_tibble(out)
  out
}

make_lavaan_factor_scores_unified <- function(dat, items, factor_name = "F") {
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

construct_scores_unified <- function(dat, 
                                     x_items, 
                                     y_items, 
                                     score_types = "sum") {
  score_types <- .validate_score_types(score_types)

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
    x_fs <- make_lavaan_factor_scores_unified(dat, x_items, factor_name = "FX")
    y_fs <- make_lavaan_factor_scores_unified(dat, y_items, factor_name = "FY")

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

latent_from_common_source <- function(z_x,
                                      z_eps,
                                      beta = 0.40,
                                      x_skew = c("none", "positive", "negative"),
                                      y_skew = c("none", "positive", "negative"),
                                      residual_sd = NULL) {
  x_skew <- match.arg(x_skew)
  y_skew <- match.arg(y_skew)

  if (is.null(residual_sd)) {
    residual_sd <- sqrt(1 - beta^2)
  }

  eta_x <- induce_skew(z_x, kind = x_skew)
  eps <- standardize(z_eps) * residual_sd
  eta_y_linear <- beta * eta_x + eps
  eta_y <- induce_skew(eta_y_linear, kind = y_skew)

  tibble::tibble(
    id = seq_along(z_x),
    eta_x = eta_x,
    eta_y = eta_y
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
    x_skew_obs = ref$x_skew_obs[[1]],
    y_skew_obs = ref$y_skew_obs[[1]],
    nca_model = ref$nca_model[[1]],
    bottleneck_table = ref$bottleneck_table[[1]]
  )
}

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
                                   shift_magnitude = 0.35,
                                   central_sd = 1.2,
                                   ceiling = "ce_fdh",
                                   test_rep = 0,
                                   steps = 10,
                                   cutoff = 0,
                                   base_reference,
                                   condition_full_reference,
                                   condition_selected_reference) {
  x_items <- simulate_scale_items_unified(
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

  y_items <- simulate_scale_items_unified(
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
  scored <- construct_scores_unified(
    dat = dat_obs,
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
        spacing = spacing,
        shift_direction = shift_direction,
        shift_magnitude = shift_magnitude,
        central_sd = central_sd,
        n_cat = n_cat,
        loading = loading,
        latent_n_full = nrow(latent_full),
        latent_n_selected = nrow(latent_selected)
      )
  })
}

# -----------------------------
# 2. One unified paired replication
# -----------------------------

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
                                           baseline_x_skew = "none",
                                           baseline_y_skew = "none",
                                           baseline_restriction = "full",
                                           score_types_1A = c("sum", "factor"),
                                           loading_levels_1B = c(0.50, 0.70, 0.85),
                                           n_cat_levels_1C = c(5, 7),
                                           shift_magnitudes_1C = c(0.25, 0.50),
                                           shift_directions_1C = c("lower", "upper"),
                                           central_sd_1C = 1.2,
                                           skew_conditions_1D = NULL,
                                           restrictions_1E = c("full", "x_upper", "x_middle", "y_upper"),
                                           ceiling = "ce_fdh",
                                           test_rep = 0,
                                           steps = 10,
                                           cutoff = 0,
                                           seed = NULL) {
  score_types_1A <- .validate_score_types(score_types_1A)

  if (!is.null(seed)) set.seed(seed + rep_id)

  if (is.null(skew_conditions_1D)) {
    skew_conditions_1D <- tibble::tibble(
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

  # Common source randomness for the entire replication.
  z_x <- stats::rnorm(n)
  z_eps <- stats::rnorm(n)

  latent_base_full <- latent_from_common_source(
    z_x = z_x,
    z_eps = z_eps,
    beta = beta,
    x_skew = baseline_x_skew,
    y_skew = baseline_y_skew
  )

  latent_base_selected <- apply_range_restriction(
    latent_base_full,
    restriction = baseline_restriction
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

  # 1A: score construction, same observed dataset, different scoring.
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

  # 1B: reliability/loading changes, same latent base.
  res_1B <- purrr::map_dfr(loading_levels_1B, function(ld) {
    run_observed_condition(
      latent_full = latent_base_full,
      latent_selected = latent_base_selected,
      block = "1B_reliability",
      condition_label = paste0("loading_", ld),
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

  # 1C: ordinal scaling changes, same latent base.
  design_1C <- tibble::tibble(
    n_cat = baseline_n_cat,
    spacing = "central",
    shift_direction = "none",
    shift_magnitude = 0,
    central_sd = central_sd_1C
  )

  if (length(n_cat_levels_1C) > 0) {
    design_1C <- dplyr::bind_rows(
      design_1C,
      tidyr::crossing(
        n_cat = n_cat_levels_1C,
        spacing = "central",
        shift_direction = "none",
        shift_magnitude = 0,
        central_sd = central_sd_1C
      )
    )
  }

  if (length(shift_magnitudes_1C) > 0 && length(shift_directions_1C) > 0) {
    design_1C <- dplyr::bind_rows(
      design_1C,
      tidyr::crossing(
        n_cat = n_cat_levels_1C,
        spacing = "shifted",
        shift_direction = shift_directions_1C,
        shift_magnitude = shift_magnitudes_1C,
        central_sd = central_sd_1C
      )
    )
  }

  design_1C <- dplyr::distinct(design_1C)

  res_1C <- purrr::pmap_dfr(design_1C, function(n_cat, spacing, shift_direction, shift_magnitude, central_sd) {
    label <- paste0(
      "cat", n_cat, "_",
      spacing, "_",
      shift_direction, "_",
      formatC(shift_magnitude, format = "f", digits = 2)
    )

    run_observed_condition(
      latent_full = latent_base_full,
      latent_selected = latent_base_selected,
      block = "1C_ordinal_scaling",
      condition_label = label,
      condition_value = label,
      score_types = "sum",
      n_items = n_items,
      loading = baseline_loading,
      n_cat = n_cat,
      spacing = spacing,
      shift_direction = shift_direction,
      shift_magnitude = shift_magnitude,
      central_sd = central_sd,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff,
      base_reference = ref_base_full,
      condition_full_reference = ref_base_full,
      condition_selected_reference = ref_base_selected
    )
  })

  # 1D: latent skewness changes from the SAME common random source.
  res_1D <- purrr::pmap_dfr(skew_conditions_1D, function(skew_label, x_skew, y_skew) {
    latent_skew_full <- latent_from_common_source(
      z_x = z_x,
      z_eps = z_eps,
      beta = beta,
      x_skew = x_skew,
      y_skew = y_skew
    )

    latent_skew_selected <- apply_range_restriction(
      latent_skew_full,
      restriction = baseline_restriction
    )

    ref_skew_full <- get_reference_metrics(
      dat = latent_skew_full,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff
    )

    ref_skew_selected <- get_reference_metrics(
      dat = latent_skew_selected,
      ceiling = ceiling,
      test_rep = test_rep,
      steps = steps,
      cutoff = cutoff
    )

    run_observed_condition(
      latent_full = latent_skew_full,
      latent_selected = latent_skew_selected,
      block = "1D_skewness",
      condition_label = skew_label,
      condition_value = skew_label,
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
      condition_full_reference = ref_skew_full,
      condition_selected_reference = ref_skew_selected
    ) %>%
      dplyr::mutate(
        x_skew_condition = x_skew,
        y_skew_condition = y_skew
      )
  })

  # 1E: range restriction changes from the same latent base.
  res_1E <- purrr::map_dfr(restrictions_1E, function(rr) {
    latent_rr_selected <- apply_range_restriction(latent_base_full, restriction = rr)

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

  dplyr::bind_rows(res_1A, res_1B, res_1C, res_1D, res_1E) %>%
    dplyr::mutate(
      rep_id = rep_id,
      n_input = n,
      beta = beta,
      baseline_loading = baseline_loading,
      baseline_n_cat = baseline_n_cat,
      baseline_spacing = baseline_spacing,
      baseline_shift_direction = baseline_shift_direction,
      baseline_shift_magnitude = baseline_shift_magnitude,
      baseline_central_sd = baseline_central_sd,
      baseline_x_skew = baseline_x_skew,
      baseline_y_skew = baseline_y_skew,
      baseline_restriction = baseline_restriction,
      ceiling = ceiling,
      test_rep = test_rep
    )
}

# -----------------------------
# 3. Replicated unified run
# -----------------------------

simulate_all_blocks_unified_paired <- function(reps = 100,
                                               progress = TRUE,
                                               ...) {
  out <- vector("list", length = reps)

  for (r in seq_len(reps)) {
    if (progress && (r %% max(1, floor(reps / 10)) == 0 || r == 1 || r == reps)) {
      message("Running unified paired replication ", r, " / ", reps)
    }

    out[[r]] <- run_unified_paired_replication(rep_id = r, ...)
  }

  dplyr::bind_rows(out)
}

# -----------------------------
# 4. Summaries
# -----------------------------

summarise_unified_paired_results <- function(results) {
  results %>%
    dplyr::group_by(block, condition_label, score_type) %>%
    dplyr::summarise(
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
      mean_n_used = mean(n_used, na.rm = TRUE),
      .groups = "drop"
    )
}

# -----------------------------
# 5. Example run (commented)
# -----------------------------

# source("simulations/nca_simulation_1A_1E.R")
# source("nca_all_blocks_unified_paired.R")
#
res_all <- simulate_all_blocks_unified_paired(
  reps = 500,
  n = 1000,
  score_types_1A = c("sum", "factor"),
  shift_magnitudes_1C = c(0.25, 0.50),
  test_rep = 0,
  progress = TRUE
)

#
# sum_all <- summarise_unified_paired_results(res_all)
# print(sum_all)

# ---------------------------- 
# 6. across sample size
# ----------------------------

simulate_all_blocks_across_n <- function(
    n_values = c(250, 500, 1000),
    reps = 100,
    ...
) {
  purrr::map_dfr(n_values, function(n_current) {
    simulate_all_blocks_unified_paired(
      reps = reps,
      n = n_current,
      ...
    ) %>%
      dplyr::mutate(sample_size = n_current)
  })
}

summarise_unified_paired_results_by_n <- function(results) {
  n_var <- if ("sample_size" %in% names(results)) "sample_size" else "n_input"
  
  results %>%
    dplyr::group_by(
      .data[[n_var]],
      block,
      condition_label,
      score_type
    ) %>%
    dplyr::summarise(
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
      mean_n_used = mean(n_used, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(sample_size = .data[[n_var]]) %>%
    dplyr::arrange(sample_size, block, condition_label, score_type)
}


res_all_n <- simulate_all_blocks_across_n(
  n_values = c(250, 500, 1000),
  reps = 500,
  test_rep = 0,
  progress = TRUE
)

tab_all_n <- summarise_unified_paired_results_by_n(res_all_n)
print(tab_all_n)

writexl::write_xlsx(tab_all_n, "tabfullsim.xlsx")

tab_all_n %>%
  dplyr::filter(sample_size == 500)

tab_all_n <- tab_all_n %>%
  dplyr::arrange(sample_size, block, score_type)


ggplot(tab_all_n, aes(y=mean_effect_size, x = condition_label, group = sample_size, colour = as.factor(sample_size))) +
  geom_point() +
  facet_wrap(~block)
