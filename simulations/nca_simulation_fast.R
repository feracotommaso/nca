# ================================================================
# Fast paired NCA simulation for psychological observed scores
# ---------------------------------------------------------------
# This version avoids NCA::nca_analysis() during simulation.
# It implements only what is needed here:
#   - CE-FDH effect size, upper-left corner, empirical scope
#   - permutation p-value for CE-FDH
#   - paired latent-reference / observed-score simulation
#   - optional checkpoint CSV files overwritten after each design cell
# ================================================================

required_packages <- c("dplyr", "purrr", "tibble", "tidyr", "lavaan", "furrr", "future")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) stop("Please install: ", paste(missing_packages, collapse = ", "))

library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(furrr)
library(future)

standardize <- function(x) as.numeric(scale(x))

sample_skewness <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(NA_real_)
  s <- stats::sd(x)
  if (is.na(s) || s == 0) return(NA_real_)
  mean(((x - mean(x)) / s)^3)
}

sample_kurtosis <- function(x, excess = FALSE) {
  x <- x[!is.na(x)]
  if (length(x) < 4) return(NA_real_)
  s <- stats::sd(x)
  if (is.na(s) || s == 0) return(NA_real_)
  k <- mean(((x - mean(x)) / s)^4)
  if (excess) k - 3 else k
}

safe_cor <- function(x, y) {
  if (length(x) < 2 || stats::sd(x, na.rm = TRUE) == 0 || stats::sd(y, na.rm = TRUE) == 0) return(NA_real_)
  stats::cor(x, y, use = "complete.obs")
}

fmt_num <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

# ----------------------------- #
# Fast CE-FDH NCA ------------- #
# ----------------------------- #

ce_fdh_fast <- function(data) {
  d <- as.data.frame(data[, 1:2])
  names(d) <- c("x", "y")
  d <- stats::na.omit(d)
  
  n_used <- nrow(d)
  x_unique <- dplyr::n_distinct(d$x)
  y_unique <- dplyr::n_distinct(d$y)
  
  if (n_used < 2 || x_unique < 2 || y_unique < 2) {
    return(list(peers = tibble(x = numeric(0), y = numeric(0)), ceiling_area = NA_real_, scope_area = NA_real_, effect = NA_real_, n_used = n_used, x_unique = x_unique, y_unique = y_unique))
  }
  
  x_min <- min(d$x)
  x_max <- max(d$x)
  y_min <- min(d$y)
  y_max <- max(d$y)
  scope_area <- (x_max - x_min) * (y_max - y_min)
  
  if (!is.finite(scope_area) || scope_area <= 0) {
    return(list(peers = tibble(x = numeric(0), y = numeric(0)), ceiling_area = NA_real_, scope_area = scope_area, effect = NA_real_, n_used = n_used, x_unique = x_unique, y_unique = y_unique))
  }
  
  d_unique <- d %>%
    group_by(x) %>%
    summarise(y = max(y), .groups = "drop") %>%
    arrange(x)
  
  keep <- logical(nrow(d_unique))
  current_y <- -Inf
  
  for (i in seq_len(nrow(d_unique))) {
    if (d_unique$y[i] > current_y) {
      keep[i] <- TRUE
      current_y <- d_unique$y[i]
    }
  }
  
  peers <- d_unique[keep, , drop = FALSE]
  
  ceiling_area <- 0
  if (nrow(peers) >= 2) {
    dx <- peers$x[-1] - x_min
    dy <- diff(peers$y)
    ceiling_area <- sum(dx * dy)
  }
  
  list(peers = peers, ceiling_area = ceiling_area, scope_area = scope_area, effect = ceiling_area / scope_area, n_used = n_used, x_unique = x_unique, y_unique = y_unique)
}

ce_fdh_p_value <- function(data, observed_effect = NULL, B = 0) {
  if (B < 1) return(NA_real_)
  
  d <- as.data.frame(data[, 1:2])
  names(d) <- c("x", "y")
  d <- stats::na.omit(d)
  if (nrow(d) < 2) return(NA_real_)
  
  if (is.null(observed_effect)) observed_effect <- ce_fdh_fast(d)$effect
  if (!is.finite(observed_effect)) return(NA_real_)
  
  sim <- numeric(B)
  for (b in seq_len(B)) {
    d_perm <- d
    d_perm$y <- sample(d$y, size = nrow(d), replace = FALSE)
    sim[b] <- ce_fdh_fast(d_perm)$effect
  }
  
  sim <- sim[is.finite(sim)]
  if (length(sim) == 0) return(NA_real_)
  
  p <- (sum(sim >= observed_effect) + 1) / (length(sim) + 1)
  max(p, 1 / B)
}

run_nca_fast <- function(dat, x_var, y_var, test_rep = 0) {
  tmp <- dat %>% select(all_of(c(x_var, y_var))) %>% tidyr::drop_na()
  names(tmp) <- c("x", "y")
  
  n_used <- nrow(tmp)
  x_unique <- if (n_used > 0) dplyr::n_distinct(tmp$x) else 0
  y_unique <- if (n_used > 0) dplyr::n_distinct(tmp$y) else 0
  
  descriptives <- tibble(
    n_used = n_used,
    x_unique = x_unique,
    y_unique = y_unique,
    correlation = if (n_used >= 2) safe_cor(tmp$x, tmp$y) else NA_real_,
    x_mean_obs = if (n_used > 0) mean(tmp$x) else NA_real_,
    y_mean_obs = if (n_used > 0) mean(tmp$y) else NA_real_,
    x_sd_obs = if (n_used > 1) stats::sd(tmp$x) else NA_real_,
    y_sd_obs = if (n_used > 1) stats::sd(tmp$y) else NA_real_,
    x_skew_obs = if (n_used > 2) sample_skewness(tmp$x) else NA_real_,
    y_skew_obs = if (n_used > 2) sample_skewness(tmp$y) else NA_real_,
    x_kurtosis_obs = if (n_used > 3) sample_kurtosis(tmp$x) else NA_real_,
    y_kurtosis_obs = if (n_used > 3) sample_kurtosis(tmp$y) else NA_real_
  )
  
  if (n_used < 10 || x_unique < 2 || y_unique < 2) {
    return(bind_cols(tibble(effect_size = NA_real_, ceiling_zone = NA_real_, scope = NA_real_, p_value = NA_real_, significant = NA, test_failed = TRUE, test_error = "Degenerate data: too few rows or unique values."), descriptives))
  }
  
  out <- tryCatch(ce_fdh_fast(tmp), error = function(e) list(effect = NA_real_, ceiling_area = NA_real_, scope_area = NA_real_, error = conditionMessage(e)))
  failed <- !is.finite(out$effect)
  err <- if (!is.null(out$error)) out$error else if (failed) "Non-finite CE-FDH effect." else NA_character_
  
  p <- if (!failed && test_rep > 0) {
    tryCatch(ce_fdh_p_value(tmp, observed_effect = out$effect, B = test_rep), error = function(e) NA_real_)
  } else NA_real_
  
  bind_cols(tibble(effect_size = out$effect, ceiling_zone = out$ceiling_area, scope = out$scope_area, p_value = p, significant = ifelse(is.na(p), NA, p < 0.05), test_failed = failed, test_error = err), descriptives)
}

reference_nca <- function(dat, x_var = "eta_x", y_var = "eta_y", test_rep = 0) {
  ref <- run_nca_fast(dat, x_var, y_var, test_rep = test_rep)
  list(effect_size = ref$effect_size[[1]], correlation = ref$correlation[[1]], n_used = ref$n_used[[1]], x_skew_obs = ref$x_skew_obs[[1]], y_skew_obs = ref$y_skew_obs[[1]])
}

# ----------------------------- #
# Latent and item generation --- #
# ----------------------------- #

simulate_latent <- function(n = 1000, beta = 0.40) {
  eta_x <- standardize(stats::rnorm(n))
  eps <- standardize(stats::rnorm(n)) * sqrt(1 - beta^2)
  eta_y <- beta * eta_x + eps
  tibble(id = seq_len(n), eta_x = eta_x, eta_y = eta_y)
}

make_thresholds <- function(n_cat = 5, spacing = c("central", "shifted"), shift_direction = c("none", "lower", "upper"), shift_magnitude = 0, central_sd = 1.2) {
  spacing <- match.arg(spacing)
  shift_direction <- match.arg(shift_direction)
  centers <- seq(-(n_cat - 1) / 2, (n_cat - 1) / 2, length.out = n_cat)
  probs <- stats::dnorm(centers, 0, central_sd)
  probs <- probs / sum(probs)
  thresholds <- stats::qnorm(cumsum(probs)[seq_len(n_cat - 1)])
  if (spacing == "shifted" && shift_direction == "lower") thresholds <- thresholds - shift_magnitude
  if (spacing == "shifted" && shift_direction == "upper") thresholds <- thresholds + shift_magnitude
  thresholds
}

latent_to_ordinal <- function(y_star, thresholds) {
  as.integer(cut(y_star, breaks = c(-Inf, thresholds, Inf), labels = FALSE, right = TRUE))
}

simulate_items <- function(eta, prefix, n_items = 6, loading = 0.70, n_cat = 5, spacing = "central", shift_direction = "none", shift_magnitude = 0, central_sd = 1.2) {
  thresholds <- make_thresholds(n_cat = n_cat, spacing = spacing, shift_direction = shift_direction, shift_magnitude = shift_magnitude, central_sd = central_sd)
  items <- vector("list", n_items)
  for (j in seq_len(n_items)) {
    err_sd <- sqrt(max(1 - loading^2, 1e-8))
    y_star <- loading * eta + stats::rnorm(length(eta), 0, err_sd)
    items[[j]] <- latent_to_ordinal(y_star, thresholds)
  }
  names(items) <- paste0(prefix, seq_len(n_items))
  as_tibble(items)
}

factor_scores <- function(dat, items, factor_name) {
  model <- paste0(factor_name, " =~ ", paste(items, collapse = " + "))
  dat_ord <- dat
  for (nm in items) dat_ord[[nm]] <- ordered(dat_ord[[nm]], levels = sort(unique(dat_ord[[nm]])))
  fit <- tryCatch(lavaan::cfa(model = model, data = dat_ord, ordered = items, std.lv = TRUE, estimator = "WLSMV", parameterization = "theta"), error = function(e) NULL, warning = function(w) invokeRestart("muffleWarning"))
  if (is.null(fit)) return(rep(NA_real_, nrow(dat)))
  tryCatch(as.numeric(lavaan::lavPredict(fit)[, 1]), error = function(e) rep(NA_real_, nrow(dat)))
}

score_items <- function(dat, x_items, y_items, score_types = "sum") {
  score_types <- unique(score_types)
  out <- dat %>% mutate(x_sum = rowSums(across(all_of(x_items))), y_sum = rowSums(across(all_of(y_items))))
  if ("factor" %in% score_types) {
    out$x_factor <- factor_scores(out, x_items, "FX")
    out$y_factor <- factor_scores(out, y_items, "FY")
  }
  out
}

restrict_latent <- function(dat, restriction = c("full", "x_lower", "x_upper", "y_lower", "y_upper"), lower_q = 0.30, upper_q = 0.70) {
  restriction <- match.arg(restriction)
  if (restriction == "full") return(dat)
  if (restriction == "x_upper") return(filter(dat, eta_x > quantile(eta_x, lower_q)))
  if (restriction == "x_lower") return(filter(dat, eta_x < quantile(eta_x, upper_q)))
  if (restriction == "y_upper") return(filter(dat, eta_y > quantile(eta_y, lower_q)))
  if (restriction == "y_lower") return(filter(dat, eta_y < quantile(eta_y, upper_q)))
}

run_condition <- function(latent_full, latent_selected, block, condition_label, condition_value = condition_label, score_types = "sum", n_items = 6, loading = 0.70, n_cat = 5, spacing = "central", shift_direction = "none", shift_magnitude = 0, central_sd = 1.2, test_rep = 0, ref_base_full, ref_condition_full, ref_condition_selected) {
  x_items <- simulate_items(latent_selected$eta_x, "x", n_items, loading, n_cat, spacing, shift_direction, shift_magnitude, central_sd)
  y_items <- simulate_items(latent_selected$eta_y, "y", n_items, loading, n_cat, spacing, shift_direction, shift_magnitude, central_sd)
  dat_scored <- score_items(bind_cols(latent_selected, x_items, y_items), names(x_items), names(y_items), score_types)
  
  map_dfr(score_types, function(st) {
    vars <- if (st == "factor") c("x_factor", "y_factor") else c("x_sum", "y_sum")
    run_nca_fast(dat_scored, vars[1], vars[2], test_rep = test_rep) %>%
      mutate(
        block = block,
        condition_label = condition_label,
        condition_value = as.character(condition_value),
        score_type = st,
        latent_effect_base_full = ref_base_full$effect_size,
        latent_correlation_base_full = ref_base_full$correlation,
        latent_effect_condition_full = ref_condition_full$effect_size,
        latent_correlation_condition_full = ref_condition_full$correlation,
        latent_effect_condition_selected = ref_condition_selected$effect_size,
        latent_correlation_condition_selected = ref_condition_selected$correlation,
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

run_replication <- function(rep_id, n = 1000, beta = 0.40, n_items = 6, baseline_loading = 0.70, baseline_n_cat = 5, baseline_central_sd = 1.2, score_types_1A = c("sum", "factor"), loading_levels_1B = c(0.50, 0.70, 0.85), n_cat_levels_1C = c(4, 5, 7), shift_directions_1D = c("lower", "upper"), shift_magnitudes_1D = c(0.50, 0.70, 1.00), restrictions_1E = c("x_lower", "x_upper", "y_lower", "y_upper"), restriction_lower_q = 0.30, restriction_upper_q = 0.70, test_rep = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed + rep_id)
  latent_full <- simulate_latent(n = n, beta = beta)
  latent_selected <- restrict_latent(latent_full, "full", restriction_lower_q, restriction_upper_q)
  ref_full <- reference_nca(latent_full, test_rep = 0)
  ref_selected <- reference_nca(latent_selected, test_rep = 0)
  
  res_1A <- run_condition(latent_full, latent_selected, "1A_score_construction", "baseline_scores", "baseline", score_types_1A, n_items, baseline_loading, baseline_n_cat, central_sd = baseline_central_sd, test_rep = test_rep, ref_base_full = ref_full, ref_condition_full = ref_full, ref_condition_selected = ref_selected)
  
  res_1B <- map_dfr(loading_levels_1B, function(ld) {
    run_condition(latent_full, latent_selected, "1B_reliability", paste0("loading_", fmt_num(ld)), ld, "sum", n_items, ld, baseline_n_cat, central_sd = baseline_central_sd, test_rep = test_rep, ref_base_full = ref_full, ref_condition_full = ref_full, ref_condition_selected = ref_selected)
  })
  
  res_1C <- map_dfr(n_cat_levels_1C, function(k) {
    run_condition(latent_full, latent_selected, "1C_response_categories", paste0("n_cat_", k), k, "sum", n_items, baseline_loading, k, central_sd = baseline_central_sd, test_rep = test_rep, ref_base_full = ref_full, ref_condition_full = ref_full, ref_condition_selected = ref_selected)
  })
  
  res_1D <- pmap_dfr(tidyr::crossing(shift_direction = shift_directions_1D, shift_magnitude = shift_magnitudes_1D), function(shift_direction, shift_magnitude) {
    label <- paste0(shift_direction, "_", fmt_num(shift_magnitude))
    run_condition(latent_full, latent_selected, "1D_threshold_shift", label, label, "sum", n_items, baseline_loading, baseline_n_cat, spacing = "shifted", shift_direction = shift_direction, shift_magnitude = shift_magnitude, central_sd = baseline_central_sd, test_rep = test_rep, ref_base_full = ref_full, ref_condition_full = ref_full, ref_condition_selected = ref_selected)
  })
  
  res_1E <- map_dfr(restrictions_1E, function(rr) {
    latent_rr <- restrict_latent(latent_full, rr, restriction_lower_q, restriction_upper_q)
    ref_rr <- reference_nca(latent_rr, test_rep = 0)
    run_condition(latent_full, latent_rr, "1E_range_restriction", rr, rr, "sum", n_items, baseline_loading, baseline_n_cat, central_sd = baseline_central_sd, test_rep = test_rep, ref_base_full = ref_full, ref_condition_full = ref_full, ref_condition_selected = ref_rr)
  })
  
  bind_rows(res_1A, res_1B, res_1C, res_1D, res_1E) %>%
    mutate(rep_id = rep_id, n_input = n, beta = beta, test_rep = test_rep, restriction_lower_q = restriction_lower_q, restriction_upper_q = restriction_upper_q)
}

simulate_cell <- function(reps = 100, n = 1000, beta = 0.40, seed = NULL, parallel = FALSE, progress = TRUE, ...) {
  ids <- seq_len(reps)
  if (!parallel) {
    out <- vector("list", reps)
    for (r in ids) {
      if (progress && (r == 1 || r == reps || r %% max(1, floor(reps / 10)) == 0)) message("Replication ", r, " / ", reps, " | n = ", n, " | beta = ", beta)
      out[[r]] <- run_replication(r, n = n, beta = beta, seed = seed, ...)
    }
    return(bind_rows(out))
  }
  furrr::future_map_dfr(ids, function(r) run_replication(r, n = n, beta = beta, seed = seed, ...), .options = furrr::furrr_options(seed = TRUE))
}

simulate_grid <- function(n_values = c(250, 500, 1000), beta_values = c(0.20, 0.40, 0.60), reps = 100, seed = NULL, parallel_inner = FALSE, progress = TRUE, results_csv = NULL, summary_csv = NULL, ...) {
  design <- tidyr::crossing(n_input = n_values, beta = beta_values) %>% mutate(design_id = row_number())
  all_results <- vector("list", nrow(design))
  for (i in seq_len(nrow(design))) {
    n_i <- design$n_input[i]
    beta_i <- design$beta[i]
    cell_seed <- if (is.null(seed)) NULL else seed + i * 100000
    if (progress) message("Design cell ", i, " / ", nrow(design), ": n = ", n_i, ", beta = ", beta_i)
    all_results[[i]] <- simulate_cell(reps = reps, n = n_i, beta = beta_i, seed = cell_seed, parallel = parallel_inner, progress = progress, ...) %>% mutate(design_id = i)
    current_results <- bind_rows(all_results[seq_len(i)])
    if (!is.null(results_csv)) write.csv(current_results, results_csv, row.names = FALSE)
    if (!is.null(summary_csv)) write.csv(summarise_results(current_results), summary_csv, row.names = FALSE)
  }
  bind_rows(all_results)
}

summarise_results <- function(results) {
  results %>%
    group_by(n_input, beta, block, condition_label, score_type) %>%
    summarise(
      n_replications = dplyr::n(),
      n_failed = sum(test_failed %in% TRUE, na.rm = TRUE),
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

save_results <- function(results, results_csv = "simulation_results_fast.csv", summary_csv = "results_table_fast.csv", results_rds = NULL) {
  write.csv(results, results_csv, row.names = FALSE)
  write.csv(summarise_results(results), summary_csv, row.names = FALSE)
  if (!is.null(results_rds)) saveRDS(results, results_rds)
  invisible(results)
}

run_latent_n_replication <- function(rep_id, n = 1000, beta = 0.40, test_rep = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed + rep_id)
  dat <- simulate_latent(n, beta)
  out <- run_nca_fast(dat, "eta_x", "eta_y", test_rep = test_rep)
  out %>% mutate(rep_id = rep_id, block = "latent_n_only", condition_label = paste0("n_", n), condition_value = as.character(n), score_type = "latent", n_input = n, beta = beta, latent_effect_base_full = effect_size, latent_correlation_base_full = correlation, latent_effect_condition_full = effect_size, latent_correlation_condition_full = correlation, latent_effect_condition_selected = effect_size, latent_correlation_condition_selected = correlation, delta_effect_from_base_full = 0, delta_effect_from_condition_full = 0, delta_effect_from_condition_selected = 0, selection_effect_on_latent_d = 0, latent_n_full = n_used, latent_n_selected = n_used, test_rep = test_rep)
}

simulate_latent_n_only <- function(n_values = c(250, 500, 1000, 2000, 4000, 8000), reps = 100, beta = 0.40, seed = NULL, test_rep = 0, parallel_inner = FALSE, progress = TRUE) {
  design <- tibble(n_input = n_values, design_id = seq_along(n_values))
  map_dfr(seq_len(nrow(design)), function(i) {
    n_i <- design$n_input[i]
    cell_seed <- if (is.null(seed)) NULL else seed + i * 100000
    if (progress) message("Latent-only cell ", i, " / ", nrow(design), ": n = ", n_i)
    ids <- seq_len(reps)
    if (!parallel_inner) {
      map_dfr(ids, function(r) run_latent_n_replication(r, n_i, beta, test_rep, cell_seed))
    } else {
      furrr::future_map_dfr(ids, function(r) run_latent_n_replication(r, n_i, beta, test_rep, cell_seed), .options = furrr::furrr_options(seed = TRUE))
    }
  })
}

summarise_latent_n_only <- function(results) {
  results %>%
    group_by(n_input, beta) %>%
    summarise(n_replications = n(), mean_effect_size = mean(effect_size, na.rm = TRUE), sd_effect_size = sd(effect_size, na.rm = TRUE), se_effect_size = sd_effect_size / sqrt(n_replications), ci_low = mean_effect_size - 1.96 * se_effect_size, ci_high = mean_effect_size + 1.96 * se_effect_size, mean_correlation = mean(correlation, na.rm = TRUE), mean_x_skew = mean(x_skew_obs, na.rm = TRUE), mean_y_skew = mean(y_skew_obs, na.rm = TRUE), mean_x_kurtosis = mean(x_kurtosis_obs, na.rm = TRUE), mean_y_kurtosis = mean(y_kurtosis_obs, na.rm = TRUE), prop_significant = mean(significant %in% TRUE, na.rm = TRUE), .groups = "drop") %>%
    arrange(n_input)
}

# ----------------------------- #
# Example server run ----------- #
# ----------------------------- #

future::plan(future::multisession, workers = 36)

res_grid <- simulate_grid(
  n_values = c(250, 500, 1000),
  beta_values = c(0.20, 0.40, 0.60),
  reps = 10,
  test_rep = 1000,
  seed = 20260428,
  parallel_inner = TRUE,
  progress = TRUE,
  results_csv = "simulation_results_fast.csv",
  summary_csv = "results_table_fast.csv"
)

# asd<-as.data.frame(summarise_results(res_grid))
# asd$prop_significant


# save_results(
#   res_grid,
#   results_csv = "simulations/simulation_results_fast.csv",
#   summary_csv = "simulations/results_table_fast.csv"
# )
#
# res_latent_n <- simulate_latent_n_only(
#   n_values = c(250, 500, 1000, 2000, 4000, 6000, 8000, 10000),
#   reps = 1000,
#   beta = 0.40,
#   test_rep = 0,
#   seed = 20260428,
#   parallel_inner = TRUE
# )
#
# write.csv(res_latent_n, "simulations/simulation_results_n_fast.csv", row.names = FALSE)
# write.csv(summarise_latent_n_only(res_latent_n), "simulations/results_table_n_fast.csv", row.names = FALSE)
