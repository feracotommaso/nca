library(dplyr)
library(tibble)
library(tidyr)
library(furrr)
library(future)
library(NCA)

ce_fdh_fast <- function(data) {
  d <- as.data.frame(data[, 1:2])
  names(d) <- c("x", "y")
  d <- stats::na.omit(d)
  
  n_used <- nrow(d)
  x_unique <- dplyr::n_distinct(d$x)
  y_unique <- dplyr::n_distinct(d$y)
  
  if (n_used < 2 || x_unique < 2 || y_unique < 2) {
    return(list(
      peers = tibble(x = numeric(0), y = numeric(0)),
      ceiling_area = NA_real_,
      scope_area = NA_real_,
      effect = NA_real_,
      n_used = n_used,
      x_unique = x_unique,
      y_unique = y_unique
    ))
  }
  
  x_min <- min(d$x)
  x_max <- max(d$x)
  y_min <- min(d$y)
  y_max <- max(d$y)
  scope_area <- (x_max - x_min) * (y_max - y_min)
  
  if (!is.finite(scope_area) || scope_area <= 0) {
    return(list(
      peers = tibble(x = numeric(0), y = numeric(0)),
      ceiling_area = NA_real_,
      scope_area = scope_area,
      effect = NA_real_,
      n_used = n_used,
      x_unique = x_unique,
      y_unique = y_unique
    ))
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
  
  list(
    peers = peers,
    ceiling_area = ceiling_area,
    scope_area = scope_area,
    effect = ceiling_area / scope_area,
    n_used = n_used,
    x_unique = x_unique,
    y_unique = y_unique
  )
}

ce_fdh_p_value <- function(data, observed_effect = NULL, B = 0) {
  if (B < 1) return(NA_real_)
  
  d <- as.data.frame(data[, 1:2])
  names(d) <- c("x", "y")
  d <- stats::na.omit(d)
  
  if (nrow(d) < 2) return(NA_real_)
  
  if (is.null(observed_effect)) {
    observed_effect <- ce_fdh_fast(d)$effect
  }
  
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

make_thresholds <- function(n_cat = 5,
                            spacing = c("central", "shifted"),
                            shift_direction = c("none", "lower", "upper"),
                            shift_magnitude = 0,
                            central_sd = 1.2) {
  spacing <- match.arg(spacing)
  shift_direction <- match.arg(shift_direction)
  
  centers <- seq(-(n_cat - 1) / 2, (n_cat - 1) / 2, length.out = n_cat)
  probs <- stats::dnorm(centers, 0, central_sd)
  probs <- probs / sum(probs)
  thresholds <- stats::qnorm(cumsum(probs)[seq_len(n_cat - 1)])
  
  if (spacing == "shifted" && shift_direction == "lower") {
    thresholds <- thresholds - shift_magnitude
  }
  if (spacing == "shifted" && shift_direction == "upper") {
    thresholds <- thresholds + shift_magnitude
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

run_simple_replication <- function(rep_id,
                                   n = 500,
                                   beta = 0.40,
                                   n_cat = 5,
                                   central_sd = 1.2,
                                   ceiling = "ce_fdh",
                                   test_rep_nca = 1000,
                                   test_rep_fast = 1000,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed + rep_id)
  
  eta_x <- as.numeric(scale(stats::rnorm(n)))
  eps <- as.numeric(scale(stats::rnorm(n))) * sqrt(1 - beta^2)
  eta_y <- beta * eta_x + eps
  
  thresholds <- make_thresholds(
    n_cat = n_cat,
    spacing = "central",
    shift_direction = "none",
    shift_magnitude = 0,
    central_sd = central_sd
  )
  
  X <- latent_to_ordinal(eta_x, thresholds)
  Y <- latent_to_ordinal(eta_y, thresholds)
  
  dat <- tibble(X = X, Y = Y)
  
  # fast/custom method
  fast_fit <- ce_fdh_fast(dat)
  effect_size_fast <- fast_fit$effect
  p_value_fast <- ce_fdh_p_value(
    dat,
    observed_effect = effect_size_fast,
    B = test_rep_fast
  )
  
  # NCA package
  nca_fit <- NCA::nca_analysis(
    data = as.data.frame(dat),
    x = "X",
    y = "Y",
    ceilings = ceiling,
    test.rep = test_rep_nca
  )
  
  effect_size_nca <- as.numeric(
    NCA::nca_extract(
      nca_fit,
      x = "X",
      ceiling = ceiling,
      param = "Effect size"
    )
  )
  
  p_value_nca <- as.numeric(nca_fit$tests$X[[ceiling]]$p_value)
  
  tibble(
    rep_id = rep_id,
    n_input = n,
    beta = beta,
    n_cat = n_cat,
    effect_size_fast = effect_size_fast,
    p_value_fast = p_value_fast,
    effect_size_nca = effect_size_nca,
    p_value_nca = p_value_nca
  )
}

simulate_simple <- function(reps = 100,
                            n = 500,
                            beta = 0.40,
                            n_cat = 5,
                            central_sd = 1.2,
                            ceiling = "ce_fdh",
                            test_rep_nca = 1000,
                            test_rep_fast = 1000,
                            seed = NULL,
                            parallel = FALSE,
                            progress = TRUE) {
  ids <- seq_len(reps)
  
  if (!parallel) {
    out <- vector("list", reps)
    
    for (r in ids) {
      if (progress && (r == 1 || r == reps || r %% max(1, floor(reps / 10)) == 0)) {
        message("Replication ", r, " / ", reps,
                " | n = ", n,
                " | beta = ", beta)
      }
      
      out[[r]] <- run_simple_replication(
        rep_id = r,
        n = n,
        beta = beta,
        n_cat = n_cat,
        central_sd = central_sd,
        ceiling = ceiling,
        test_rep_nca = test_rep_nca,
        test_rep_fast = test_rep_fast,
        seed = seed
      )
    }
    
    return(bind_rows(out))
  }
  
  furrr::future_map_dfr(
    ids,
    function(r) {
      run_simple_replication(
        rep_id = r,
        n = n,
        beta = beta,
        n_cat = n_cat,
        central_sd = central_sd,
        ceiling = ceiling,
        test_rep_nca = test_rep_nca,
        test_rep_fast = test_rep_fast,
        seed = seed
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
}

run_simple_replication <- function(rep_id,
                                   n = 500,
                                   beta = 0.40,
                                   n_cat = 5,
                                   central_sd = 1.2,
                                   ceiling = "ce_fdh",
                                   test_rep_nca = 1000,
                                   test_rep_fast = 1000,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed + rep_id)
  
  eta_x <- as.numeric(scale(stats::rnorm(n)))
  eps <- as.numeric(scale(stats::rnorm(n))) * sqrt(1 - beta^2)
  eta_y <- beta * eta_x + eps
  
  thresholds <- make_thresholds(
    n_cat = n_cat,
    spacing = "central",
    shift_direction = "none",
    shift_magnitude = 0,
    central_sd = central_sd
  )
  
  X <- latent_to_ordinal(eta_x, thresholds)
  Y <- latent_to_ordinal(eta_y, thresholds)
  
  dat <- tibble(X = X, Y = Y)
  
  # fast/custom method
  fast_fit <- ce_fdh_fast(dat)
  effect_size_fast <- fast_fit$effect
  p_value_fast <- ce_fdh_p_value(
    dat,
    observed_effect = effect_size_fast,
    B = test_rep_fast
  )
  
  # NCA package with fallback
  nca_fit <- tryCatch(
    NCA::nca_analysis(
      data = as.data.frame(dat),
      x = "X",
      y = "Y",
      ceilings = ceiling,
      test.rep = test_rep_nca
    ),
    error = function(e) e
  )
  
  nca_error <- NA_character_
  
  if (inherits(nca_fit, "error")) {
    nca_error <- conditionMessage(nca_fit)
    
    nca_fit <- tryCatch(
      NCA::nca_analysis(
        data = as.data.frame(dat),
        x = "X",
        y = "Y",
        ceilings = ceiling,
        test.rep = 0
      ),
      error = function(e) e
    )
  }
  
  if (inherits(nca_fit, "error")) {
    effect_size_nca <- NA_real_
    p_value_nca <- NA_real_
    nca_error <- paste0(nca_error, " | fallback failed: ", conditionMessage(nca_fit))
  } else {
    effect_size_nca <- tryCatch(
      as.numeric(
        NCA::nca_extract(
          nca_fit,
          x = "X",
          ceiling = ceiling,
          param = "Effect size"
        )
      ),
      error = function(e) NA_real_
    )
    
    p_value_nca <- tryCatch(
      as.numeric(nca_fit$tests$X[[ceiling]]$p_value),
      error = function(e) NA_real_
    )
  }
  
  tibble(
    rep_id = rep_id,
    n_input = n,
    beta = beta,
    n_cat = n_cat,
    effect_size_fast = effect_size_fast,
    p_value_fast = p_value_fast,
    effect_size_nca = effect_size_nca,
    p_value_nca = p_value_nca,
    nca_error = nca_error
  )
}

future::plan(future::multisession, workers = 36)

res_simple <- simulate_simple(
  reps = 1000,
  n = 500,
  beta = 0.40,
  n_cat = 5,
  central_sd = 1.2,
  ceiling = "ce_fdh",
  test_rep_nca = 1000,
  test_rep_fast = 1000,
  seed = 20260505,
  parallel = TRUE,
  progress = TRUE
)

res_simple

write.csv(res_simple, "simulation_results_concordance.csv", row.names = FALSE)
