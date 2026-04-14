library(dplyr)
library(ggplot2)

res_grid %>%
  {if ("sample_size" %in% names(.)) mutate(., sample_size_plot = sample_size) else mutate(., sample_size_plot = n_input)} %>%
  {if ("design_beta" %in% names(.)) mutate(., beta_plot = design_beta) else mutate(., beta_plot = beta)} %>%
  mutate(
    condition_plot = case_when(
      block == "1A_score_construction" ~ score_type,
      block == "1B_reliability" ~ paste0("loading=", condition_value),
      block == "1C_response_categories" ~ paste0("cats=", condition_value),
      block == "1D_threshold_shift" ~ as.character(condition_value),
      block == "1E_range_restriction" ~ as.character(condition_value),
      TRUE ~ as.character(condition_value)
    ),
    block = factor(
      block,
      levels = c(
        "1A_score_construction",
        "1B_reliability",
        "1C_response_categories",
        "1D_threshold_shift",
        "1E_range_restriction"
      ),
      labels = c(
        "1A Score construction",
        "1B Reliability",
        "1C Response categories",
        "1D Threshold shift",
        "1E Range restriction"
      )
    )
  ) %>%
  group_by(block, condition_plot, sample_size_plot, beta_plot) %>%
  summarise(
    mean_effect = mean(effect_size, na.rm = TRUE),
    se_effect = sd(effect_size, na.rm = TRUE) / sqrt(sum(!is.na(effect_size))),
    ci_low = mean_effect - 1.96 * se_effect,
    ci_high = mean_effect + 1.96 * se_effect,
    .groups = "drop"
  ) %>%
  ggplot(
    aes(
      x = condition_plot,
      y = mean_effect,
      color = factor(beta_plot),
      shape = factor(sample_size_plot)
    )
  ) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 2.8
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    position = position_dodge(width = 0.5),
    width = 0.15
  ) +
  facet_wrap(~ block, scales = "free_x", ncol = 1) +
  labs(
    x = NULL,
    y = "Mean NCA effect size",
    color = expression(beta),
    shape = "Sample size"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )



library(dplyr)
library(ggplot2)

res_grid %>%
  {if ("sample_size" %in% names(.)) mutate(., sample_size_plot = sample_size) else mutate(., sample_size_plot = n_input)} %>%
  {if ("design_beta" %in% names(.)) mutate(., beta_plot = design_beta) else mutate(., beta_plot = beta)} %>%
  mutate(
    condition_plot = case_when(
      block == "1A_score_construction" ~ score_type,
      block == "1B_reliability" ~ paste0("loading=", condition_value),
      block == "1C_response_categories" ~ paste0("cats=", condition_value),
      block == "1D_threshold_shift" ~ as.character(condition_value),
      block == "1E_range_restriction" ~ as.character(condition_value),
      TRUE ~ as.character(condition_value)
    ),
    block = factor(
      block,
      levels = c(
        "1A_score_construction",
        "1B_reliability",
        "1C_response_categories",
        "1D_threshold_shift",
        "1E_range_restriction"
      ),
      labels = c(
        "1A Score construction",
        "1B Reliability",
        "1C Response categories",
        "1D Threshold shift",
        "1E Range restriction"
      )
    ),
    beta_plot = factor(beta_plot, levels = sort(unique(beta_plot))),
    sample_size_plot = factor(sample_size_plot, levels = c(250, 500, 1000))
  ) %>%
  group_by(block, beta_plot, condition_plot, sample_size_plot) %>%
  summarise(
    mean_effect = mean(effect_size, na.rm = TRUE),
    se_effect = sd(effect_size, na.rm = TRUE) / sqrt(sum(!is.na(effect_size))),
    ci_low = mean_effect - 1.96 * se_effect,
    ci_high = mean_effect + 1.96 * se_effect,
    .groups = "drop"
  ) %>%
  ggplot(
    aes(
      x = condition_plot,
      y = mean_effect,
      color = sample_size_plot,
      group = sample_size_plot
    )
  ) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 2.8
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    position = position_dodge(width = 0.5),
    width = 0.15
  ) +
  facet_grid(block ~ beta_plot, scales = "free_x", labeller = label_both) +
  labs(
    x = NULL,
    y = "Mean NCA effect size",
    color = "Sample size"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )



library(dplyr)
library(ggplot2)
library(patchwork)

dat_plot <- res_grid %>%
  {if ("sample_size" %in% names(.)) mutate(., sample_size_plot = sample_size) else mutate(., sample_size_plot = n_input)} %>%
  {if ("design_beta" %in% names(.)) mutate(., beta_plot = design_beta) else mutate(., beta_plot = beta)} %>%
  mutate(
    condition_plot = case_when(
      block == "1A_score_construction" ~ score_type,
      block == "1B_reliability" ~ paste0("loading=", condition_value),
      block == "1C_response_categories" ~ paste0("cats=", condition_value),
      block == "1D_threshold_shift" ~ as.character(condition_value),
      block == "1E_range_restriction" ~ as.character(condition_value),
      TRUE ~ as.character(condition_value)
    ),
    block = factor(
      block,
      levels = c(
        "1A_score_construction",
        "1B_reliability",
        "1C_response_categories",
        "1D_threshold_shift",
        "1E_range_restriction"
      ),
      labels = c(
        "1A Score construction",
        "1B Reliability",
        "1C Response categories",
        "1D Threshold shift",
        "1E Range restriction"
      )
    ),
    beta_plot = factor(beta_plot, levels = sort(unique(beta_plot))),
    sample_size_plot = factor(sample_size_plot, levels = c(250, 500, 1000))
  ) %>%
  group_by(block, condition_plot, sample_size_plot, beta_plot) %>%
  summarise(
    mean_effect = mean(effect_size, na.rm = TRUE),
    se_effect = sd(effect_size, na.rm = TRUE) / sqrt(sum(!is.na(effect_size))),
    ci_low = mean_effect - 1.96 * se_effect,
    ci_high = mean_effect + 1.96 * se_effect,
    .groups = "drop"
  )

p1 <- dat_plot %>%
  filter(block == "1A Score construction") %>%
  ggplot(aes(x = condition_plot, y = mean_effect,
             color = beta_plot, shape = sample_size_plot,
             group = interaction(beta_plot, sample_size_plot))) +
  geom_point(position = position_dodge(width = 0.5), size = 2.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.5), width = 0.15) +
  labs(title = "1A Score construction", x = NULL, y = "Mean NCA effect size",
       color = expression(beta), shape = "Sample size") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

p2 <- dat_plot %>%
  filter(block == "1B Reliability") %>%
  ggplot(aes(x = condition_plot, y = mean_effect,
             color = beta_plot, shape = sample_size_plot,
             group = interaction(beta_plot, sample_size_plot))) +
  geom_point(position = position_dodge(width = 0.5), size = 2.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.5), width = 0.15) +
  labs(title = "1B Reliability", x = NULL, y = "Mean NCA effect size",
       color = expression(beta), shape = "Sample size") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

p3 <- dat_plot %>%
  filter(block == "1C Response categories") %>%
  ggplot(aes(x = condition_plot, y = mean_effect,
             color = beta_plot, shape = sample_size_plot,
             group = interaction(beta_plot, sample_size_plot))) +
  geom_point(position = position_dodge(width = 0.5), size = 2.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.5), width = 0.15) +
  labs(title = "1C Response categories", x = NULL, y = "Mean NCA effect size",
       color = expression(beta), shape = "Sample size") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

p4 <- dat_plot %>%
  filter(block == "1D Threshold shift") %>%
  ggplot(aes(x = condition_plot, y = mean_effect,
             color = beta_plot, shape = sample_size_plot,
             group = interaction(beta_plot, sample_size_plot))) +
  geom_point(position = position_dodge(width = 0.5), size = 2.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.5), width = 0.15) +
  labs(title = "1D Threshold shift", x = NULL, y = "Mean NCA effect size",
       color = expression(beta), shape = "Sample size") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

p5 <- dat_plot %>%
  filter(block == "1E Range restriction") %>%
  ggplot(aes(x = condition_plot, y = mean_effect,
             color = beta_plot, shape = sample_size_plot,
             group = interaction(beta_plot, sample_size_plot))) +
  geom_point(position = position_dodge(width = 0.5), size = 2.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.5), width = 0.15) +
  labs(title = "1E Range restriction", x = NULL, y = "Mean NCA effect size",
       color = expression(beta), shape = "Sample size") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

p1 / p2 / p3 / p4 / p5
