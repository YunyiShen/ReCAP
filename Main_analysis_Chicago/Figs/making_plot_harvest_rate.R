require(ggplot2)

output_dir <- "./monograph_figs/harvest_vs_lambda"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

period <- 17
n_draws <- nrow(Chicago_RES$mcmc.objs$living.mcmc)
age_code <- c(
  "Female.Fawn", "Female.Yearling", "Female.Adult",
  "Male.Fawn", "Male.Yearling", "Male.Adult"
)
age_display_order <- c(
  "Female.Fawn", "Male.Fawn",
  "Female.Yearling", "Male.Yearling",
  "Female.Adult", "Male.Adult"
)

## Age-specific harvest rates
harvest_rate_draws <- plogis(as.matrix(Chicago_RES$mcmc.objs$H.mcmc))
median_harvest <- matrix(
  apply(harvest_rate_draws, 2, median),
  ncol = period,
  byrow = TRUE
)
low_harvest <- matrix(
  apply(harvest_rate_draws, 2, quantile, probs = 0.025),
  ncol = period,
  byrow = TRUE
)
high_harvest <- matrix(
  apply(harvest_rate_draws, 2, quantile, probs = 0.975),
  ncol = period,
  byrow = TRUE
)

all_summary <- do.call(rbind, lapply(seq_along(age_code), function(i) {
  data.frame(
    harvest.rate = median_harvest[i, ],
    low = low_harvest[i, ],
    high = high_harvest[i, ],
    Year = 1992:2008,
    group = age_code[i]
  )
}))
all_summary$group <- factor(all_summary$group, levels = age_display_order)
all_summary <- all_summary[order(all_summary$group, all_summary$Year), ]

harvest_plot <- ggplot(
  all_summary,
  aes(x = Year, y = harvest.rate, colour = group)
) +
  geom_line() +
  geom_errorbar(aes(ymin = low, ymax = high), linewidth = 0.2) +
  labs(x = "Year", y = "Harvest rate", colour = "Age group") +
  theme_classic()

ggsave(
  file.path(output_dir, "harv_porp.pdf"),
  plot = harvest_plot,
  width = 10,
  height = 6,
  scale = 0.8
)
write.csv(
  all_summary,
  file.path(output_dir, "harvest_porp.csv"),
  row.names = FALSE
)

## Annual population growth rate
lambda_mcmc <- matrix(NA_real_, nrow = n_draws, ncol = period - 1)
for (i in seq_len(n_draws)) {
  living <- matrix(
    Chicago_RES$mcmc.objs$living.mcmc[i, ],
    ncol = period,
    byrow = FALSE
  )
  total_living <- colSums(living)
  lambda_mcmc[i, ] <- total_living[-1] / total_living[-period]
}

lambda_val <- data.frame(
  Year = 1993:2008,
  lambda = apply(lambda_mcmc, 2, median),
  lambda.low = apply(lambda_mcmc, 2, quantile, probs = 0.025),
  lambda.high = apply(lambda_mcmc, 2, quantile, probs = 0.975)
)
write.csv(
  lambda_val,
  file.path(output_dir, "lambda_summary.csv"),
  row.names = FALSE
)

rate_lambda_data <- merge(
  all_summary[all_summary$Year != 1992, ],
  lambda_val,
  by = "Year",
  sort = TRUE
)
rate_lambda_data$group <- factor(
  rate_lambda_data$group,
  levels = age_display_order
)

fit_stability_threshold <- function(data) {
  trend_model <- lm(lambda ~ harvest.rate, data = data)
  model_coefficients <- coef(trend_model)
  intercept <- unname(model_coefficients[["(Intercept)"]])
  slope <- unname(model_coefficients[["harvest.rate"]])
  threshold <- if (is.finite(slope) && abs(slope) > .Machine$double.eps) {
    (1 - intercept) / slope
  } else {
    NA_real_
  }
  observed_min <- min(data$harvest.rate)
  observed_max <- max(data$harvest.rate)

  data.frame(
    intercept = intercept,
    slope = slope,
    harvest_rate_at_lambda_1 = threshold,
    observed_harvest_rate_min = observed_min,
    observed_harvest_rate_max = observed_max,
    crossing_within_observed_range = is.finite(threshold) &&
      threshold >= observed_min && threshold <= observed_max,
    r_squared = summary(trend_model)$r.squared,
    n_years = nrow(data)
  )
}

age_stability_thresholds <- do.call(
  rbind,
  lapply(age_display_order, function(age_group) {
    threshold <- fit_stability_threshold(
      rate_lambda_data[rate_lambda_data$group == age_group, ]
    )
    data.frame(age_group = age_group, threshold, check.names = FALSE)
  })
)
write.csv(
  age_stability_thresholds,
  file.path(output_dir, "rate_vs_lambda_stability_thresholds.csv"),
  row.names = FALSE
)

rate_lambda_plot <- ggplot(
  rate_lambda_data,
  aes(x = harvest.rate, y = lambda)
) +
  geom_smooth(
    method = lm,
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.6,
    colour = "gray10"
  ) +
  geom_point() +
  geom_errorbar(
    aes(ymin = lambda.low, ymax = lambda.high),
    width = 0.01,
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(xmin = low, xmax = high),
    orientation = "y",
    width = 0.01,
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 1, colour = "gray30", linetype = 2) +
  facet_wrap(~group, ncol = 2) +
  coord_cartesian(xlim = c(0, 0.9), ylim = c(0.3, 1.5)) +
  labs(x = "Harvest rate", y = "Annual growth rate (lambda)") +
  theme_classic()

ggsave(
  file.path(output_dir, "rate_vs_lambda.pdf"),
  plot = rate_lambda_plot,
  width = 12,
  height = 8,
  scale = 0.8
)

## Overall harvest rate
overall_harvest_mcmc <- matrix(NA_real_, nrow = n_draws, ncol = period)
for (i in seq_len(n_draws)) {
  living <- matrix(
    Chicago_RES$mcmc.objs$living.mcmc[i, ],
    ncol = period,
    byrow = FALSE
  )
  harvested <- matrix(
    Chicago_RES$mcmc.objs$harvest.mcmc[i, ],
    ncol = period,
    byrow = FALSE
  )
  total_living <- colSums(living)
  total_harvested <- colSums(harvested)
  overall_harvest_mcmc[i, ] <-
    total_harvested / (total_living + total_harvested)
}

overall_harvest <- data.frame(
  Year = 1992:2008,
  harvest.rate = apply(overall_harvest_mcmc, 2, median),
  Hall.low = apply(overall_harvest_mcmc, 2, quantile, probs = 0.025),
  Hall.high = apply(overall_harvest_mcmc, 2, quantile, probs = 0.975)
)
write.csv(
  overall_harvest,
  file.path(output_dir, "harvest_rate_overall.csv"),
  row.names = FALSE
)

overall_lambda_data <- merge(
  overall_harvest[overall_harvest$Year != 1992, ],
  lambda_val,
  by = "Year",
  sort = TRUE
)
overall_lambda_data$Year <- factor(overall_lambda_data$Year)

overall_stability_threshold <- fit_stability_threshold(overall_lambda_data)
write.csv(
  overall_stability_threshold,
  file.path(output_dir, "rate_overall_vs_lambda1_stability_threshold.csv"),
  row.names = FALSE
)

overall_lambda_plot <- ggplot(
  overall_lambda_data,
  aes(x = harvest.rate, y = lambda)
) +
  geom_smooth(
    method = lm,
    formula = y ~ x,
    linewidth = 0.6,
    colour = "gray10",
    se = FALSE
  ) +
  geom_point(aes(colour = Year)) +
  geom_errorbar(
    aes(ymin = lambda.low, ymax = lambda.high),
    width = 0.01,
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(xmin = Hall.low, xmax = Hall.high),
    orientation = "y",
    width = 0.01,
    linewidth = 0.3
  ) +
  geom_hline(yintercept = 1, colour = "gray30", linetype = 2) +
  coord_cartesian(xlim = c(0.12, 0.8), ylim = c(0.3, 1.7)) +
  labs(
    x = "Harvest rate",
    y = "Annual growth rate (lambda)",
    colour = "Year"
  ) +
  theme_classic()

ggsave(
  file.path(output_dir, "rate_overall_vs_lambda1.pdf"),
  plot = overall_lambda_plot,
  width = 8,
  height = 5.5,
  scale = 0.8
)
