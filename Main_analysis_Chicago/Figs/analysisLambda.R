require(ggplot2)

output_dir <- "./monograph_figs/one_year_growth"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summarize_scenarios <- function(results) {
  scenario_rows <- c(
    "observed (w/ culling)" = 6,
    "uniform age structure" = 2,
    "stable age structure" = 3
  )
  years <- seq_len(ncol(results[[1]])) + 1992

  do.call(rbind, lapply(seq_along(scenario_rows), function(i) {
    samples <- vapply(
      results,
      function(result) result[scenario_rows[[i]], ],
      numeric(length(years))
    )

    data.frame(
      point = names(scenario_rows)[[i]],
      lambda = apply(samples, 1, median),
      low = apply(samples, 1, quantile, probs = 0.025),
      high = apply(samples, 1, quantile, probs = 0.975),
      year = years
    )
  }))
}

make_growth_plot <- function(plot_data, y_label, reference_value) {
  ggplot(plot_data, aes(x = year, y = lambda, shape = point, linetype = point)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.1) +
    geom_hline(yintercept = reference_value, colour = "gray30", linetype = 2) +
    labs(y = y_label, x = "Year") +
    theme_classic()
}

# Exclude the "skip culling" scenario from both one-year-growth outputs.
lambda_results <- ReCAP::analysisLambda(
  Chicago_RES$mcmc.objs, Assumptions, nage, 16
)
lambda_data <- summarize_scenarios(lambda_results)
lambda_plot <- make_growth_plot(lambda_data, "Lambda", 1)

write.csv(
  lambda_data,
  file.path(output_dir, "one_year_lambda.csv"),
  row.names = FALSE
)
ggsave(
  file.path(output_dir, "one_year_lambda.pdf"),
  plot = lambda_plot,
  width = 6,
  height = 3.5,
  scale = 1
)

change_results <- ReCAP::analysisRecruitment(
  Chicago_RES$mcmc.objs, Assumptions, nage, 16
)
change_data <- summarize_scenarios(change_results)
change_plot <- make_growth_plot(change_data, "Net population change", 0)

write.csv(
  change_data,
  file.path(output_dir, "one_year_change.csv"),
  row.names = FALSE
)
ggsave(
  file.path(output_dir, "one_year_change.pdf"),
  plot = change_plot,
  width = 6,
  height = 3.5,
  scale = 1
)
