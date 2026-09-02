require(ggplot2)

output_dir <- "./monograph_figs/harvest_vs_change"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

summarize_change <- function(samples, group) {
  data.frame(
    group = group,
    Year = 1993:2008,
    pop_median = apply(samples, 1, median),
    pop_CI_low = apply(samples, 1, quantile, probs = 0.025),
    pop_CI_high = apply(samples, 1, quantile, probs = 0.975)
  )
}

n_draws <- nrow(Chicago_RES$mcmc.objs$living.mcmc)
living_matrices <- lapply(seq_len(n_draws), function(i) {
  matrix(Chicago_RES$mcmc.objs$living.mcmc[i, ], nrow = sum(nage))
})

annual_change <- function(rows) {
  sapply(living_matrices, function(living) {
    totals <- colSums(living[rows, , drop = FALSE])
    totals[-1] - totals[-length(totals)]
  })
}

fawns <- summarize_change(annual_change(c(1, 9)), "fawn")
females <- summarize_change(annual_change(2:8), "female")
males <- summarize_change(annual_change(10:11), "male")
overall <- summarize_change(annual_change(seq_len(sum(nage))), "overall")

changes <- list(
  fawns = fawns,
  females = females,
  males = males,
  overall = overall
)

for (name in names(changes)) {
  write.csv(
    changes[[name]],
    file.path(output_dir, paste0(sub("s$", "", name), "_changes.csv")),
    row.names = FALSE
  )
}

make_change_time_plot <- function(data, y_label) {
  ggplot(data, aes(x = Year, y = pop_median)) +
    geom_line() +
    geom_point() +
    geom_errorbar(
      aes(ymin = pop_CI_low, ymax = pop_CI_high),
      width = 0.3,
      linewidth = 0.4
    ) +
    geom_hline(yintercept = 0, colour = "gray30", linetype = 2) +
    labs(x = "Year", y = y_label) +
    theme_classic()
}

female_change_plot <- make_change_time_plot(
  females,
  "Net change of females"
)
overall_change_plot <- make_change_time_plot(
  overall,
  "Net change of overall population"
)

ggsave(
  file.path(output_dir, "female_changes.pdf"),
  plot = female_change_plot,
  width = 8,
  height = 5.5,
  scale = 0.8
)
ggsave(
  file.path(output_dir, "overall_changes.pdf"),
  plot = overall_change_plot,
  width = 8,
  height = 5.5,
  scale = 0.8
)

harvest_overall <- read.csv(
  "./monograph_figs/harvest_vs_lambda/harvest_rate_overall.csv",
  check.names = FALSE
)

change_labels <- c(
  fawns = "fawns",
  females = "females",
  males = "males",
  overall = "overall population"
)
harvest_change_data <- do.call(rbind, lapply(names(changes), function(name) {
  result <- merge(
    harvest_overall[harvest_overall$Year != 1992, ],
    changes[[name]],
    by = "Year",
    sort = TRUE
  )
  result$population_group <- change_labels[[name]]
  result
}))
harvest_change_data$population_group <- factor(
  harvest_change_data$population_group,
  levels = unname(change_labels)
)
harvest_change_data$Year <- factor(harvest_change_data$Year)

harvest_change_plot <- ggplot(
  harvest_change_data,
  aes(x = harvest.rate, y = pop_median)
) +
  geom_smooth(method = lm, linewidth = 0.5, colour = "gray10", se = TRUE) +
  geom_point(aes(colour = Year)) +
  geom_errorbar(
    aes(ymin = pop_CI_low, ymax = pop_CI_high),
    width = 0.01,
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(xmin = Hall.low, xmax = Hall.high),
    orientation = "y",
    width = 0.01,
    linewidth = 0.3
  ) +
  geom_hline(yintercept = 0, colour = "gray30", linetype = 2) +
  facet_wrap(~population_group, nrow = 2, ncol = 2) +
  coord_cartesian(xlim = c(0.12, 0.8)) +
  labs(
    x = "Harvest rate",
    y = "Net population change",
    colour = "Year"
  ) +
  theme_classic()

ggsave(
  file.path(output_dir, "harvest_rate_vs_change.pdf"),
  plot = harvest_change_plot,
  width = 15,
  height = 9,
  scale = 0.8
)

overall_harvest_change <- harvest_change_data[
  harvest_change_data$population_group == "overall population",
]
overall_harvest_change$population_group <- NULL

overall_harvest_change_plot <- ggplot(
  overall_harvest_change,
  aes(x = harvest.rate, y = pop_median)
) +
  geom_smooth(method = lm, linewidth = 0.5, colour = "gray10", se = TRUE) +
  geom_point(aes(colour = Year)) +
  geom_errorbar(
    aes(ymin = pop_CI_low, ymax = pop_CI_high),
    width = 0.01,
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(xmin = Hall.low, xmax = Hall.high),
    orientation = "y",
    width = 0.01,
    linewidth = 0.3
  ) +
  geom_hline(yintercept = 0, colour = "gray30", linetype = 2) +
  coord_cartesian(xlim = c(0.12, 0.8)) +
  labs(
    x = "Harvest rate",
    y = "Net change of overall population",
    colour = "Year"
  ) +
  theme_classic()

ggsave(
  file.path(output_dir, "harvest_rate_vs_overallchange.pdf"),
  plot = overall_harvest_change_plot,
  width = 8,
  height = 5.5,
  scale = 0.8
)
