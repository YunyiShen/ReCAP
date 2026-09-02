generate_plot_data <- function(w, group, year = 1:17+1991){
  CI_low <- apply(w, 1, quantile, 0.025)
  CI_high <- apply(w, 1, quantile, 0.975)
  
  data.frame(
    group = group, year = year,
    mean = rowMeans(w),
    CI_low = CI_low,
    CI_high = CI_high
  )
}

sigmoid <- function(x){
  1/(exp(-x)+1)
}



n_draws <- nrow(Chicago_RES$mcmc.objs$living.mcmc)

living_matrix <- lapply(seq_len(n_draws), function(i, living_mcmc){
  matrix(living_mcmc[i,], nrow = 11)
}, Chicago_RES$mcmc.objs$living.mcmc)


all_inid <- sapply(living_matrix, function(w){
  t(colSums(w))
})

males <- sapply(living_matrix, function(w){
  t(colSums(w[10:11,]))
})

females <- sapply(living_matrix, function(w){
  t(colSums(w[2:8,]))
})

fawns <- sapply(living_matrix, function(w){
  t(colSums(w[c(1,9),]))
})

fawn_ratio <- fawns / all_inid
female_ratio <- females / all_inid
male_ratio <- males / all_inid

fawns <- generate_plot_data(fawns, "fawn")
females <- generate_plot_data(females, "female")
males <- generate_plot_data(males, "male")
total <- generate_plot_data(all_inid, "all")
ratio_data <- rbind(
  generate_plot_data(fawn_ratio, "fawn"),
  generate_plot_data(female_ratio, "female"),
  generate_plot_data(male_ratio, "male")
)

aerial_det <- generate_plot_data(t(sigmoid(Chicago_RES$mcmc.objs$aerial.detection.mcmc)),"aerial_det", year = 1:17+1991)


plot_data <- rbind(fawns, females, males)
write.csv(plot_data, "./monograph_figs/postcull_population/post_cull_population.csv", row.names = F)
write.csv(total,"./monograph_figs/postcull_population/total_post_cull.csv", row.names = F)
write.csv(ratio_data,"./monograph_figs/postcull_population/post_cull_population_ratio.csv", row.names = F)
write.csv(aerial_det,"./monograph_figs/vital_rates/aerial_detection.csv", row.names = F)

library(ggplot2)
postcull_population_plot <- ggplot(plot_data, aes(x = year, y = mean, shape = group, lty = group))+
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), linewidth = .5, width = 0.3)+
  labs(x = "Year", y = "Post-cull population", shape = "Age group", linetype = "Age group")+
  theme_classic() 
ggsave("./monograph_figs/postcull_population/postcull_population.pdf", plot = postcull_population_plot, width = 6, height = 3.5, scale = .9)

library(ggplot2)
postcull_total_plot <- ggplot(total, aes(x = year, y = mean))+
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), linewidth = .5, width = 0.3)+
  xlab("Year") + 
  ylab("Post-cull population")+
  theme_classic() 
ggsave("./monograph_figs/postcull_population/postcull_all.pdf", plot = postcull_total_plot, width = 6, height = 3.5, scale = .9)

postcull_ratio_plot <- ggplot(ratio_data, aes(x = year, y = mean, shape = group, lty = group))+
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), linewidth = .5, width = 0.3)+
  labs(x = "Year", y = "Post-cull sex-age structure", shape = "Age group", linetype = "Age group")+
  theme_classic()
ggsave("./monograph_figs/postcull_population/postcull_population_ratio.pdf", plot = postcull_ratio_plot, width = 6, height = 3.5, scale = .9)


ggplot(aerial_det, aes(x = year, y = mean))+
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), linewidth = .5, width = 0.3)+
  xlab("Year") + 
  ylab("Estimated aerial detection rate")+
  theme_classic() 
ggsave("./monograph_figs/vital_rates/aerial_det.pdf", width = 6, height = 3.5, scale = .9)


# survival
survival_matrix <- lapply(seq_len(n_draws), function(i, surv_mcmc){
  t(matrix(surv_mcmc[i,], nrow = 16)) # we stored by row, check the initial value
}, Chicago_RES$mcmc.objs$survival.mcmc)

age_names <- (expand.grid(c("fawn","yearling","adult"),c("female","male")))


plot_data <- lapply(1:4, function(i, thenames,mats){
  name_temp <- paste0(as.character(thenames[i,]),collapse = "-")
  temp <- sapply(mats, function(w){
    t(sigmoid(w[i,]))
  })|>
  generate_plot_data(name_temp,year = 1:16+1991)
  temp$age <- thenames[i,1]
  temp$sex <- thenames[i,2]
  return(temp[,-1])
}, age_names, survival_matrix) |>
  Reduce(f = rbind)


write.csv(plot_data,"./monograph_figs/vital_rates/non-harvest-survival.csv", row.names = F)

ggplot(plot_data, aes(x = year, y = mean))+
  geom_point(size = 2) + 
  geom_line() + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), linewidth = .5, width = 0.3)+
  xlab("Year") + 
  ylab("Reconstructed survival")+
  facet_grid(age~sex)+
  theme_classic() 
ggsave("./monograph_figs/vital_rates/survival.pdf", width = 6, height = 3.5, scale = .9)
