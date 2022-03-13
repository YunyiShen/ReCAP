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



living_matrix <- lapply(1:2000, function(i, living_mcmc){
  matrix(living_mcmc[i,], nrow = 11)
}, Chicago_RES$mcmc.objs$living.mcmc)


survival_matrix <- lapply(1:2000, function(i, surv_mcmc){
  matrix(surv_mcmc[i,], nrow = 6, byrow = T)
}, Chicago_RES$mcmc.objs$survival.mcmc)

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

fawns <- generate_plot_data(fawns, "fawn")
females <- generate_plot_data(females, "female")
males <- generate_plot_data(males, "male")
total <- generate_plot_data(all_inid, "all")

aerial_det <- generate_plot_data(t(sigmoid(Chicago_RES$mcmc.objs$aerial.detection.mcmc)),"aerial_det", year = 1:17+1991)


plot_data <- rbind(fawns, females, males)
write.csv(plot_data, "post_cull_population.csv", row.names = F)
write.csv(total,"total_post_cull.csv", row.names = F)
write.csv(aerial_det,"aerial_detection.csv", row.names = F)

library(ggplot2)
ggplot(plot_data, aes(x = year, y = mean, shape = group, lty = group))+
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),size = .5, width = 0.3)+
  xlab("Year") + 
  ylab("Reconstructed post-cull population")+
  theme_classic() 
ggsave("postcull_population.pdf", width = 6, height = 3.5, scale = .9)

library(ggplot2)
ggplot(total, aes(x = year, y = mean))+
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),size = .5, width = 0.3)+
  xlab("Year") + 
  ylab("Reconstructed post-cull population")+
  theme_classic() 
ggsave("postcull_all.pdf", width = 6, height = 3.5, scale = .9)


ggplot(aerial_det, aes(x = year, y = mean))+
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),size = .5, width = 0.3)+
  xlab("Year") + 
  ylab("Estimated aerial detection rate")+
  theme_classic() 
ggsave("aerial_det.pdf", width = 6, height = 3.5, scale = .9)


# survival
survival_matrix <- lapply(1:2000, function(i, surv_mcmc){
  t(matrix(surv_mcmc[i,], nrow = 16)) # we stored by row, check the initial value
}, Chicago_RES$mcmc.objs$survival.mcmc)

age_names <- (expand.grid(c("fawn","yearling","adult"),c("female","male")))


plot_data <- lapply(1:6, function(i, thenames,mats){
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


write.csv(plot_data,"non-harvest-survival.csv", row.names = F)

ggplot(plot_data, aes(x = year, y = mean))+
  geom_point(size = 2) + 
  geom_line() + 
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),size = .5, width = 0.3)+
  xlab("Year") + 
  ylab("Reconstructed survival")+
  facet_grid(age~sex)+
  theme_classic() 
ggsave("survival.pdf", width = 6, height = 3.5, scale = .9)


