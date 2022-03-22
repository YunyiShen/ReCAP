generate_plot_data <- function(w, group, year = 1:17+1991){
  CI_low <- apply(w, 1, quantile, 0.025)
  CI_high <- apply(w, 1, quantile, 0.975)

  data.frame(
    group = group, year = year,
    pop_mean = rowMeans(w),
    pop_CI_low = CI_low,
    pop_CI_high = CI_high
  )
}

living_matrix <- lapply(1:2000, function(i, living_mcmc){
  matrix(living_mcmc[i,], nrow = 11)
}, Chicago_RES$mcmc.objs$living.mcmc)

males <- sapply(living_matrix, function(w){
  temp <- t(colSums(w[10:11,]))
  temp[-1]-temp[-length(temp)]
})

females <- sapply(living_matrix, function(w){
  temp <- t(colSums(w[2:8,]))
  temp[-1]-temp[-length(temp)]
})

fawns <- sapply(living_matrix, function(w){
  temp <- t(colSums(w[c(1,9),]))
  temp[-1]-temp[-length(temp)]
})

overall <- sapply(living_matrix, function(w){
  temp <- t(colSums(w))
  temp[-1]-temp[-length(temp)]
})


fawns <- generate_plot_data(fawns, "fawn",year = 2:17+1991)
females <- generate_plot_data(females, "female",year = 2:17+1991)
males <- generate_plot_data(males, "male",year = 2:17+1991)
overall <- generate_plot_data(overall, "overall",year = 2:17+1991)

changes <- list(fawns, females, males, overall)
names(changes) <- c("fawns", "females", "males","all")
write.csv(fawns,"fawn_changes.csv", row.names = F)
write.csv(females,"female_changes.csv", row.names = F)
write.csv(males,"male_changes.csv", row.names = F)
write.csv(overall,"overall_changes.csv", row.names = F)

harvest_overall <- read.csv("./monograph_figs/harvest_vs_lambda/harvest_rate_overall.csv")[-1,]

library(ggplot2)
library(ggpubr)
plotsss <- lapply(1:4, function(i,changes, harvest_overall){
  popchange <- changes[[i]]
  thename <- names(changes)[i]
  hallvschange <- cbind(harvest_overall, popchange[,-c(1:2)])
  hallvschange$year <- as.factor(hallvschange$year)
  ggplot(data = hallvschange, aes(x=harvest.rate,y=pop_mean))+
    geom_smooth(method=lm, size = .5, col = "gray10") +
    geom_point(aes(col = year))+
    #scale_color_grey()+
    #geom_point()+
    geom_errorbar(aes(xmin=Hall.low,xmax = Hall.high, ymin = pop_CI_low,ymax = pop_CI_high),size = 0.2)+
    geom_errorbarh(aes(xmin=Hall.low,xmax = Hall.high),size = 0.3)+
    xlim(0.12,.8)+
    #ylim(0.3,1.7)+
    xlab("Harvest rate")+
    ylab(paste("Net change of",thename))+
    geom_hline(yintercept = 1,col = "gray30", linetype = 2) +
    theme_classic()+
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
}, changes, harvest_overall)

ggarrange(plotlist = plotsss, nrow = 2,ncol = 2, labels = "AUTO", legend = "right", common.legend = T)
ggsave("./harvest_rate_vs_change.pdf", width = 15, height = 9, scale = 0.8)

plotsss[[4]]
ggsave("./harvest_rate_vs_overallchange.pdf", width = 8, height = 5.5, scale = 0.8)
