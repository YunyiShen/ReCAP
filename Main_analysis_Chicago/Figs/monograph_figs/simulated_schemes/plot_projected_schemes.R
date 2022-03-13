library(ggplot2)
library(ggpubr)

all_pop_est <- list.files(pattern = "_est_")
all_goal_rate <- list.files(pattern = "goal")

pop_fig <- list()
goal_fig <- list()
for(i in 1:4){
  pop_data <- read.csv(all_pop_est[i])
  pop_fig[[i]] <- ggplot(data = pop_data, aes(x = year, y = pop_mean))+
    geom_line()+
    geom_point() +
    geom_errorbar(aes(ymin = CI_low,ymax = CI_high),width = .2) +
    xlab("Year")+
    ylab("Projected population")+
    ylim(0,1000)+
    geom_hline(yintercept = 150,col = "gray40", linetype = 2) +
    theme_classic()

  goal_data <- read.csv(all_goal_rate[i])
  goal_fig[[i]] <- ggplot(goal_data, aes(x = year, y = p_ctrl))+
    geom_line()+
    geom_point() +
    xlab("Year")+
    ylab("Posterior probability \n reaching the goal")+
    ylim(0,1)+
    theme_classic()


}

ggarrange(plotlist = pop_fig, ncol = 2, nrow = 2, labels = "AUTO")
ggsave("4_schemes_population.pdf", width = 12,height=8, scale = .7)

ggarrange(plotlist = goal_fig, ncol = 2, nrow = 2, labels = "AUTO")
ggsave("4_schemes_goal.pdf", width = 12,height=8, scale = .7)
