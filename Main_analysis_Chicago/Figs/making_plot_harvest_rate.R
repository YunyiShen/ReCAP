require(ggplot2)
# Just for plotting things
period = 17

## Harvest prediction
mean.harv = apply(Chicago_RES$mcmc.objs$H.mcmc,2,mean)
mean.harv.matrix = matrix(mean.harv,ncol = period, byrow=T)
mean.harv.mean = data.frame(age = c(paste0("F",1:3),paste0("M",1:3)), mean.harv.matrix)

BI.low.harv = apply(Chicago_RES$mcmc.objs$H.mcmc,2,quantile,probs = .025)
BI.low.harv.matrix = matrix(BI.low.harv,ncol = period, byrow = T)
BI_harv_low = data.frame(age = c(paste0("F",1:3),paste0("M",1:3)),BI.low.harv.matrix)


BI.high.harv = apply(Chicago_RES$mcmc.objs$H.mcmc,2,quantile,probs = .975)
BI.high.harv.matrix = matrix(BI.high.harv,ncol = period, byrow = T)
BI_harv_high = data.frame(age = c(paste0("F",1:3),paste0("M",1:3)),BI.high.harv.matrix)

har_data = data.frame(matrix(nrow = 1,ncol = 5))
colnames(har_data) = c("age","mean","low","high","time")
har_data = har_data[-1,]

age_code <- c("Female.Fawn","Female.Yearling","Female.Adult","Male.Fawn","Male.Yearling","Male.Adult")

require(ggplot2)

all_summary = data.frame(harvest.rate = mean.harv.matrix[1,],low = BI.low.harv.matrix[1,],high = BI.high.harv.matrix[1,],time = 1992:2008, group = age_code[1])


for(i in 2:6){
  temp1 = data.frame(harvest.rate  = mean.harv.matrix[i,],low = BI.low.harv.matrix[i,],high = BI.high.harv.matrix[i,],time = 1992:2008, group = age_code[i])
  all_summary <- rbind(temp1,all_summary)
}

all_summary[,1:3] <- 1/(1+exp(-all_summary[,1:3]))

ggplot(data = all_summary, aes(x=time,y=harvest.rate,col = group))+
  geom_line()+
  geom_errorbar(aes(ymin = low, ymax = high),size = .2)

ggsave("./harv_porp.pdf",width = 10,height = 6, scale = .8)

write.csv(all_summary,"harvest_porp.csv",row.names = F)


## now calculate overall lambda:

lambda_MCMC <- matrix(NA, nrow = 2000, ncol = period-1)
for(i in 1:2000){
  temp <- matrix(Chicago_RES$mcmc.objs$living.mcmc[i,], ncol = period, byrow = F)
  temp <- colSums(temp)
  lambda_MCMC[i,] <- temp[-1]/temp[-period]
}


mean.lambda = apply(lambda_MCMC,2,mean)
BI.low.lambda = apply(lambda_MCMC,2,quantile,probs = .025)
BI.high.lambda = apply(lambda_MCMC,2,quantile,probs = .975)

lambda_val <- data.frame(year = 1993:2008, lambda = mean.lambda, lambda.low = BI.low.lambda, lambda.high = BI.high.lambda)
write.csv(lambda_val, "lambda_summary.csv",row.names = F)
lambda_val <- read.csv("lambda_summary.csv")

pp <- list()

for(i in 1:6){
  temp <- all_summary[all_summary$group==age_code[i],]
  temp <- temp[temp$time!=1992,]
  temp <- cbind(temp,lambda_val)
  pp_temp <- ggplot(data = temp, aes(x=harvest.rate,y=lambda))+
    geom_point()+
    geom_errorbar(aes(xmin=low,xmax = high, ymin = lambda.low,ymax = lambda.high),size = 0.2)+
    geom_errorbarh(aes(xmin=low,xmax = high),size = 0.2)+
    xlim(0,.9)+
    ylim(0.3,1.5)+
    #geom_smooth(method=lm) +
    geom_hline(yintercept = 1,col = "red", linetype = 2)#+
    #xlab(age_code[i])
  pp[[i]] <- pp_temp
}

ggpubr::ggarrange(plotlist = pp, labels =age_code, label.x = 0.1,align = "hv")

ggsave("./rate_vs_lambda.pdf",width = 12, height = 8, scale = .8)

H_all_MCMC <- matrix(NA, nrow = 2000, ncol = period)
for(i in 1:2000){
  temp_living <- matrix(Chicago_RES$mcmc.objs$living.mcmc[i,], ncol = period, byrow = F)
  temp_harving <- matrix(Chicago_RES$mcmc.objs$harvest.mcmc[i,], ncol = period, byrow = F)
  temp_living <- colSums(temp_living)
  temp_harving <- colSums(temp_harving)
  H_all_MCMC[i,] <- temp_harving/(temp_living+temp_harving)
}

mean.Hall = apply(H_all_MCMC ,2,mean)
BI.low.Hall = apply(H_all_MCMC ,2,quantile,probs = .025)
BI.high.Hall = apply(H_all_MCMC ,2,quantile,probs = .975)

Hall_val <- data.frame(year = 1992:2008, harvest.rate = mean.Hall, Hall.low = BI.low.Hall, Hall.high = BI.high.Hall)

write.csv(Hall_val,"harvest_rate_overall.csv",row.names = F)
Hall_val <- read.csv("harvest_rate_overall.csv")

hallvslambda <- cbind(Hall_val[Hall_val$year!=1992,-1], lambda_val)

hallvslambda$year <- as.factor(hallvslambda$year)
ggplot(data = hallvslambda, aes(x=harvest.rate,y=lambda))+
  geom_smooth(method=lm, size = .5, col = "gray10") +
  geom_point(aes(col = year))+
  #scale_color_grey()+
  #geom_point()+
  geom_errorbar(aes(xmin=Hall.low,xmax = Hall.high, ymin = lambda.low,ymax = lambda.high),size = 0.2)+
  geom_errorbarh(aes(xmin=Hall.low,xmax = Hall.high),size = 0.3)+
  xlim(0.12,.8)+
  ylim(0.3,1.7)+
  xlab("Harvest rate")+
  ylab("Annual growth rate (lambda)")+
  geom_hline(yintercept = 1,col = "gray30", linetype = 2) +
  theme_classic()+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
ggsave("./rate_overall_vs_lambda1.pdf", width = 8, height = 5.5, scale = 0.8)
