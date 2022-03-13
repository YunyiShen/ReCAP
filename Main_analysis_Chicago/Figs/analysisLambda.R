analyisoflambda =
  ReCAP::analysisLambda(Chicago_RES$mcmc.objs,Assumptions,nage,16)
  #ReCAP::analysisRecruitment(Chicago_RES$mcmc.objs,Assumptions,nage,16)
mean_lambda = Reduce("+",analyisoflambda)/nrow(Chicago_RES$mcmc.objs$survival.mcmc)

list_each_lambda = lapply(1:length(mean_lambda),function(i,analyisoflambda){
  temp = lapply(analyisoflambda,function(ana,i){ana[i]},i)
  Reduce(rbind,temp)
},analyisoflambda)

lower_025 = sapply(list_each_lambda,quantile,probs = .025)
lower_025 = matrix(lower_025,ncol=16)

higher_975 = sapply(list_each_lambda,quantile,probs = .975)
higher_975 = matrix(higher_975,ncol=16)

year = 1:16+1992

observed = data.frame(point = "observed (w/ culling)"
                      ,lambda=mean_lambda[6,]
                      ,low = lower_025[6,]
                      ,high = higher_975[6,]
                      ,year = year)
even = data.frame(point = "uniform age structure"
                      ,lambda=mean_lambda[2,]
                      ,low = lower_025[2,]
                      ,high = higher_975[2,]
                      ,year = year)
nocull = data.frame(point = "skip culling"
                      ,lambda=mean_lambda[4,]
                      ,low = lower_025[4,]
                      ,high = higher_975[4,]
                      ,year = year)
stable = data.frame(point = "stable age structure"
                    ,lambda=mean_lambda[3,]
                    ,low = lower_025[3,]
                    ,high = higher_975[3,]
                    ,year = year)

plot_data = rbind(observed,even,nocull,stable)

require(ggplot2)

ggplot(data = plot_data,aes(x=year,y=lambda,shape=point, lty = point))+
  geom_line()+
  geom_point() +
  geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
  labs(y = "Lambda", x = "Year")+
  #labs(y = "Net population change", x = "Year")+
  geom_hline(yintercept = 1, col = "gray30", lty = 2)+
  #theme(legend.position = "top")+
  theme_classic()

write.csv(plot_data,"one_year_lambda.csv", row.names = F)
ggsave("./one_year_lambda.pdf", plot = last_plot(), width = 6, height = 3.5, scale = 1)
