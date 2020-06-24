require(Rcpp)
require(RcppArmadillo)
require(coda)
sourceCpp("./src/Projections.cpp")
source("./R/misc.R")
source("./R/ReCAP.R")
nage = matrix( c(8,3),2,1) # nage is female first and then male, a vector with lenght usually 2
period = 16


mean.s = read.csv("./_data_/Survival_mean_Etter.csv",row.names = 1)#[c(1:3,9:11),]
mean.f = read.csv("./_data_/Fecundity_mean.csv",row.names = 1)[1:3,]
mean.SRB = read.csv("./_data_/SRB_mean.csv",row.names = 1)
Harv.data = read.csv("./_data_/Culling_1992_2008.csv",row.names = 1)
Aeri.data = read.csv("./_data_/Aerial_count_1992_2008.csv",row.names = 1)
#mean.H = read.csv("./_data_/Harvest_rate_prior.csv",row.names = 1)
Harv_assump = read.csv("./_data_/Harv_assump.csv",header = F)
Harv_assump = as.matrix(Harv_assump) # this is the assumption matrix for specific harvest!
#Surv_assump = read.csv("./_data_/Surv_assump_age.csv",header = F)

Assumptions = list()
Assumptions$Harv = list(time = eyes(period+1),age = Harv_assump) # tons of assumptions on vital rates


Assumptions$err = list(time = eyes(period))


prior.mean = list(Fec = c(0.5,1,2),Surv = 0.6, SRB = 0.5, Harv = 0.5, AerialDet = 0.7)
prior.var = list(Fec = c(.2,.2,.2),Surv = 2, SRB = .1, Harv = .8, AerialDet = .8)

FYA = matrix(0,8,3)
FYA[1,1] = 1
FYA[2,2] = 1
FYA[3:8,3] = 1
prior.ageclass = list(Fec = FYA)
Observations = list(Fec = FYA)

FYA_sex = matrix(0,11,6)
FYA_sex[9:11,4:6] = eyes(3)
FYA_sex[1:3,1:3] = eyes(3)
FYA_sex[3:8,3] = 1


Assumptions$Fec = list(time = eyes(period),age = eyes(8))
Assumptions$err  =list(time = matrix(1,1,period))
Assumptions$Surv = list(time = eyes(period),age = FYA_sex)
#Assumptions$Surv = list(time = eyes(period),age = eyes(11))


# full matrix for e.g. Harvest will be:
#  Assumptions$Harv$age %*% as.matrix(mean.H) %*% Assumptions$Harv$time
#  It is a good idea to try the command above to see how to use assumption matrices.

mean.A = matrix(0.7,1,period+1)
mean.aK0 = list(matrix(0,nage[1],1),matrix(0,sum(nage),1),matrix(10,1,1))
prop.vars = list(Fec = matrix(1,nrow = nage[1],ncol = period),
                 Surv = matrix(1,nrow = sum(nage), ncol = period),
                 SRB = matrix(0.1,nage[1],period), # vital rates has period cols
                 AerialDet = matrix(1,1,period+1),
                 Harv = matrix(1,nrow = 4,ncol = period+1),
                 aK0=list(5e-8,5e-8,50),
                 baseline.pop.count = matrix(.1,nrow = sum(nage),ncol = 1))

prior.measurement.error = list(Alpha = list(Fec = 3,Surv = 3,SRB = 3),Beta = list(Fec=1,Surv = 1,SRB=.5))

ProjectAllCpp(as.matrix(mean.s),mean.H.full,as.matrix(mean.f), as.matrix(mean.SRB),mean.aK0,T,F,as.matrix(mean.b),period,c(8,3))


set.seed(42)

Chicago_RES = ReCAP_sampler( Harv.data = as.matrix(Harv.data)
                            , Aerial.data = as.matrix( Aeri.data)
							, nage = nage
							, measure.Fec = as.matrix(mean.f)
							, measure.Surv = as.matrix(NA + mean.s)
							, measure.SRB =  mean.SRB
							, prior.mean = prior.mean
							, prior.var = prior.var
							, prior.ageclass = prior.ageclass
							, n.iter = 100000, burn.in = 50000,thin.by = 50
							, prior.measurement.err = prior.measurement.error
                            , min.aK0 = list(matrix(-.001,nage[1],1),matrix(-.001,sum(nage),1),100)
                            , max.aK0 = list(matrix(.001,nage[1],1),matrix(.001,sum(nage),1),1500)
							, Assumptions = Assumptions
							, Observations = Observations
                            , prop.vars = prop.vars, estFec = T,estaK0 = T,null = F,global = T)

save.image("4harv_8fec_6surv_equal.RData")
