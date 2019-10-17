
Check_assumptions = function(Assumptions, nage, proj.period){
	if(is.null(Assumptions$Fec)) Assumptions$Fec = list(time = eyes(proj.period),age = as.matrix(eyes(nage[1])))
	if(is.null(Assumptions$Surv)) Assumptions$Surv = list(time = eyes(proj.period),age = as.matrix(eyes(sum(nage))))
	if(is.null(Assumptions$SRB)) Assumptions$SRB = list(time = eyes(proj.period),age = eyes(1))
	if(is.null(Assumptions$AerialDet)) Assumptions$AerialDet  = list(time = eyes(proj.period+1),age = eyes(1))
	if(is.null(Assumptions$Harv)) Assumptions$Harv = list(time = eyes(proj.period+1),age = eyes(sum(nage)))
	if(is.null(Assumptions$aK0)) Assumptions$aK0 = list(eyes(nage[1]),eyes(sum(nage)),eyes(1))

	return(Assumptions)

}

Check_prop_var = function(prop.var,nage,proj.period){
	if(is.null(prop.var$fert.rate)) prop.var$fert.rate = matrix(1,nrow = nage[1],ncol = proj.period)
	if(is.null(prop.var$surv.prop)) prop.var$surv.prop = matrix(1,nrow = sum(nage),ncol = proj.period)
	if(is.null(prop.var$SRB)) prop.var$SRB = matrix(.1,nage[1],proj.period)
	if(is.null(prop.var$A)) prop.var$A = matrix(1,1,proj.period+1)
	if(is.null(prop.var$H)) prop.var$H = matrix(1,sum(nage),proj.period+1)
	if(is.null(prop.var$aK0)) prop.var$aK0 = list(5e-8,5e-8,50)
	if(is.null(prop.var$baseline.pop.count)) prop.var$baseline.pop.count = matrix(1,nrow = sum(nage),ncol = 1)
	return(prop.var)
}

Check_dimensions = function(mean_vital,Assumption,nage,period_used){
	errs = 0
	assu_time = Assumption$time
	assu_age = Assumption$age
	if(nrow(assu_age)!=nage){
		cat("    Age assumption matrix must have row number same to number of age classes considered\n")
		errs = errs + 1
	}

	if(ncol(assu_time) != period_used){
		cat("    Time assumption matrix must have column number same to number of periods considered\n")
		errs = errs + 1
	}

	if(ncol(assu_age) != nrow(mean_vital)){
		cat("    Incompatable age assumption matrix, if you did not specify the assumption matrix, check the prior mean matrix\n")
		errs = errs + 1
	}

	if(nrow(assu_time) != ncol(mean_vital)){
		cat("    Incompatable age assumption matrix, if you did not specify the assumption matrix, check the prior mean matrix\n")
		errs = errs + 1
	}
	if(errs==0) cat("    all clear\n")

	return(errs)
}

Check_data = function(data_in,nage,nperiod){
	errs = 0
	if(nrow(data_in)!=nage){
		cat("    Row number of data does not match age classes\n")
		errs = errs + 1
	}
	if(ncol(data_in)!=nperiod){
		cat("    Column number of data does not match period of measurements\n")
		errs = errs + 1
	}
	if(errs==0) cat("    all clear\n")
	return(errs)
}

ProjectHarvest = function(Surv, Harvpar, Fec, SRB, bl, period, nage, aK0 = list(matrix(0,nage[1],1),matrix(0,sum(nage),1),matrix(0,1,1)), global = T, null = T){
	ProjectAllCpp(Surv, Harvpar, Fec, SRB, aK0, global, null, bl, period, nage)
}

########
# Coming functions are
# modified from popReconstruct package

## ........... Misc Functions .......... ##
## ..................................... ##

logitf = function(p){
    log(p / (1 - p))
 } # checked 10/24/2018


## logit function
invlogit = function(x){
        if(any(is.infinite(exp(x)))) {
                y = x
                y[is.infinite(exp(x))] = 1
                y[!is.infinite(exp(x))] =
                        invlogit(y[!is.infinite(exp(x))])
                return(y)
        }
        else return(exp(x) / (1 + exp(x)))
} # checked 10/24/2018

##--- Generates random draws from inverse gamma ---##
rinvGamma = function(n, shape, scale){
        return(1/(rgamma(n, shape = shape, rate = scale)))
} # checked 10/24/2018 # debug mode with print

##--- Returns value of inverse gamma pdf ---##
dinvGamma = function(x, shape, scale, log = FALSE){
        if(log) d =
                shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
        else d = scale^shape / gamma(shape) * (1/x)^(shape + 1) * exp(-scale/x)
        return(d)
} # checked 10/24/2018


## ............. Likelihood ............ ##
## ..................................... ##

#log likelihood function of gaussian distributed
log.lhood =function(n.census, n.hat){
    ##-- value of log likelihoods --##
    density = dpois(x=n.census,lambda = n.hat,log = TRUE) # Use Poisson instead in Poisson_likelihood
    ##-- joint log likelihood --##

        return(sum(density))
} # checked 10/24/2018

## .............. Posterior ............ ##
## ..................................... ##    keep it, change in sampler function, but no migration here, should add estaK0

log.post = function(## estimated vitals
                    f, s, SRB,baseline.n, aK0, A, H # A for Aerial count detection probability
                    , estFec, estaK0
                    ## fixed prior means on vitals
                    , prior.mean.f, prior.mean.s, prior.mean.SRB
                    , prior.mean.b #, prior.mean.aK0
                    , prior.mean.A, prior.mean.H
                    ## fixed prior parameters on variance distns
                    , alpha.f, beta.f, alpha.s, beta.s, alpha.SRB, beta.SRB
                    , min.aK0, max.aK0
                    , alpha.A, beta.A, alpha.H, beta.H
                    ## updated variances on prior distns
                    , sigmasq.f, sigmasq.s, sigmasq.SRB,sigmasq.n#, sigmasq.aK0
                    , sigmasq.A ,sigmasq.H
                    ## value of the log likelihood
                    , log.like
                    ## non zero rows of fertility matrix
                    , non.zero.fert){

        ##-- Values of prior densities for vitals --##

        ##.. Note that log densities are calculated for numerical stability.
        ##         f, baseline.n, prior.mean.f, prior.mean.b are logged coming
        ##         in, s, prior.mean.s is logit transformed coming in, g and
        ##         prior.mean.g are not transformed coming in.
        ##-- prior for f and baseline K0 if needed to be estimatedd --##

        if(estFec){
            log.f.prior = dnorm(as.vector(f[non.zero.fert,])
                                , mean = as.vector(prior.mean.f[non.zero.fert,])
                                , sd = sqrt(sigmasq.f)
                                , log = TRUE)
            log.sigmasq.f.prior =
                log(dinvGamma(sigmasq.f, alpha.f, beta.f))
        }
        else {
            log.f.prior = 0
            log.sigmasq.f.prior = 0
        }

        if(estaK0){
            log.aK0.prior =#sum( dunif(aK0[[1]],min.aK0[[1]],max.aK0[[1]],T) , dunif(aK0[[2]],min.aK0[[2]],max.aK0[[2]],T))
                Reduce(sum,
                       lapply(1:length(aK0),
                              function(i,aK0,al.aK0,be.aK0){
                                  dunif(aK0[[i]],al.aK0[[i]],be.aK0[[i]],T)}
                              ,aK0,min.aK0,max.aK0))

        }
        else {
            log.aK0.prior = 0
            log.sigmasq.aK0.prior = 0
        }

        ##-- prior for s and Sex Ratio at Birth (SRB), Aerial count detection rate, and hunting rate H --##
        log.s.prior = dnorm(s, mean = prior.mean.s, sd = sqrt(sigmasq.s)
                            ,log = TRUE)
        log.SRB.prior = dnorm(SRB, mean = prior.mean.SRB, sd = sqrt(sigmasq.SRB)
                              ,log = TRUE)
        log.H.prior = dnorm(H, mean = prior.mean.H
                            ,sd = sqrt(sigmasq.H)
                            ,log = TRUE)
        log.A.prior = dnorm(A, mean = prior.mean.A
                            ,sd = sqrt(sigmasq.A)
                            ,log = TRUE)
        log.sigmasq.s.prior =
                log(dinvGamma(sigmasq.s, alpha.s, beta.s))
        log.sigmasq.SRB.prior =
                log(dinvGamma(sigmasq.SRB, alpha.SRB, beta.SRB))
        log.sigmasq.A.prior =
                log(dinvGamma(sigmasq.A, alpha.A, beta.A))
        log.sigmasq.H.prior =
                log(dinvGamma(sigmasq.H, alpha.H, beta.H))

        ##-- The log posterior is the SUM of these with the log.like --##

    return(sum(log.f.prior, log.s.prior, log.SRB.prior#, log.b.prior
    , log.aK0.prior, log.H.prior,log.A.prior,
                             log.sigmasq.f.prior
                             ,log.sigmasq.s.prior
                             ,log.sigmasq.SRB.prior
                             #,log.sigmasq.aK0.prior
                             ,log.sigmasq.H.prior
                             ,log.sigmasq.A.prior
                             ,log.like,na.rm = T)) # keep NA for missing data, remove NAs

}

## ......... Acceptance Ratio .......... ##
## ..................................... ##
acc.ra = function(log.prop, log.current){
        min(1, exp(log.prop - log.current))
}

acc.ra.var = function(log.prop.post, log.curr.post, log.prop.var, log.curr.var){
        min(1, exp(log.curr.var + log.prop.post - log.prop.var - log.curr.post))
}

## sensitivity analysis
HarvestSen = function(Fec,Surv,SRB,Harvpar,nage,Harv_assump){
        L = getLeslie(Surv,Fec,SRB) # intrinsic growth
        H = matrix(0,sum(nage),sum(nage))
        diag(H) = Harv_assump %*% Harvpar # propotional harvest

        EL = H %*% L # effactive growth

        ev = eigen(EL)
        lmax = which(Re(ev$values) == max(Re(ev$values)))
        lambda = Re(ev$values[lmax])
        W = ev$vectors
        w = abs(Re(W[, lmax]))
        V = Conj(solve(W))
        v = abs(Re(V[lmax, ]))
        E_Sen = v %o% w # sensitivity analysis of the effective growth

        # apply chain rule:
        H_Sen = rowSums( t(Harv_assump) %*% (E_Sen*L)) # sensitivity of harvest
        return(H_Sen)
}



getListmcmc_full = function(mcmc_obj,Assumptions = list(),nage,n_proj){
  nsample = nrow(mcmc_obj[[1]])
  Assumptions = Check_assumptions(Assumptions,nage,n_proj)
  lapply(1:nsample,function(i,mcmc_obj,nage_female,nage_total,Assumptions){
    temp = lapply(mcmc_obj,function(obj,i){obj[i,]},i=i)
    temp$survival.mcmc = Assumptions$Surv$age %*% matrix(temp$survival.mcmc,nrow = ncol(Assumptions$Surv$age)) %*% Assumptions$Surv$time
    temp$SRB.mcmc = Assumptions$SRB$age %*% matrix(temp$SRB.mcmc,nrow = ncol(Assumptions$SRB$age)) %*% Assumptions$SRB$time
    temp$aerial.detection.mcmc = Assumptions$AerialDet$age %*% matrix(temp$aerial.detection.mcmc,nrow = ncol(Assumptions$AerialDet$age)) %*% Assumptions$AerialDet$time
    temp$H.mcmc = Assumptions$Harv$age %*% matrix(temp$H.mcmc,nrow = ncol(Assumptions$Harv$age)) %*% Assumptions$Harv$time
    temp$fecundity.mcmc = Assumptions$Fec$age %*% matrix(temp$fecundity.mcmc,nrow = ncol(Assumptions$Fec$age)) %*% Assumptions$Fec$time
    temp$harvest.mcmc =  matrix(temp$harvest.mcmc,ncol = n_proj+1)
	temp$living.mcmc =  matrix(temp$living.mcmc,ncol = n_proj+1)
    temp$aerial.count.mcmc = matrix(temp$aerial.count.mcmc,ncol = n_proj+1)
    return(temp)
  } , mcmc_obj,nage[1],sum(nage),Assumptions)
}

## analysis Lambda, only work for no aK0 settings.
analysisLambda = function(mcmc_obj,Assumptions = list(),nage,n_proj){
  mcmc_list = getListmcmc_full(mcmc_obj,Assumptions,nage,n_proj)
  res = lapply(1:length(mcmc_list),function(i,mcmc_list){
    temp = mcmc_list[[i]]
    hypo_lambdas = get_hypo_Lambdas(temp$harvest.mcmc,temp$living.mcmc,temp$H.mcmc,temp$survival.mcmc,temp$fecundity.mcmc,temp$SRB.mcmc)
    obs_lambda = get_obs_LambdasA(temp$living.mcmc)
    lambdas = rbind(hypo_lambdas,obs_lambda)
	  row.names(lambdas) = c("maximum","uniform_age","stable_age","no_harvest","minimum","observed")
	  return(lambdas)
  },mcmc_list)
  class(res) = "ReCAP_lambda"
  return(res)
}


Lambda_plot = function(ReCAP_lambda_obj,start_year=1,alpha=.05){
	mean_lambda = Reduce("+",ReCAP_lambda_obj)/(length(ReCAP_lambda_obj))
	nyear = ncol(mean_lambda)
	list_each_lambda = lapply(1:length(mean_lambda),function(i,ReCAP_lambda_obj1){
		temp = lapply(ReCAP_lambda_obj1,function(ana,i){ana[i]},i)
		Reduce(rbind,temp)
	},ReCAP_lambda_obj)
	lower_025 = sapply(list_each_lambda,quantile,probs = alpha/2)
	lower_025 = matrix(lower_025,ncol=nyear)
	higher_975 = sapply(list_each_lambda,quantile,probs = 1-alpha/2)
	higher_975 = matrix(higher_975,ncol=nyear)

	year = 1:nyear + start_year
	observed = data.frame(point = "observed (w/ harvest)"
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
	ggplot2::ggplot(data = plot_data,aes(x=year,y=lambda,color=point))+
		geom_line()+
		geom_point() +
		geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
		labs(y = "Lambda   X(t+1)/X(t)")

}

