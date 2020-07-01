
Check_assumptions = function(Assumptions, nage, proj.period){
	if(is.null(Assumptions$Fec)) Assumptions$Fec = list(time = eyes(proj.period),age = as.matrix(eyes(nage[1])))
	if(is.null(Assumptions$Surv)) Assumptions$Surv = list(time = eyes(proj.period),age = as.matrix(eyes(sum(nage))))
	if(is.null(Assumptions$SRB)) Assumptions$SRB = list(time = eyes(proj.period),age = matrix(1,nage[1],1))
	if(is.null(Assumptions$AerialDet)) Assumptions$AerialDet  = list(time = eyes(proj.period+1),age = matrix(1,sum(nage),1))
	if(is.null(Assumptions$Harv)) Assumptions$Harv = list(time = eyes(proj.period+1),age = eyes(sum(nage)))
	if(is.null(Assumptions$aK0)) Assumptions$aK0 = list(eyes(nage[1]),eyes(sum(nage)),eyes(1))
	if(is.null(Assumptions$err)) Assumptions$err = list(time = matrix(1,1,proj.period)) # assume no age structure
	return(Assumptions)

}

Check_prior_measurement_err = function(obj){
	if(is.null(obj$Alpah$Fec)) obj$Alpah$Fec = 3
	if(is.null(obj$Beta$Fec)) obj$Beta$Fec = 1
	if(is.null(obj$Alpah$Surv)) obj$Alpah$Surv = 3
	if(is.null(obj$Beta$Surv)) obj$Beta$Surv = 1
	if(is.null(obj$Alpah$SRB)) obj$Alpah$SRB = 3
	if(is.null(obj$Beta$SRB)) obj$Beta$SRB = .5
	return(obj)

}

Check_start_measurement_err = function(obj,Ass_var){
	n = nrow(Ass_var)
	if(is.null(obj$Fec)) obj$Fec = matrix(.05,1,n)
	if(is.null(obj$Surv)) obj$Surv = matrix(.05,1,n)
	if(is.null(obj$SRB)) obj$SRB = matrix(.05,1,n)

	return(obj)
}

Check_observations = function(Observations, nage){
	if(is.null(Observations$Fec)) Observations$Fec = eyes(nage[1])
	if(is.null(Observations$Surv)) Observations$Surv = eyes(sum(nage))
	if(is.null(Observations$AerialCount)) Observations$AerialCount = matrix(1,sum(nage),1)

	return(Observations)

}

Check_Designs = function(Designs,nage,proj.period){

     if(is.null(Designs$Fec))  Designs$Fec = eyes(proj.period)
     if(is.null(Designs$Surv))  Designs$Surv = eyes(proj.period)
     if(is.null(Designs$SRB))  Designs$SRB = eyes(proj.period)
     if(is.null(Designs$Harv))  Designs$Harv = eyes(proj.period + 1)
     if(is.null(Designs$AerialDet))  Designs$AerialDet = eyes(proj.period + 1)
     return(Designs)

}

Check_prop_var = function(prop.var,nage,proj.period){
	if(is.null(prop.var$Fec)) prop.var$Fec = matrix(1,nrow = nage[1],ncol = proj.period)
	if(is.null(prop.var$Surv)) prop.var$Surv = matrix(1,nrow = sum(nage),ncol = proj.period)
	if(is.null(prop.var$SRB)) prop.var$SRB = matrix(.1,nage[1],proj.period)
	if(is.null(prop.var$AerialDet)) prop.var$AerialDet = matrix(1,1,proj.period+1)
	if(is.null(prop.var$Harv)) prop.var$Harv = matrix(1,sum(nage),proj.period+1)
	if(is.null(prop.var$aK0)) prop.var$aK0 = list(5e-8,5e-8,50)
	if(is.null(prop.var$baseline.pop.count)) prop.var$baseline.pop.count = matrix(1,nrow = sum(nage),ncol = 1)
	return(prop.var)
}




random_beta = function(Design, Assumption, nage,period_used){
	npara = ncol(Design)
	assu_time = Assumption$time
	assu_age = Assumption$age

	if(nrow(assu_age)!=nage){
		stop("    Age assumption matrix must have row number same to number of age classes considered\n")
	}

	if(ncol(assu_time) != period_used){
		stop("    Time assumption matrix must have column number same to number of periods considered\n")
	}



	if(nrow(Design)!= nrow(assu_time)) {
		stop("    Incompatable age assumption matrix with design matrix, if you did not specify the assumption matrix, check the Design matrix\n")
	}



	return(matrix( runif(ncol( assu_age ) * npara,-1,1),ncol = ncol( assu_age )))

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
		cat("    Incompatable time assumption matrix, if you did not specify the assumption matrix, check the prior mean matrix\n")
		errs = errs + 1
	}
	if(errs==0) cat("    all clear\n")

	return(errs)
}

Check_data = function(data_in,nage,nperiod){
	if(nrow(data_in)!=nage){
		stop("    Row number of data does not match age classes\n")
	}
	if(ncol(data_in)!=nperiod){
		stop("    Column number of data does not match period of measurements\n")
	}
	cat("    all clear\n")
	return(0)
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
log.lhood_popu =function(n.census, n.hat){
    ##-- value of log likelihoods --##
    density = dpois(x=n.census,lambda = n.hat,log = TRUE) # Use Poisson instead in Poisson_likelihood
    ##-- joint log likelihood --##

        return(sum(density,na.rm=T))
} # checked 10/24/2018

log.lhood_vital = function(f, s, SRB,estFec,measure.f, measure.s, measure.SRB
						   ,sigmasq.f, sigmasq.s, sigmasq.SRB
						   ){
	nperiod = ncol(f)
	if(estFec)

	log.f.lhood = sapply(1:nperiod,function(i,v,mv,sigsq){
		sum(dnorm(mv[,i],v[,i],sqrt(sigsq[i]),log=T),na.rm=T)
	},v=f,mv=measure.f,sigsq = sigmasq.f)


	else log.f.lhood = 0

	log.s.lhood = sapply(1:nperiod,function(i,v,mv,sigsq){
		sum(dnorm(mv[,i],v[,i],sqrt(sigsq[i]),log=T),na.rm=T)
	},v=s,mv=measure.s,sigsq = sigmasq.s)

	log.SRB.lhood = sapply(1:nperiod,function(i,v,mv,sigsq){
		sum(dnorm(mv[i],v[i],sqrt(sigsq[i]),log=T),na.rm=T)
	},v=SRB,mv=measure.SRB,sigsq = sigmasq.SRB)


	return(sum(log.f.lhood,log.SRB.lhood,log.s.lhood,na.rm = T)) # need to deal with no measurement case
}

## .............. Posterior ............ ##
## ..................................... ##    keep it, change in sampler function, but no migration here, should add estaK0

log.post = function(## estimated vitals
					f, s, SRB
                    , aK0, A, H # A for Aerial count detection probability
                    , estaK0,estFec
                    ## fixed prior means on detections
					, prior.mean.f, prior.mean.s
					, prior.var.f, prior.var.s
					, prior.mean.SRB, prior.var.SRB
                    , prior.mean.A, prior.mean.H
					, prior.var.A, prior.var.H
					, min.aK0, max.aK0
                    ## fixed prior parameters on variance distns
                    , alpha.f, beta.f, alpha.s, beta.s, alpha.SRB, beta.SRB # prior for measurement error
					, sigmasq.f, sigmasq.s, sigmasq.SRB # measurement errors
                    ## value of the log likelihood
                    , log.like
                    ## non zero rows of fertility matrix
                    ){

        ##-- Values of prior densities for vitals --##

        ##.. Note that log densities are calculated for numerical stability.
        ##         f, baseline.n, prior.mean.f, prior.mean.b are logged coming
        ##         in, s, prior.mean.s is logit transformed coming in, g and
        ##         prior.mean.g are not transformed coming in.
        ##-- prior for f and baseline K0 if needed to be estimatedd --##


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

		if(estFec){
			log.f.prior = dnorm(f, mean = prior.mean.f, sd = sqrt(prior.var.f)
                            ,log = TRUE)
			log.sigmasq.f.prior =
                log(dinvGamma(sigmasq.f, alpha.f, beta.f))
		}
		else{
			log.f.prior = 0
			log.sigmasq.f.prior = 0
		}

        ##-- prior for s and Sex Ratio at Birth (SRB), Aerial count detection rate, and hunting rate H --##
        log.s.prior = dnorm(s, mean = prior.mean.s, sd = sqrt(prior.var.s)
                            ,log = TRUE)
        log.SRB.prior = dnorm(SRB, mean = prior.mean.SRB, sd = sqrt(prior.var.SRB)
                              ,log = TRUE)
        log.H.prior = dnorm(H, mean = prior.mean.H
                            ,sd = sqrt(prior.var.H)
                            ,log = TRUE)
        log.A.prior = dnorm(A, mean = prior.mean.A
                            ,sd = sqrt(prior.var.A)
                            ,log = TRUE)
        log.sigmasq.s.prior =
                log(dinvGamma(sigmasq.s, alpha.s, beta.s))
        log.sigmasq.SRB.prior =
                log(dinvGamma(sigmasq.SRB, alpha.SRB, beta.SRB))
        ##-- The log posterior is the SUM of these with the log.like --##

    return(sum(log.f.prior, log.s.prior, log.SRB.prior
    						 , log.aK0.prior, log.H.prior,log.A.prior
                             , log.sigmasq.f.prior
                             , log.sigmasq.s.prior
                             , log.sigmasq.SRB.prior
                             , log.like,na.rm = T)) # keep NA for missing data, remove NAs

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


# this helper make mcmc object a large list, each entry has a sample
getListmcmc_full = function(mcmc_obj,Assumptions = list(),nage,n_proj){
  nsample = nrow(mcmc_obj[[1]])
  Assumptions = Check_assumptions(Assumptions,nage,n_proj)
  lapply(1:nsample,function(i,mcmc_obj,nage_female,nage_total,Assumptions){
    temp = lapply(mcmc_obj,function(obj,i){obj[i,]},i=i)
    temp$survival.mcmc = invlogit( Assumptions$Surv$age %*% matrix(temp$survival.mcmc,nrow = ncol(Assumptions$Surv$age),byrow = T) %*% Assumptions$Surv$time)
    temp$SRB.mcmc = invlogit(  Assumptions$SRB$age %*% matrix(temp$SRB.mcmc,nrow = ncol(Assumptions$SRB$age)) %*% Assumptions$SRB$time)
    temp$aerial.detection.mcmc = Assumptions$AerialDet$age %*% matrix(temp$aerial.detection.mcmc,nrow = ncol(Assumptions$AerialDet$age)) %*% Assumptions$AerialDet$time
    temp$H.mcmc = invlogit( Assumptions$Harv$age %*% matrix(temp$H.mcmc,nrow = ncol(Assumptions$Harv$age)) %*% Assumptions$Harv$time)
    temp$fecundity.mcmc = exp( Assumptions$Fec$age %*% matrix(temp$fecundity.mcmc,nrow = ncol(Assumptions$Fec$age),byrow = T) %*% Assumptions$Fec$time)
    temp$harvest.mcmc = matrix(temp$harvest.mcmc,ncol = n_proj+1)
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
    living = matrix(temp$living.mcmc,ncol = n_proj+1)
    living_t = matrix( colSums(living),nrow = 1)
    obs_lambda = get_obs_LambdasA(living_t)
    lambdas = rbind(hypo_lambdas,obs_lambda)
	  row.names(lambdas) = c("maximum","uniform_age","stable_age","no_harvest","minimum","observed")
	  return(lambdas)
  },mcmc_list)
  class(res) = "ReCAP_lambda"
  return(res)
}

analysisRecruitment = function(mcmc_obj,Assumptions = list(),nage,n_proj){
  mcmc_list = getListmcmc_full(mcmc_obj,Assumptions,nage,n_proj)
  res = lapply(1:length(mcmc_list),function(i,mcmc_list){
    temp = mcmc_list[[i]]
    hypo_lambdas = get_hypo_Lambdas(temp$harvest.mcmc,temp$living.mcmc,temp$H.mcmc,temp$survival.mcmc,temp$fecundity.mcmc,temp$SRB.mcmc)
    living = matrix(temp$living.mcmc,ncol = n_proj+1)
    living_t = matrix( colSums(living),nrow = 1)
	hypo_recu = apply(hypo_lambdas,1,"*",living_t[1:n_proj]) - living_t[1:n_proj]
    obs_recu = living_t[1:n_proj+1] - living_t[1:n_proj]
    lambdas = rbind(t(hypo_recu),obs_recu)
	  row.names(lambdas) = c("maximum","uniform_age","stable_age","no_harvest","minimum","observed")
	  return(lambdas)
  },mcmc_list)
  class(res) = "ReCAP_recruitment"
  return(res)
}

analysisquotaScheme = function(mcmc_obj,Assumptions = list(),nage,n_proj,harv_weight,skip = rep(0,n_proj+1)){
  mcmc_list = getListmcmc_full(mcmc_obj,Assumptions,nage,n_proj)
  harv_weight = apply(harv_weight,2,function(k){k/sum(k)})
  res = lapply(1:length(mcmc_list),function(i,mcmc_list,harv_weight,skip){
	  temp = mcmc_list[[i]]
	  get_hypo_harvest_quotaCpp(matrix(temp$living.mcmc[,1]+temp$harvest.mcmc[,1]),
						  temp$harvest.mcmc*(1-skip)+1e-5,
						  temp$survival.mcmc,
						  temp$fecundity.mcmc,
						  temp$SRB.mcmc,
						  harv_weight,
						  TRUE,
						  list(matrix(0,nage[1],1),matrix(0,sum(nage),1),matrix(10,1,1)), # no DD for now
						  TRUE)
  },mcmc_list,harv_weight,skip)

  return(res)
}

analysisportionScheme = function(mcmc_obj,Assumptions = list(),nage,n_proj,harv_weight,skip = rep(0,n_proj+1)){
    mcmc_list = getListmcmc_full(mcmc_obj,Assumptions,nage,n_proj)
    harv_weight = apply(harv_weight,2,function(k){k/sum(k)})
    res = lapply(1:length(mcmc_list),function(i,mcmc_list,harv_weight,skip){
        temp = mcmc_list[[i]]
        nage = c(nrow(temp$fecundity.mcmc),nrow(temp$survival.mcmc)-nrow(temp$fecundity.mcmc))
        period = ncol(temp$fecundity.mcmc)
        harvest_rate = matrix( colSums(temp$harvest.mcmc)/(colSums(temp$living.mcmc)+colSums(temp$harvest.mcmc)))
        harvest_rate = (harvest_rate * matrix(1-skip))+1e-5
        get_hypo_harvest_portionCpp(matrix(temp$living.mcmc[,1]+temp$harvest.mcmc[,1]),
                                  harvest_rate,
                                  temp$survival.mcmc,
                                  temp$fecundity.mcmc,
                                  temp$SRB.mcmc,
                                  harv_weight,
                                  TRUE,
                                  list(matrix(0,nage[1],1),matrix(0,sum(nage),1),matrix(10,1,1)), # no DD for now
                                  TRUE,
                                  period,nage)
    },mcmc_list,harv_weight,skip)
    class(res) = "ReCAP_Scheme"
    return(res)
}


analysisScheme = function(mcmc_obj,Assumptions = list(),nage,n_proj,harv_weight,quota=F,skip = rep(0,n_proj+1)){
    if(quota){
        res = analysisquotaScheme(mcmc_obj,Assumptions,nage,n_proj,harv_weight,skip)
    }
    else {
        res = analysisportionScheme(mcmc_obj,Assumptions,nage,n_proj,harv_weight,skip)


    }
    class(res) = "ReCAP_Scheme"
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

	if(class(ReCAP_lambda_obj) == "ReCAP_lambda") label_y = "Lambda   X(t+1)/X(t)"
	if(class(ReCAP_lambda_obj) == "ReCAP_recruitment") label_y = "Recruitment  X(t+1)-X(t)"
	res_plot = ggplot2::ggplot(data = plot_data,aes(x=year,y=lambda,color=point))+
		geom_line()+
		geom_point() +
		geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
		labs(y = label_y)
	return(list(plot = res_plot,data = plot_data))
}

Scheme_plot = function(ReCAP_Scheme_obj,start_year=1,alpha=.05){
	Sum_living = lapply(ReCAP_Scheme_obj,colSums)
    Sum_living = lapply(Sum_living, function(w){w[is.nan(w)]=0;return(w)})

	nyear = length(Sum_living[[1]])
	Sum_living_m = Reduce(rbind,Sum_living)
	mean_living = colMeans(Sum_living_m,na.rm = T)
	lower_025 = apply(Sum_living_m,2,quantile,probs = alpha/2,na.rm = T)

	higher_975 = apply(Sum_living_m,2,quantile,probs = 1-alpha/2,na.rm = T)


	year = 1:nyear + start_year
	plot_data = data.frame(point = "Post-cull population"
                      ,liv=mean_living
                      ,low = lower_025
                      ,high = higher_975
                      ,year = year)


	label_y = "Post-cull population (# individuals)"

	ggplot2::ggplot(data = plot_data,aes(x=year,y=liv))+
		geom_line()+
		geom_point() +
		geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
		labs(y = label_y)
}

plotthings = function(YD_obj,pathsave="./figs/temp/age",nage,period,years,ppt=F,ylabs = "individuals"){ # YD_obj should be a mcmc object, with vital rates in it
    if(ppt){require(export)}
    mean.harv = apply(YD_obj,2,mean)
    mean.harv.matrix = matrix(mean.harv,nrow = nage,ncol = period)

    BI.low.harv = apply(YD_obj,2,quantile,probs = .025)
    BI.low.harv.matrix = matrix(BI.low.harv,nrow = nage,ncol = period)
    BI_harv_low = data.frame(age = 1:nage,BI.low.harv.matrix)


    BI.high.harv = apply(YD_obj,2,quantile,probs = .975)
    BI.high.harv.matrix = matrix(BI.high.harv,nrow = nage,ncol = period)
    BI_harv_high = data.frame(age = 1:nage,BI.high.harv.matrix)

    har_data = data.frame(matrix(nrow = 1,ncol = 5))
    colnames(har_data) = c("age","mean","low","high","time")
    har_data = har_data[-1,]

    for(i in 1:nage){
        temp = data.frame(age = i,mean = mean.harv.matrix[i,],low = BI.low.harv.matrix[i,],high = BI.high.harv.matrix[i,],time = years)
        har_data = rbind(har_data,temp)
    }
    require(ggplot2)
    res = list(nage)
    for(i in 1:nage){
        temp = data.frame(point = "model predict (95% CI)",mean = mean.harv.matrix[i,],low = BI.low.harv.matrix[i,],high = BI.high.harv.matrix[i,],time = years)
        write.csv(temp,paste0(pathsave,i,".csv"))
        filename = paste0(pathsave,i,".jpg")
        plot_temp = ggplot(data.frame(temp),aes(x=time, y=mean, colour = point)) +
            geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
            geom_point() +
            geom_line() +
            ylab(ylabs) +
            theme(text = element_text(size=16))

        res[[i]] = plot_temp
        plot_temp
        if(ppt)    graph2ppt(file=sub("jpg","pptx",filename))
        else ggsave(filename, plot = last_plot())
    }
    return(res)
}


