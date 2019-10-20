### --------------------------- SAMPLER --------------------------- ###
### --------------------------------------------------------------- ###

ReCAP_sampler =
        function( Harv.data
                         , Aerial.data
                         , nage
                         #.. fixed prior means
                         , measure.f, measure.s, measure.SRB
                         , prior.mean.f, prior.mean.s, prior.mean.SRB
                         , prior.var.f, prior.var.s, prior.var.SRB              
                         , prior.mean.H, prior.mean.A
                         , prior.var.H = 2, prior.var.A = 2
                         , measure.b=Harv.data[,1]
                         , Aerialcount_time = "post"
                         , n.iter=50000, burn.in = 5000, thin.by = 10

                         #.. fixed variance hyper-parameters
                         ,al.f = 1, be.f = 0.05, al.s = 1, be.s = 0.1,al.SRB = 1,be.SRB = 0.05
                         , min.aK0 = list(matrix(-.001,nage[1],1),matrix(-.001,sum(nage),1),0)
                         , max.aK0 = list(matrix(.001,nage[1],1),matrix(.001,sum(nage),1),500)
                         
                         #.. census data
                         #     *not transformed coming in*
                         , Assumptions = list() # Assumption matrices
                         , Designs = list() # Design matrices 
                         , Observations = list() # observation matrices, row as raw most detailed age, columns as age classes (e.g. fawn, yearling, adults)
                          
                         #.. inital values for vitals and variances
                         #     *vitals not transformed coming in* all not transfer, will transfer later before sample and transfer back when saving
                         , start.sigmasq.f = .05, start.sigmasq.s = .05, start.sigmasq.SRB = .05

                         #.. **variances** for proposal distributions used in M-H
                         #     steps which update vital rates.
                         ,prop.vars = list() # col names should be as follow:
                                     # "fert.rate", "surv.prop", "SRB","H", "A","aK0"
                                     # ,"baseline.count"

                         #.. number of periods to project forward over (e.g.,
                         #         number of projection steps, usually a year
                         ,proj.periods = (ncol(Harv.data)-1)

                         #.. age group, if multiple sex, the one reproduce should be at first.


                         ,estFec=T, estaK0 = F
                         ,aK0 = list(matrix(0,nage[1],1),matrix(0,sum(nage),1),matrix(0,1,1)), global = T, null = T
                         # control parameters for the model, global is whether density dependency is global rather than age specific, null is whether exist density dependency (True of not ).
                         ,point.est = mean
                         #.. print algorithm progress
                         ,verb = FALSE
                         #.. tolerance defining allowable survival probabilities, should be >0, or things is gonna die out.
                         ,s.tol = 10^(-10)
                         )
{
        ## .............. Sampler .............. ##
        ## ..................................... ##
        ## -------- Checking input dimensions ------- ##
        measure.f = as.matrix( measure.f)
        measure.s = as.matrix( measure.s)
        measure.SRB = as.matrix( measure.SRB)
        measure.b = as.matrix( measure.b)
        Harv.data = as.matrix(Harv.data)
        Aerial.data = as.matrix(Aerial.data)

        missing_harv = which(colSums(Harv.data)==0 |
                             colSums(is.na(Harv.data))==nrow(Harv.data)) # for missing harvest/skipped year

        missing_aerial = which((Aerial.data)==0 |
                               is.na(Aerial.data)) # for missing aerial count years
		Aerial.data[missing_aerial] = 0 # set it to be 0, since aerial count is non-invasive, we just set detection to be 0 at that year.

        if(Aerialcount_time=="pre") getAerialCount=getAerialCountPre
        else getAerialCount=getAerialCountPost
        cat("Checking input dimensions...\n")
        
        
        Assumptions = Check_assumptions(Assumptions, nage, proj.periods)
        prop.vars = Check_prop_var(prop.vars,nage,proj.periods)
                
        Observations = Check_observations(Observations, nage)
        errs_dim = 0
        cat("  Check Harvest data:\n")
        Check_data(Harv.data,sum(nage),proj.periods + 1 )
        cat("\n  Check Aerial count data:\n")
        Check_data(Aerial.data,1,proj.periods + 1)
        
       
        
                
        start.f = random_beta(Designs$Fec, Assumption$Fec, nage,proj.periods)
        srart.s = random_beta(Designs$Surv, Assumption$Surv,nage,proj.periods)
        srart.SRB = random_beta(Designs$SRB, Assumption$SRB,nage,proj.periods)
        srart.H = random_beta(Designs$Harv, Assumption$Harv,nage,proj.periods+1)
        srart.A = random_beta(Designs$AerialDet, Assumption$AerialDet,nage,proj.periods+1)        
                

        start.b = measure.b
        start.b[is.na(start.b)] = 5
        start.aK0 = min.aK0
                


        ## -------- Begin timing ------- ##
        cat("\n")
        ptm = proc.time()


        # ## ------- Determine fert.rows --------- ##

        zero.elements = measure.f == 0
        fert.rows = as.logical(apply(zero.elements, 1, function(z) !all(z)))


        ## ---------- Storage ---------- ##

        #.. MCMC objects for posterior samples
        # Samples are stored as 2D arrays for compatibility with coda's
        # mcmc format with iterations as rows, year*age.group as columns.
        # Age.group cycles fastest across columns, e.g.,
        # _______________________________________________________________
        #     1992    | 1992    | 1992    | ... | 1993    | 1993    | ...
        #     15.19   | 20.24   | 25.29   | ... | 15.19   | 20.24   | ...
        # 1    --     |    --   |    --   | ... |    --   |    --   | ...
        # 2    --     |    --   |    --   | ... |    --   |    --   | ...
        #     etc.
        # _______________________________________________________________

        ## How many (samples) stored?
        cat("Preparing...")

        n.stored = ceiling(n.iter / thin.by)
            # Fertility

        if(estFec){
            fert.rate.mcmc =
                    mcmc(matrix(nrow = n.stored
                             ,ncol = length(start.f))
                             ,start = burn.in + 1
                             ,thin = thin.by
                             )
            colnames(fert.rate.mcmc) = NULL

        }
        else{fert.rate.mcmc = NULL}
            # Survival proportions
        surv.prop.mcmc =
                    mcmc(matrix(nrow = n.stored
                             ,ncol = length(start.s))
                             ,start = burn.in + 1
                             ,thin = thin.by
                             )
        colnames(surv.prop.mcmc) = NULL
             # Sex Ratio at Birth
        SRB.mcmc =
                    mcmc(matrix(nrow = n.stored
                             ,ncol = length(start.SRB))
                             ,start = burn.in + 1
                             ,thin = thin.by
                             )
        colnames(surv.prop.mcmc) = NULL

        log.like.mcmc =
                mcmc(matrix(nrow = n.stored
                           ,ncol = 1)
                                 ,start = burn.in + 1
                       ,thin = thin.by)
        colnames(log.like.mcmc) = NULL

            # lx
            # this is current culling
        living.mcmc = mcmc(matrix(nrow = n.stored
                                  ,ncol = nrow(start.b) * (proj.periods+1))
                           ,start = burn.in + 1
                           ,thin = thin.by
        )
        colnames(living.mcmc) = NULL
            lx.mcmc =
                    mcmc(matrix(nrow = n.stored
                             ,ncol = nrow(start.b) * (proj.periods+1))
                             ,start = burn.in + 1
                             ,thin = thin.by
                             )
            colnames(lx.mcmc) = NULL

            # this is for aerial count
            ae.mcmc =
                    mcmc(matrix(nrow = n.stored
                             ,ncol = proj.periods+1)
                             ,start = burn.in + 1
                             ,thin = thin.by)

            # carrying capacity assumed to be time homogeneous
            if(estaK0){
                aK0.Fec.mcmc =
                    mcmc(matrix(nrow = n.stored
                             ,ncol = length(start.aK0[[1]]))
                             ,start = burn.in + 1
                             ,thin = thin.by)
                colnames(aK0.Fec.mcmc) = NULL
                aK0.Surv.mcmc =
                    mcmc(matrix(nrow = n.stored
                             ,ncol = length(start.aK0[[2]]))
                             ,start = burn.in + 1
                             ,thin = thin.by)
                colnames(aK0.Surv.mcmc) = NULL
                aK0.midPopulation.mcmc =
                    mcmc(matrix(nrow = n.stored
                             ,ncol = length(start.aK0[[3]]))
                             ,start = burn.in + 1
                             ,thin = thin.by)
                colnames(aK0.midPopulation.mcmc) = NULL
            }
            else{
                aK0.Fec.mcmc = NULL
                aK0.Surv.mcmc = NULL
                aK0.midPopulation.mcmc = NULL
            }

            # Harvest proportion, can be either time homo or not
            H.mcmc =
                    mcmc(matrix(nrow = n.stored
                             ,ncol = length(start.H))
                             ,start = burn.in + 1
                             ,thin = thin.by)
            colnames(H.mcmc) = NULL
            # Aerial counts
            A.mcmc =
                    mcmc(matrix(nrow = n.stored
                             ,ncol = length(start.A))
                             ,start = burn.in + 1
                             ,thin = thin.by)
            colnames(H.mcmc) = NULL

            # baseline counts
            baseline.count.mcmc =
                    mcmc(matrix(nrow = n.stored, ncol = length(start.b))
                             ,start = burn.in + 1
                             ,thin = thin.by)
            colnames(baseline.count.mcmc) = NULL

            # variances
            variances.mcmc =
                    mcmc(matrix(nrow = n.stored, ncol = 3)
                             ,start = burn.in + 1
                             ,thin = thin.by)
            colnames(variances.mcmc) =
                    c("fert.rate.var", "surv.prop.var", "SRB.var") 

        #.. Record acceptance rate

        acc.count =
                list(fert.rate = matrix(0, nrow = nrow(as.matrix( measure.f[fert.rows,]))
                         ,ncol = ncol(as.matrix( measure.f[fert.rows,]))
                         ,dimnames = dimnames(( measure.f[fert.rows,]))
                         )
                         ,surv.prop = matrix(0, nrow = nrow(as.matrix( measure.s))
                            ,ncol = ncol(as.matrix( measure.s))
                            ,dimnames = dimnames(measure.s)
                            )
                         ,SRB = matrix(0, nrow = nrow(as.matrix( measure.SRB))
                            ,ncol = ncol(as.matrix( measure.SRB))
                            ,dimnames = dimnames(measure.SRB))
                         ,A = matrix(0, nrow = nrow(as.matrix(prior.mean.A)), ncol = ncol(as.matrix( prior.mean.A))
                            ,dimnames = dimnames(prior.mean.A)
                            )
                         ,H = matrix(0, nrow = nrow(as.matrix(prior.mean.H)), ncol = ncol(as.matrix( prior.mean.H))
                            ,dimnames = dimnames(prior.mean.H)
                            )
                         ,aK0 = matrix(0, nrow = nrow(as.matrix( start.aK0)), ncol = ncol( as.matrix(start.aK0))
                            ,dimnames = dimnames(start.aK0)
                            )
                         ,baseline.count = matrix(0, nrow = nrow(as.matrix( measure.b))
                            ,dimnames = dimnames( measure.b)
                            )
                         ,sigmasq.f = 0
                         ,sigmasq.s = 0
                         ,sigmasq.SRB = 0

                         )


        #.. Count how often acceptance ratio missing or na

        ar.na = acc.count


        #.. Count how often projection gives negative population

        pop.negative =
                list(fert.rate = matrix(0, nrow = nrow(as.matrix( measure.f[fert.rows,]))
                         ,ncol = ncol(as.matrix( measure.f[fert.rows,]))
                         ,dimnames = dimnames(( measure.f[fert.rows,]))
                         )
                         ,surv.prop = matrix(0, nrow = nrow(as.matrix( measure.s))
                            ,ncol = ncol(as.matrix( measure.s))
                            ,dimnames = dimnames(measure.s)
                            )
                         ,SRB = matrix(0, nrow = nrow(as.matrix( measure.SRB))
                            ,ncol = ncol(as.matrix( measure.SRB)))
                            ,dimnames = dimnames(measure.SRB)
                         ,A = matrix(0, nrow = nrow(as.matrix(prior.mean.A)), ncol = ncol(as.matrix( prior.mean.A))
                            ,dimnames = dimnames(prior.mean.A)
                            )
                         ,H = matrix(0, nrow = nrow(as.matrix(prior.mean.H)), ncol = ncol(as.matrix( prior.mean.H))
                            ,dimnames = dimnames(prior.mean.H)
                            )
                         ,aK0 = matrix(0, nrow = nrow(as.matrix( start.aK0)), ncol = ncol( as.matrix(start.aK0))
                            ,dimnames = dimnames(start.aK0)
                            )
                         ,baseline.count = matrix(0, nrow = nrow(as.matrix( measure.b))
                            ,dimnames = dimnames( measure.b)
                            )
                         )


        #.. Count how often surv probs are outside tolerance

        s.out.tol = matrix(0, nrow = nrow(measure.s), ncol = ncol(measure.s)
                                                ,dimnames = dimnames(measure.s))

        cat("done\n")
        ## -------- Initialize -------- ## Restart here in 10/19/2018
        cat("Initializing...")
        #.. Set current vitals and variances to inital values
        #     Take logs/logits here where required
        if(estFec){
            log.curr.f = log(start.f)    #<- log(0) stored as "-Inf". Gets
            log.prop.f = log(start.f)    #        converted to 0 under exponentiation

        }
        else{
            log.curr.f =    (!estFec)*log(start.f) #<- log(0) stored as "-Inf". Gets
            log.prop.f =    (!estFec)*log(start.f) #        converted to 0 under exponentiation
        }
        logit.curr.s = logitf(start.s)
        logit.curr.SRB = logitf(start.SRB)
        logit.curr.H = logitf(start.H)
        logit.curr.A = logitf(start.A)
        curr.aK0=(start.aK0)
        #curr.aK0=(aK0)}
        #curr.aK0 = estaK0 * log(start.aK0) + (!estaK0) * log(K0)
        log.curr.b = log(start.b)

        curr.sigmasq.f = start.sigmasq.f
        curr.sigmasq.s = start.sigmasq.s
        curr.sigmasq.SRB = start.sigmasq.SRB
        curr.sigmasq.A = start.sigmasq.A
        curr.sigmasq.H = start.sigmasq.H
        #curr.sigmasq.aK0 = start.sigmasq.aK0



        #.. Fixed means for vitals and baseline
        #     Set these to inputs, take logs where required.

        log.measure.f = log(measure.f)
        logit.measure.s = logitf(measure.s)
        logit.measure.SRB = logitf(measure.SRB)
                                     
        log.prior.mean.f = log(prior.mean.f)
        logit.prior.mean.s = logitf(prior.mean.s)
        logit.prior.mean.SRB = logitf(prior.mean.SRB)
        
        logit.prior.mean.A = logitf(prior.mean.A)
        logit.prior.mean.H = logitf(prior.mean.H)
        log.measure.b = log(measure.b)


        #.. Fixed Harvest data
        #     Take logs here

        log.Harv.mat = log(Harv.data)
        log.Aeri.mat = log(Aerial.data)

        Fec_assump = Assumptions$Fec
        Surv_assump = Assumptions$Surv
        SRB_assump = Assumptions$SRB
        A_assump = Assumptions$AerialDet
        Harv_assump = Assumptions$Harv
        aK0_assump = Assumptions$aK0

        #.. Set current projection: base on initial values # homo or not is important, determin it use homo = T

        # make full vital matrix, not good for ram saving
        logit.curr.s.full = Surv_assump$age %*% logit.curr.s %*%Surv_assump$time
        log.curr.f.full = Fec_assump$age %*% log.curr.f %*% Fec_assump$time
        logit.curr.H.full = Harv_assump$age %*% logit.curr.H %*%Harv_assump$time
        logit.curr.SRB.full = SRB_assump$age %*% logit.curr.SRB %*%SRB_assump$time
        logit.curr.A.full = A_assump$age %*% logit.curr.A %*%A_assump$time
        curr.aK0.full = lapply(1:length(aK0_assump),function(i,aK0,assump){
            assump[[i]] %*% aK0[[i]]
        },aK0 = curr.aK0,assump = aK0_assump)


        curr.proj =
                (ProjectHarvest(Surv = invlogit(logit.curr.s.full), Harvpar = invlogit(logit.curr.H.full),Fec=exp(log.curr.f.full), SRB = invlogit(logit.curr.SRB.full), aK0 = (curr.aK0.full), global = global, null = null, bl = exp(log.curr.b)    , period = proj.periods, nage = nage))

        curr.aeri = ( getAerialCount( curr.proj,A = invlogit(logit.curr.A.full),obsMat = Observations$AerialCount))

        log.curr.obs_f = log(getobsVitals(curr.proj$Fec_obs,curr.proj$living,Observations$Fec))
        logit.curr.obs_s = logitf(getobsVitals(curr.proj$Surv_obs,curr.proj$living,Observations$Surv))
#.. Current log posterior


        log.curr.posterior =
                log.post(f = log.curr.f
                                             ,s = logit.curr.s
                                             ,SRB = logit.curr.SRB
                                             ,A = logit.curr.A
                                             ,H = logit.curr.H
                                             ,aK0 = curr.aK0
                                             
                                             ,estFec=estFec, estaK0=estaK0
                         
                                             ,prior.mean.f = log.prior.mean.f
                                             ,prior.mean.s = log.prior.mean.s
                                             ,prior.mean.SRB = log.prior.mean.SRB
                                             ,prior.mean.A = logit.prior.mean.A
                                             ,prior.mean.H = logit.prior.mean.H

                                             ,prior.var.f = prior.var.f
                                             ,prior.var.s = prior.var.s
                                             ,prior.var.SRB = prior.var.SRB
                                             ,prior.var.A = prior.var.A
                                             ,prior.var.H = prior.var.H
                                             
                                             ,alpha.f = al.f, beta.f = be.f
                                             ,alpha.s = al.s, beta.s = be.s
                                             ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                                             ,alpha.A = al.A, beta.A = be.A
                                             ,alpha.H = al.H, beta.H = be.H
                                             ,min.aK0 = min.aK0, max.aK0 = max.aK0


                                             ,sigmasq.f = curr.sigmasq.f
                                             ,sigmasq.s = curr.sigmasq.s
                                             ,sigmasq.SRB = curr.sigmasq.SRB
                                             



                                             ,log.like =
                                                        log.lhood_popu(
                                                                n.census = Harv.data
                                                                ,n.hat = curr.proj$Harvest
                                                                ) +
                                                        log.lhood_popu(
                                                                n.census = Aerial.data
                                                                ,n.hat = curr.aeri
                                                                ) +
                                                        log.lhood_vital(f = log.curr.obs_f
                                                                        ,s = logit.curr.obs_s # observed survival and fecundity after adding DD.
                                                                        ,SRB = logit.curr.SRB
                                                                        ,estFec = estFec
                                                                        ,measure.f = log.measure.f
                                                                        ,measure.s = log.measure.s
                                                                        ,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = curr.sigmasq.f
                                                                        ,sigmasq.s = curr.sigmasq.s
                                                                        ,sigmasq.SRB = curr.sigmasq.SRB
                                                                        ,non.zero.fert = fert.rows)
                                                        #　both harvest and aerial count, as well as vital rates
                                             ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                                             )

        ## -------- Begin loop ------- ##
        #...............................#

        if(verb) {
                cat("\n\ntotal iterations = ", n.iter+burn.in
                        ,"\nburn in = ", burn.in
                        ,"\nthin = ", thin.by
                        ,",\nnumber stored = ", n.stored, sep = "")
                cat("\n\nfert.rows = ", which(fert.rows)
                        )
                cat("\n\n"
                        ,"iter ", " quantity\n", "---- ", " --------"
                        ,sep = "")
        }
        cat("done\n")
        cat("Start sampling...\n")
        for(i in 1:(n.iter + burn.in)) {
            svMisc::progress(((i-1)/(n.iter + burn.in))*100,progress.bar = T)
            # k is the index into the storage objects
            k = (i - burn.in - 1) / thin.by + 1


            ## -------- Vital Rate M-H Steps ------- ##
        if(estFec){
            ##...... Fertility .....##

            if(verb && identical(i%%1000, 0)) cat("\n", i, " Fertility")

            # - Proposal

            #.. cycle through components
            for(j in 1:length(log.curr.f[fert.rows,])) {

                #.. make a matrix conformable w fertility rate matrix
                log.prop.f.mat =
                        matrix(0, nrow = nrow(log.curr.f), ncol = ncol(log.curr.f))
                log.prop.f.mat[fert.rows,][j] =
                        rnorm(1, 0, sqrt(prop.vars$fert.rate[j])) #pop vars col names

                #.. make proposal
                log.prop.f = log.curr.f + log.prop.f.mat
                # - Run CCMP (project on the original scale)
                #     ** Don't allow negative population
                logit.curr.s.full = Surv_assump$age %*% logit.curr.s %*%Surv_assump$time
                log.prop.f.full = Fec_assump$age %*% log.prop.f %*% Fec_assump$time
                logit.curr.H.full = Harv_assump$age %*% logit.curr.H %*%Harv_assump$time
                logit.curr.SRB.full = SRB_assump$age %*% logit.curr.SRB %*%SRB_assump$time
                logit.curr.A.full = A_assump$age %*% logit.curr.A %*%A_assump$time
                curr.aK0.full = lapply(1:length(aK0_assump),function(i,aK0,assump){
                            assump[[i]] %*% aK0[[i]]
                },aK0 = curr.aK0,assump = aK0_assump)

                full.proj =(ProjectHarvest(Surv = invlogit(logit.curr.s.full), Harvpar = invlogit(logit.curr.H.full),Fec=exp(log.prop.f.full)#<- use proposal
                                , SRB = invlogit(logit.curr.SRB.full)
                                , aK0 = (curr.aK0.full), global = global, null = null, bl = exp(log.curr.b)    , period = proj.periods, nage = nage))



                if(sum(full.proj$Harvest < 0) > 0 || is.na(sum(full.proj$Harvest))
                     || is.nan(sum(full.proj$Harvest))) {
                        if(i > burn.in) {
                                pop.negative$fert.rate[j] =
                                        pop.negative$fert.rate[j] + 1/n.iter
                        }
                } else {

                    prop.aeri = ( getAerialCount( full.proj,A = invlogit(logit.curr.A.full),obsMat = Observations$AerialCount))
                    log.full.obs_f = log(getobsVitals(full.proj$Fec_obs,curr.proj$living,Observations$Fec))
                    logit.full.obs_s = logitf(getobsVitals(full.proj$Surv_obs,curr.proj$living,Observations$Surv))


                    # - Calculate log posterior of proposed vital under projection

        log.prop.posterior =
                log.post(f = log.prop.f #<- use proposal
                                             ,s = logit.curr.s
                                             ,SRB = logit.curr.SRB
                                             ,A = logit.curr.A
                                             ,H = logit.curr.H
                                             ,aK0 = curr.aK0
                                             
                                             ,estFec=estFec, estaK0=estaK0
                         
                                             ,prior.mean.f = log.prior.mean.f
                                             ,prior.mean.s = log.prior.mean.s
                                             ,prior.mean.SRB = log.prior.mean.SRB
                                             ,prior.mean.A = logit.prior.mean.A
                                             ,prior.mean.H = logit.prior.mean.H

                                             ,prior.var.f = prior.var.f
                                             ,prior.var.s = prior.var.s
                                             ,prior.var.SRB = prior.var.SRB
                                             ,prior.var.A = prior.var.A
                                             ,prior.var.H = prior.var.H
                                             
                                             ,alpha.f = al.f, beta.f = be.f
                                             ,alpha.s = al.s, beta.s = be.s
                                             ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                                             ,alpha.A = al.A, beta.A = be.A
                                             ,alpha.H = al.H, beta.H = be.H
                                             ,min.aK0 = min.aK0, max.aK0 = max.aK0


                                             ,sigmasq.f = curr.sigmasq.f
                                             ,sigmasq.s = curr.sigmasq.s
                                             ,sigmasq.SRB = curr.sigmasq.SRB
                                             



                                             ,log.like =
                                                        log.lhood_popu(
                                                                 n.census = Harv.data
                                                                ,n.hat = full.proj$Harvest#<- use proposal
                                                                ) +
                                                        log.lhood_popu(
                                                                 n.census = Aerial.data
                                                                ,n.hat = prop.aeri#<- use proposal
                                                                ) +
                                                        log.lhood_vital(f = log.full.obs_f
                                                                        ,s = logit.full.obs_s # observed survival and fecundity after adding DD.
                                                                        ,SRB = logit.curr.SRB
                                                                        ,estFec = estFec
                                                                        ,measure.f = log.measure.f
                                                                        ,measure.s = log.measure.s
                                                                        ,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = curr.sigmasq.f
                                                                        ,sigmasq.s = curr.sigmasq.s
                                                                        ,sigmasq.SRB = curr.sigmasq.SRB
                                                                        ,non.zero.fert = fert.rows)
                                                        #　both harvest and aerial count, as well as vital rates
                                             ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                                             )

                        
                        
                    #- Acceptance ratio
                    ar = acc.ra(log.prop = log.prop.posterior,
                                                         log.current = log.curr.posterior)

                    # - Move or stay
                    #.. stay if acceptance ratio 0, missing, infinity, etc.
                    if(is.na(ar) || is.nan(ar) || ar < 0) {
                        if(i > burn.in) ar.na$fert.rate[j] =
                                ar.na$fert.rate[j] + 1/n.iter
                    } else {
                        #.. if accept, update current fert rates, store proposed
                        #     rate, update current projection and count acceptance
                        if(runif(1) <= ar) {
                            if(i > burn.in) acc.count$fert.rate[j] =
                                    acc.count$fert.rate[j] + 1/n.iter
                            log.curr.f = log.prop.f
                            curr.proj = full.proj
                            curr.aeri = (prop.aeri)
                            log.curr.posterior = log.prop.posterior # change log curr
                        }
                        #.. if reject, leave current fert rates and projections
                        #     alone

                    } # close else after checking for ar=0, missing, inf

                } # close else after checking neg or zero population

            } # close loop over all age-spec fertility rates

            #.. Store proposed fertility rate matrix
            if(k %% 1 == 0 && k > 0) fert.rate.mcmc[k,] =
                    as.vector(exp(log.curr.f[fert.rows,]))
        }
        # pause 0519
            ##...... Survival ......##

            if(verb && identical(i%%1000, 0)) cat("\n", i, " Survival")

            # - Proposal

            #.. cycle through components
            for(j in 1:length(logit.curr.s)) {

                #.. make a matrix conformable w rate matrix
                ##     this is a strange structure inherite from popReconstruct
                ##     the logit.prop.s.mat renew every loop, do not know why they do this.
                logit.prop.s.mat =
                        matrix(0, nrow = nrow(logit.curr.s)
                                     ,ncol = ncol(logit.curr.s)) # this result depends on whether time-homo assumed.
                logit.prop.s.mat[j] = rnorm(1, 0, sqrt(prop.vars$surv.prop[j]))

                #.. make proposal
                logit.prop.s = logit.curr.s + logit.prop.s.mat

                #.. If proposal resulted in back-transformed s = 0 or 1, do
                #     nothing
                if(invlogit(logit.prop.s[j]) > 1 - s.tol ||
                     invlogit(logit.prop.s[j]) < s.tol) {
                    #.. leave current surv rates and projections
                    #     alone (simply do not propose
                    #     extreme survival probabilities)
                    s.out.tol[j] = s.out.tol[j] + 1/n.iter
                } else {

                    # - Run CCMP (project on the original scale)
                    #     ** Don't allow negative population; again, simply treat
                    #            this as if the proposal were never made
                        logit.prop.s.full = Surv_assump$age %*% logit.prop.s %*%Surv_assump$time
                        log.curr.f.full = Fec_assump$age %*% log.curr.f %*% Fec_assump$time
                        logit.curr.H.full = Harv_assump$age %*% logit.curr.H %*%Harv_assump$time
                        logit.curr.SRB.full = SRB_assump$age %*% logit.curr.SRB %*%SRB_assump$time
                        logit.curr.A.full = A_assump$age %*% logit.curr.A %*%A_assump$time
                        curr.aK0.full = lapply(1:length(aK0_assump),function(i,aK0,assump){
                            assump[[i]] %*% aK0[[i]]
                        },aK0 = curr.aK0,assump = aK0_assump)



                        full.proj =
                                (ProjectHarvest(Surv = invlogit(logit.prop.s.full)#<- use proposal
                                , Harvpar = invlogit(logit.curr.H.full),Fec=exp(log.curr.f.full), SRB = invlogit(logit.curr.SRB.full), aK0 = (curr.aK0.full), global = global, null = null, bl = exp(log.curr.b)    , period = proj.periods, nage = nage))


                        if(sum(full.proj$Harvest < 0) > 0 || is.na(sum(full.proj$Harvest))
                             || is.nan(sum(full.proj$Harvest))) {
                                if(i > burn.in) {
                                        pop.negative$surv.prop[j] =
                                                pop.negative$surv.prop[j] + 1/n.iter
                                }
                        } else {

                    prop.aeri = ( getAerialCount( full.proj,A = invlogit(logit.curr.A.full),obsMat = Observations$AerialCount))
                    log.full.obs_f = log(getobsVitals(full.proj$Fec_obs,curr.proj$living,Observations$Fec))
                    logit.full.obs_s = logitf(getobsVitals(full.proj$Surv_obs,curr.proj$living,Observations$Surv))


                        # - Calculate log posterior of proposed vital under projection
                                
        log.prop.posterior =
                log.post(f = log.curr.f
                                             ,s = logit.prop.s #<- use proposal
                                             ,SRB = logit.curr.SRB
                                             ,A = logit.curr.A
                                             ,H = logit.curr.H
                                             ,aK0 = curr.aK0
                                             
                                             ,estFec=estFec, estaK0=estaK0
                         
                                             ,prior.mean.f = log.prior.mean.f
                                             ,prior.mean.s = log.prior.mean.s
                                             ,prior.mean.SRB = log.prior.mean.SRB
                                             ,prior.mean.A = logit.prior.mean.A
                                             ,prior.mean.H = logit.prior.mean.H

                                             ,prior.var.f = prior.var.f
                                             ,prior.var.s = prior.var.s
                                             ,prior.var.SRB = prior.var.SRB
                                             ,prior.var.A = prior.var.A
                                             ,prior.var.H = prior.var.H
                                             
                                             ,alpha.f = al.f, beta.f = be.f
                                             ,alpha.s = al.s, beta.s = be.s
                                             ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                                             ,alpha.A = al.A, beta.A = be.A
                                             ,alpha.H = al.H, beta.H = be.H
                                             ,min.aK0 = min.aK0, max.aK0 = max.aK0


                                             ,sigmasq.f = curr.sigmasq.f
                                             ,sigmasq.s = curr.sigmasq.s
                                             ,sigmasq.SRB = curr.sigmasq.SRB
                                             



                                             ,log.like =
                                                        log.lhood_popu(
                                                                n.census = Harv.data
                                                                ,n.hat = full.proj$Harvest#<- use proposal
                                                                ) +
                                                        log.lhood_popu(
                                                                n.census = Aerial.data
                                                                ,n.hat = prop.aeri#<- use proposal
                                                                ) +
                                                        log.lhood_vital(f = log.full.obs_f
                                                                        ,s = logit.full.obs_s # observed survival and fecundity after adding DD.
                                                                        ,SRB = logit.curr.SRB
                                                                        ,estFec = estFec
                                                                        ,measure.f = log.measure.f
                                                                        ,measure.s = log.measure.s
                                                                        ,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = curr.sigmasq.f
                                                                        ,sigmasq.s = curr.sigmasq.s
                                                                        ,sigmasq.SRB = curr.sigmasq.SRB
                                                                        ,non.zero.fert = fert.rows)
                                                        #　both harvest and aerial count, as well as vital rates
                                             ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                                             )


                        #- Acceptance ratio
                        ar = acc.ra(log.prop = log.prop.posterior,
                                                            log.current = log.curr.posterior)

                        # - Move or stay
                        #.. stay if acceptance ratio 0, missing, infinity, etc.
                        if(is.na(ar) || is.nan(ar) || ar < 0) {
                            if(i > burn.in) ar.na$surv.prop[j] =
                                    ar.na$surv.prop[j] + 1/n.iter
                        } else {
                            #.. if accept, update current surv rates,
                            #     update current projection and count acceptance
                            if(runif(1) <= ar) {
                                if(i > burn.in) acc.count$surv.prop[j] =
                                        acc.count$surv.prop[j] + 1/n.iter
                                logit.curr.s = logit.prop.s
                                curr.proj = full.proj
                                curr.aeri = (prop.aeri)
                                log.curr.posterior = log.prop.posterior
                            }

                        } # close else{ after checking for undefined ar

                    } # close else{ after checking for negative pop

                } # close else{ after checking for s outside tol

            } # close loop over all age-spec survival probabilities

            #.. Store proposed survival probability matrix
            if(k %% 1 == 0 && k > 0) surv.prop.mcmc[k,] =
                as.vector(invlogit(logit.curr.s))





                    ##...... SRB ......##

            if(verb && identical(i%%1000, 0)) cat("\n", i, " SRB")

            # - Proposal

            #.. cycle through components
            for(j in 1:length(logit.curr.SRB)) {
                logit.prop.SRB.mat =
                        matrix(0, nrow = nrow(logit.curr.SRB)
                                     ,ncol = ncol(logit.curr.SRB)) # this result depends on whether time-homo assumed.
                logit.prop.SRB.mat[j] = rnorm(1, 0, sqrt(prop.vars$SRB[j]))

                #.. make proposal
                logit.prop.SRB = logit.curr.SRB + logit.prop.SRB.mat

                #.. If proposal resulted in back-transformed s = 0 or 1, do
                #     nothing
                if(invlogit(logit.prop.SRB[j]) > 1 - s.tol ||
                     invlogit(logit.prop.SRB[j]) < s.tol) {
                    #.. leave current surv rates and projections
                    #     alone (simply do not propose
                    #     extreme survival probabilities)
                    s.out.tol[j] = s.out.tol[j] + 1/n.iter
                } else {

                    #  (project on the original scale)
                    #     ** Don't allow negative population; again, simply treat
                    #            this as if the proposal were never made
                        logit.curr.s.full = Surv_assump$age %*% logit.curr.s %*%Surv_assump$time
                        log.curr.f.full = Fec_assump$age %*% log.curr.f %*% Fec_assump$time
                        logit.curr.H.full = Harv_assump$age %*% logit.curr.H %*%Harv_assump$time
                        logit.prop.SRB.full = SRB_assump$age %*% logit.prop.SRB %*%SRB_assump$time
                        logit.curr.A.full = A_assump$age %*% logit.curr.A %*%A_assump$time
                        curr.aK0.full = lapply(1:length(aK0_assump),function(i,aK0,assump){
                            assump[[i]] %*% aK0[[i]]
                        },aK0 = curr.aK0,assump = aK0_assump)




                        full.proj =
                                (ProjectHarvest(Surv = invlogit(logit.curr.s.full)
                                , Harvpar = invlogit(logit.curr.H.full),Fec=exp(log.curr.f.full), SRB = invlogit(logit.prop.SRB.full)#<- use proposal
                                , aK0 = (curr.aK0.full), global = global, null = null, bl = exp(log.curr.b)    , period = proj.periods, nage = nage))


                        if(sum(full.proj$Harvest < 0) > 0 || is.na(sum(full.proj$Harvest ))
                             || is.nan(sum(full.proj$Harvest ))) {
                                if(i > burn.in) {
                                        pop.negative$surv.prop[j] =
                                                pop.negative$surv.prop[j] + 1/n.iter
                                }
                        } else {

                                prop.aeri = ( getAerialCount( full.proj,A = invlogit(logit.curr.A.full),obsMat = Observations$AerialCount))
                                log.full.obs_f = log(getobsVitals(full.proj$Fec_obs,curr.proj$living,Observations$Fec))
                                logit.full.obs_s = logitf(getobsVitals(full.proj$Surv_obs,curr.proj$living,Observations$Surv))


                        # - Calculate log posterior of proposed vital under projection
                           
                                
      log.prop.posterior =
                log.post(f = log.curr.f
                                             ,s = logit.curr.s
                                             ,SRB = logit.prop.SRB #<- use proposal
                                             ,A = logit.curr.A
                                             ,H = logit.curr.H
                                             ,aK0 = curr.aK0
                                             
                                             ,estFec=estFec, estaK0=estaK0
                         
                                             ,prior.mean.f = log.prior.mean.f
                                             ,prior.mean.s = log.prior.mean.s
                                             ,prior.mean.SRB = log.prior.mean.SRB
                                             ,prior.mean.A = logit.prior.mean.A
                                             ,prior.mean.H = logit.prior.mean.H

                                             ,prior.var.f = prior.var.f
                                             ,prior.var.s = prior.var.s
                                             ,prior.var.SRB = prior.var.SRB
                                             ,prior.var.A = prior.var.A
                                             ,prior.var.H = prior.var.H
                                             
                                             ,alpha.f = al.f, beta.f = be.f
                                             ,alpha.s = al.s, beta.s = be.s
                                             ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                                             ,alpha.A = al.A, beta.A = be.A
                                             ,alpha.H = al.H, beta.H = be.H
                                             ,min.aK0 = min.aK0, max.aK0 = max.aK0


                                             ,sigmasq.f = curr.sigmasq.f
                                             ,sigmasq.s = curr.sigmasq.s
                                             ,sigmasq.SRB = curr.sigmasq.SRB
                                             



                                             ,log.like =
                                                        log.lhood_popu(
                                                                n.census = Harv.data
                                                                ,n.hat = full.proj$Harvest#<- use proposal
                                                                ) +
                                                        log.lhood_popu(
                                                                n.census = Aerial.data
                                                                ,n.hat = prop.aeri#<- use proposal
                                                                ) +
                                                        log.lhood_vital(f = log.full.obs_f
                                                                        ,s = logit.full.obs_s # observed survival and fecundity after adding DD.
                                                                        ,SRB = logit.curr.SRB
                                                                        ,estFec = estFec
                                                                        ,measure.f = log.measure.f
                                                                        ,measure.s = log.measure.s
                                                                        ,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = curr.sigmasq.f
                                                                        ,sigmasq.s = curr.sigmasq.s
                                                                        ,sigmasq.SRB = curr.sigmasq.SRB
                                                                        ,non.zero.fert = fert.rows)
                                                        #　both harvest and aerial count, as well as vital rates
                                             ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                                             )



                        #- Acceptance ratio
                        ar = acc.ra(log.prop = log.prop.posterior,
                                                            log.current = log.curr.posterior)

                        # - Move or stay
                        #.. stay if acceptance ratio 0, missing, infinity, etc.
                        if(is.na(ar) || is.nan(ar) || ar < 0) {
                            if(i > burn.in) ar.na$SRB[j] =
                                    ar.na$SRB[j] + 1/n.iter
                        } else {
                            #.. if accept, update current surv rates,
                            #     update current projection and count acceptance
                            if(runif(1) <= ar) {
                                if(i > burn.in) acc.count$SRB[j] =
                                        acc.count$SRB[j] + 1/n.iter
                                        logit.curr.SRB = logit.prop.SRB
                                        curr.proj = full.proj
                                        curr.aeri = (prop.aeri)
                                        log.curr.posterior = log.prop.posterior
                            }

                        } # close else{ after checking for undefined ar

                    } # close else{ after checking for negative pop

                } # close else{ after checking for s outside tol

            } # close loop over all age-spec survival probabilities

            #.. Store proposed survival probability matrix
            if(k %% 1 == 0 && k > 0) SRB.mcmc[k,] =
                as.vector(invlogit(logit.curr.SRB))



            ##...... Harvesting ......##

            if(verb && identical(i%%1000, 0)) cat("\n", i, " Harvesting")

            # - Proposal

            #.. cycle through components
            for(j in 1:length(logit.curr.H)) {

                #.. make a matrix conformable w rate matrix
                prop.H.mat =
                        0*logit.curr.H
                prop.H.mat[j] = rnorm(1, 0, sqrt(prop.vars$H[j])) # if need age-imspecific harvest, simply give a 1 by 1 start.H

                #.. make proposal
                logit.prop.H = logit.curr.H + prop.H.mat

            # - Run CCMP (project on the original scale)
            #     ** Don't allow negative population
                # - Run CCMP (project on the original scale)
                    #     ** Don't allow negative population; again, simply treat
                    #            this as if the proposal were never made
                        logit.curr.s.full = Surv_assump$age %*% logit.curr.s %*%Surv_assump$time
                        log.curr.f.full = Fec_assump$age %*% log.curr.f %*% Fec_assump$time
                        logit.prop.H.full = Harv_assump$age %*% logit.prop.H %*%Harv_assump$time
                        logit.curr.SRB.full = SRB_assump$age %*% logit.curr.SRB %*%SRB_assump$time
                        logit.curr.A.full = A_assump$age %*% logit.curr.A %*%A_assump$time
                        curr.aK0.full = lapply(1:length(aK0_assump),function(i,aK0,assump){
                            assump[[i]] %*% aK0[[i]]
                        },aK0 = curr.aK0,assump = aK0_assump)


                        full.proj =
                                (ProjectHarvest(Surv = invlogit(logit.curr.s.full), Harvpar =( invlogit(logit.prop.H.full))#<- use proposal
                                ,Fec=exp(log.curr.f.full), SRB = invlogit(logit.curr.SRB.full), aK0 = (curr.aK0.full), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))


                        if(sum(full.proj$Harvest  < 0) > 0 || is.na(sum(full.proj$Harvest ))
                             || is.nan(sum(full.proj$Harvest ))) {
                                if(i > burn.in) {
                                        pop.negative$surv.prop[j] =
                                                pop.negative$surv.prop[j] + 1/n.iter
                                }
                        } else {

                    prop.aeri = ( getAerialCount( full.proj,A = invlogit(logit.curr.A.full),obsMat = Observations$AerialCount))

                    log.full.obs_f = log(getobsVitals(full.proj$Fec_obs,curr.proj$living,Observations$Fec))
                    logit.full.obs_s = logitf(getobsVitals(full.proj$Surv_obs,curr.proj$living,Observations$Surv))

                
                                
        log.prop.posterior =
                log.post(f = log.curr.f
                                             ,s = logit.curr.s
                                             ,SRB = logit.curr.SRB
                                             ,A = logit.curr.A
                                             ,H = logit.prop.H #<- use proposal
                                             ,aK0 = curr.aK0
                                             
                                             ,estFec=estFec, estaK0=estaK0
                         
                                             ,prior.mean.f = log.prior.mean.f
                                             ,prior.mean.s = log.prior.mean.s
                                             ,prior.mean.SRB = log.prior.mean.SRB
                                             ,prior.mean.A = logit.prior.mean.A
                                             ,prior.mean.H = logit.prior.mean.H

                                             ,prior.var.f = prior.var.f
                                             ,prior.var.s = prior.var.s
                                             ,prior.var.SRB = prior.var.SRB
                                             ,prior.var.A = prior.var.A
                                             ,prior.var.H = prior.var.H
                                             
                                             ,alpha.f = al.f, beta.f = be.f
                                             ,alpha.s = al.s, beta.s = be.s
                                             ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                                             ,alpha.A = al.A, beta.A = be.A
                                             ,alpha.H = al.H, beta.H = be.H
                                             ,min.aK0 = min.aK0, max.aK0 = max.aK0


                                             ,sigmasq.f = curr.sigmasq.f
                                             ,sigmasq.s = curr.sigmasq.s
                                             ,sigmasq.SRB = curr.sigmasq.SRB
                                             



                                             ,log.like =
                                                        log.lhood_popu(
                                                                n.census = Harv.data
                                                                ,n.hat = full.proj$Harvest#<- use proposal
                                                                ) +
                                                        log.lhood_popu(
                                                                n.census = Aerial.data
                                                                ,n.hat = prop.aeri#<- use proposal
                                                                ) +
                                                        log.lhood_vital(f = log.full.obs_f
                                                                        ,s = logit.full.obs_s # observed survival and fecundity after adding DD.
                                                                        ,SRB = logit.curr.SRB
                                                                        ,estFec = estFec
                                                                        ,measure.f = log.measure.f
                                                                        ,measure.s = log.measure.s
                                                                        ,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = curr.sigmasq.f
                                                                        ,sigmasq.s = curr.sigmasq.s
                                                                        ,sigmasq.SRB = curr.sigmasq.SRB
                                                                        ,non.zero.fert = fert.rows)
                                                        #　both harvest and aerial count, as well as vital rates
                                             ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                                             )


                #- Acceptance ratio
                ar = acc.ra(log.prop = log.prop.posterior,
                                                     log.current = log.curr.posterior)

                # - Move or stay
                #.. stay if acceptance ratio 0, missing, infinity, etc.
                if(is.na(ar) || is.nan(ar) || ar < 0) {
                        if(i > burn.in) ar.na$H[j] =
                                ar.na$H[j] + 1/n.iter
                } else {
                        #.. if accept, update current vital rates, store proposed
                        #     rate, update current projection and count acceptance
                        if(runif(1) <= ar) {
                                if(i > burn.in) acc.count$H[j] =
                                        acc.count$H[j] + 1/n.iter
                                logit.curr.H = logit.prop.H
                                curr.proj = full.proj
                                curr.aeri = (prop.aeri)
                                log.curr.posterior = log.prop.posterior
                        }

                } # close else after checking for ar=na, nan, zero

            } # close else after checking for negative population

        } # close loop over all age-specific Harvest proportions

            #.. Store proposed Harvest proportion matrix
            if(k %% 1 == 0 && k > 0) H.mcmc[k,] = as.vector(invlogit(logit.curr.H))

            if(verb && identical(i%%1000, 0)) cat("\n", i, " Aerial Counts Detection")

            # - Proposal

            #.. cycle through components
            for(j in 1:length(logit.curr.A)) {

                #.. make a matrix conformable w rate matrix
                prop.A.mat =
                        0*logit.curr.A
                prop.A.mat[j] = rnorm(1, 0, sqrt(prop.vars$A[j])) # if need age-imspecific harvest, simply give a 1 by 1 start.H

                #.. make proposal
                logit.prop.A = logit.curr.A + prop.A.mat

            # - Run CCMP (project on the original scale)
            #     ** Don't allow negative population
                # - Run CCMP (project on the original scale)
                    #     ** Don't allow negative population; again, simply treat
                    #            this as if the proposal were never made
                logit.curr.s.full = Surv_assump$age %*% logit.curr.s %*%Surv_assump$time
                log.curr.f.full = Fec_assump$age %*% log.curr.f %*% Fec_assump$time
                logit.curr.H.full = Harv_assump$age %*% logit.curr.H %*%Harv_assump$time
                logit.curr.SRB.full = SRB_assump$age %*% logit.curr.SRB %*%SRB_assump$time
                logit.prop.A.full = A_assump$age %*% logit.prop.A %*%A_assump$time
                curr.aK0.full = lapply(1:length(aK0_assump),function(i,aK0,assump){
                    assump[[i]] %*% aK0[[i]]
                },aK0 = curr.aK0,assump = aK0_assump)


                full.proj =
                                (ProjectHarvest(Surv = invlogit(logit.curr.s.full), Harvpar = invlogit(logit.curr.H.full)
                                ,Fec=exp(log.curr.f.full), SRB = invlogit(logit.curr.SRB.full), aK0 = (curr.aK0.full), global = global, null = null, bl = exp(log.curr.b) , period = proj.periods, nage = nage))


                        if(sum(full.proj$Harvest  < 0) > 0 || is.na(sum(full.proj$Harvest ))
                             || is.nan(sum(full.proj$Harvest ))) {
                                if(i > burn.in) {
                                        pop.negative$A[j] =
                                                pop.negative$A[j] + 1/n.iter
                                }
                        } else {

                    prop.aeri = ( getAerialCount( full.proj,A = invlogit(logit.prop.A.full),obsMat = Observations$AerialCount)) #<- use proposal
                    log.full.obs_f = log(getobsVitals(full.proj$Fec_obs,curr.proj$living,Observations$Fec))
                    logit.full.obs_s = logitf(getobsVitals(full.proj$Surv_obs,curr.proj$living,Observations$Surv))


                                
        log.prop.posterior =
                log.post(f = log.curr.f
                                             ,s = logit.curr.s
                                             ,SRB = logit.curr.SRB
                                             ,A = logit.prop.A #<- use proposal
                                             ,H = logit.curr.H
                                             ,aK0 = curr.aK0
                                             
                                             ,estFec=estFec, estaK0=estaK0
                         
                                             ,prior.mean.f = log.prior.mean.f
                                             ,prior.mean.s = log.prior.mean.s
                                             ,prior.mean.SRB = log.prior.mean.SRB
                                             ,prior.mean.A = logit.prior.mean.A
                                             ,prior.mean.H = logit.prior.mean.H

                                             ,prior.var.f = prior.var.f
                                             ,prior.var.s = prior.var.s
                                             ,prior.var.SRB = prior.var.SRB
                                             ,prior.var.A = prior.var.A
                                             ,prior.var.H = prior.var.H
                                             
                                             ,alpha.f = al.f, beta.f = be.f
                                             ,alpha.s = al.s, beta.s = be.s
                                             ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                                             ,alpha.A = al.A, beta.A = be.A
                                             ,alpha.H = al.H, beta.H = be.H
                                             ,min.aK0 = min.aK0, max.aK0 = max.aK0


                                             ,sigmasq.f = curr.sigmasq.f
                                             ,sigmasq.s = curr.sigmasq.s
                                             ,sigmasq.SRB = curr.sigmasq.SRB
                                             



                                             ,log.like =
                                                        log.lhood_popu(
                                                                n.census = Harv.data
                                                                ,n.hat = full.proj$Harvest#<- use proposal
                                                                ) +
                                                        log.lhood_popu(
                                                                n.census = Aerial.data
                                                                ,n.hat = prop.aeri#<- use proposal
                                                                ) +
                                                        log.lhood_vital(f = log.full.obs_f
                                                                        ,s = logit.full.obs_s # observed survival and fecundity after adding DD.
                                                                        ,SRB = logit.curr.SRB
                                                                        ,estFec = estFec
                                                                        ,measure.f = log.measure.f
                                                                        ,measure.s = log.measure.s
                                                                        ,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = curr.sigmasq.f
                                                                        ,sigmasq.s = curr.sigmasq.s
                                                                        ,sigmasq.SRB = curr.sigmasq.SRB
                                                                        ,non.zero.fert = fert.rows)
                                                        #　both harvest and aerial count, as well as vital rates
                                             ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                                             )



                                
                                
                #- Acceptance ratio
                ar = acc.ra(log.prop = log.prop.posterior,
                                                     log.current = log.curr.posterior)

                # - Move or stay
                #.. stay if acceptance ratio 0, missing, infinity, etc.
                if(is.na(ar) || is.nan(ar) || ar < 0) {
                        if(i > burn.in) ar.na$A[j] =
                                ar.na$A[j] + 1/n.iter
                } else {
                        #.. if accept, update current vital rates, store proposed
                        #     rate, update current projection and count acceptance
                        if(runif(1) <= ar) {
                                if(i > burn.in) acc.count$A[j] =
                                        acc.count$A[j] + 1/n.iter
                                logit.curr.A = logit.prop.A
                                curr.proj = full.proj
                                curr.aeri = (prop.aeri)
                                log.curr.posterior = log.prop.posterior
                        }

                } # close else after checking for ar=na, nan, zero

            } # close else after checking for negative population

        } # close loop over all age-specific Harvest proportions

            #.. Store proposed Harvest proportion matrix
            if(k %% 1 == 0 && k > 0) A.mcmc[k,] = as.vector(invlogit(logit.curr.A))


            ##...... Carrying Capacity ......##
        if(estaK0){

            for(j in 1:length(start.aK0)){
                for(w in 1:length(curr.aK0[[j]])){
                    prop.aK0 = curr.aK0
                    prop.aK0[[j]][w] = curr.aK0[[j]][w] + rnorm(1, 0, sqrt(prop.vars$aK0[[j]]))

                    logit.curr.s.full = Surv_assump$age %*% logit.curr.s %*%Surv_assump$time
                    log.curr.f.full = Fec_assump$age %*% log.curr.f %*% Fec_assump$time
                    logit.curr.H.full = Harv_assump$age %*% logit.curr.H %*%Harv_assump$time
                    logit.curr.SRB.full = SRB_assump$age %*% logit.curr.SRB %*%SRB_assump$time
                    logit.curr.A.full = A_assump$age %*% logit.curr.A %*%A_assump$time
                    prop.aK0.full = lapply(1:length(aK0_assump),function(i,aK0,assump){
                        assump[[i]] %*% aK0[[i]]
                    },aK0 = prop.aK0,assump = aK0_assump)



                     full.proj =
                                (ProjectHarvest(Surv = invlogit(logit.curr.s.full), Harvpar = invlogit(logit.curr.H.full),Fec=exp(log.curr.f.full), SRB = invlogit(logit.curr.SRB.full), aK0 = prop.aK0.full#<- use proposal
                                , global = global, null = null, bl = exp(log.curr.b)    , period = proj.periods, nage = nage))


                     if(sum(full.proj$Harvest  < 0) > 0 || is.na(sum(full.proj$Harvest ))
                             || is.nan(sum(full.proj$Harvest ))) {
                                if(i > burn.in) {
                                        pop.negative$surv.prop[j] =
                                                pop.negative$surv.prop[j] + 1/n.iter
                                }
                        }
                     else {

                    prop.aeri = ( getAerialCount( full.proj,A = invlogit(logit.curr.A.full),obsMat = Observations$AerialCount))
                    log.full.obs_f = log(getobsVitals(full.proj$Fec_obs,curr.proj$living,Observations$Fec))
                    logit.full.obs_s = logitf(getobsVitals(full.proj$Surv_obs,curr.proj$living,Observations$Surv))


                             
         log.prop.posterior =
                log.post(f = log.curr.f
                                             ,s = logit.curr.s
                                             ,SRB = logit.curr.SRB
                                             ,A = logit.curr.A
                                             ,H = logit.curr.H
                                             ,aK0 = prop.aK0 #<- use proposal
                                             
                                             ,estFec=estFec, estaK0=estaK0
                         
                                             ,prior.mean.f = log.prior.mean.f
                                             ,prior.mean.s = log.prior.mean.s
                                             ,prior.mean.SRB = log.prior.mean.SRB
                                             ,prior.mean.A = logit.prior.mean.A
                                             ,prior.mean.H = logit.prior.mean.H

                                             ,prior.var.f = prior.var.f
                                             ,prior.var.s = prior.var.s
                                             ,prior.var.SRB = prior.var.SRB
                                             ,prior.var.A = prior.var.A
                                             ,prior.var.H = prior.var.H
                                             
                                             ,alpha.f = al.f, beta.f = be.f
                                             ,alpha.s = al.s, beta.s = be.s
                                             ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                                             ,alpha.A = al.A, beta.A = be.A
                                             ,alpha.H = al.H, beta.H = be.H
                                             ,min.aK0 = min.aK0, max.aK0 = max.aK0


                                             ,sigmasq.f = curr.sigmasq.f
                                             ,sigmasq.s = curr.sigmasq.s
                                             ,sigmasq.SRB = curr.sigmasq.SRB
                                             



                                             ,log.like =
                                                        log.lhood_popu(
                                                                n.census = Harv.data
                                                                ,n.hat = full.proj$Harvest#<- use proposal
                                                                ) +
                                                        log.lhood_popu(
                                                                n.census = Aerial.data
                                                                ,n.hat = prop.aeri#<- use proposal
                                                                ) +
                                                        log.lhood_vital(f = log.full.obs_f
                                                                        ,s = logit.full.obs_s # observed survival and fecundity after adding DD.
                                                                        ,SRB = logit.curr.SRB
									,estFec = estFec
									,measure.f = log.measure.f
									,measure.s = log.measure.s
									,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = curr.sigmasq.f
                                                                        ,sigmasq.s = curr.sigmasq.s
                                                                        ,sigmasq.SRB = curr.sigmasq.SRB
                                                                        ,non.zero.fert = fert.rows)
                                                        #　both harvest and aerial count, as well as vital rates
                                             ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                                             )



                            
                             
                    
                #- Acceptance ratio
                ar = acc.ra(log.prop = log.prop.posterior, log.current = log.curr.posterior)

                # - Move or stay
                #.. stay if acceptance ratio 0, missing, infinity, etc.
                if(is.na(ar) || is.nan(ar) || ar < 0) {
                        if(i > burn.in) ar.na$aK0[j] =
                                ar.na$aK0[j] + 1/n.iter
                } else {
                        #.. if accept, update current vital rates, store proposed
                        #     rate, update current projection and count acceptance
                        if(runif(1) <= ar) {
                                if(i > burn.in) acc.count$K0[j] =
                                        acc.count$aK0[j] + 1/n.iter
                                curr.aK0 = prop.aK0
                                curr.proj = full.proj
                                curr.aeri=prop.aeri
                                log.curr.posterior = log.prop.posterior
                        }

                } # close else after checking for ar=na, nan, zero

        } # close else after checking for negative population
        }
                }
            }
            #.. Store proposed aK0 matrix
            if(k %% 1 == 0 && k > 0 && estaK0){
                aK0.Fec.mcmc[k,] = as.vector((curr.aK0[[1]]))
                aK0.Surv.mcmc[k,] = as.vector((curr.aK0[[2]]))
                aK0.midPopulation.mcmc[k,] = as.vector((curr.aK0[[3]]))

            }


            ##...... Baseline population ......##

            if(verb && identical(i%%1000, 0)) cat("\n", i, " Baseline")

            # - Proposal

            #.. cycle through components (never update last
            #     value as this affects years beyond the estimation period)
            for(j in 1:length(log.curr.b)) {

            #.. make a matrix conformable w rate matrix
            #log.prop.b.mat = matrix(0, nrow = nrow(log.curr.b), ncol = 1)
                log.prop.b.mat = 0 * log.curr.b
                log.prop.b.mat[j] = rnorm(1, 0, sqrt(prop.vars$baseline.pop.count[j]))

            #.. make proposal
                log.prop.b = log.curr.b + log.prop.b.mat

            # - Run CCMP (project on the original scale)
            #     ** Don't allow negative population
                logit.curr.s.full = Surv_assump$age %*% logit.curr.s %*%Surv_assump$time
                log.curr.f.full = Fec_assump$age %*% log.curr.f %*% Fec_assump$time
                logit.curr.H.full = Harv_assump$age %*% logit.curr.H %*%Harv_assump$time
                logit.curr.SRB.full = SRB_assump$age %*% logit.curr.SRB %*%SRB_assump$time
                logit.curr.A.full = A_assump$age %*% logit.curr.A %*%A_assump$time
                curr.aK0.full = lapply(1:length(aK0_assump),function(i,aK0,assump){
                     assump[[i]] %*% aK0[[i]]
                },aK0 = curr.aK0,assump = aK0_assump)




                full.proj =(ProjectHarvest(Surv = invlogit(logit.curr.s.full), Harvpar = invlogit(logit.curr.H.full), SRB = invlogit(logit.curr.SRB.full),Fec=exp(log.curr.f.full), aK0 = (curr.aK0.full), global = global, null = null, bl = exp(log.prop.b)    #<- use proposal
                                , period = proj.periods, nage = nage))


            if(sum(full.proj$Harvest  < 0) > 0 || is.na(sum(full.proj$Harvest ))
                 || is.nan(sum(full.proj$Harvest ))) {
                if(i > burn.in) {
                    pop.negative$baseline.count[j] =
                            pop.negative$baseline.count[j] + 1/n.iter
                }
            } else {

                    prop.aeri = ( getAerialCount( full.proj,A = invlogit(logit.curr.A.full),obsMat = Observations$AerialCount))
                    log.full.obs_f = log(getobsVitals(full.proj$Fec_obs,curr.proj$living,Observations$Fec))
                    logit.full.obs_s = logitf(getobsVitals(full.proj$Surv_obs,curr.proj$living,Observations$Surv))


                # - Calculate log posterior of proposed vital under projection

        log.prop.posterior =
                log.post(f = log.curr.f
                                             ,s = logit.curr.s
                                             ,SRB = logit.curr.SRB
                                             ,A = logit.curr.A
                                             ,H = logit.curr.H
                                             ,aK0 = curr.aK0
                                             
                                             ,estFec=estFec, estaK0=estaK0
                         
                                             ,prior.mean.f = log.prior.mean.f
                                             ,prior.mean.s = log.prior.mean.s
                                             ,prior.mean.SRB = log.prior.mean.SRB
                                             ,prior.mean.A = logit.prior.mean.A
                                             ,prior.mean.H = logit.prior.mean.H

                                             ,prior.var.f = prior.var.f
                                             ,prior.var.s = prior.var.s
                                             ,prior.var.SRB = prior.var.SRB
                                             ,prior.var.A = prior.var.A
                                             ,prior.var.H = prior.var.H
                                             
                                             ,alpha.f = al.f, beta.f = be.f
                                             ,alpha.s = al.s, beta.s = be.s
                                             ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                                             ,alpha.A = al.A, beta.A = be.A
                                             ,alpha.H = al.H, beta.H = be.H
                                             ,min.aK0 = min.aK0, max.aK0 = max.aK0


                                             ,sigmasq.f = curr.sigmasq.f
                                             ,sigmasq.s = curr.sigmasq.s
                                             ,sigmasq.SRB = curr.sigmasq.SRB
                                             



                                             ,log.like =
                                                        log.lhood_popu(
                                                                n.census = Harv.data
                                                                ,n.hat = full.proj$Harvest#<- use proposal
                                                                ) +
                                                        log.lhood_popu(
                                                                n.census = Aerial.data
                                                                ,n.hat = prop.aeri#<- use proposal
                                                                ) +
                                                        log.lhood_vital(f = log.full.obs_f
                                                                        ,s = logit.full.obs_s # observed survival and fecundity after adding DD.
                                                                        ,SRB = logit.curr.SRB
                                                                        ,estFec = estFec
                                                                        ,measure.f = log.measure.f
                                                                        ,measure.s = log.measure.s
                                                                        ,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = curr.sigmasq.f
                                                                        ,sigmasq.s = curr.sigmasq.s
                                                                        ,sigmasq.SRB = curr.sigmasq.SRB
                                                                        ,non.zero.fert = fert.rows)
                                                        #　both harvest and aerial count, as well as vital rates
                                             ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                                             )



                    
                    
                    
                 #- Acceptance ratio
                ar = acc.ra(log.prop = log.prop.posterior,
                                                     log.current = log.curr.posterior)

                # - Move or stay
                #.. stay if acceptance ratio 0, missing, infinity, etc.
                if(is.na(ar) || is.nan(ar) || ar < 0) {
                        if(i > burn.in) ar.na$baseline.count[j] =
                                ar.na$baseline.count[j] + 1/n.iter
                } else {
                        #.. if accept, update current mig rates, store proposed
                        #     rate, update current projection and count acceptance
                        if(runif(1) <= ar) {
                                if(i > burn.in) acc.count$baseline.count[j] =
                                        acc.count$baseline.count[j] + 1/n.iter
                                        log.curr.b = log.prop.b
                                        curr.proj = full.proj
                                        curr.aeri = prop.aeri
                                        log.curr.posterior = log.prop.posterior
                        } #.. if reject, leave current fert rates and projections
                        #     alone, store current rate

                } # close else after checking for ar=na, nan, zero

        } # close else after checking for negative population

    } # close loop over all age-specific baseline counts

            #.. Store proposed baseline count matrix
            if(k %% 1 == 0 && k > 0) baseline.count.mcmc[k,] =
                    as.vector(exp(log.curr.b))

# TEST stop here 10/26/2018
            ## ------- Variance Updates ------- ##
            log.curr.obs_f = log(getobsVitals(curr.proj$Fec_obs,curr.proj$living,Observations$Fec))
            logit.curr.obs_s = logitf(getobsVitals(curr.proj$Surv_obs,curr.proj$living,Observations$Surv))


            if(verb && identical(i%%1000, 0)) cat("\n", i, " Variances")

            ##...... Fertility rate ......##
            if(estFec){ # if not est Fer, this is not needed
                prop.sigmasq.f = rinvGamma(1, al.f + length(measure.f[fert.rows,])/2-sum(is.na(measure.f))/2,be.f + 0.5*sum((log.curr.f[fert.rows,] -log.measure.f[fert.rows,])^2,na.rm = T))

                # - Calculate log posterior of proposed vital under projection
        log.prop.posterior =
                log.post(f = log.curr.f
                                             ,s = logit.curr.s
                                             ,SRB = logit.curr.SRB
                                             ,A = logit.curr.A
                                             ,H = logit.curr.H
                                             ,aK0 = curr.aK0
                                             
                                             ,estFec=estFec, estaK0=estaK0
                         
                                             ,prior.mean.f = log.prior.mean.f
                                             ,prior.mean.s = log.prior.mean.s
                                             ,prior.mean.SRB = log.prior.mean.SRB
                                             ,prior.mean.A = logit.prior.mean.A
                                             ,prior.mean.H = logit.prior.mean.H

                                             ,prior.var.f = prior.var.f
                                             ,prior.var.s = prior.var.s
                                             ,prior.var.SRB = prior.var.SRB
                                             ,prior.var.A = prior.var.A
                                             ,prior.var.H = prior.var.H
                                             
                                             ,alpha.f = al.f, beta.f = be.f
                                             ,alpha.s = al.s, beta.s = be.s
                                             ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                                             ,alpha.A = al.A, beta.A = be.A
                                             ,alpha.H = al.H, beta.H = be.H
                                             ,min.aK0 = min.aK0, max.aK0 = max.aK0


                                             ,sigmasq.f = prop.sigmasq.f #<- use proposal
                                             ,sigmasq.s = curr.sigmasq.s
                                             ,sigmasq.SRB = curr.sigmasq.SRB
                                             



                                             ,log.like =
                                                        log.lhood_popu(
                                                                n.census = Harv.data
                                                                ,n.hat = curr.proj$Harvest#<- use current
                                                                ) +
                                                        log.lhood_popu(
                                                                n.census = Aerial.data
                                                                ,n.hat = curr.aeri#<- use current
                                                                ) +
                                                        log.lhood_vital(f = log.curr.obs_f
                                                                        ,s = log.curr.obs_s
                                                                        ,SRB = logit.curr.SRB
                                                                        ,estFec = estFec
                                                                        ,measure.f = log.measure.f
                                                                        ,measure.s = log.measure.s
                                                                        ,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = prop.sigmasq.f #<- use proposal
                                                                        ,sigmasq.s = curr.sigmasq.s
                                                                        ,sigmasq.SRB = curr.sigmasq.SRB
                                                                        ,non.zero.fert = fert.rows)
                                                        #　both harvest and aerial count, as well as vital rates
                                             ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                                             )




                    

            #- Acceptance ratio
            ar = acc.ra.var(log.prop.post = log.prop.posterior
                                                         ,log.curr.post = log.curr.posterior
                                                         ,log.prop.var = dinvGamma(prop.sigmasq.f
                                                            ,al.f + length(measure.f[fert.rows,])/2
                                                            ,be.f + 0.5*sum((log.curr.f[fert.rows,] - log.measure.f[fert.rows,])^2)
                                                            ,log = TRUE)
                                                         ,log.curr.var = dinvGamma(curr.sigmasq.f
                                                            ,al.f + length(measure.f[fert.rows,])/2
                                                            ,be.f + 0.5*sum((log.curr.f[fert.rows,] - log.measure.f[fert.rows,])^2)
                                                            ,log = TRUE)
                                                         )

                # - Move or stay
                #.. stay if acceptance ratio 0, missing, infinity, etc.
                if(is.na(ar) || is.nan(ar) || ar < 0) {
                        if(i > burn.in) ar.na$sigmasq.f =
                                ar.na$sigmasq.f + 1/n.iter
                } else {
                        #.. if accept, update current, store proposed
                        #     and count acceptance
                        if(runif(1) <= ar) {
                                if(i > burn.in) acc.count$sigmasq.f =
                                        acc.count$sigmasq.f + 1/n.iter
                                curr.sigmasq.f = prop.sigmasq.f
                                log.curr.posterior = log.prop.posterior
                        } #.. if reject, leave current and posterior
                } # close else after checking for ar=na, nan, zero

            if(k %% 1 == 0 && k > 0) variances.mcmc[k,"fert.rate.var"] = curr.sigmasq.f
        }
            ##...... Survival Proportion ......##

            prop.sigmasq.s =
                rinvGamma(1, al.s + length(measure.s)/2-sum(is.na(measure.s))/2,
                                    be.s +
                                        0.5*sum((logit.curr.s - logit.measure.s)^2,na.rm = T))
        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior =
                log.post(f = log.curr.f
                                             ,s = logit.curr.s
                                             ,SRB = logit.curr.SRB
                                             ,A = logit.curr.A
                                             ,H = logit.curr.H
                                             ,aK0 = curr.aK0
                                             
                                             ,estFec=estFec, estaK0=estaK0
                         
                                             ,prior.mean.f = log.prior.mean.f
                                             ,prior.mean.s = log.prior.mean.s
                                             ,prior.mean.SRB = log.prior.mean.SRB
                                             ,prior.mean.A = logit.prior.mean.A
                                             ,prior.mean.H = logit.prior.mean.H

                                             ,prior.var.f = prior.var.f
                                             ,prior.var.s = prior.var.s
                                             ,prior.var.SRB = prior.var.SRB
                                             ,prior.var.A = prior.var.A
                                             ,prior.var.H = prior.var.H
                                             
                                             ,alpha.f = al.f, beta.f = be.f
                                             ,alpha.s = al.s, beta.s = be.s
                                             ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                                             ,alpha.A = al.A, beta.A = be.A
                                             ,alpha.H = al.H, beta.H = be.H
                                             ,min.aK0 = min.aK0, max.aK0 = max.aK0


                                             ,sigmasq.f = curr.sigmasq.f
                                             ,sigmasq.s = prop.sigmasq.s #<- use proposal
                                             ,sigmasq.SRB = curr.sigmasq.SRB
                                             



                                             ,log.like =
                                                        log.lhood_popu(
                                                                n.census = Harv.data
                                                                ,n.hat = curr.proj$Harvest#<- use current
                                                                ) +
                                                        log.lhood_popu(
                                                                n.census = Aerial.data
                                                                ,n.hat = curr.aeri#<- use current
                                                                ) +
                                                        log.lhood_vital(f = log.curr.obs_f
                                                                        ,s = log.curr.obs_s
                                                                        ,SRB = logit.curr.SRB
                                                                        ,estFec = estFec
                                                                        ,measure.f = log.measure.f
                                                                        ,measure.s = log.measure.s
                                                                        ,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = curr.sigmasq.f
                                                                        ,sigmasq.s = prop.sigmasq.s #<- use proposal
                                                                        ,sigmasq.SRB = curr.sigmasq.SRB
                                                                        ,non.zero.fert = fert.rows)
                                                        #　both harvest and aerial count, as well as vital rates
                                             ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                                             )

            #- Acceptance ratio
            ar = acc.ra.var(log.prop.post = log.prop.posterior
                                                         ,log.curr.post = log.curr.posterior
                                                         ,log.prop.var = dinvGamma(prop.sigmasq.s
                                                            ,al.s + length(measure.s)/2
                                                            ,be.s + 0.5*sum((logit.curr.s - logit.measure.s)^2)
                                                            ,log = TRUE)
                                                         ,log.curr.var = dinvGamma(curr.sigmasq.s
                                                            ,al.s + length(measure.s)/2
                                                            ,be.s + 0.5*sum((logit.curr.s - logit.measure.s)^2)
                                                            ,log = TRUE)
                                                         )

                # - Move or stay
                #.. stay if acceptance ratio 0, missing, infinity, etc.
                if(is.na(ar) || is.nan(ar) || ar < 0) {
                        if(i > burn.in) ar.na$sigmasq.s =
                                ar.na$sigmasq.s + 1/n.iter
                } else {
                        #.. if accept, update current, store proposed
                        #     and count acceptance
                        if(runif(1) <= ar) {
                                if(i > burn.in) acc.count$sigmasq.s =
                                        acc.count$sigmasq.s + 1/n.iter
                                curr.sigmasq.s = prop.sigmasq.s
                                log.curr.posterior = log.prop.posterior
                        } #.. if reject, leave current and posterior
                } # close else after checking for ar=na, nan, zero

            if(k %% 1 == 0 && k > 0) variances.mcmc[k,"surv.prop.var"] = curr.sigmasq.s


            ##...... Sex Ratio at Birth ......##
            prop.sigmasq.SRB =
                rinvGamma(1, al.SRB + length(measure.SRB)/2-sum(is.na(measure.SRB))/2,be.SRB + 0.5*sum((logit.curr.SRB - logit.measure.SRB)^2,na.rm = T))

                # - Calculate log posterior of proposed vital under projection
        log.prop.posterior =
                log.post(f = log.curr.f
                                             ,s = logit.curr.s
                                             ,SRB = logit.curr.SRB
                                             ,A = logit.curr.A
                                             ,H = logit.curr.H
                                             ,aK0 = curr.aK0
                                             
                                             ,estFec=estFec, estaK0=estaK0
                         
                                             ,prior.mean.f = log.prior.mean.f
                                             ,prior.mean.s = log.prior.mean.s
                                             ,prior.mean.SRB = log.prior.mean.SRB
                                             ,prior.mean.A = logit.prior.mean.A
                                             ,prior.mean.H = logit.prior.mean.H

                                             ,prior.var.f = prior.var.f
                                             ,prior.var.s = prior.var.s
                                             ,prior.var.SRB = prior.var.SRB
                                             ,prior.var.A = prior.var.A
                                             ,prior.var.H = prior.var.H
                                             
                                             ,alpha.f = al.f, beta.f = be.f
                                             ,alpha.s = al.s, beta.s = be.s
                                             ,alpha.SRB = al.SRB, beta.SRB = be.SRB
                                             ,alpha.A = al.A, beta.A = be.A
                                             ,alpha.H = al.H, beta.H = be.H
                                             ,min.aK0 = min.aK0, max.aK0 = max.aK0


                                             ,sigmasq.f = curr.sigmasq.f
                                             ,sigmasq.s = curr.sigmasq.s
                                             ,sigmasq.SRB = prop.sigmasq.SRB #<- use proposal
                                             



                                             ,log.like =
                                                        log.lhood_popu(
                                                                n.census = Harv.data
                                                                ,n.hat = curr.proj$Harvest
                                                                ) +
                                                        log.lhood_popu(
                                                                n.census = Aerial.data
                                                                ,n.hat = curr.aeri
                                                                ) +
                                                        log.lhood_vital(f = log.curr.obs_f
                                                                        ,s = log.curr.obs_s
                                                                        ,SRB = logit.curr.SRB
                                                                        ,estFec = estFec
                                                                        ,measure.f = log.measure.f
                                                                        ,measure.s = log.measure.s
                                                                        ,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = curr.sigmasq.f
                                                                        ,sigmasq.s = curr.sigmasq.s
                                                                        ,sigmasq.SRB = prop.sigmasq.SRB #<- use proposal
                                                                        ,non.zero.fert = fert.rows)
                                                        #　both harvest and aerial count, as well as vital rates
                                             ,non.zero.fert = fert.rows # tell algorithm where the fert has to be 0
                                             )





            #- Acceptance ratio
            ar = acc.ra.var(log.prop.post = log.prop.posterior
                                                         ,log.curr.post = log.curr.posterior
                                                         ,log.prop.var = dinvGamma(prop.sigmasq.SRB
                                                            ,al.SRB + length(measure.SRB)/2
                                                            ,be.SRB + 0.5*sum((logit.curr.SRB - logit.measure.SRB)^2)
                                                            ,log = TRUE)
                                                         ,log.curr.var = dinvGamma(curr.sigmasq.SRB
                                                            ,al.SRB + length(measure.SRB)/2
                                                            ,be.SRB + 0.5*sum((logit.curr.SRB - logit.measure.SRB)^2)
                                                            ,log = TRUE)
                                                         )

                # - Move or stay
                #.. stay if acceptance ratio 0, missing, infinity, etc.
                if(is.na(ar) || is.nan(ar) || ar < 0) {
                        if(i > burn.in) ar.na$sigmasq.SRB =
                                ar.na$sigmasq.SRB + 1/n.iter
                } else {
                        #.. if accept, update current, store proposed
                        #     and count acceptance
                        if(runif(1) <= ar) {
                                if(i > burn.in) acc.count$sigmasq.SRB =
                                        acc.count$sigmasq.s + 1/n.iter
                                curr.sigmasq.SRB = prop.sigmasq.SRB
                                log.curr.posterior = log.prop.posterior
                        } #.. if reject, leave current and posterior
                } # close else after checking for ar=na, nan, zero

            if(k %% 1 == 0 && k > 0) variances.mcmc[k,"SRB.var"] = curr.sigmasq.SRB



            ## ------- Store current population ------- ##
                logit.curr.s.full = Surv_assump$age %*% logit.curr.s %*%Surv_assump$time
                log.curr.f.full = Fec_assump$age %*% log.curr.f %*% Fec_assump$time
                logit.curr.H.full = Harv_assump$age %*% logit.curr.H %*%Harv_assump$time
                logit.curr.SRB.full = SRB_assump$age %*% logit.curr.SRB %*%SRB_assump$time
                logit.curr.A.full = A_assump$age %*% logit.curr.A %*%A_assump$time
                curr.aK0.full = lapply(1:length(aK0_assump),function(i,aK0,assump){
                            assump[[i]] %*% aK0[[i]]
                },aK0 = curr.aK0,assump = aK0_assump)


                full.proj = (ProjectHarvest(Surv = invlogit(logit.curr.s.full), Harvpar = invlogit(logit.curr.H.full),Fec=exp(log.curr.f.full), SRB = invlogit(logit.curr.SRB.full), aK0 = (curr.aK0.full), global = global, null = null, bl = exp(log.curr.b)    , period = proj.periods, nage = nage))

                full.aeri = getAerialCount( full.proj,A = invlogit(logit.curr.A.full))
                
                log.full.obs_f = log(getobsVitals(full.proj$Fec_obs,curr.proj$living,Observations$Fec))
                logit.full.obs_s = logitf(getobsVitals(full.proj$Surv_obs,curr.proj$living,Observations$Surv))


                #full.living = (1-invlogit(logit.curr.H.full))*full.proj/(invlogit(logit.curr.H.full))

                if(k %% 1 == 0 && k > 0){
                    lx.mcmc[k,] =
                    as.vector(full.proj$Harvest) # to delete all age class' baseline count, because of as.vector,thus need to do like this
                    ae.mcmc[k,] = as.vector(full.aeri)
                    living.mcmc[k,] = as.vector(full.proj$Living)
                    log.like.mcmc[k,] = log.lhood_popu(
                                                                n.census = Harv.data
                                                                ,n.hat = curr.proj$Harvest
                                                                ) +
                                                        log.lhood_popu(
                                                                n.census = Aerial.data
                                                                ,n.hat = curr.aeri
                                                                ) +
                                                        log.lhood_vital(f = log.full.obs_f
                                                                        ,s = logit.full.obs_s # observed survival and fecundity after adding DD.
                                                                        ,SRB = logit.curr.SRB
                                                                        ,estFec = estFec
                                                                        ,measure.f = log.measure.f
                                                                        ,measure.s = log.measure.s
                                                                        ,measure.SRB = log.measure.SRB
                                                                        ,sigmasq.f = curr.sigmasq.f
                                                                        ,sigmasq.s = curr.sigmasq.s
                                                                        ,sigmasq.SRB = curr.sigmasq.SRB
                                                                        ,non.zero.fert = fert.rows)

                  }

                if(verb && identical(i%%1000, 0)) cat("\n\n")

            } # Ends outer-most loop

        ## ......... End Loop ........ ##
        #...............................#
    cat("\n","done","\n","model checking...")
    # calculate DIC
    mcmc.objs = list(survival.mcmc = surv.prop.mcmc
                     ,SRB.mcmc = SRB.mcmc
                     ,aerial.detection.mcmc = A.mcmc
                     ,H.mcmc = H.mcmc
                     ,living.mcmc = living.mcmc

                     ,baseline.count.mcmc = baseline.count.mcmc
                     ,harvest.mcmc = lx.mcmc
                     ,aerial.count.mcmc = ae.mcmc
                     ,variances.mcmc = variances.mcmc)
    mean.vital = lapply(mcmc.objs,function(kk){
        apply(as.matrix(kk),2,point.est)
    })
    if(estaK0){
        mcmc.objs$invK0.Fec = aK0.Fec.mcmc
        mean.vital$invK0.Fec = apply( as.matrix( aK0.Fec.mcmc),2,point.est)
        mcmc.objs$invK0.Surv = aK0.Surv.mcmc
        mean.vital$invK0.Surv = apply( as.matrix( aK0.Surv.mcmc),2,point.est)
        mcmc.objs$invK0.midPopulation.mcmc = aK0.midPopulation.mcmc
        mean.vital$invK0.midPopulation.mcmc = apply( as.matrix( aK0.midPopulation.mcmc),2,point.est)

    }
    else {mean.vital$invK0.mcmc = c(0,0)}
    if(estFec){
        mcmc.objs$fecundity.mcmc = fert.rate.mcmc
        mean.vital$fecundity.mcmc = apply( as.matrix( fert.rate.mcmc),2,point.est)
    }
    else{mean.vital$fecundity.mcmc = start.f}

    pD_Gelman04 = 2 * var(log.like.mcmc)
    DIC_Gelman04 = -2* mean(log.like.mcmc) + (pD_Gelman04)
    DIC = list(pD_Gelman04,DIC_Gelman04)
    names(DIC) = c("pD_Gelman04","DIC_Gelman04")

    ## abs_dif
    abs_dif = abs(mean.vital$harvest.mcmc-as.vector(Harv.data))
    mean_abs_dif_harv = mean(abs_dif,na.rm = T)
    se_abs_dif_harv = sd(abs_dif,na.rm = T)/(sqrt(proj.periods))

    abs_dif_aerial = abs(mean.vital$aerial.count.mcmc-as.vector(Aerial.data))
    mean_abs_dif_ae = mean(abs_dif_aerial,na.rm = T)
    se_abs_dif_ae = sd(abs_dif_aerial,na.rm = T)/sqrt(proj.periods)

    #obs_f = (getobsVitals(full.proj$Fec_obs,curr.proj$living,Observations$Fec))
    #obs_s = (getobsVitals(full.proj$Surv_obs,curr.proj$living,Observations$Surv))

                                         
    ## sd
    sd_counts = lapply(list(harvest = lx.mcmc,aerial.count = ae.mcmc)
                       ,function(kk){
                            apply(kk,2,sd)
                       })
    mean_sd_counts = lapply(sd_counts,mean,na.rm = T)

    ## precision

    model.checking = list(
        DIC=DIC,
        absolute.difference = list(MAD.harv = mean_abs_dif_harv
                                   ,SEAD.harv = se_abs_dif_harv
                                   ,MAD.aerial = mean_abs_dif_ae
                                   ,SEAD.aerial = se_abs_dif_ae)
                                   ,sd = mean_sd_counts
    )

        ## ---------- Output --------- ##

        #cat("inital values", "\n\n")
        #.. initial values
        start.vals = list(fecundity= start.f
                          ,survrvival = start.s
                          ,SRB = start.SRB
                          ,H = start.H
                          ,Aerial.detection = start.A
                          ,K0 = start.aK0
                          ,baseline.count = start.b
                          ,start.sigmasq.f = start.sigmasq.f
                          ,start.sigmasq.s = start.sigmasq.s
                          ,start.sigmasq.SRB = start.sigmasq.SRB
                          ,start.sigmasq.A = start.sigmasq.A
                          ,start.sigmasq.H = start.sigmasq.H
                         #,start.sigmasq.aK0 = start.sigmasq.aK0

                          ,Harv.data = Harv.data
                          ,Aerial.data = Aerial.data)

        #.. fixed parameters
        fixed.params = list(alpha.fecundity = al.f
                           ,beta.fecundity = be.f
                           ,alpha.survival = al.s
                           ,beta.survival = be.s
                           ,alpha.SRB = al.SRB
                           ,beta.SRB = be.SRB
                           ,alpha.aerial.det = al.A
                           ,beta.aerial.det = be.A
                           ,alpha.Harvest = al.H
                           ,beta.Hervest = be.H
                           ,alpha.1overK = min.aK0
                           ,beta.1overK = max.aK0

                           ,measure.fert.rate = measure.f
                           ,measure.surv.prop = measure.s
                           ,measure.sRB = measure.sRB
                           ,prior.mean.Aerial.detection = prior.mean.A
                           ,prior.mean.Harvest.proportion = prior.mean.H

                           ,measure.baseline.count = measure.b
                           ,prior.mean.Harv.data = Harv.data
                           ,Aerial.data = Aerial.data
                           ,Assumptions = Assumptions
                           ,point.est = point.est)


        #.. algorithm statistics
        alg.stats =
                list(acceptance.proportions = acc.count
                         ,pop.went.neg = pop.negative
                         ,acc.rat.na = ar.na
                         ,surv.outside.tol = s.out.tol
                         ,run.time = proc.time() - ptm
                         )

        #.. algorithm parameters
        alg.params = list(prop.vars = prop.vars
                                       ,vital.transformations = list(fecundity = "log"
                                       ,survival = "logit",SRB = "logit",aerial.detection="logit", H = "logit", invK0="identical"
                                       ,baseline.count = "log"
                                       ,harvest.count = "id"
                                       ,aerial.count = "id")
                                       ,projection.periods = proj.periods
                                       ,age.gp.size = nage
                                       ,non.zero.fert.rows = fert.rows
                                       ,surv.tolerance = s.tol
                                       ,burn.in = burn.in
                                       ,iters = n.iter)


        #.. results
        ret.list = list(mcmc.objs = mcmc.objs
                         ,log.like.mcmc = log.like.mcmc
                         ,alg.stats = alg.stats
                         ,model.checking = model.checking
                         ,fixed.params = fixed.params
                         ,start.vals = start.vals
                         ,alg.params = alg.params)
        class(ret.list) = "ReCAP_sample"
        cat("done \n","all done \n")
        return(ret.list)
}
