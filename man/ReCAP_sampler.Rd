\name{ReCAP_sampler}
\alias{ReCAP_sampler}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrated Population Model for Harvesting and Aerial Counting data under Bayesian Framework
}
\description{
Main sampler of the posterior distribution of vital rates and harvest count in \code{ReCAP}.
}
\usage{
ReCAP_sampler(Harv.data, Aerial.data
, nage, mean.f, mean.s
, mean.SRB, mean.H
, mean.A, mean.b = Harv.data[, 1]
, n.iter = 50000, burn.in = 5000, thin.by = 10
, al.f = 1, be.f = 0.05, al.s = 1, be.s = 0.1
, al.SRB = 1, be.SRB = 0.05
, min.aK0 = list(matrix(-0.001, nage[1], 1), matrix(-0.001, sum(nage), 1), 0)
, max.aK0 = list(matrix(0.001, nage[1], 1), matrix(0.001, sum(nage), 1), 500)
, al.H = 1, be.H = 0.05, al.A = 1, be.A = 0.05
, Assumptions = list()
, start.sigmasq.f = 0.05, start.sigmasq.s = 0.05
, start.sigmasq.SRB = 0.05, start.sigmasq.H = 0.05
, start.sigmasq.A = 0.05
, prop.vars = list()
, proj.periods = (ncol(Harv.data) - 1), estFec = T
, estaK0 = F
, aK0 = list(matrix(0, nage[1], 1), matrix(0, sum(nage), 1), matrix(0, 1, 1))
, global = T, null = T, point.est = mean, verb = FALSE, s.tol = 10^(-10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Harv.data}{
Harvest data, rows to be age classes and columns for years.
}
  \item{Aerial.data}{
Aerial count data, should be only one row and same number of columns with harvest.
}
  \item{nage}{
Number of age classes, female first.
}
  \item{mean.f,mean.s}{
Best estimation of fecundity/survival. Number of rows and columns should be less than harvest, after multiply by Assumptions matrices (see \code{Assumptions}), a full fecundity matrix (i.e. assign fecundity for each age class in each projection year) should have one column less and same number of rows with harvest data(one less year interval). Default input (i.e. no specific Assumptions matrix) should be a full fecundity matrix by itself, i.e.:

\code{nrow(mean.f) == nrow(Harv.data) & ncol(mean.f)==ncol(Harv.data)-1}.
\code{nrow(mean.s) == nrow(Harv.data) & ncol(mean.s)==ncol(Harv.data)-1}.
}

  \item{mean.SRB}{
Similar to mean.f, sex ratio at birth, but note SRB is assumed having no age structure. Thus a full SRB matrix should only have one row.
}
  \item{mean.H}{
Similar to fecundity, for harvest rate. Difference is since we assume a pre-repoduction harvest, a full harvest rate matrix should have same size of harvest data.
}
  \item{mean.A}{
Similar to H, Aerial count detection rate. By default, a full aerial count matrix should have same size as aerial count data.
}
  \item{mean.b}{
Best estimation of baseline harvest, usually just firest year's harvest count.
}
  \item{n.iter}{
Number of posterior sample saved.
}
  \item{burn.in}{
Burn in iteration.
}
  \item{thin.by}{
Sample thined by?
}
  \item{al.f,al.s,al.SRB,al.H,al.A}{
Hyperparameter for prior distribution of fecundity/survival/SRB/Harvest rate/Aerial Detection Rate's variance. Variance of vital rates v beside aK0 follow inverse Gamma distribution with parameter al.v and be.v, v can be fecundity, survival, harvest rate, aerial count detection rate and SRB.
}
  \item{be.f,be.s,be.SRB,be.H,be.A}{
Hyperparameter for prior distribution of fecundity/survival/SRB/Harvest rate/Aerial Detection Rate's variance. Variance of vital rates v beside aK0 follow inverse Gamma distribution with parameter al.v and be.v, v can be fecundity, survival, harvest rate, aerial count detection rate and SRB.
}

  \item{min.aK0,max.aK0}{
Prior distribution of aK0s. If fit density dependent Leslie matrix model:

\code{fecundity(Population) = (1 - aK0 [[1]] * (Population - aK0[[3]])) * Fec}

\code{survival(Population) = (1 - aK0[[2]] * (Population - aK0[[3]])) * Surv}
In which fecundity is a function of population and Fec and Surv is fecundity and survival at population aK0[[3]]

aK0's component were assumed to be uniformly distributed with parameter min.aK0 and max.aK0.
}

  \item{Assumptions}{
Should be a list, each element is a list and named Fec, Surv, SRB, Harv, AerialDet, aK0 for fecundity, Survival, SRB, Harvest rate, Aerial count detection rate, and linear coeffient of Density Dependency.

In each list, there should be two matrices, naming \code{age} and \code{time}. In sampler function, a full vital rate matrix (i.e. assign vital rates for each age class in each year needed, for repoduction related, will be one year less than harvest and aerial count related viral rates) used in projection will be produced by matrix multiplication, e.g.

\code{Fec_full <- Assumptions$Fec$age \%*\% mean.f \%*\% Assumptions$Fec$time}

By default, all Assumptions are identity matrices, i.e. prior mean vital rates matrices should be a full one by itself.
}
  \item{start.sigmasq.f}{
Starting value of fecundity's variance.
}
  \item{start.sigmasq.s}{
Starting value of survival's variance.
}
  \item{start.sigmasq.SRB}{
Starting value of SRB's variance.
}
  \item{start.sigmasq.H}{
Starting value of harvest rate's variance.
}
  \item{start.sigmasq.A}{
Starting value of aerial count detection rate's variance.
}
  \item{prop.vars}{
Proposal variance, should be a list named fert.rate, surv.prop, SRB, A,H, aK0, and baseline.pop.count, same dimension with mean.f, mean.s, mean.SRB, mean.A, mean.H, min.aK0, mean.b.
}
  \item{proj.periods}{
How many years need to be projected? Usually ncol(Harv.data)-1
}
  \item{estFec}{
Estimate fecundity? Bool. If not, mean.f will be used.
}
  \item{estaK0}{
Estimate density dependency parameters? Bool. Note this will be masked if set null=TRUE.
}
  \item{aK0}{
If not estimate aK0, what aK0 is.
}
  \item{global}{
Whether aK0[[3]] is global or age specific (i.e. will there be age specific carrying capcacity). Bool.
}
  \item{null}{
Null model for density dependency? Bool. If null=TRUE, no density dependency will be considered, but model can fit time inhomogeneous vital rates.
}
  \item{point.est}{
How to do point estimation? Can be mean or median.
}
  \item{verb}{
Want a verbose run? Bool
}
  \item{s.tol}{
Tolerant for lowest survival rate (below will be considered as 0).
}
}
\details{
Do not use Assumptions unless you understand it.
}
\value{
%%  ~Describe the value returned
A list with length 7:
\item{mcmc.objs }{A list of mcmc objects returned by the sampling algorithm}
\item{log.like.mcmc }{mcmc object for log likelihood }
\item{alg.stats}{algorithm statistics, e.g. acceptance rate}
\item{model.checking}{Model checking, including DIC, absolute difference calculated at point estimation and standard deviation of fitted values}
\item{fixed.params}{fixed parameters}
\item{start.vals}{starting value of the algorithm}
\item{alg.params}{parameters for algorithm setting}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Yunyi Shen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##Use our Chicago Deer data.
set.seed(42)
ReCAP_sampler(ChicagoDeerdata$Harv.data,ChicagoDeerdata$Aeri.data
    ,ChicagoDeerdata$nage,ChicagoDeerdata$mean.f,ChicagoDeerdata$mean.s
    ,ChicagoDeerdata$mean.SRB,ChicagoDeerdata$mean.H,ChicagoDeerdata$mean.A
    ,n.iter = 50,burn.in = 50,thin.by = 1
    ,Assumptions = ChicagoDeerdata$Assumptions)


## another example
### This makes a full prior harvest matrix and the model will run a fully age
### and time specific reconstruction of harvest rate together with other vital rates.


mean.H.full <- ChicagoDeerdata$Assumptions$Harv$age \%*\%
    as.matrix( ChicagoDeerdata$mean.H) \%*\%
    ChicagoDeerdata$Assumptions$Harv$time
set.seed(42)
ReCAP_sampler(ChicagoDeerdata$Harv.data,ChicagoDeerdata$Aeri.data
    ,ChicagoDeerdata$nage,ChicagoDeerdata$mean.f,ChicagoDeerdata$mean.s
    ,ChicagoDeerdata$mean.SRB,mean.H.full,ChicagoDeerdata$mean.A
    ,n.iter = 50,burn.in = 50,thin.by = 1
    )
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
