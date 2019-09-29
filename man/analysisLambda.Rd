\name{analysisLambda}
\alias{analysisLambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reconstruct Posterior Annual Growth Rate in Different Scenario
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
analysisLambda(mcmc_obj, Assumptions = list(), nage, n_proj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc_obj}{
  A list, mcmc.obj from ReCAP_sampler result
}
  \item{Assumptions}{
Assumption matrices used in ReCAP_sampler
}
  \item{nage}{
Number of age classes, female first, same as ReCAP_sampler
}
  \item{n_proj}{
Number of years did projection, usually ncol(Harvest_data)-1
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of class "ReCAP_lambda" has same length as number of MCMC samples taken in ReCAP_sampler. Each element is one posterior sample which is a matrix, with each column as years. Rows are anaual growth rate of: maximum possible, uniform age structure, stable age structure, observed age structure w/o harvest, minimum possible, observed growth rate (w/ harvest).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
set.seed(42)
ReCAP_sample = ReCAP_sampler(ChicagoDeerdata$Harv.data
    ,ChicagoDeerdata$Aeri.data
    ,ChicagoDeerdata$nage
    ,ChicagoDeerdata$mean.f
    ,ChicagoDeerdata$mean.s
    ,ChicagoDeerdata$mean.SRB
    ,ChicagoDeerdata$mean.H
    ,ChicagoDeerdata$mean.A
    ,n.iter = 50,burn.in = 5,thin.by = 1
    ,Assumptions = ChicagoDeerdata$Assumptions)

analysisoflambda = analysisLambda(ReCAP_sample$mcmc.objs
    ,ChicagoDeerdata$Assumptions
    ,ChicagoDeerdata$nage
    ,ChicagoDeerdata$proj.periods)

Lambda_plot(analysisoflambda,1992,0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Annual Growth Rate }% use one of  RShowDoc("KEYWORDS")
