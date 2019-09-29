\name{ProjectHarvest}
\alias{ProjectHarvest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Project Harvest Using Vital Rates
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ProjectHarvest(Surv, Harvpar, Fec, SRB
        , bl, period, nage
        , aK0 = list(matrix(0, nage[1], 1)
                , matrix(0, sum(nage), 1)
                , matrix(0, 1, 1))
        , global = T, null = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Surv,Harvpar,Fec,SRB,aK0}{
Full sruvival, Harvest, Fecundity, SRB and aK0 matrices.
}
  \item{bl}{
Baseline harvest to start with
}
  \item{period}{
Number of year to project
}
  \item{nage}{
Number of age classes, female first.
}

  \item{global}{
Bool. Whether DD is global
}
  \item{null}{
Bool. DD or not
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
mean.H.full <- ChicagoDeerdata$Assumptions$Harv$age \%*\%
    as.matrix( ChicagoDeerdata$mean.H) \%*\%
    ChicagoDeerdata$Assumptions$Harv$time

ProjectHarvest(as.matrix(ChicagoDeerdata$mean.s)
    ,mean.H.full,as.matrix(ChicagoDeerdata$mean.f)
    ,as.matrix(ChicagoDeerdata$mean.SRB)
    ,matrix(ChicagoDeerdata$Harv.data[,1])
    ,ChicagoDeerdata$proj.periods
    ,ChicagoDeerdata$nage)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

