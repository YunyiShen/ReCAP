\name{ProjectHarvest}
\alias{ProjectHarvest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Project Harvest Matrix Using Vital Rates
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ProjectHarvest(Surv, Harvpar, Fec, SRB, bl, period, nage, aK0 = list(matrix(0, nage[1], 1), matrix(0, sum(nage), 1), matrix(0, 1, 1)), global = T, null = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Surv}{
%%     ~~Describe \code{Surv} here~~
}
  \item{Harvpar}{
%%     ~~Describe \code{Harvpar} here~~
}
  \item{Fec}{
%%     ~~Describe \code{Fec} here~~
}
  \item{SRB}{
%%     ~~Describe \code{SRB} here~~
}
  \item{bl}{
%%     ~~Describe \code{bl} here~~
}
  \item{period}{
%%     ~~Describe \code{period} here~~
}
  \item{nage}{
%%     ~~Describe \code{nage} here~~
}
  \item{aK0}{
%%     ~~Describe \code{aK0} here~~
}
  \item{global}{
%%     ~~Describe \code{global} here~~
}
  \item{null}{
%%     ~~Describe \code{null} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Surv, Harvpar, Fec, SRB, bl, period, nage, aK0 = list(matrix(0,
    nage[1], 1), matrix(0, sum(nage), 1), matrix(0, 1, 1)), global = T,
    null = T)
{
    ProjectHarvestCpp(Surv, Harvpar, Fec, SRB, aK0, global, null,
        bl, period, nage)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

