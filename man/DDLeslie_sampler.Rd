\name{DDLeslie_sampler}
\alias{DDLeslie_sampler}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrated Population Model for Harvesting and Aerial Counting data under Baysian Framework
}
\description{
Main sampler of the posterior distribution of vital rates and harvest count in \code{ReDDLeslie}.
}
\usage{
DDLeslie_sampler(Harv.data, Aerial.data
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
%%     ~~Describe \code{Harv.data} here~~
}
  \item{Aerial.data}{
%%     ~~Describe \code{Aerial.data} here~~
}
  \item{nage}{
%%     ~~Describe \code{nage} here~~
}
  \item{mean.f}{
%%     ~~Describe \code{mean.f} here~~
}
  \item{mean.s}{
%%     ~~Describe \code{mean.s} here~~
}
  \item{mean.SRB}{
%%     ~~Describe \code{mean.SRB} here~~
}
  \item{mean.H}{
%%     ~~Describe \code{mean.H} here~~
}
  \item{mean.A}{
%%     ~~Describe \code{mean.A} here~~
}
  \item{mean.b}{
%%     ~~Describe \code{mean.b} here~~
}
  \item{n.iter}{
%%     ~~Describe \code{n.iter} here~~
}
  \item{burn.in}{
%%     ~~Describe \code{burn.in} here~~
}
  \item{thin.by}{
%%     ~~Describe \code{thin.by} here~~
}
  \item{al.f}{
%%     ~~Describe \code{al.f} here~~
}
  \item{be.f}{
%%     ~~Describe \code{be.f} here~~
}
  \item{al.s}{
%%     ~~Describe \code{al.s} here~~
}
  \item{be.s}{
%%     ~~Describe \code{be.s} here~~
}
  \item{al.SRB}{
%%     ~~Describe \code{al.SRB} here~~
}
  \item{be.SRB}{
%%     ~~Describe \code{be.SRB} here~~
}
  \item{min.aK0}{
%%     ~~Describe \code{min.aK0} here~~
}
  \item{max.aK0}{
%%     ~~Describe \code{max.aK0} here~~
}
  \item{al.H}{
%%     ~~Describe \code{al.H} here~~
}
  \item{be.H}{
%%     ~~Describe \code{be.H} here~~
}
  \item{al.A}{
%%     ~~Describe \code{al.A} here~~
}
  \item{be.A}{
%%     ~~Describe \code{be.A} here~~
}
  \item{Assumptions}{
%%     ~~Describe \code{Assumptions} here~~
}
  \item{start.sigmasq.f}{
%%     ~~Describe \code{start.sigmasq.f} here~~
}
  \item{start.sigmasq.s}{
%%     ~~Describe \code{start.sigmasq.s} here~~
}
  \item{start.sigmasq.SRB}{
%%     ~~Describe \code{start.sigmasq.SRB} here~~
}
  \item{start.sigmasq.H}{
%%     ~~Describe \code{start.sigmasq.H} here~~
}
  \item{start.sigmasq.A}{
%%     ~~Describe \code{start.sigmasq.A} here~~
}
  \item{prop.vars}{
%%     ~~Describe \code{prop.vars} here~~
}
  \item{proj.periods}{
%%     ~~Describe \code{proj.periods} here~~
}
  \item{estFec}{
%%     ~~Describe \code{estFer} here~~
}
  \item{estaK0}{
%%     ~~Describe \code{estaK0} here~~
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
  \item{point.est}{
%%     ~~Describe \code{point.est} here~~
}
  \item{verb}{
%%     ~~Describe \code{verb} here~~
}
  \item{s.tol}{
%%     ~~Describe \code{s.tol} here~~
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
