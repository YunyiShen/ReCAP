\name{getAerialCount}
\alias{getAerialCount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Expacted post-Harvest Aerial Count from Harvest Count, Harvest Rate and Detection rate
}
\description{

}
\usage{
getAerialCount(Harv, H, A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Harv}{
  Harvest count, a matrix.
}
  \item{H}{
  Harvest rate, a matrix with same dimension with \code{Harv}.
}
  \item{A}{
  Aerial Detection rate, a 1 by p matrix with p=period=ncol(Harv).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A 1 by p matrix with p=period=ncol(Harv), contains expacted post-harvest aerial count for each harvest year.
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
Harv <- matrix(rpois(15,10),3,5) # 3 age and 5 years
H <- matrix(runif(15),3,5)
A <- matrix(runif(5),1,5) # 5 years
getAerialCount(Harv,H,A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{AerialCount }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
