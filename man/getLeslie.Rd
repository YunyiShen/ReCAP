\name{getLeslie}
\alias{getLeslie}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Leslie Matrix from Fecundity Survival and Sex Ratio at Birth (SRB)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getLeslie(Surv,Fec,SRB)
}
%- maybe also 'usage' for other objects documented here.


\arguments{
  \item{Surv}{
%%     ~~Describe \code{x} here~~
Survival rate, should be a column matrix.
}
  \item{Fec}{

%%     ~~Describe \code{x} here~~
Fecundity, should be a column matrix.}
  \item{SRB}{
%%     ~~Describe \code{x} here~~
Sex Ratio at Birth, should be a 1x1 matrix.}

}


\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
A Leslie matrix
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
##---- NOT RUN ----##
Surv <- matrix(runif(6),6,1) # 6 age classes in total
Fec <- matrix(2*runif(3),3,1) # 3 female age classes
SRB <- matrix(runif(1),1,1)
getLeslie(Surv,Fec,SRB)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Leslie }% use one of  RShowDoc("KEYWORDS")
\keyword{ DDLeslie }% __ONLY ONE__ keyword per line
