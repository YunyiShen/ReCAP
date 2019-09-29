\name{Lambda_plot}
\alias{Lambda_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Annual Growth Rate at Each Year under Different Scanerio
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Lambda_plot(ReCAP_lambda_obj,start_year=1,alpha=.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ReCAP_lambda_obj}{
An object returned by \code{analysisLambda}.
}

\item{start_year}{
The start year.
}

\item{alpha}{
Alpha level for CI.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A ggplot object will be returned, with annual growth rates of 4 scanerios (uniform age structure, stable age structure, skip harvest and observed with harvest).
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
## see
?analysisLambda
## for example
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Annual growth rate}% use one of  RShowDoc("KEYWORDS")

