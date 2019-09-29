\name{ChicagoDeerdata}
\alias{ChicagoDeerdata}
\docType{data}
\title{
Population data of Chicago deer under intensive culling
}
\description{
Population data of Chicago deer under intensive culling, including culling count, aerial count, best estimation of fecundity, survival, etc.
}
\usage{data("ChicagoDeerdata")}
\format{
  The format is:

List of 10

 $ Harv.data   :'data.frame':	11 obs. of  15 variables

 $ Aeri.data   :'data.frame':	1 obs. of  15 variables

 $ Assumptions :List of 1

 $ mean.f      :'data.frame':	8 obs. of  14 variables

 $ mean.s      :'data.frame':	11 obs. of  14 variables

 $ mean.H      :'data.frame':	4 obs. of  15 variables

 $ mean.SRB    :'data.frame':	1 obs. of  14 variables

 $ mean.A      : num [1, 1:15] 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 ...

 $ nage        : num [1:2, 1] 8 3

 $ proj.periods: num 14
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Etter, Dwayne Robert. Ecology and management of overabundant white-tailed deer from suburban Chicago, Illinois. Diss. University of Illinois at Urbana-Champaign, 2001.
}
\examples{
plot(1992:2006,ChicagoDeerdata$Harv.data[1,]) # fawn harvest
}
\keyword{datasets}
