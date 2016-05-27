\name{hcrabs}
\alias{hcrabs}
\docType{data}
\title{
Satellites of Female Horseshoe Crabs
}
\description{
Recorded are the numbers of male satellites, and other characteristics of 173 female horseshoe crabs. 
}
\usage{data("hcrabs")}
\format{
  A data frame with 173 observations on the following 5 variables.
  \describe{
    \item{\code{color}}{ the female crab's color, coded 1: light; 2: medium light; 3: medium; 4: medium dark; 5: dark. Not all of these colors appear.}
    \item{\code{spine}}{ the female crab's spine condition, coded 1: both good; 2: one worn or broken; 3: both worn or broker.}
    \item{\code{width}}{ the female crab's carapace width (cm).}
    \item{\code{num.satellites}}{ the number of satellite males.}
    \item{\code{weight}}{ the female crab's weight (kg).}
  }
}
\details{
A nesting female horseshoe crab may have male crabs residing nearby, called satellites, besides the male crab residing in her nest. Brockmann (1996) investigated factors (including the female crab's color, spine condition, weight, and carapace width) which may influence the presence/obsence of satellite males. This data set has been discussed by Agresti (2002).
}
\source{
Agresti, A. (2012). \emph{An Introduction to Categorical Data Analysis}, 3rd edition. Wiley: New Jersey.
}
\references{
Brockmann, H. J. (1996). Satellite male groups in horseshoe crabs. \emph{Limulus polyphemus. Ethology}, 102: 1-21.
}
\examples{
data(hcrabs)
summary(hcrabs)
head(hcrabs)

attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
rsq(bnfit)
rsq(bnfit,adj=TRUE)
rsq.partial(bnfit)
rsq.partial(bnfit,adj=TRUE)

quasips <- glm(num.satellites~color+spine+width+weight,family=quasipoisson)
rsq(quasips)
rsq(quasips,adj=TRUE)
rsq.partial(quasips)
rsq.partial(quasips,adj=TRUE)
}
\keyword{datasets}
