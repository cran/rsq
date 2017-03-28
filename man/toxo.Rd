\name{toxo}
\alias{toxo}
\docType{data}
\title{
Toxoplasmosis Test in El Salvador
}
\description{
Recorded are the numbers of subjects testing positive for toxoplasmosis in 34 cities of El Salvador. 
}
\usage{data("toxo")}
\format{
  A data frame with the test results in 34 cities of El Salvador, includingthe following 4 variables.
  \describe{
    \item{\code{city}}{ index of each city.}
    \item{\code{positive}}{ the number of subjects testing positive for toxoplasmosis.}
    \item{\code{nsubs}}{ the total number of subjects tested.}
    \item{\code{rainfall}}{ annual rainfall (mm) in home city of subject.}
  }
}
\details{
All subjects are between 11 and 15 year old. The data set was abstracted from a larger data set in Rmington et al. (1970).
}
\source{
Efron, B. (1978). Regression and ANOVA with zero-one data: measures of residual variation. \emph{JASA}, 73: 113-121.
}
\references{
Remington, J.S., Efron, B., Cavanaugh, E., Simon,
H.J., and Trejos, A. (1970). Studies on toxoplasmosis
in El Salvador, prevalence and incidence of toxoplasmosis as
measured by the Sabin-Feldman Dye test. \emph{Transactions of the
Royal Society of Tropical Medicine and Hygiene}, 64: 252-267.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq}, \link{rsq.partial}, \link{pcor}, \link{simglm}}.
}
\examples{
data(toxo)
summary(toxo)
attach(toxo)

toxofit<-glm(cbind(positive,nsubs-positive)~rainfall+I(rainfall^2)+I(rainfall^3),family=binomial)

rsq(toxofit)
rsq(toxofit,adj=TRUE)
rsq.partial(toxofit)

detach(toxo)
}
\keyword{datasets}
