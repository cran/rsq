\name{lifetime}
\alias{lifetime}
\docType{data}
\title{
Lifetimes in Two Different Environments.
}
\description{
There are 27 tests in each of the two environments.
}
\usage{data("lifetime")}
\format{
  A data frame with 54 observations on the following 2 variables.
  \describe{
    \item{\code{time}}{ the lifetime (x10).}
    \item{\code{env}}{ the environment of each test (kg/mm^2).}
  }
}
\details{
This data set is discussed by Wang et al. (1992).
}
\source{
Wang, H., Ma, B., and Shi, J. (1992). Estimation of environmental factors for the inverse gaussian distribution. \emph{Microelectron. Reliab.}, 32: 931-934.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq}, \link{rsq.partial}, \link{pcor}, \link{simglm}}.
}\examples{
data(lifetime)
summary(lifetime)
head(lifetime)

attach(lifetime)
igfit <- glm(time~env,family=inverse.gaussian)
rsq(igfit)
rsq(igfit,adj=TRUE)
}
\keyword{datasets}
