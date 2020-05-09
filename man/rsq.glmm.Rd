\name{rsq.glmm}
\alias{rsq.glmm}
\title{R-Squared for Generalized Linear Mixed Models}
\description{Calculate the variance-function-based R-squared for generalized linear mixed models.}
\usage{rsq.glmm(fitObj,adj=FALSE,DFUN=NULL)}
\arguments{
    \item{fitObj}{ an object of class "glmerMod", usually, a result of a call to \link{glmer} or \link{glmer.nb} in \link{lme4}.}
    \item{adj}{ logical; if TRUE, calculate the adjusted R^2.}
    \item{DFUN}{  function which calculates the directed length along the variance function, and should be specified for quasi models (with \link{family}=quasi()).}
}
\details{
There are three types of R^2 calculated on the basis of observed response values, estimates of fixed effects, and variance components, i.e., model-based R_M^2 (proportion of variation explained by the model in total, including both fixed-effects and random-efffects factors), fixed-effects R_F^2 (proportion of variation explained by the fixed-effects factors), and random-effects R_R^2 (proportion of variation explained by the random-effects factors).
}
\value{
  \item{R_M^2}{proportion of variation explained by the model in total, including both fixed-effects and random-efffects factors.}
  \item{R_F^2}{proportion of variation explained by the fixed-effects factors.}
  \item{R_R^2}{proportion of variation explained by the random-effects factors.}
}
\references{
Zhang, D. (2017). A coefficient of determination for generalized linear models. \emph{The American Statistician}, 71(4): 310-316.

Zhang, D. (2020). Coefficients of determination for generalized linear mixed models. \emph{Technical Report}, 20-01, Department of Statistics, Purdue University.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{vresidual}, \link{rsq}, \link{rsq.v}}.
}
\examples{
require(lme4)
data(cbpp)
glmm1 <- glmer(cbind(incidence,size-incidence)~period+(1|herd),data=cbpp,family=binomial)
rsq.glmm(glmm1)
rsq(glmm1)
}
