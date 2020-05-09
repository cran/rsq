\name{rsq.lmm}
\alias{rsq.lmm}
\title{R-Squared for Linear Mixed Models}
\description{Calculate the R-squared for linear mixed models.}
\usage{rsq.lmm(fitObj,adj=FALSE)}
\arguments{
    \item{fitObj}{ an object of class "merMod" or "lmerMod" or "lme", usually, a result of a call to \link{lmer} in \link{lme4}, or \link{lme} in \link{nlme}.}
    \item{adj}{ logical; if TRUE, calculate the adjusted R^2.}
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
Zhang, D. (2020). Coefficients of determination for generalized linear mixed models. \emph{Technical Report}, 20-01, Department of Statistics, Purdue University.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq}, \link{rsq.v}}.
}
\examples{
# lmer in lme4
require(lme4)
lmm1 <- lmer(Reaction~Days+(Days|Subject),data=sleepstudy)
rsq(lmm1)
rsq.lmm(lmm1)

# lme in nlme
require(nlme)
lmm2 <- lme(Reaction~Days,data=sleepstudy,random=~Days|Subject)
rsq(lmm2)
rsq.lmm(lmm2)
}
