\name{rsq.v}
\alias{rsq.v}
\title{Variance-Function-Based R-Squared}
\description{Calculate the variance-function-based R-squared for generalized linear (mixed) models.}
\usage{rsq.v(fitObj,adj=FALSE,DFUN=NULL)}
\arguments{
    \item{fitObj}{ an object of class "lm", "glm", "lme", or "glmerMod", usually, a result of a call to \link{lm}, \link{glm}, \link{glm.nb}, \link{glmer}, or \link{glmer.nb}.}
    \item{adj}{ logical; if TRUE, calculate the adjusted R^2.}
    \item{DFUN}{ function which calculates the directed length along the variance function, and should be specified for quasi models (with \link{family}=quasi()).}
}
\details{
The R^2 relies on the variance function, and is well-defined for quasi models. It reduces to the classical R^2 when the variance function is constant or linear. For quasi models specified with family=quasi(), DFUN must be specified when the variance function is not constant or linear. For (generalized) linear mixed models, there are three types of R^2 calculated on the basis of observed response values, estimates of fixed effects, and variance components, i.e., model-based R_M^2 (proportion of variation explained by the model in total, including both fixed-effects and random-efffects factors), fixed-effects R_F^2 (proportion of variation explained by the fixed-effects factors), and random-effects R_R^2 (proportion of variation explained by the random-effects factors).
}
\value{The R^2 or adjusted R^2. For (generalized) linear mixed models,
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
\code{\link{vresidual}, \link{rsq}, \link{rsq.glmm}, \link{rsq.partial}, \link{pcor}}.
}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
rsq.v(bnfit)
rsq.v(bnfit,adj=TRUE)

quasibn <- glm(y~color+spine+width+weight,family=quasibinomial)
rsq.v(quasibn)
rsq.v(quasibn,adj=TRUE)

# Generalized linear mixed models
require(lme4)
data(cbpp)
glmm1 <- glmer(cbind(incidence,size-incidence)~period+(1|herd),data=cbpp,family=binomial)
rsq.v(glmm1)
}
