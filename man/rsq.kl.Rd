\name{rsq.kl}
\alias{rsq.kl}
\title{KL-Divergence-Based R-Squared}
\description{The Kullback-Leibler-divergence-based R^2 for generalized linear models.}
\usage{rsq.kl(fitObj,adj=FALSE)}
\arguments{
    \item{fitObj}{ an object of class "lm" or "glm", usually, a result of a call to \link{lm}, \link{glm}, or \link{glm.nb}.}
    \item{adj}{ logical; if TRUE, calculate the adjusted R^2.}
}
\details{
This version of R^2 was proposed by Cameron and Windmeijer (1997). It is extended to quasi models (Zhang, 2017) based on the quasi-likelihood function (McCullagh, 1983).}
\value{The R^2 or adjusted R^2.}
\references{
Cameron, A. C. and Windmeijer, A. G. (1997) An R-squared measure of goodness of fit for some common nonlinear regression models. \emph{Journal of Econometrics}, 77: 329-342.

McCullagh, P. (1983) Quasi-likelihood functions. \emph{Annals of Statistics}, 11: 59-67.

Zhang, D. (2017) rsq: An R package to calculate the R-squared for generalized linear models. In preparation.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq}, \link{rsq.partial}, \link{pcor}}.
}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
rsq.kl(bnfit)
rsq.kl(bnfit,adj=TRUE)

psfit <- glm(num.satellites~color+spine+width+weight,family=poisson)
rsq.kl(psfit)
rsq.kl(psfit,adj=TRUE)

# Effectiveness of Bycycle Safety Helmets in Thompson et al. (1989)
y <- matrix(c(17,218,233,758),2,2)
x <- factor(c("yes","no"))
tbn <- glm(y~x,family=binomial)
rsq.kl(tbn)
rsq.kl(tbn,adj=TRUE)
}
