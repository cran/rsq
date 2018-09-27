\name{rsq.v}
\alias{rsq.v}
\title{Variance-Function-Based R-Squared}
\description{Calculate the variance-function-based R-squared for generalized linear models.}
\usage{rsq.v(fitObj,adj=FALSE,data=NULL)}
\arguments{
    \item{fitObj}{ an object of class "lm" or "glm", usually, a result of a call to \link{lm}, \link{glm}, or \link{glm.nb}.}
    \item{adj}{ logical; if TRUE, calculate the adjusted R^2.}
    \item{data}{an optional data frame used by fitObj.}
}
\details{
The R^2 relies on the variance function, and is well-defined for quasi models. It reduces to the classical R^2 when the variance function is constant or linear.
}
\value{The R^2 or adjusted R^2.}
\references{
Zhang, D. (2016). A coefficient of determination for generalized linear models. \emph{The American Statistician}, accepted. <http://dx.doi.org/10.1080/00031305.2016.1256839>.
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
rsq.v(bnfit)
rsq.v(bnfit,adj=TRUE)

quasibn <- glm(y~color+spine+width+weight,family=quasibinomial)
rsq.v(quasibn)
rsq.v(quasibn,adj=TRUE)

psfit <- glm(num.satellites~color+spine+width+weight,family=poisson)
rsq.v(psfit)
rsq.v(psfit,adj=TRUE)

quasips <- glm(num.satellites~color+spine+width+weight,family=quasipoisson)
rsq.v(quasips)
rsq.v(quasips,adj=TRUE)

# Effectiveness of Bycycle Safety Helmets in Thompson et al. (1989)
y <- matrix(c(17,218,233,758),2,2)
x <- factor(c("yes","no"))
tbn <- glm(y~x,family=binomial)
rsq.v(tbn)
rsq.v(tbn,adj=TRUE)
}
