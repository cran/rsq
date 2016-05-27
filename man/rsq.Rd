\name{rsq}
\alias{rsq}
\title{Coefficient of Determination}
\description{Calculate the coefficient of determination, aka R^2, for both linear regression models, and generalized linear regression models.}
\usage{rsq(fitObj,adj=FALSE)}
\arguments{
    \item{fitObj}{ an object of class "lm" or "glm", usually, a result of a call to \link{lm}, \link{glm}, or \link{glm.nb}.}
    \item{adj}{ logical; if TRUE, calculate the adjusted R^2.}
}
\details{
The R^2 relies on the variance function, and reduces to the traditional R^2 when the variance function is constant or linear.
}
\value{The R^2 or adjusted R^2.}
\references{
Zhang, D. (2006). A coefficient of determination for generalized linear models. \emph{The American Statistician}, tentatively accepted.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq.partial}}.
}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
rsq(bnfit)
rsq(bnfit,adj=TRUE)

quasibn <- glm(y~color+spine+width+weight,family=quasibinomial)
rsq(quasibn)
rsq(quasibn,adj=TRUE)

psfit <- glm(num.satellites~color+spine+width+weight,family=poisson)
rsq(psfit)
rsq(psfit,adj=TRUE)

quasips <- glm(num.satellites~color+spine+width+weight,family=quasipoisson)
rsq(quasips)
rsq(quasips,adj=TRUE)

# Effectiveness of Bycycle Safety Helmets in Thompson et al. (1989)
y <- matrix(c(17,218,233,758),2,2)
x <- factor(c("yes","no"))
tbn <- glm(y~x,family=binomial)
rsq(tbn)
rsq(tbn,adj=TRUE)
}
