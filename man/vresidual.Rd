\name{vresidual}
\alias{vresidual}
\title{Variance-Function-Based Residuals}
\description{Calculate the variance-function-based residuals for generalized linear models, which are used to calculate the variance-function-based R-squared.}
\usage{vresidual(y,yfit,family="binomial",DFUN=NULL,theta=1)}
\arguments{
    \item{y}{ a vector of observed values.}
    \item{yfit}{ a vector of fitted values.}
    \item{family}{ family of the distribution.}
    \item{DFUN}{ function which calculates the directed distance along the variance function (specified by family by default).}
    \item{theta}{ extra parameter for ceterain family like negative binomial.}
}
\details{
The calcualted residual relies on the variance function, and is well-defined for quasi models. It reduces to the classical residual when the variance function is constant or linear. Note that only one of DFUN and family (when the family is not "quasi") is required to specify. When the family is "quasi" and the variance is not constant, DFUN should be specified.
}
\value{Variance-function-based residuals.}
\references{
Zhang, D. (2017). A coefficient of determination for generalized linear models. \emph{The American Statistician}, 71(4): 310-316.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq.v}, \link{rsq}}.
}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family="binomial")
vresidual(y,bnfit$fitted.values,family="binomial")

# Effectiveness of Bycycle Safety Helmets in Thompson et al. (1989)
y <- matrix(c(17,218,233,758),2,2)
x <- factor(c("yes","no"))
tbn <- glm(y~x,family="binomial")
vresidual(y,tbn$fitted.values,family="binomial")
}
