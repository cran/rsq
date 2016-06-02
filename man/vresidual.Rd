\name{vresidual}
\alias{vresidual}
\alias{vresidual.nb}
\title{
Variance-Function-Based Residual
}
\description{
Calculate the residual which is defined on the basis of the variance function.}
\usage{
vresidual(y,yfit,family="binomial")
vresidual.nb(y,yfit,theta=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ the response values;}
  \item{yfit}{ the fitted values;}
  \item{theta}{ \eqn{\theta} of the negative binomial distribution;}
  \item{family}{ a character string naming the family of the generalized linear model. (See \link{family} for details.)}
}
\details{
They are internal functions called by \link{vsse}.
}
\value{The variance-function-based residual.}
\references{
Zhang, D. (2016). A coefficient of determination for generalized linear models. The American Statistician, tentatively accepted.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
ypred <- bnfit$fitted.values

resid<-vresidual(y,ypred,family="binomial")
summary(resid)
head(resid)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
