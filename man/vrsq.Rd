\name{vrsq}
\alias{vrsq}
\title{
Variance-Function-Based Coefficient of Determination
}
\description{
Calculate the coefficient of partial determination, aka partial R^2, which is defined on the basis of variance function.}
\usage{
vrsq(y,yfit0,yfit1,family="binomial",adj=FALSE,df.int=1,p=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ the response values.}
  \item{yfit0}{ the fitted values based on the null model.}
  \item{yfit1}{ the fitted values based on the model.}
  \item{family}{ a character string naming the family of the generalized linear model. (See \link{family} for details.)}
  \item{adj}{ logical; if TRUE, calculate the adjusted partial R^2.}
  \item{df.int}{ = 1, if the intercept is included in the model; = 0, otherwise.}
  \item{p}{ number of predictors, including the intercept if it is modeled.}
}
\details{
This is an internal function called by \link{rsq} and \link{rsq.partial}. Note that, when calculating the adjusted R^2, it is necessary to specify p as the total number of regression coefficients in the model.
}
\value{The R^2 or adjusted R^2.}
\references{
Zhang, D. (2006). A coefficient of determination for generalized linear models. The American Statistician, tentatively accepted.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
ypred <- bnfit$fitted.values
bnfit0 <- update(bnfit,.~1)
ypred0 <- bnfit0$fitted.values
 
vrsq(y,ypred0,ypred,family="binomial")
p <- bnfit$rank
vrsq(y,ypred0,ypred,family="binomial",adj=TRUE,p=p)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
