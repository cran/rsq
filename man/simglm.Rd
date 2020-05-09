\name{simglm}
\alias{simglm}
\title{Simulate Data from Generalized Linear Models}
\description{Simulate data from linear and generalized linear models. Only the first covariate truely affects the response variable with coefficient equal to \code{lambda}.}
\usage{simglm(family=c("binomial", "gaussian", "poisson","Gamma"),lambda=3,n=50,p=3)}
\arguments{
    \item{family}{ the family of the distribution.}
    \item{lambda}{ size of the coefficient of the first covariate.}
    \item{n}{ the sample size.}
    \item{p}{ the number of covarites.}
}
\details{
The first covariate takes 1 in half of the observations, and 0 or -1 in the other half. When \code{lambda} gets larger, it is supposed to easier to predict the response variable.
}
\value{
Returned values include \code{yx} and \code{beta}.
\item{yx}{ a data frame including the response \code{y} and covariates \code{x.1}, \code{x.2}, and so on.}
\item{beta}{ true values of the regression coefficients.}
}
\references{
Zhang, D. (2017). A coefficient of determination for generalized linear models. \emph{The American Statistician}, 71(4): 310-316.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq}, \link{rsq.partial}, \link{pcor}}.
}
\examples{
# Poisson Models
sdata <- simglm(family="poisson",lambda=4)
fitf <- glm(y~x.1+x.2+x.3,family=poisson,data=sdata$yx)
rsq(fitf)  # type='v'

fitr <- glm(y~x.2+x.3,family=poisson,data=sdata$yx)
rsq(fitr)  # type='v'
rsq(fitr,type='kl')
rsq(fitr,type='lr')
rsq(fitr,type='n')

pcor(fitr)  # type='v'
pcor(fitr,type='kl')
pcor(fitr,type='lr')
pcor(fitr,type='n')

# Gamma models with shape=100
n <- 50
sdata <- simglm(family="Gamma",lambda=4,n=n)
fitf <- glm(y~x.1+x.2+x.3,family=Gamma,data=sdata$yx)
rsq(fitf)  # type='v'
rsq.partial(fitf)  # type='v'

fitr <- glm(y~x.2,family=Gamma,data=sdata$yx)
rsq(fitr)  # type='v'
rsq(fitr,type='kl')
rsq(fitr,type='lr')
rsq(fitr,type='n')

# Likelihood-ratio-based R-squared
y <- sdata$yx$y
yhatr <- fitr$fitted.values
fit0 <- update(fitr,.~1)
yhat0 <- fit0$fitted.values
llr <- sum(log(dgamma(y,shape=100,scale=yhatr/100)))
ll0 <- sum(log(dgamma(y,shape=100,scale=yhat0/100)))

# Likelihood-ratio-based R-squared
1-exp(-2*(llr-ll0)/n)

# Corrected likelihood-ratio-based R-squared
(1-exp(-2*(llr-ll0)/n))/(1-exp(2*ll0/n))
}
