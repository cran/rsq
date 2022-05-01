\name{simglmm}
\alias{simglmm}
\title{Simulate Data from Generalized Linear Mixed Models}
\description{Simulate data from linear and generalized linear mixed models. The coefficients of the two covariate are specified by \code{beta}.}
\usage{simglmm(family=c("binomial","gaussian","poisson","negative.binomial"),
beta=c(2,0),tau=1,n=200,m=10,balance=TRUE)}
\arguments{
    \item{family}{ the family of the distribution.}
    \item{beta}{ regression coefficients (excluding the intercept which is set as zero).}
    \item{tau}{ the variance of the random intercept.}
    \item{n}{ the sample size.}
    \item{m}{ the number of groups.}
    \item{balance}{ simulate balanced data if TRUE, unbalanced data otherwise.}
}
\details{
The first covariate takes 1 in half of the observations, and 0 or -1 in the other half. When \code{beta} gets larger, it is supposed to easier to predict the response variable.
}
\value{
Returned values include \code{yx}, \code{beta}, and \code{u}.
\item{yx}{ a data frame including the response \code{y} and covariates \code{x1}, \code{x2}, and so on.}
\item{beta}{ true values of the regression coefficients.}
\item{u}{ the random intercepts.}
}
\references{
Zhang, D. (2020). Coefficients of determination for generalized linear mixed models. \emph{Technical Report}, 20-01, Department of Statistics, Purdue University.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq}, \link{rsq.lmm}, \link{rsq.glmm}, \link{simglm}},
}
\examples{
require(lme4)

# Linear mixed models
gdata <- simglmm(family="gaussian")
lmm1 <- lmer(y~x1+x2+(1|subject),data=gdata$yx)
rsq(lmm1)

# Generalized linear mixed models
bdata <- simglmm(family="binomial",n=400,m=20)
glmm1 <- glmer(y~x1+x2+(1|subject),family="binomial",data=bdata$yx)
rsq(glmm1)
}
