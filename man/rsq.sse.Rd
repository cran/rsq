\name{rsq.sse}
\alias{rsq.sse}
\title{SSE-Based R-Squared}
\description{The sum-of-squared-errors-based R^2 for generalized linear models.}
\usage{rsq.sse(fitObj,adj=FALSE,data=NULL)}
\arguments{
    \item{fitObj}{ an object of class "lm" or "glm", usually, a result of a call to \link{lm}, \link{glm}, or \link{glm.nb}.}
    \item{adj}{ logical; if TRUE, calculate the adjusted R^2.}
    \item{data}{an optional data frame used by fitObj.}
}
\details{
This version of R^2 was proposed by Efron (1978). It is calculated on the basis of the formula of the classical R^2. }
\value{The R^2 or adjusted R^2.}
\references{
Efron, B. (1978) Regression and ANOVA with zero-one data: measures of residual variation. \emph{Journal of the American Statistical Association}, 73: 113-121.
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
rsq.sse(bnfit)
rsq.sse(bnfit,adj=TRUE)

psfit <- glm(num.satellites~color+spine+width+weight,family=poisson)
rsq.sse(psfit)
rsq.sse(psfit,adj=TRUE)

# Effectiveness of Bycycle Safety Helmets in Thompson et al. (1989)
y <- matrix(c(17,218,233,758),2,2)
x <- factor(c("yes","no"))
tbn <- glm(y~x,family=binomial)
rsq.sse(tbn)
rsq.sse(tbn,adj=TRUE)
}
