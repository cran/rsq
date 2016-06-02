\name{rsq.partial}
\alias{rsq.partial}
\title{
Coefficient of Partial Determination
}
\description{
Calculate the coefficient of partial determination, aka partial R^2, for both linear regression models, and generalized linear regression models.
}
\usage{
rsq.partial(objF,objR=NULL,adj=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{objF}{  an object of class "lm" or "glm", a result of a call to \link{lm}, \link{glm}, or \link{glm.nb} to fit the full model.}
  \item{objR}{  an object of class "lm" or "glm", a result of a call to \link{lm}, \link{glm}, or \link{glm.nb} to fit the reduced model.}
  \item{adj}{ logical; if TRUE, calculate the adjusted partial R^2.}
}
\details{
When the fitting object of the reduced model is not specified, the partial R^2 of each term in the model will be calculated.
}
\value{
Returned values include \code{adjustment} and \code{partial.rsq}. When \code{objR} is not \code{NULL}, \code{variable.full} and \code{variable.reduced} are returned; otherwise \code{variable} is returned.
\item{adjustment}{ logical; if TRUE, calculate the adjusted partial R^2.}
\item{variable.full}{ all covariates in the full model.}
\item{variable.reduced}{ all covariates in the reduced model.}
\item{variable}{ all covariates in the full model.}
\item{partial.rsq}{ partial R^2 or the adjusted partial R^2.}
}
\references{
Zhang, D. (2016). A coefficient of determination for generalized linear models. \emph{The American Statistician}, tentatively accepted.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq}}.
}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
rsq.partial(bnfit)
rsq.partial(bnfit,adj=TRUE)

bnfitr <- glm(y~color+weight,family=binomial)
rsq.partial(bnfit,bnfitr)
rsq.partial(bnfit,bnfitr,adj=TRUE)

quasibn <- glm(y~color+spine+width+weight,family=quasibinomial)
rsq.partial(quasibn)
rsq.partial(quasibn,adj=TRUE)

quasibnr <- glm(y~color+weight,family=binomial)
rsq.partial(quasibn,quasibnr)
rsq.partial(quasibn,quasibnr,adj=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
