\name{qvresidual}
\alias{qvresidual}
\title{
Residual Defined on the Basis of a Quardratic Variance Function
}
\description{
Calculate the residual for given observed y and its fitted value yfit, i.e. the distance between them along the quardratic variance function:
 \eqn{V(\mu) = v_2*\mu^2+v_1*\mu+v_0}.}
\usage{qvresidual(y,yfit,v2,v1)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ observed value;}
  \item{yfit}{ fitted value;}
  \item{v2}{ \eqn{v_2};}
  \item{v1}{ \eqn{v_1}.}
}
\details{
This is an internal function called by another internal function \link{vresidual}and \link{vresidual.nb}.
}
\value{The quardratic-variance-function-based residual.}
\references{
Zhang, D. (2016). A coefficient of determination for generalized linear models. \emph{The American Statistician}, tentatively accepted.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\examples{
# Variance function of binomial distribution: V(mu)=mu*(1-mu)
v2 <- -1
v1 <- 1
y <- 0
yfit <- 0.5
vr <- qvresidual(y,yfit,v2,v1)
vr
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
