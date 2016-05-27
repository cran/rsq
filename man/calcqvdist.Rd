\name{calcqvdist}
\alias{calcqvdist}
\title{
Distance Between Two Points Along a Quardratic Variance Function
}
\description{
Calculate the distance between two points along a quadratic variance function: \eqn{V(\mu) = v_2*\mu^2+v_1*\mu+v_0}.}
\usage{
calcqvdist(a,b,v2,v1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ starting point.}
  \item{b}{ end point.}
  \item{v2}{ \eqn{v_2}.}
  \item{v1}{ \eqn{v_1}.}
}
\details{
This is an internal function called by another internal function \link{vrsq}.
}
\value{The distance between \code{a} and \code{b} along the variance function.}
\references{
Zhang, D. (2006). A coefficient of determination for generalized linear models. \emph{The American Statistician}, tentatively accepted.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\examples{
# Variance function of binomial distribution: V(mu)=mu*(1-mu)
v2 <- -1
v1 <- 1
a <- 0
b <- 1
vd <- calcqvdist(a,b,v2,v1)
vd
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
