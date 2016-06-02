\name{vsse}
\alias{vsse}
\title{
Variance-Function-Based Error Sum of Squares
}
\description{
Calculate the error sum of squares which is defined on the basis of the variance function.}
\usage{vsse(fitObj)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fitObj}{ an object of class "lm" or "glm", usually, a result of a call to \link{lm}, \link{glm}, or \link{glm.nb}.}
}
\details{
This is an internal function called by \link{rsq}.
}
\value{Variance-function-based error sum of squares.}
\references{
Zhang, D. (2016). A coefficient of determination for generalized linear models. The American Statistician, tentatively accepted.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
vsse(bnfit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
