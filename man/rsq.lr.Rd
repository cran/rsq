\name{rsq.lr}
\alias{rsq.lr}
\title{Likelihood-Ratio-Based R-Squared}
\description{Calculate the likelihood-ratio-based R^2 for generalized linear models.}
\usage{rsq.lr(fitObj,adj=FALSE)}
\arguments{
    \item{fitObj}{ an object of class "lm" or "glm", usually, a result of a call to \link{lm}, \link{glm}, or \link{glm.nb}.}
    \item{adj}{ logical; if TRUE, calculate the adjusted R^2.}
}
\details{
Proposed by Maddala (1983), Cox and Snell (1989), and Magee (1990), this version of R^2 is defined with the likelihood ratio statistics, so it is not defined for quasi models. It reduces to the classical R^2 when the variance function is constant or linear.
}
\value{The R^2 or adjusted R^2.}
\references{
Cox, D. R. and Snell, E. J. (1989) \emph{The Analysis of Binary Data}, 2nd ed. London: Chapman and Hall.

Maddala, G. S. (1983) \emph{Limited-Dependent and Qualitative Variables in Econometrics}. Cambridge University.

Magee, L. (1990) R^2 measures based on Wald and likelihood ratio joint significance tests. \emph{The American Statistician}, 44: 250-253.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq}, \link{rsq.partial}, \link{pcor}, \link{rsq.n}}.
}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
rsq.lr(bnfit)
rsq.lr(bnfit,adj=TRUE)

psfit <- glm(num.satellites~color+spine+width+weight,family=poisson)
rsq.lr(psfit)
rsq.lr(psfit,adj=TRUE)

# Effectiveness of Bycycle Safety Helmets in Thompson et al. (1989)
y <- matrix(c(17,218,233,758),2,2)
x <- factor(c("yes","no"))
tbn <- glm(y~x,family=binomial)
rsq.lr(tbn)
rsq.lr(tbn,adj=TRUE)
}
