\name{rsq.n}
\alias{rsq.n}
\title{Corrected Likelihood-Ratio-Based R-Squared}
\description{Corrected likelihood-ratio-based R^2 for generalized linear models.}
\usage{rsq.n(fitObj,adj=FALSE,data=NULL)}
\arguments{
    \item{fitObj}{ an object of class "lm" or "glm", usually, a result of a call to \link{lm}, \link{glm}, or \link{glm.nb}.}
    \item{adj}{ logical; if TRUE, calculate the adjusted R^2.}
    \item{data}{an optional data frame used by fitObj.}
}
\details{
Nagelkerke (1991) proposed this version of R^2 to correct the likelihood-ratio-statistic-based one which was proposed by Maddala (1983), Cox and Snell (1989), and Magee (1990). This corrected generalization of R^2 cannot reduce to the classical R^2 in case of linear models. It is not defined for quasi models.
}
\value{The R^2 or adjusted R^2.}
\references{
Cox, D. R. and Snell, E. J. (1989) \emph{The Analysis of Binary Data}, 2nd ed. London: Chapman and Hall.

Maddala, G. S. (1983) \emph{Limited-Dependent and Qualitative Variables in Econometrics}. Cambridge University.

Magee, L. (1990) R^2 measures based on Wald and likelihood ratio joint significance tests. \emph{The American Statistician}, 44: 250-253.

Nagelkerke, N. J. D. (1991) A note on a general definition of the coefficient of determination. \emph{Biometrika}, 78: 691-692.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq}, \link{rsq.partial}, \link{pcor}, \link{rsq.lr}}.
}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
rsq.n(bnfit)
rsq.n(bnfit,adj=TRUE)

psfit <- glm(num.satellites~color+spine+width+weight,family=poisson)
rsq.n(psfit)
rsq.n(psfit,adj=TRUE)

# Effectiveness of Bycycle Safety Helmets in Thompson et al. (1989)
y <- matrix(c(17,218,233,758),2,2)
x <- factor(c("yes","no"))
tbn <- glm(y~x,family=binomial)
rsq.n(tbn)
rsq.n(tbn,adj=TRUE)
}
