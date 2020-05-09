\name{pcor}
\alias{pcor}
\title{
Partial Correlation for Generalized Linear Models
}
\description{
Calculate the partial correlation for both linear and generalized linear models.
}
\usage{
pcor(objF,objR=NULL,adj=FALSE,type=c('v','kl','sse','lr','n'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{objF}{  an object of class "lm" or "glm", a result of a call to \link{lm}, \link{glm}, or \link{glm.nb} to fit the full model.}
  \item{objR}{  an object of class "lm" or "glm", a result of a call to \link{lm}, \link{glm}, or \link{glm.nb} to fit the reduced model.}
  \item{adj}{ logical; if TRUE, calculate the adjusted partial R^2.}
  \item{type}{  the type of R-squared used:
  
  'v' (default) -- variance-function-based (Zhang, 2016), calling \link{rsq.v};
  
  'kl' -- KL-divergence-based (Cameron and Windmeijer, 1997), calling \link{rsq.kl}; 
  
  'sse' -- SSE-based (Efron, 1978), calling \link{rsq.sse}; 
  
  'lr' -- likelihood-ratio-based (Maddala, 1983; Cox and Snell, 1989; Magee, 1990), calling \link{rsq.lr}; 
  
  'n' -- corrected version of 'lr' (Nagelkerke, 1991), calling \link{rsq.n}.}
}
\details{
When the fitting object of the reduced model is not specified, the partial correlation of each covariate (excluding factor covariates with more than two levels) in the model will be calculated.
}
\value{
The partial correlation coefficient is returned.
}
\references{
Cameron, A. C. and Windmeijer, A. G. (1997) An R-squared measure of goodness of fit for some common nonlinear regression models. \emph{Journal of Econometrics}, 77: 329-342.

Cox, D. R. and Snell, E. J. (1989) \emph{The Analysis of Binary Data}, 2nd ed. London: Chapman and Hall.

Efron, B. (1978) Regression and ANOVA with zero-one data: measures of residual variation. \emph{Journal of the American Statistical Association}, 73: 113-121.

Maddala, G. S. (1983) \emph{Limited-Dependent and Qualitative Variables in Econometrics}. Cambridge University.

Magee, L. (1990) R^2 measures based on Wald and likelihood ratio joint significance tests. \emph{The American Statistician}, 44: 250-253.

Nagelkerke, N. J. D. (1991) A note on a general definition of the coefficient of determination. \emph{Biometrika}, 78: 691-692.

Zhang, D. (2017). A coefficient of determination for generalized linear models. \emph{The American Statistician}, 71(4): 310-316.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq}, \link{rsq.partial}}.
}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
rsq.partial(bnfit)

bnfitr <- glm(y~color+weight,family=binomial)
rsq.partial(bnfit,bnfitr)

quasibn <- glm(y~color+spine+width+weight,family=quasibinomial)
rsq.partial(quasibn)

quasibnr <- glm(y~color+weight,family=binomial)
rsq.partial(quasibn,quasibnr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
