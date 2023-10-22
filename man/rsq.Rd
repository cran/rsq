\name{rsq}
\alias{rsq}
\title{R-Squared for Generalized Linear (Mixed) Models}
\description{Calculate the coefficient of determination, aka R^2, for  both linear and generalized linear (mixed) models.}
\usage{rsq(fitObj,adj=FALSE,type=c('v','kl','sse','lr','n'))}
\arguments{
    \item{fitObj}{ an object of class "lm", "glm", "merMod", "lmerMod", or "lme"; usually a result of a call to \link{lm}, \link{glm}, \link{glm.nb}, \link{lmer} or \link{glmer} or \link{glmer.nb} in \link{lme4}, or \link{lme} in \link{nlme}.}
    \item{adj}{ logical; if TRUE, calculate the adjusted R^2.}
    \item{type}{  the type of R-squared (only applicable for generalized linear models):
    
        'v' (default) -- variance-function-based (Zhang, 2017), calling \link{rsq.v}; 
        
        'kl' -- KL-divergence-based (Cameron and Windmeijer, 1997), calling \link{rsq.kl}; 
        
        'sse' -- SSE-based (Efron, 1978), calling \link{rsq.sse}; 
        
        'lr' -- likelihood-ratio-based (Maddala, 1983; Cox and Snell, 1989; Magee, 1990), calling \link{rsq.lr}; 
        
        'n' -- corrected version of 'lr' (Nagelkerke, 1991), calling \link{rsq.n}.}
}
\details{
Calculate the R-squared for (generalized) linear models. For (generalized) linear mixed models, there are three types of R^2 calculated on the basis of observed response values, estimates of fixed effects, and variance components, i.e., model-based R_M^2 (proportion of variation explained by the model in total, including both fixed-effects and random-efffects factors), fixed-effects R_F^2 (proportion of variation explained by the fixed-effects factors), and random-effects R_R^2 (proportion of variation explained by the random-effects factors).
}
\value{The R^2 or adjusted R^2. For (generalized) linear mixed models,
  \item{R_M^2}{proportion of variation explained by the model in total, including both fixed-effects and random-efffects factors.}
  \item{R_F^2}{proportion of variation explained by the fixed-effects factors.}
  \item{R_R^2}{proportion of variation explained by the random-effects factors.}
}
\references{
Cameron, A. C. and Windmeijer, A. G. (1997) An R-squared measure of goodness of fit for some common nonlinear regression models. \emph{Journal of Econometrics}, 77: 329-342.

Cox, D. R. and Snell, E. J. (1989) \emph{The Analysis of Binary Data}, 2nd ed. London: Chapman and Hall.

Efron, B. (1978) Regression and ANOVA with zero-one data: measures of residual variation. \emph{Journal of the American Statistical Association}, 73: 113-121.

Maddala, G. S. (1983) \emph{Limited-Dependent and Qualitative Variables in Econometrics}. Cambridge University.

Magee, L. (1990) R^2 measures based on Wald and likelihood ratio joint significance tests. \emph{The American Statistician}, 44: 250-253.

Nagelkerke, N. J. D. (1991) A note on a general definition of the coefficient of determination. \emph{Biometrika}, 78: 691-692.

Zhang, D. (2017). A coefficient of determination for generalized linear models. \emph{The American Statistician}, 71(4): 310-316.

Zhang, D. (2020). Coefficients of determination for generalized linear mixed models. \emph{Technical Report}, 20-01, Department of Statistics, Purdue University.
}
\author{Dabao Zhang, Department of Statistics, Purdue University}
\seealso{
\code{\link{rsq.partial},\link{pcor},\link{simglm}}.
}
\examples{
data(hcrabs)
attach(hcrabs)
y <- ifelse(num.satellites>0,1,0)
bnfit <- glm(y~color+spine+width+weight,family=binomial)
rsq(bnfit)
rsq(bnfit,adj=TRUE)

quasibn <- glm(y~color+spine+width+weight,family=quasibinomial)
rsq(quasibn)
rsq(quasibn,adj=TRUE)

psfit <- glm(num.satellites~color+spine+width+weight,family=poisson)
rsq(psfit)
rsq(psfit,adj=TRUE)

quasips <- glm(num.satellites~color+spine+width+weight,family=quasipoisson)
rsq(quasips)
rsq(quasips,adj=TRUE)

# Linear mixed models
require(lme4)
lmm1 <- lmer(Reaction~Days+(Days|Subject),data=sleepstudy)
rsq(lmm1)
rsq.lmm(lmm1)

# Generalized linear mixed models
data(cbpp)
glmm1 <- glmer(cbind(incidence,size-incidence)~period+(1|herd),data=cbpp,family=binomial)
rsq(glmm1)
}
