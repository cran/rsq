# fitObj: object from lm() or glm() or glm.nb();
# adj: whether to calculate the adjusted R-squared.
# type: specifying the extension of R-squared to generalized linear models,
#   v -- variance-function-based (Zhang, 2017),
#   kl -- Kullback-Leibler-divergence-based (Cameron and Windmeijer, 1997),
#   sse -- sum-of-squared-errors-based (Efron, 1978),
#   lr -- likelihood-ratio-statistic-based (Maddala, 1983; Cox & Snell, 1989; Magee, 1990),
#   n -- corrected generalization of 'lr' (Nagelkerke, 1991)
#
rsq<-function(fitObj,adj=FALSE,type=c('v','kl','sse','lr','n'))
{
  if( is(fitObj,"glm")|is(fitObj,"glmerMod") ) # glm in stats, glmer & glmer.nb in lme4
  {
    type <- type[[1]]
    rsq <- switch(type,
        v = rsq.v(fitObj,adj=adj), 
        kl = rsq.kl(fitObj,adj=adj),
        sse = rsq.sse(fitObj,adj=adj),
        lr = rsq.lr(fitObj,adj=adj),
        n = rsq.n(fitObj,adj=adj))
  }
  else if( is(fitObj,"glmmPQL") ) # glmmPQL in MASS, which is also "belongs to "lme" 
    warning("Unsupported object!")
  else if( is(fitObj,"lmerMod")|is(fitObj,"lme") )  # lmer in lme4, lme in nlme
    rsq <- rsq.lmm(fitObj,adj=adj)
  else if( is(fitObj,"lm") )
    rsq <- ifelse(adj,summary(fitObj)$adj.r.squared,summary(fitObj)$r.squared)
  else if( is(fitObj,"deming")) {
    if( is.null(fitObj$model) ) warning("Set model=TRUE in deming(...)!")
    rsq <- 1-ifelse(adj,(var(fitObj$residuals)/var(fitObj$model[,1]))*(1+1/(fitObj$n-2)),
                    var(fitObj$residuals)/var(fitObj$model[,1]))
  }
  else
    warning("Unsupported object!")

  rsq  
}



# objF: object from fitting the full model;
# objR: object from fitting the reduced model.
rsq.partial<-function(objF,objR=NULL,adj=FALSE,type=c('v','kl','sse','lr','n'))
{
  if( !is(objF,"glm") && !is(objF,"lm") ) # only supports lm & glm in stats
    stop("Unsupported object!")
  
  if( is.null(objF$model) )
    stop('The model frame of objF is unavailable! Turn on model=TRUE in fitting!')
  
  type <- type[[1]]
  n <- ifelse(is.null(weights(objF)),nobs(objF),sum(weights(objF)))
  pF <- length(objF$coefficients)

  rsqF <- rsq(objF,adj=FALSE,type=type)
  if( !is.null(objR) )
  {
    rsqR <- rsq(objR,adj=FALSE,type=type)

    pR <- length(objR$coefficients)
    prsq <- 1-((1-rsqF)/(1-rsqR))*ifelse(adj,(n-pR)/(n-pF),1)
    
    list(adjustment=adj,variables.full=attr(terms(objF),"term.labels"),
         variables.reduced=attr(terms(objR),"term.labels"),partial.rsq=prsq)
  }
  else # Return the partial correlation coefficient for each term in the model
  {
    nterms <- length(attr(terms(objF),"term.labels"))
    prsq <- double(length=nterms)
    for(j in 1:nterms)
    {
      if( pmatch("Negative Binomial",family(objF)$family,nomatch=F) )
      {
        theta <- ifelse(is.null(objF$theta),
                        as.numeric(gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.","",
                                        family(objF)$family, perl=T)),objF$theta)
        yy <- with(attributes(terms(objF)),as.character(variables[response+1]))
        if( nterms==1 )
          m1 <- as.formula(paste(yy,"~1"))
        else
          m1<-as.formula(paste(yy,"~",paste(attr(terms(objF),"term.labels")[-j],collapse="+")))
        
        tfit <- glm(m1,family=negative.binomial(theta=theta,link=family(objF)$link),data=objF$model)
      }
      else
      {
        if( nterms==1 )
          m1 <- as.formula(".~1")
        else
          m1<-as.formula(paste(".~",paste(attr(terms(objF),"term.labels")[-j],collapse="+")))
        
        tfit <- update(objF,m1,family=family(objF))
      }
        
      rsqR <- rsq(tfit,adj=FALSE,type=type)
      if( rsqR>1-1e-12 )
      {
        rsqR <- 1-1e-12
        rsqF <- 1-1e-12
      }

      pR <- length(tfit$coefficients)
      prsq[j] <- 1-((1-rsqF)/(1-rsqR))*ifelse(adj,(n-pR)/(n-pF),1)
      
    }
    list(adjustment=adj,variable=attr(terms(objF),"term.labels"), partial.rsq=prsq)
  }
}



# objF: object from fitting the full model;
# objR: object from fitting the reduced model.
pcor<-function(objF,objR=NULL,adj=FALSE,type=c('v','kl','sse','lr','n'))
{
  
  if( !is(objF,"glm") && !is(objF,"lm") ) # only supports lm & glm in stats
    stop("Unsupported object!")
  
  type <- type[[1]]
  nterms <- length(attr(terms(objF),"term.labels"))
  
  if( !is.null(objR) )
  {
    if( objF$rank-objR$rank==1 )
    {
      tmp <- rsq.partial(objF,objR=objR,adj=adj,type=type)
      tmp$partial.rsq <- tmp$partial.rsq*(tmp$partial.rsq>0)
      
      idVar <- !(names(objF$coefficients)%in%names(objR$coefficients))
      pcorr <- sign(objF$coefficients[idVar])*sqrt(tmp$partial.rsq)
      variable <- names(objF$coefficients)[idVar]
    }
    else
      stop("Too many covariates!")
  }
  else if( nterms==1 )
  {
    tmp <- rsq.partial(objF,objR=objR,adj=adj,type=type)
    pcorr <- sqrt(tmp$partial.rsq*(tmp$partial.rsq>0))
    variable <- tmp$variable
  }
  else
  {
    tmp <- rsq.partial(objF,objR=objR,adj=adj,type=type)
    pcorr <- sqrt(tmp$partial.rsq*(tmp$partial.rsq>0))

    # For a categorical predictor with more than two categories, output "NA"
    #     "factor" or "ordered"
    nterms <- length(attr(terms(objF),"term.labels"))
    nLevels <- rep(1,nterms) # No. of categories in each term

    varcls <- attr(terms(objF),"dataClasses")[-1]  # Only for main effects!
    tfactors <- attr(terms(objF),"factors")[-1,]
    trnames <- row.names(tfactors)
    for( j in 1:nterms )
    {
      tcovs <- trnames[tfactors[,j]>0] # covariates defined current term
      if( (length(tcovs)==1)&&(any(varcls[tcovs]=="factor")|any(varcls[tcovs]=="ordered")) )
      { # main effects with factors
        nLevels[j] <- nlevels(objF$model[1,tcovs])-1
      }
      else if( (length(tcovs)>1)&&(any(varcls[tcovs]=="factor")|any(varcls[tcovs]=="ordered")) )
      {  # interactive effects with factors
        tlevels <- rep(1,length(tcovs))
        for( kk in 1:length(tcovs) )
        {
          tlevels[kk] <- ifelse((varcls[tcovs[kk]]=="factor")|(varcls[tcovs[kk]]=="ordered"),
                                nlevels(objF$model[1,tcovs[kk]])-1,1)
        }
        nLevels[j] <- prod(tlevels)
      }
    }

    pcorr[nLevels>1] <- NA
    outIdx <- cumsum(c(1,nLevels))[-1]
    outIdx <- outIdx[nLevels==1]
    
    pcorr[nLevels==1] <- sign(objF$coefficients[outIdx])*pcorr[nLevels==1]
    variable <- tmp$variable
  }

  list(adjustment=adj,variable=variable, partial.cor=as.numeric(pcorr))
}



# Cameron and Windmeijer (1997): Kullback-Leibler-divergence-based R-squared.
#    = 1-Deviance(Current Model)/Deviance(Null Model)
# 
# Make sure that all models have the same scale parameters!!!!!! (glm.nb?)
# Q: How to extend it to quasi models (using quasi-likelihood?)
rsq.kl<-function(fitObj,adj=FALSE)
{
  if( is(fitObj,"glm") )
  {
    sse1 <- fitObj$deviance
    sse0 <- fitObj$null.deviance
    
    n <- ifelse(is.null(weights(fitObj)),nobs(fitObj),sum(weights(fitObj)))
    pM <- fitObj$df.null-fitObj$df.residual+1
    #rsq <- 1-(sse1/sse0)*ifelse(adj,fitObj$df.null/fitObj$df.residual,1)
    rsq <- 1-(sse1/sse0)*ifelse(adj,(n-1)/(n-pM),1)
  }
  else if( is(fitObj,"lm") )
    rsq <- ifelse(adj,summary(fitObj)$adj.r.squared,summary(fitObj)$r.squared)
  else
    stop("Unsupported object!")
  
  rsq
}



# Efron (1978): sums of squared errors.
rsq.sse<-function(fitObj,adj=FALSE)
{
  # Calculate the sum of squared errors  
  sse<-function(fitObj)
  {
    y <- fitObj$y
    tw <- weights(fitObj)
    if( is.null(tw) )
      tw <- y*0+1  # tw <- rep(0,length(y))
    
    if( is.null(y) ) stop("The glm object does not include response values!")
    yfit <- fitObj$fitted.values
    sse <- sum(tw*(y-yfit)^2)
    
    sse
  }

  if( is(fitObj,"glm") )
  {
    sse1 <- sse(fitObj)
    
    #fitObj0 <- update(fitObj,.~1,data=fitObj$model)
    fitObj0 <- update(fitObj,.~1,family=family(fitObj))
    sse0 <- sse(fitObj0)

    n <- ifelse(is.null(weights(fitObj)),nobs(fitObj),sum(weights(fitObj)))
    pM <- fitObj$df.null-fitObj$df.residual+1
    #rsq <- 1-(sse1/sse0)*ifelse(adj,fitObj0$df.residual/fitObj$df.residual,1)
    rsq <- 1-(sse1/sse0)*ifelse(adj,(n-1)/(n-pM),1)
  }
  else if( is(fitObj,"lm") )
    rsq <- ifelse(adj,summary(fitObj)$adj.r.squared,summary(fitObj)$r.squared)
  else
    stop("Unsupported object!")
  
  rsq
}



# Maddala (1983), Cox and Snell (1989), and Magee (1990): R-squared based on the
# likelihood ratio statistics.
rsq.lr<-function(fitObj,adj=FALSE)
{
  if( is(fitObj,"glm") )
  {
    n <- fitObj$df.null+1
    #k <- fitObj$rank
    logLik1 <- as.numeric(logLik(fitObj))
    
    #fitObj0 <- update(fitObj,.~1,data=fitObj$model)
    fitObj0 <- update(fitObj,.~1,family=family(fitObj))
    logLik0 <- as.numeric(logLik(fitObj0))

    rsq <- 1-exp(-2*(logLik1-logLik0)/n)

    if( adj )
      warning("No adjusted version of rsq.lr. Reset 'adj = FALSE' ...")
  }
  else if( is(fitObj,"lm") )
    rsq <- ifelse(adj,summary(fitObj)$adj.r.squared,summary(fitObj)$r.squared)
  else
    stop("Unsupported object!")

  rsq
}



# Nagelkerke (1991):  corrected generalization of 'lr'.
rsq.n<-function(fitObj,adj=FALSE)
{
  if( is(fitObj,"glm") )
  {
    n <- fitObj$df.null+1
    #k <- fitObj$rank
    logLik1 <- as.numeric(logLik(fitObj))
    
    #fitObj0 <- update(fitObj,.~1,data=fitObj$model)
    fitObj0 <- update(fitObj,.~1,family=family(fitObj))
    logLik0 <- as.numeric(logLik(fitObj0))
    
    rsq <- (1-exp(-2*(logLik1-logLik0)/n))/(1-exp(logLik0*2/n))
  }
  else
    stop("Unsupported object!")
  
  if( adj )
    warning("No adjusted version of rsq.n! Reset 'adj=FALSE' ...")

  rsq
}


# Simulate generalized linear models as in Zhang (2017)
simglm<-function(family=c("binomial","gaussian","poisson","Gamma"),lambda=3,n=50,p=3) 
{
  family <- family[[1]]

  beta <- matrix(rep(0,p),nrow=p,ncol=1) # No intercept
  beta[1] <- lambda
  x <- matrix(rnorm(n*p,mean=0,sd=1),nrow=n,ncol=p)
  x[,1] <- rbind(rep(1,n/2),rep(-1,n/2))
  
  switch(family,
    binomial={
      mu <- 1/(1+exp(-x%*%beta))
      y <- rbinom(n,1,mu)},
    gaussian={
      mu <- x%*%beta
      y <- rnorm(n,mean=mu,sd=1)},
    poisson={ 
      mu <- exp(x%*%beta)
      y <- rpois(n,lambda=mu)},
    #negative.binomial={   # Need to specify a better model here!
    #  size <- 2
    #  mu <- size/(exp(-x%*%beta)-1)
    #  y <- rnbinom(n,size=2,mu=mu)},
    Gamma={
      tscale <- 0.01/(lambda+1+x%*%beta)
      y <- rgamma(n,shape=100,scale=tscale)})

  list(yx=data.frame(y=y,x=x),beta=beta)
}



# Linear Mixed Models (LMM):
#fitObj: object from lmer() in lme4 or lme() in nlme;
#adj: whether to calculate the adjusted R-squared.
rsq.lmm<-function(fitObj,adj=FALSE)
{
  if( adj ) n <- nobs(fitObj)

  if( is(fitObj,"lmerMod") ) # lme4 & methods(class="merMod")
  { 
    if( adj ) ndf <- n-length(lme4::fixef(fitObj))

    y <- getME(fitObj,"y")
    mresidual <- y-getME(fitObj,"X")%*%lme4::fixef(fitObj)
    
    ZLambda <- getME(fitObj,"Z")%*%getME(fitObj,"Lambda")
    tau2 <- Matrix::diag(ZLambda%*%Matrix::t(ZLambda))
  }
  else if( is(fitObj,"lme") ) # nlme & methods(class="lme")
  {
    if( adj ) ndf <- n-length(nlme::fixef(fitObj))

    y <- getResponse(fitObj)
    mresidual <- resid(fitObj,type="response",level=0) # marginal
    
    Z <- model.matrix(formula(fitObj$modelStruct$reStr)[[1]],data=fitObj$data)
    tau2 <- apply(Z,1,function(x) matrix(x,nrow=1)%*%getVarCov(fitObj,
                  type="random.effects")%*%t(matrix(x,nrow=1)))
    tau2 <- tau2/(sigma(fitObj)^2)
  }
  else
    stop("Unsupported object!")
  
  SST <- sum((y-mean(y))^2)
  rsqFix <- 1-sum(mresidual^2)/SST
  if( adj )
    rsqFix <- 1-(1-rsqFix)*(n-1)/ndf

  tmpRatio <- 1/(1+tau2)
  if( adj )
    rsqMod <- 1-(sum(tmpRatio*(tau2*(sigma(fitObj)^2)+
                tmpRatio*(mresidual^2)))/ndf)/(SST/(n-1))
  else
    rsqMod <- 1-sum(tmpRatio*(tau2*(sigma(fitObj)^2)+tmpRatio*(mresidual^2)))/SST
  
  rsqRnd <- rsqMod-rsqFix
  list(model=rsqMod,fixed=rsqFix,random=rsqRnd)
}



# Calculate the variance-function-based residual for given observed y and
# its fitted value yfit, with prespecified family or variance function.
vresidual<-function(y,yfit,family=binomial(),variance=NULL)
{
  # Calculate the residual for given observed y and its fitted value yfit: 
  # the length between y and yfit along the quardratic variance function:
  #           V(mu) = v2*mu^2+v1*mu+v0
  qvresidual<-function(y,yfit,v2,v1)
  {
    vpa <- 2*v2*yfit+v1
    svpa2 <- sqrt(1+vpa*vpa)
    
    vpb <- 2*v2*y+v1
    svpb2 <- sqrt(1+vpb*vpb)
    
    vr <- (log((vpb+svpb2)/(vpa+svpa2))+vpb*svpb2-vpa*svpa2)/(4*v2)
    vr
  }
  
  if( is.character(family) ) {
    cf <- family
    family <- get(family, mode="function", envir=parent.frame())
  } else
    cf <- family$family

  if( pmatch("Negative Binomial",cf,nomatch=F) )
  {
    theta <- as.numeric(gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.","", cf, perl=T))
    cf <- "negative.binomial"
  } else if( pmatch("Tweedie",cf,nomatch=F) )
  {
    dv <- Deriv(family$variance,"mu")
    theta <- dv(1)
  }
  
  if( is.null(variance) )
  {
    switch(cf,
      binomial={DFUN<-function(x) qvresidual(x[1],x[2],-1,1)}, # Need modify for Y~Bin(n,p)
      gaussian={DFUN<-function(x) x[1]-x[2]},
      Gamma={DFUN<-function(x) qvresidual(x[1],x[2],1,0)},
      negative.binomial={DFUN<-function(x) qvresidual(x[1],x[2],1/theta,1)},
      poisson={DFUN<-function(x) x[1]-x[2]},
      quasibinomial={DFUN<-function(x) qvresidual(x[1],x[2],-1,1)},
      quasipoisson={DFUN<-function(x) x[1]-x[2]},
      inverse.gaussian={
        DFUN<-function(x) integrate(function(mu){sqrt(1+9*mu^4)},x[1],x[2])$value},
      Tweedie={ # var.power: 0, 1, (1,2), 2, >2
        if( (theta==0)|(theta==1) )
          DFUN<-function(x) x[1]-x[2]
        else if( theta==2 ) 
          DFUN<-function(x) qvresidual(x[1],x[2],1,0)
        else
          DFUN<-function(x) integrate(function(mu){sqrt(1+theta^2*mu^(2*theta-2))},x[1],x[2])$value},
      quasi={  # variance for quasi: "constant","mu(1-mu)","mu","mu^2","mu^3", or other
        if( (family$varfun=="constant")|(family$varfun=="mu") )
          DFUN <- function(x) x[1]-x[2]
        else if( family$varfun=="mu(1-mu)" )
          DFUN<-function(x) qvresidual(x[1],x[2],-1,1)
        else if( family$varfun=="mu^2" )
          DFUN<-function(x) qvresidual(x[1],x[2],1,0)
        else
          DFUN<-function(x) integrate(function(mu){sqrt(1+Deriv(family$variance,"mu")^2)},x[1],x[2])$value})
  }
  else
    DFUN<-function(x) integrate(function(mu){sqrt(1+Deriv(variance,"mu")^2)},x[1],x[2])$value
  
  vresidual <- apply(cbind(y,yfit),1,DFUN)
}



# Zhang (2017): variance-function-based R-squared.
rsq.v<-function(fitObj,adj=FALSE)
{
  if( is(fitObj,"glmerMod") ) # glmer or glmer.nb in lme4
  { # methods(class="merMod")
    rsq <- rsq.glmm(fitObj,adj=adj)
  }
  else if( is(fitObj,"glm") ) # glm in stats
  {
    y <- fitObj$y
    wt <- weights(fitObj)
    if( is.null(wt) )
      wt <- y*0+1
    
    n <- sum(wt)
    if(is.null(y)) stop("The glm object does not include response values!")
    yfit <- fitObj$fitted.values
    
    if(pmatch("Negative Binomial",family(fitObj)$family,nomatch=F))
    {
      theta <- ifelse(is.null(fitObj$theta),
                      as.numeric(gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.","",
                                      family(fitObj)$family, perl=T)), fitObj$theta)
      sse1 <- sum(wt*vresidual(y,yfit,family=negative.binomial(theta))^2)
      
      f0 <- glm(y~1,family=negative.binomial(theta))
      yf0 <- f0$fitted.values
      sse0 <- sum(wt*vresidual(y,yf0,family=negative.binomial(theta))^2)
    }
    else if(family(fitObj)$family=="binomial")
    {
      nSuc <- wt*y
      nFai <- wt-nSuc
      tone <- rep(1,length(nSuc))
      sse1 <- sum(nSuc*vresidual(tone,yfit,family=family(fitObj))^2)+
        sum(nFai*vresidual(1-tone,yfit,family=family(fitObj))^2)
      
      f0 <- update(fitObj,.~1,family=family(fitObj))
      yf0 <- f0$fitted.values
      sse0 <- sum(nSuc*vresidual(tone,yf0,family=family(f0))^2)+
                sum(nFai*vresidual(1-tone,yf0,family=family(f0))^2)
    }
    else
    {
      sse1 <- sum(wt*vresidual(y,yfit,family=family(fitObj))^2)
      
      f0 <- update(fitObj,.~1,family=family(fitObj))
      sse0 <- sum(wt*vresidual(y,f0$fitted.values,family=family(f0))^2)
    }
    
    #rsq <- 1-(sse1/sse0)*ifelse(adj,f0$df.residual/fitObj$df.residual,1)
    pM <- fitObj$df.null-fitObj$df.residual+1
    rsq <- 1-(sse1/sse0)*ifelse(adj,(n-1)/(n-pM),1)
    rsq <- max(0,rsq)
  }
  
  rsq
}



#fitObj: object from glmer() or glmer.nb in lme4;
#adj: whether to calculate the adjusted R-squared.
rsq.glmm<-function(fitObj,adj=FALSE)
{
  y <- getME(fitObj,"y")
  yhat <- getME(fitObj,"mu")

  wt <- weights(fitObj)
  if( is.null(wt) )
    wt <- y*0+1
  
  glmf <- formula(paste(as.character(terms(fitObj))[2],
            as.character(terms(fitObj))[1],as.character(terms(fitObj))[3]))
  if( pmatch("Negative Binomial",family(fitObj)$family,nomatch=F) )
  {
    theta <- getME(fitObj,"glmer.nb.theta")
    theta <- ifelse(is.null(theta),
                    as.numeric(gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.","",
                                    family(fitObj)$family, perl=T)), theta)
    
    # SST
    f0 <- glm(y~1,family=negative.binomial(theta=theta,link=family(fitObj)$link),weights=wt)
    yf0 <- f0$fitted.values
    SST <- sum(wt*vresidual(y,yf0,family=negative.binomial(theta))^2)
    
    # SSE for GLM only
    f1 <- glm(glmf,family=negative.binomial(theta=theta,link=family(fitObj)$link),
              data=getData(fitObj),weights=wt)
    yf1 <- f1$fitted.values
    SSE <- sum(wt*vresidual(y,yf1,family=negative.binomial(theta))^2)

    rsqFix <- 1-SSE/SST
    rsqMod <- 1-sum(wt*vresidual(y,yhat,family=negative.binomial(theta))^2)/SST
  }
  else if( family(fitObj)$family=="binomial" )
  {
    #SST
    nSuc <- wt*y
    nFai <- wt-nSuc
    tone <- rep(1,length(nSuc))
    ymean <- sum(nSuc)/sum(wt)
    SST <- sum(nSuc*vresidual(tone,ymean,family=family(fitObj))^2)+
      sum(nFai*vresidual(1-tone,ymean,family=family(fitObj))^2)
    
    #SSE for GLM only
    f1 <- glm(glmf,family=binomial(link=family(fitObj)$link),
              data=getData(fitObj),weights=wt)
    yf1 <- f1$fitted.values
    SSE <- sum(nSuc*vresidual(tone,yf1,family=family(fitObj))^2+
               nFai*vresidual(1-tone,yf1,family=family(fitObj))^2)

    rsqFix <- 1-SSE/SST
    rsqMod <- 1-sum(nSuc*vresidual(tone,yhat,family=family(fitObj))^2+
                      nFai*vresidual(1-tone,yhat,family=family(fitObj))^2)/SST
  }
  else
  {
    SST <- sum(wt*vresidual(y,mean(y),family=family(fitObj))^2)
    
    # SSE for GLM only
    f1 <- glm(glmf,family=family(fitObj),data=getData(fitObj),weights=wt)
    yf1 <- f1$fitted.values
    SSE <- sum(wt*vresidual(y,yf1,family=family(fitObj))^2)
    
    rsqFix <- 1-SSE/SST
    rsqMod <- 1-sum(wt*vresidual(y,yhat,family=family(fitObj))^2)/SST
  }
  
  if( adj ) {
    n <- nobs(fitObj)
    rsqFix <- 1-(1-rsqFix)*(n-1)/df.residual(f1)
    rsqMod <- 1-(1-rsqMod)*(n-1)/df.residual(fitObj)
  }
  rsqRnd <- rsqMod-rsqFix
  
  list(model=max(0,rsqMod),fixed=max(0,rsqFix),random=max(0,rsqRnd))
}



# Simulate generalized linear models as in Zhang (2019+)
#    g(mu_{ij}) = beta*X_{ij} + u_i, j=1, 2, ..., [n/m]; i=1, 2, ..., m.
#    u_i ~ N(0,tau^2)
# So there are m groups with K=[n/m] subjects in each panel.
#
simglmm<-function(family=c("binomial","gaussian","poisson","negative.binomial"),
                  beta=c(2,0),tau=1,n=200,m=10,balance=TRUE)
{
  family <- family[[1]]
  
  p <- 3
  if( balance ) {
    pSize <- rep(round(n/m),m)
    pSize[m] <- n-sum(pSize[-m])
  } else
    pSize <- rmultinom(1,n-2*m,rep(1,m)/m)+2
  
  # Group ID
  grpID <- NULL
  for( j in 1:m )
    grpID <- rbind(grpID,matrix(j,nrow=pSize[j],ncol=1))

  # Random Effects
  rndEff <- rep(mvrnorm(n=1,mu=rep(0,m),Sigma=diag(m)*(tau^2)),pSize) # vector
  
  # X & eta
  x <- matrix(rnorm(n*p,mean=0,sd=1),nrow=n,ncol=p)
  x[,1] <- rep(c(1,-1),n/2)
  eta <- x[,1]*beta[1]+x[,2]*beta[2]+rndEff
  
  # Y
  switch(family,
         binomial={y<-rbinom(n,1,1/(1+exp(-eta)))},
         gaussian={y<-rnorm(n,mean=eta,sd=1)},
         poisson={y<-rpois(n,lambda=exp(eta))},
         negative.binomial={y<-rnbinom(n,size=10,prob=1/(1+exp(-eta)))})
  
  list(yx=data.frame(y=y,x1=x[,1],x2=x[,2],x3=x[,3],subject=grpID),beta=beta,u=rndEff)
}