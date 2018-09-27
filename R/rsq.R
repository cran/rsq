#fitObj: object from lm() or glm() or glm.nb();
#adj: whether to calculate the adjusted R-squared.
#type: specifying the extension of R-squared to generalized linear models,
#      v -- variance-function-based (Zhang, 2016),
#      kl -- Kullback-Leibler-divergence-based (Cameron and Windmeijer, 1997),
#      sse -- sum-of-squared-errors-based (Efron, 1978),
#      lr -- likelihood-ratio-statistic-based (Maddala, 1983; Cox and Snell, 1989; Magee, 1990),
#      n -- corrected generalization of 'lr' (Nagelkerke, 1991)
rsq<-function(fitObj,adj=FALSE,type=c('v','kl','sse','lr','n'),data=NULL)
{
  if( is.null(data) )
  {
    if( is.null(fitObj$model) )
      stop('The model frame in fitObj is unavailable! Turn on model=TRUE in fitting!')

    # Missing values have been removed from model:
    data <- fitObj$model
  }

  type <- type[[1]]
  
  switch(type,
    v = rsq.v(fitObj,adj=adj,data=data), 
    kl = rsq.kl(fitObj,adj=adj),
    sse = rsq.sse(fitObj,adj=adj,data=data),
    lr = rsq.lr(fitObj,adj=adj,data=data),
    n = rsq.n(fitObj,adj=adj,data=data))
}

# objF: object from fitting the full model;
# objR: object from fitting the reduced model.
rsq.partial<-function(objF,objR=NULL,adj=FALSE,type=c('v','kl','sse','lr','n'))
{
  if( is.null(objF$model) )
    stop('The model frame of objF is unavailable! Turn on model=TRUE in fitting!')
  
  # Missinge values are removed from model:
  cdata <- objF$model
  if( !is.null(objF$na.action) )
  {
    # Make sure same missing values are removed from both objF & objR
    if( !is.null(objR) )
      objR <- update(objR,data=cdata)
  }

  type <- type[[1]]

  rsqF <- rsq(objF,adj=FALSE,type=type,data=cdata)
  if( !is.null(objR) )
  {
    rsqR <- rsq(objR,adj=FALSE,type=type,data=cdata)
    prsq <- 1-((1-rsqF)/(1-rsqR))*ifelse(adj,objR$df.residual/objF$df.residual,1)
    
    list(adjustment=adj,variables.full=attr(terms(objF),"term.labels"),
         variables.reduced=attr(terms(objR),"term.labels"),partial.rsq=prsq)
  }
  else # Return the partial correlation coefficient for each term in the model
  {
    nterms <- length(attr(terms(objF),"term.labels"))
    prsq <- double(length=nterms)
    for(j in 1:nterms)
    {
      if( (type=='kl')&is(objF,"negbin") )
      {
        theta <- ifelse(is.null(objF$theta),
                        as.numeric(gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.","",
                                        family(objF)$family, perl=T)),
                        objF$theta)
        y <- model.response(objF$model)
        tfit <- glm(y~terms(objF)[-j],family=negative.binomial(theta=theta,link=family(objF)$link),data=cdata)
      }
      else
        tfit <- update(objF,terms(objF)[-j],data=cdata)
      
      rsqR <- rsq(tfit,adj=FALSE,type=type,data=cdata)
      if( rsqR>1-1e-12 )
      {
        rsqR <- 1-1e-12
        rsqF <- 1-1e-12
      }

      prsq[j] <- 1-((1-rsqF)/(1-rsqR))*ifelse(adj,tfit$df.residual/objF$df.residual,1)
    }
    list(adjustment=adj,variable=attr(terms(objF),"term.labels"), partial.rsq=prsq)
  }
}

# objF: object from fitting the full model;
# objR: object from fitting the reduced model.
pcor<-function(objF,objR=NULL,adj=FALSE,type=c('v','kl','sse','lr','n'))
{
  type <- type[[1]]

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
      warning("Too many covariates!")
  }
  else
  {
    tmp <- rsq.partial(objF,objR=objR,adj=adj,type=type)
    pcorr <- sqrt(tmp$partial.rsq*(tmp$partial.rsq>0))

    nterms <- length(attr(terms(objF),"term.labels"))
    nLevels <- rep(1,nterms)
    varcls <- attr(terms(objF),"dataClasses")[-1]
    nLevels[varcls=="factor"] <- as.numeric(sapply(objF$xlevels,length))-1

    naIdx <- (nLevels>1)&(varcls=="factor")
    pcorr[naIdx] <- NA
    
    outIdx <- cumsum(c(1,nLevels))[-1]
    outIdx <- outIdx[!naIdx]
    
    pcorr[!naIdx] <- sign(objF$coefficients[outIdx])*pcorr[!naIdx]
    variable <- tmp$variable
  }

  list(adjustment=adj,variable=variable, partial.cor=as.numeric(pcorr))
}

# Zhang (2016): variance-function-based R-squared.
rsq.v<-function(fitObj,adj=FALSE,data=NULL)
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
  
  # Calculate the residual for given observed y and its fitted value yfit
  # theta is the additional parameter, e.g., in negative binomial
  vresidual<-function(y,yfit,family="binomial",theta=1)
  {
    vresidual <- 0
    switch(family,
      binomial={       # V(mu) = mu*(1-mu)
        v2 <- -1
        v1 <- 1
        vresidual <- ifelse(ncol(as.matrix(y))>1,
                            c(qvresidual(1,yfit,v2,v1),qvresidual(0,yfit,v2,v1)),
                            qvresidual(y,yfit,v2,v1)) },
      gaussian={       # V(mu) = 1
        vresidual <- y-yfit },
      Gamma={          # V(mu) = mu*mu
        v2 <- 1
        v1 <- 0
        vresidual <- qvresidual(y,yfit,v2,v1) },
      inverse.gaussian={    # V(mu) = mu*mu*mu
        vresidual <- integrate(function(mu){sqrt(1+9*mu^4)},yfit,y)$value },
      negative.binomial={   # V(mu)=mu+mu*mu/theta
        v2 <- 1/theta
        v1 <- 1
        vresidual <- qvresidual(y,yfit,v2,v1) },
      poisson={        # V(mu) = mu
        vresidual <- y-yfit },
      quasibinomial={  # V(mu) = mu*(1-mu)
        v2 <- -1
        v1 <- 1
        vresidual <- qvresidual(y,yfit,v2,v1) },
      quasipoisson={   # V(mu) = mu
        vresidual <- y-yfit },
      quasi={          # V(mu) specified by variance including varfun?
        warning("Not Implemented!") })
    
    vresidual
  }
  
  # Calculate the sum of squared errors
  vsse<-function(fitObj)
  {
    y <- fitObj$y
    if(is.null(y)) warning("The glm object does not include response values!")
    yfit <- fitObj$fitted.values
    
    if(pmatch("Negative Binomial",family(fitObj)$family,nomatch=F))
    {
      theta <- ifelse(is.null(fitObj$theta),
                      as.numeric(gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.","",
                                      family(fitObj)$family, perl=T)),
                      fitObj$theta)
      tresid <- function(x){vresidual(x[1],x[2],family="negative.binomial",theta=theta)}
    }
    else
    {
      tresid <- function(x){vresidual(x[1],x[2],family=family(fitObj)$family)}
    }
    vsse <- sum(weights(fitObj)*apply(cbind(y,yfit),1,tresid)^2)
    
    vsse
  }
  
  if( is.null(data) )
    data <- fitObj$model
  
  #if( is(fitObj,"negbin") )
  if(pmatch("Negative Binomial",family(fitObj)$family,nomatch=F))
  {
    sse1 <- vsse(fitObj)
    theta <- ifelse(is.null(fitObj$theta),
                    as.numeric(gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.","",
                                    family(fitObj)$family, perl=T)),
                    fitObj$theta)
    y <- model.response(fitObj$model)
    fitObj0 <- glm(y~1,family=negative.binomial(theta=theta,link=family(fitObj)$link))
    sse0 <- vsse(fitObj0)
    rsq <- 1-(sse1/sse0)*ifelse(adj,fitObj0$df.residual/fitObj$df.residual,1)
  }
  else if( is(fitObj,"glm") )
  {
    sse1 <- vsse(fitObj)
               
    fitObj0 <- update(fitObj,.~1,data=data)
    sse0 <- vsse(fitObj0)
    rsq <- 1-(sse1/sse0)*ifelse(adj,fitObj0$df.residual/fitObj$df.residual,1)
  }
  else if( is(fitObj,"lm") )
    rsq <- ifelse(adj,summary(fitObj)$adj.r.squared,summary(fitObj)$r.squared)
  else
    warning("Unsupported object!")

  rsq
}

# Cameron and Windmeijer (1997): Kullback-Leibler-divergence-based R-squared.
#    -- 1-Deviance(Current Model)/Deviance(Null Model)
#       Make sure that all models have the same scale parameters!!!!!! (glm.nb?)
#       Q: How to extend it to quasi models (using quasi-likelihood?)
rsq.kl<-function(fitObj,adj=FALSE)
{
  if( is(fitObj,"glm") )
  {
    sse1 <- fitObj$deviance
    sse0 <- fitObj$null.deviance
    rsq <- 1-(sse1/sse0)*ifelse(adj,fitObj$df.null/fitObj$df.residual,1)
  }
  else if( is(fitObj,"lm") )
    rsq <- ifelse(adj,summary(fitObj)$adj.r.squared,summary(fitObj)$r.squared)
  else
    warning("Unsupported object!")
  
  rsq
}

# Efron (1978): sums of squared errors.
rsq.sse<-function(fitObj,adj=FALSE,data=NULL)
{
  if( is.null(data) )
    data <- fitObj$model
  
  # Calculate the sum of squared errors  
  sse<-function(fitObj)
  {
    y <- fitObj$y
    if( is.null(y) ) warning("The glm object does not include response values!")
    yfit <- fitObj$fitted.values
    sse <- sum(weights(fitObj)*(y-yfit)^2)
    
    sse
  }

  if( is(fitObj,"glm") )
  {
    sse1 <- sse(fitObj)
    
    fitObj0 <- update(fitObj,.~1,data=data)
    sse0 <- sse(fitObj0)

    rsq <- 1-(sse1/sse0)*ifelse(adj,fitObj0$df.residual/fitObj$df.residual,1)
  }
  else if( is(fitObj,"lm") )
    rsq <- ifelse(adj,summary(fitObj)$adj.r.squared,summary(fitObj)$r.squared)
  else
    warning("Unsupported object!")
  
  rsq
}

# Maddala (1983), Cox and Snell (1989), and Magee (1990): R-squared based on the likelihood ratio statistics.
rsq.lr<-function(fitObj,adj=FALSE,data=NULL)
{
  if( is.null(data) )
    data <- fitObj$model
  
  if( is(fitObj,"glm") )
  {
    n <- fitObj$df.null+1
    k <- fitObj$rank
    logLik1 <- as.numeric(logLik(fitObj))
    
    fitObj0 <- update(fitObj,.~1,data=data)
    logLik0 <- as.numeric(logLik(fitObj0))

    rsq <- 1-exp(-2*(logLik1-logLik0)/n)
  }
  else if( is(fitObj,"lm") )
    rsq <- ifelse(adj,summary(fitObj)$adj.r.squared,summary(fitObj)$r.squared)
  else
    warning("Unsupported object!")

  if( adj )
    warning("No adjusted version of rsq.lr. Reset 'adj = FALSE' ...")
    
  rsq
}

# Nagelkerke (1991):  corrected generalization of 'lr'.
rsq.n<-function(fitObj,adj=FALSE,data=NULL)
{
  if( is.null(data) )
    data <- fitObj$model
  
  if( is(fitObj,"glm") )
  {
    n <- fitObj$df.null+1
    k <- fitObj$rank
    logLik1 <- as.numeric(logLik(fitObj))
    
    fitObj0 <- update(fitObj,.~1,data=data)
    logLik0 <- as.numeric(logLik(fitObj0))
    
    rsq <- (1-exp(-2*(logLik1-logLik0)/n))/(1-exp(logLik0*2/n))
  }
  else
    warning("Unsupported object!")
  
  if( adj )
    warning("No adjusted version of rsq.clr. Reset 'adj = FALSE' ...")

  rsq
}

simglm<-function(family=c("binomial", "gaussian", "poisson","Gamma"),lambda=3,n=50,p=3) 
{
  family <- family[[1]]
  
  switch(family,
    binomial={
      beta <- matrix(rep(0,p+1),nrow=p+1,ncol=1)
      beta[2] <- lambda
      x <- matrix(rnorm(n*(p+1),mean=0,sd=1),nrow=n,ncol=p+1)
      x[,1] <- rep(1,n)
      x[,2] <- rbind(rep(1,n/2),rep(-1,n/2))
      mu <- 1/(1+exp(-x%*%beta))
      y <- rbinom(n,1,mu)
      list(yx=data.frame(y=y,x=x[,-1]),beta=beta)},
    gaussian={ 
      beta <- matrix(rep(0,p+1),nrow=p+1,ncol=1)
      beta[2] <- lambda
      x <- matrix(rnorm(n*(p+1),mean=0,sd=1),nrow=n,ncol=p+1)
      x[,1] <- rep(1,n)
      x[,2] <- rbind(rep(1,n/2),rep(-1,n/2))
      mu <- x%*%beta
      y <- rnorm(n,mean=mu,sd=1)
      list(yx=data.frame(y=y,x=x[,-1]),beta=beta)},
    poisson={ 
      beta <- matrix(rep(0,p+1),nrow=p+1,ncol=1)
      beta[2] <- lambda
      x <- matrix(rnorm(n*(p+1),mean=0,sd=1),nrow=n,ncol=p+1)
      x[,1] <- rep(1,n)
      x[,2] <- rbind(rep(1,n/2),rep(-1,n/2))
      mu <- exp(x%*%beta)
      y <- rpois(n,lambda=mu)
      list(yx=data.frame(y=y,x=x[,-1]),beta=beta)},
    Gamma={
      beta <- matrix(rep(0,p+1),nrow=p+1,ncol=1)
      beta[2] <- lambda
      x <- matrix(rnorm(n*(p+1),mean=0,sd=1),nrow=n,ncol=p+1)
      x[,1] <- rep(1,n)
      x[,2] <- rbind(rep(1,n/2),rep(-1,n/2))
      tscale <- 0.01/(lambda+1+x%*%beta)
      y <- rgamma(n,shape=100,scale=tscale)
      list(yx=data.frame(y=y,x=x[,-1]),beta=beta)})
    #negative.binomial={   # Need to specify a better model here!
    #  size <- 2
    #  beta <- matrix(rep(0,p+1),nrow=p+1,ncol=1)
    #  beta[1] <- log(0.05/(0.05+size))
    #  beta[2] <- log(lambda/(lambda+size))-log(0.05/(0.05+size))
    #  x <- matrix(rnorm(n*(p+1),mean=0,sd=1),nrow=n,ncol=p+1)
    #  x[,1] <- rep(1,n)
    #  x[,2] <- rbind(rep(0,n/2),rep(1,n/2))
    #  mu <- size/(exp(-x%*%beta)-1)
    #  y <- rnbinom(n,size=2,mu=mu)
    #  list(yx=data.frame(y=y,x=x[,-1]),beta=beta)})
}

