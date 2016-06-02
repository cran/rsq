#fitObj: object from lm() or glm() or glm.nb();
#adj: whether to calculate the adjusted r-squared.
rsq<-function(fitObj,adj=FALSE)
{
  if( is(fitObj,"glm") )
  {
    sse1 <- vsse(fitObj)
               
    fitObj0 <- update(fitObj,.~1)
    sse0 <- vsse(fitObj0)
    
    if(adj)
      rsq <- 1-(sse1/summary(fitObj)$df.residual)/(sse0/summary(fitObj0)$df.residual)
    else
      rsq <- 1-sse1/sse0
  }
  else if( is(fitObj,"lm") )
  {
    if( adj )
      summary(fitObj)$adj.r.squared
    else
      summary(fitObj)$r.squared
  }
  else
    warning("Unsupported object!")

  rsq
}

# objF: object from fitting the full model;
# objR: object from fitting the reduced model.
rsq.partial<-function(objF,objR=NULL,adj=FALSE)
{
  rsqF <- rsq(objF,adj=adj)
  if(!is.null(objR))
  {
    rsqR <- rsq(objR,adj=adj)
    prsq <- (rsqF-rsqR)/(1-rsqR)
    list(adjustment=adj,variables.full=attr(terms(objF),"term.labels"),
         variables.reduced=attr(terms(objR),"term.labels"),partial.rsq=prsq)
  }
  else # Return the partial correlation coefficient for each term in the model
  {
    nterms <- length(attr(terms(objF),"term.labels"))
    prsq <- double(length=nterms)
    for(j in 1:nterms)
    {
      tfit <-update(objF,terms(objF)[-j])
      rsqR <- rsq(tfit,adj=adj)
      prsq[j] <- (rsqF-rsqR)/(1-rsqR)
    }
    list(adjustment=adj,variable=attr(terms(objF),"term.labels"), partial.rsq=prsq)
  }
}

# Calculate the sum of squared errors
vsse<-function(fitObj)
{
  yyfit <- cbind(model.response(fitObj$model),fitObj$fitted.values)

  if( (family(fitObj)$family=="binomial")&(ncol(as.matrix(yyfit))>2) )
  {
    rsse <- function(x){
      tres<-vresidual(x[1:2],x[-(1:2)],family=family(fitObj)$family)
      rsse<-sum((tres^2)*x)
    }
    vsse <- sum(apply(yyfit,1,rsse))
  }
  else if(pmatch("Negative Binomial",family(fitObj)$family,nomatch=F))
  {
    tresid <- function(x){vresidual.nb(x[1],x[2],fitObj$theta)}
    vsse <- sum(apply(yyfit,1,tresid)^2)
  }
  #else if( family(fitObj)$family=="quasi" )
  #{ # Note that we need to specify the first order derivative of the variance function
  #  tresid <- function(x){vresidual.quasi(x[1],x[2],variancep=family(fitObj)$variancep)}
  #  vsse <- sum(apply(yyfit,1,tresid)^2)
  #}
  else
  {
    tresid <- function(x){vresidual(x[1],x[2],family=family(fitObj)$family)}
    vsse <- sum(apply(yyfit,1,tresid)^2)
  }
  
  vsse
}

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
vresidual<-function(y,yfit,family="binomial")
{
  vresidual <- 0
  switch(family,
    binomial=         # V(mu) = mu*(1-mu)
    {
      v2 <- -1
      v1 <- 1
      if(ncol(as.matrix(y))>1)
        vresidual <- c(qvresidual(1,yfit,v2,v1), qvresidual(0,yfit,v2,v1))
      else
        vresidual <- qvresidual(y,yfit,v2,v1)
    }, 
    gaussian=         # V(mu) = 1
    {
      vresidual <- y-yfit
    },
    Gamma=            # V(mu) = mu*mu
    {
      v2 <- 1
      v1 <- 0
      vresidual <- qvresidual(y,yfit,v2,v1)
    },
    inverse.gaussian=     # V(mu) = mu*mu*mu
    {
      vresidual <- integrate(function(mu){sqrt(1+9*mu^4)},yfit,y)$value;
    },
    poisson=          # V(mu) = mu
    {
      vresidual <- y-yfit
    },
    quasibinomial=    # V(mu) = mu*(1-mu)
    {
      v2 <- -1
      v1 <- 1
      vresidual <- qvresidual(y,yfit,v2,v1)
    },
    quasipoisson=     # V(mu) = mu
    {
      vresidual <- y-yfit
    },
    quasi=            # V(mu) specified by variance including varfun?
    {
      warning("Not Implemented!")
    })

  vresidual
}

# Calculate the residual for given observed y and its fitted value yfit
# for the negative binomial model
vresidual.nb<-function(y,yfit,theta=1)
{
  # V(mu) = mu+mu*mu/theta
  v2 <- 1/theta
  v1 <- 1
  vresidual <- qvresidual(y,yfit,v2,v1)

  vresidual
}

# Calculate the residual for given observed y and its fitted value yfit for quasi models
#vresidual.quasi(y, yfit,variancep=function(mu){1},dispersion=1)
#{
#  vresidual <- integrate(function(mu){sqrt(1+variancep(mu)^2)},yfit,y)$value;
#  vresidual
#}
