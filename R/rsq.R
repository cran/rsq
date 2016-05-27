#fitObj: object from lm() or glm() or glm.nb();
#adj: whether to calculate the adjusted r-squared.
rsq<-function(fitObj,adj=FALSE)
{
  if( is(fitObj,"glm") )
  {
    y <- model.response(fitObj$model)
    yfit1 <- fitObj$fitted.values
    p <- fitObj$rank
    df.int <- if(attr(fitObj$terms,"intercept")) 1L else 0L
               
    tfit0 <- update(fitObj,.~1)
    yfit0 <- tfit0$fitted.values
 
    vrsq(y,yfit0,yfit1,family=family(fitObj)$family,adj=adj,df.int,p)
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
}


# Calculate the variance-based distance for quardratic variance functions:
#           V(mu) = v2*mu^2+v1*mu+v0
calcqvdist<-function(a,b,v2,v1)
{
  vpa <- 2*v2*a+v1
  svpa2 <- sqrt(1+vpa*vpa)
  
  vpb <- 2*v2*b+v1
  svpb2 <- sqrt(1+vpb*vpb)
  
  dv <- (log((vpb+svpb2)/(vpa+svpa2))+vpb*svpb2-vpa*svpa2)/(4*v2)
  dv <- dv*dv
}

# Calculate variance-based R^2 using observed & fitted response values
vrsq<-function(y,yfit0,yfit1,family="binomial",adj=FALSE,df.int=1,p=1)
{
  switch(family,
     binomial=         # V(mu) = mu*(1-mu)
     {
       v2 <- -1
       v1 <- 1
       if(ncol(as.matrix(y))>1)
       {
         rsq <- (calcqvdist(1,yfit1,v2,v1)%*%y[,1]+calcqvdist(0,yfit1,v2,v1)%*%y[,2])/
                (calcqvdist(1,yfit0,v2,v1)%*%y[,1]+calcqvdist(0,yfit0,v2,v1)%*%y[,2])
         n <- sum(y)
       } 
       else
       {
         rsq <- sum(calcqvdist(y,yfit1,v2,v1))/sum(calcqvdist(y,yfit0,v2,v1))
         n <- length(y)
       }
     }, 
     gaussian=         # V(mu) = 1
     {
       rsq <- sum((y-yfit1)^2)/sum((y-yfit0)^2)
       n <- length(y)
     },
     Gamma=            # V(mu) = mu*mu
     {
       v2 <- 1
       v1 <- 0
       rsq <- sum(calcqvdist(y,yfit1,v2,v1))/sum(calcqvdist(y,yfit0,v2,v1))
       n <- length(y)
     },
     inverse.gaussian=     # V(mu) = mu*mu*mu
     {
       tint <- function(x){integrate(function(mu){mu*mu*mu},x[1],x[2])$value}
       nv <- cbind(y,yfit1)
       dv <- cbind(y,yfit0)
       rsq <- sum(apply(nv,1,tint))/sum(apply(dv,1,tint))
       n <- length(y)
     },
     poisson=          # V(mu) = mu
     {
       rsq <- sum((y-yfit1)^2)/sum((y-yfit0)^2)
       n <- length(y)
     },
     quasibinomial=    # V(mu) = phi*mu*(1-mu)
     {
       v2 <- -1
       v1 <- 1
       rsq <- sum(calcqvdist(y,yfit1,v2,v1))/sum(calcqvdist(y,yfit0,v2,v1))
       n <- length(y)
     },
     quasipoisson=     # V(mu) = phi*mu
     {
       rsq <- sum((y-yfit1)^2)/sum((y-yfit0)^2)
       n <- length(y)
     },
     quasi=            # V(mu) specified by variance including varfun?
     {
       warning("Not Implemented!")
     }
  )
  
  if(pmatch("Negative Binomial",family,nomatch=F)) 
  {
    # V(mu) = mu+mu*mu
    v2 <- 1
    v1 <- 1
    rsq <- sum(calcqvdist(y,yfit1,v2,v1))/sum(calcqvdist(y,yfit0,v2,v1))
    n <- length(y)
  }
  
  if(adj)
    rsq <- 1-rsq*(n-df.int)/(n-p)
  else
    rsq <- 1-rsq
  
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
