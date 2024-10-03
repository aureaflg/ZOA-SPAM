## Functions related to ZOABR distribution
library(GoFKernel)
dBRr<-function(y, mu=0.5, alpha=0.5, sigma=1, log = FALSE){
  if(any(sigma <= 0)) stop(paste("sigma must be positive","\n",""))
  #if(any(mu < 0)|any(mu > 1)) stop(paste("mu must be between 0 and 1","\n",""))
  if(any(alpha <= 0)|any(alpha >= 1)) stop(paste("alpha must be between 0 and 1","\n",""))
  if (any(y <= 0) | any(y >= 1)) stop(paste("y must be between 0 and 1", "\n",""))
  epi<- 1-sqrt(1-4*alpha*mu*(1-mu))
  delta<- (mu-0.5*epi)/(1-epi)
  fy<-epi+(1-epi)*dbeta(y,shape1=delta*sigma, shape2=(1-delta)*sigma,ncp=0, log = log)
  if (log==T) 
    fy <- exp(fy)
  fy
}

pBRr<- function(q, mu=0.5, alpha=0.5, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if(any(sigma <= 0)) stop(paste("sigma must be positive","\n",""))
  if(any(mu <= 0)|any(mu >= 1)) stop(paste("mu must be between 0 and 1","\n",""))
  if(any(alpha <= 0)|any(alpha >= 1)) stop(paste("alpha must be between 0 and 1","\n",""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("y must be between 0 and 1", "\n",""))
  pdf<-function(y,mu,alpha,sigma){
    epi<-(1-sqrt(1-4*alpha*mu*(1-mu)))
    delta<- (mu-epi/2)/(1-epi)
    pdf<-epi+(1-epi)*dbeta(y,shape1=delta*sigma, shape2=(1-delta)*sigma) 
    pdf
  }
  #cdf<-integrate(pdf,lower=0,upper=q,mu=mu,alpha=alpha,sigma=sigma)$value
  cdfun <- function(upper, mu, alpha, sigma) {
    int <- integrate(pdf, lower = 0, upper = upper, mu, alpha, sigma)
    int$value
  }
  Vcdf <- Vectorize(cdfun)
  cdf <- Vcdf(upper = q, mu = mu,alpha=alpha, sigma = sigma)
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  
  if(any(cdf>0.99999999)){
    cdf[which(cdf>0.99999999)]=0.99999999
  }
  cdf
}



qBRr<- function(p, mu=0.5, alpha=0.5, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if(any(sigma <= 0)) stop(paste("sigma must be positive","\n",""))
  if(any(mu <= 0)|any(mu >= 1)) stop(paste("mu must be between 0 and 1","\n",""))
  if(any(alpha <= 0)|any(alpha >= 1)) stop(paste("alpha must be between 0 and 1","\n",""))
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  lp <- max(length(p), length(alpha), length(mu), length(sigma))
  p <- rep(p, length = lp)
  sigma <- rep(sigma, length = lp)
  mu <- rep(mu, length = lp)
  alpha <- rep(alpha, length = lp)
  q <- rep(0, lp)
  for(i in lp){
    h1 <- function(x) pBRr(x, mu[i], alpha[i], sigma[i]) 
    qq <- inverse(h1, lower=.Machine$double.eps, upper=1-.Machine$double.eps)
    q[i] <- qq(p[i])
  }
  q
}

rBRr<-function(n,mu,alpha,phi){
  epi<- 1-sqrt(1-4*alpha*mu*(1-mu))
  delta<- (mu-0.5*epi)/(1-epi)
  aux_unif<-runif(n)
  sim_unif<-mapply(runif,1,rep(0,n),rep(1,n))
  sim_beta<-mapply(rbeta,1,delta*phi,(1-delta)*phi)
  y<-ifelse(aux_unif<epi,sim_unif,sim_beta)
  while(any(y>0.9999999)){
    y[which(y>0.9999999)]=rbeta(length(which(y>0.9999999)),
                                delta[which(y>0.9999999)]*phi[which(y>0.9999999)],
                                (1-delta[which(y>0.9999999)])*phi[which(y>0.9999999)])
  }
  return(y)
}


## Function to create P-splines similar the pb function of GAMLSS
pbfake <- function(x, df = NULL, max.df = NULL, 
               inter = 20, degree= 3, order = 2,quantiles=F) 
{

bbase <- function(x, xl, xr, ndx, deg, quantiles=FALSE)
{
  tpower <- function(x, t, p)
    # Truncated p-th power function
    (x - t) ^ p * (x > t)
  # DS xl= min, xr=max, ndx= number of points within 
  # Construct B-spline basis
  # if quantiles=TRUE use different bases
  dx <- (xr - xl) / ndx # DS increment 
  if (quantiles) # if true use splineDesign
  {  #  this is not working and should be taken out
    knots <-  sort(c(seq(xl-deg*dx, xl, dx),quantile(x, prob=seq(0, 1, length=ndx)), seq(xr, xr+deg*dx, dx))) 
    B <- splineDesign(knots, x = x, outer.ok = TRUE, ord=deg+1)
    n <- dim(B)[2]  
    attr(B, "knots") <- knots[-c(1:(deg-1), (n-(deg-2)):n)]
    return(B)    
  }
  else # if false use Paul's
  { 
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, knots, tpower, deg)# calculate the power in the knots
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg) # 
    B <- (-1) ^ (deg + 1) * P %*% t(D)
    attr(B, "knots") <- knots[-c(1:(deg-1), (n-(deg-2)):n)]
    B 
  }
}

no.dist.val <-  length(table(x))
#if (is.matrix(x)) stop("x is a matric declare it as a vector")
#lx <- length(x)
#control$inter <- if (lx<99) 10 else control$inter # this is to prevent singularities when length(x) is small:change to 99 30-11-11 MS
inter <- if (no.dist.val<=inter)  no.dist.val else inter 
xl <- min(x)
xr <- max(x)
xmax <- xr + 0.01 * (xr - xl)
xmin <- xl - 0.01 * (xr - xl)  
##                create the basis
X <- bbase(x, xmin, xmax, inter, degree, quantiles) # 
#        getting the base out    assign("X", X, envir = .GlobalEnv) 
r <- ncol(X)
##                the penalty matrix
D <- if(order==0) diag(r) else diff(diag(r), diff=order)
## ------      if df are set                
if(!is.null(df)) # degrees of freedom
{
  if (df>(dim(X)[2]-2)) 
  {df <- 3}  
  #warning("The df's exceed the number of columns of the design matrix", "\n",  "   they are set to 3") }
  if (df < 0)  #warning("the extra df's are set to 0")   
  df <- if (df < 0)  2  else  df+2
}
## -------- check max.df   (new 7-2018 MS)  
if (is.null(max.df)) max.df <- dim(X)[2]
if (max.df>(dim(X)[2])) 
{
  max.df <- dim(X)[2]
  warning("The max.df's are set to",  dim(X)[2],  "\n")
}   
xvar <- x
attr(xvar, "D") <- D
attr(xvar, "X") <- X
attr(xvar, "df")<- df
xvar
}


## Function to estimate the smooth parameter
  gmu1<-function(mu,lig){
  x=lig
  if(x=="logit"){
    dmu=(1/mu+1/(1-mu))^(-1)
    d2mu2=-1/mu^2+1/(1-mu)^2
  }
  if(x=="cloglog"){
    dmu=(1/(log(1-mu)*(mu-1)))^(-1)
    d2mu2=-(log(1-mu)+1)/(log(1-mu)*(mu-1))^2
  }
  if(x=="loglog"){
    dmu=(-1/(mu*log(mu)))^(-1)
    d2mu2=(log(mu)+1)/(mu^2*log(mu)^2)
  }
  if(x=="cauchit"){
    dmu=(pi/(cos(pi*(mu-1/2)))^2)^(-1)
    d2mu2=(2*pi^2*sin(pi*(mu-1/2)))/(cos(pi*(mu-1/2)))^3
  }
  if(x=="probit"){
    dmu=(1/dnorm(qnorm(mu)))^(-1)
    d2mu2=-dnorm(qnorm(mu))*(-qnorm(mu))/(dnorm(qnorm(mu)))^3
  }
  ddmu<-list("dmu"=dmu,"d2mu2"=d2mu2)
}

  estlambdainf=function(lambda,y,zast,X,N,Ds,D,Bs,B,qj,bet,gamm,kapp,alpha,
                   method,k,c,statuslam,lig,lig.phi,s0,s1){
  s<-length(kapp)
  p<-length(bet)
  ntol<- nrow(B)
  eta1<-X%*%bet+B%*%gamm
  eta2<-N%*%kapp
  
  if(lig.phi=='log'){
    phi<-exp(eta2)
  }else if(lig.phi=='sqrt'){
    phi<-eta2^2
  } else if(lig.phi=='identity'){
    phi<-eta2
  }  else {print('Link function not defined')
  } 
  
  if(lig=="logit"){
    mu<-exp(eta1)/(1+exp(eta1)) 
  }
  if(lig=="cauchit"){
    mu<-(1/pi)*(atan(eta1)+0.5*pi)
  }
  if(lig=="cloglog"){
    mu<-1-exp(-exp(eta1))
  }
  if(lig=="loglog"){
    mu<-exp(-exp(-eta1))
  }
  if(lig=="probit"){
    mu<-pnorm(eta1)
  }
  epi<-1-sqrt(1-4*alpha*mu*(1-mu))
  delta<-(mu-0.5+0.5*(1-epi))/(1-epi)
  a<-delta*phi
  b<-(1-delta)*phi
  v<-epi^(1-zast)/((epi+(1-epi)*dbeta(y,a,b))^(1-zast))
  logy=ifelse(y>0 & y<1,log(y),0)
  logyc=ifelse(y>0 & y<1,log(1-y),0)
  
  f1<-(2*alpha*(1-2*mu))/(1-epi)
  f2<-(v/epi-(1-v)/(1-epi))
  ff<-(1-zast)*(f1*f2+(1-v)*(1-alpha)*(phi/((1-epi)^3))*(logy-
                 logyc-digamma(a)+digamma(b)))
  f3<-(4*(alpha^2)*(1-2*mu)^2)/((1-epi)^2)
  f4<-(-v/(epi^2)-(1-v)/((1-epi)^2))
  f5<-(-(4*alpha)/(1-epi)+(4*alpha^2*(1-2*mu)^2)/((1-epi)^3))
  f6<-(((1-alpha)*phi^2)/((1-epi)^6))*(trigamma(a)+trigamma(b))
  f7<-((6*alpha*phi*(1-2*mu))/((1-epi)^5))*(logy-logyc-
                                             digamma(a)+digamma(b)) 
  Rant<-diag(as.vector((f3*f4+f2*f5+(1-alpha)*(1-v)*(-f6+f7))*gmu1(mu=mu,lig=lig)$dmu^2-
                     ff*gmu1(mu=mu,lig=lig)$dmu^3*gmu1(mu=mu,lig=lig)$d2mu2),ntol)
  R<-diag(as.vector(1-zast),ntol)%*%Rant
  T1<-diag(as.vector(gmu1(mu=mu,lig=lig)$dmu),ntol)
  z<-X%*%bet+B%*%gamm+ginv(R)%*%T1%*%ff
  fnGAIC <- function(lambda, k , c , B , D , gamm)
  {
    Dlam<-lambda*D
    zast<-ifelse(y==0 | y==1, 1,0)
    dens<-ifelse(y>0 & y<1, dbeta(y,a,b,log=TRUE),0)
    vero <-sum((1-zast)*(v*log(epi)+(1-v)*log(1-epi)+(1-v)*dens))-
      0.5*t(gamm)%*%Dlam%*%gamm
    df <- sum(diag(B%*%ginv(t(B)%*%R%*%B-0.5*lambda*D)%*%t(B)%*%R))   
    GAIC <- -2*vero+k*(1+s0+s1+s+p+df)
    
    if(c==T){
      kstar=1+s0+s1+s+p+df
      GAIC <- GAIC+2*kstar*(kstar+1)/(ntol-kstar-1)
    }else GAIC <- GAIC
    
    GAIC   
  }
  
  fnGCV <- function(lambda, k , B , D , gamm)
  {
    df<-sum(diag(B%*%solve(t(B)%*%R%*%B+0.5*lambda*D)%*%t(B)%*%R))
    GCV <- sum(diag(R)%*%(z-X%*%bet)^2/(1-ntol^(-1)*df)^2)
    GCV
  }  
  
  if(statuslam==1){
    
    if(method=="GAIC"){
      edf=rep(0,length(qj))
      lambdaf=rep(0,sum(qj))
      for(au in 1:length(qj)){
        if(au==1){
          lambda0=lambda[1]
          aux <- nlminb(lambda0, fnGAIC,  lower = 10, upper = 1.0e7, 
                        k=k,c=c,B=Bs[[au]],D=Ds[[au]],gamm=gamm[1:qj[au]])$par
          lambdaf[1:qj[au]]<-rep(aux,qj[au])
          edf[au]<- sum(diag(Bs[[au]]%*%solve(t(Bs[[au]])%*%R%*%Bs[[au]]-
                                                0.5*aux*Ds[[au]])%*%t(Bs[[au]])%*%R))
        }else{lambda0=lambda[(sum(qj[1:(au-1)])+1)]
        aux <- nlminb(lambda0, fnGAIC,  lower = 1.0e-7, upper = 1.0e7,
                      k=k,c=c,B=Bs[[au]],D=Ds[[au]],
                      gamm=gamm[(sum(qj[1:(au-1)])+1):(sum(qj[1:(au-1)])+qj[au])])$par
        lambdaf[(sum(qj[1:(au-1)])+1):(sum(qj[1:(au-1)])+qj[au])]<-rep(aux,qj[au])
        edf[au]<- sum(diag(Bs[[au]]%*%solve(t(Bs[[au]])%*%R%*%Bs[[au]]-
                                              0.5*aux*Ds[[au]])%*%t(Bs[[au]])%*%R))
        }
        
      }
    }
    
    if(method=="GCV"){
      edf=rep(0,length(qj))
      lambdaf=rep(0,sum(qj))
      for(au in 1:length(qj)){
        if(au==1){
          lambda0=lambda[1]
          aux <- nlminb(lambda0, fnGCV,  lower = 10, upper = 1.0e7, 
                        k=k,B=Bs[[au]],D=Ds[[au]],gamm=gamm[1:qj[au]])$par
          lambdaf[1:qj[au]]<-rep(aux,qj[au])
          edf[au]<- sum(diag(Bs[[au]]%*%solve(t(Bs[[au]])%*%R%*%Bs[[au]]-
                                                0.5*aux*Ds[[au]])%*%t(Bs[[au]])%*%R))
        }else{lambda0=lambda[(sum(qj[1:(au-1)])+1)]
        aux <- nlminb(lambda0, fnGCV,  lower = 1.0e-7, upper = 1.0e7,
                      k=k,B=Bs[[au]],D=Ds[[au]],
                      gamm=gamm[(sum(qj[1:(au-1)])+1):(sum(qj[1:(au-1)])+qj[au])])$par
        lambdaf[(sum(qj[1:(au-1)])+1):(sum(qj[1:(au-1)])+qj[au])]<-rep(aux,qj[au])
        edf[au]<- sum(diag(Bs[[au]]%*%solve(t(Bs[[au]])%*%R%*%Bs[[au]]-
                                              0.5*aux*Ds[[au]])%*%t(Bs[[au]])%*%R))
        }
      }
      
    }
    
  }else{
    edf=rep(0,length(qj))
    lambdaf=lambda
    for(au in 1:length(qj)){
      if(au==1){
        lambda0=lambda[1]
        edf[au]<- sum(diag(Bs[[au]]%*%ginv(t(Bs[[au]])%*%R%*%Bs[[au]]-
                                             0.5*lambda0*Ds[[au]])%*%t(Bs[[au]])%*%R))
      }else{
        lambda0=lambda[(sum(qj[1:(au-1)])+1)]
        edf[au]<- sum(diag(Bs[[au]]%*%ginv(t(Bs[[au]])%*%R%*%Bs[[au]]-
                                             0.5*lambda0*Ds[[au]])%*%t(Bs[[au]])%*%R))
      }
    }
  }
  aic<-fnGAIC(lambda=lambda, k=2 , c=c , B=B , D=D , gamm=gamm)
  bic<-fnGAIC(lambda, k=log(ntol) , c , B , D , gamm)
  aicc <- fnGAIC(lambda, k=2 , c=TRUE , B , D , gamm)
  sabic <- fnGAIC(lambda, k=log(ntol+2)/24 , c , B , D , gamm)
  hqic <- fnGAIC(lambda, k=2*log(log(ntol)) , c , B , D , gamm)
  return(list(lambda=lambdaf,edf=edf,aic=aic,bic=bic,aicc=aicc,
              sabic=sabic,hqic=hqic))
}

