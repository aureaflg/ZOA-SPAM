##This model was implemented in GAMLSS package
library(gamlss)
library(purrr)
library(betareg)
library(Matrix)
library(matrixcalc)
library(GoFKernel)
library(pracma)
library(cubature)
library(rootSolve)
library(aod)
library(ggplot2)
library(qqplotr)
library(ggpubr)
library(ggrepel)

####################################################################################
######################### Estimation process #######################################
####################################################################################


dBET<-function(y, mu=0.5, sigma=1, log = FALSE){
  if(any(0.99999<mu & mu<1.000001)){
    mu[which(0.99999<mu & mu<1.000001)]=0.9999
  }
  if (any(mu <= 0)) 
    stop(paste("mu must be between 0 and 1", "\n",""))
  if(any(sigma <= 0)) stop(paste("sigma must be positive","\n",""))
  if (any(y <= 0) | any(y >= 1)) stop(paste("y must be between 0 and 1", "\n",""))
  fy<-dbeta(y,shape1=mu*sigma, shape2=(1-mu)*sigma,ncp=0, log = log)
  fy
}

pBET<- function(q, mu=0.5, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if(any(sigma <= 0)) stop(paste("sigma must be positive","\n",""))
  if(any(mu <= 0)|any(mu >= 1)) stop(paste("mu must be between 0 and 1","\n",""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("y must be between 0 and 1", "\n",""))
  cdf<-pbeta(q,shape1=mu*sigma, shape2=(1-mu)*sigma,ncp=0, lower.tail = lower.tail, log.p = log.p)
  cdf
}

qBET<-function(p, mu=0.5, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if(any(sigma <= 0)) stop(paste("sigma must be positive","\n",""))
  if(any(mu <= 0)|any(mu >= 1)) stop(paste("mu must be between 0 and 1","\n",""))
  if(any(p <= 0)|any(p >= 1)) stop(paste("p must be between 0 and 1","\n",""))
  q<-qbeta(p,shape1=mu*sigma, shape2=(1-mu)*sigma,ncp=0, lower.tail = lower.tail, log.p = log.p)
  q
}

dBETINF<-function (x, mu = 0.5, sigma = 1, nu = 0.1, tau = 0.1, 
                   log = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma greated than 0", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must greated than 0", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must greated than 0", "\n", ""))
  if (any(x < 0) | any(x > 1)) 
    stop(paste("x must be 0<=x<=1, i.e. 0 to 1 inclusively", 
               "\n", ""))
  lp <- pmax.int(length(x), length(mu), length(sigma),
                 length(nu),length(tau))
  x <- rep(x, length = lp)
  sigma <- rep(sigma, length = lp)
  mu <- rep(mu, length = lp)
  nu <- rep(nu, length = lp)
  tau <- rep(tau, length = lp)
  logfy <- rep(0, length = lp)
  logfy <- rep(0, length(x))
  for(i in 1:lp){
    logfy[i] <- ifelse((x[i]>0 & x[i]<1), dBET(x[i], mu=mu[i], sigma=sigma[i], log = TRUE), 0)
    logfy[i] <- ifelse((x[i]==0), log(nu[i]), logfy[i])          
    logfy[i] <- ifelse((x[i]==1), log(tau[i]) , logfy[i])
    logfy[i] <- logfy[i] - log(1+nu[i]+tau[i])   
  }
  if (log == FALSE) 
    fy <- exp(logfy)
  else fy <- logfy
  fy
}

pBETINF<-function (q, mu = 0.5, sigma = 1, nu = 0.1, tau = 0.1, lower.tail = TRUE, 
         log.p = FALSE) {
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greated than 0", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must greated than 0", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must greated than 0", "\n", ""))
  if (any(q < 0) | any(q > 1)) 
    stop(paste("y must be 0<=y<=1, i.e. 0 to 1 inclusively", 
               "\n", ""))
  lp <- pmax.int(length(q), length(mu), length(sigma),
                 length(nu),length(tau))
  q <- rep(q, length = lp)
  sigma <- rep(sigma, length = lp)
  mu <- rep(mu, length = lp)
  nu <- rep(nu, length = lp)
  tau <- rep(tau, length = lp)
  p <- rep(0, length = lp)
  for(k in 1:lp){
    if(q[k]==0){p[k]<-nu[k]}
    if(q[k]==1){p[k]<-1+nu[k]+tau[k]}
    if(q[k]>0 & q[k]<1){
      p[k]<-nu[k] + pBET(q[k], mu = mu[k], 
                         sigma=sigma[k],lower.tail=TRUE,log.p=FALSE)
    }
    p[k] <- p[k]/(1+nu[k]+tau[k])
  }
  if (lower.tail == TRUE) p <- p
  else p = 1 - p
  if (log.p == FALSE) p <- p
  else p <- log(p)
  p
}

qBETINF<-function (p, mu = 0.5, sigma = 1, nu = 0.1, tau = 0.1, lower.tail = TRUE, 
          log.p = FALSE) {
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greated than 0", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must greated than 0", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must greated than 0", "\n", ""))
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  lp <- pmax.int(length(p), length(mu), length(sigma),
                 length(nu),length(tau))
  p <- rep(p, length = lp)
  sigma <- rep(sigma, length = lp)
  mu <- rep(mu, length = lp)
  nu <- rep(nu, length = lp)
  tau <- rep(tau, length = lp)
  q <- rep(0, length = lp)
  for(k in 1:lp){
    if(p[k]<=(nu[k]/(1+nu[k]+tau[k]))){
      q[k]<-0
      next
    }
    if(p[k]>=((1+nu[k])/(1+nu[k]+tau[k]))){
      q[k]<-1
      next
    }else{
      q[k]<-qBET((p[k]- (nu[k]/(1 + nu[k] + tau[k])))/
                   (1/(1 + nu[k] + tau[k])),
                 mu = mu[k],sigma = sigma[k],
                 lower.tail=TRUE, log.p=FALSE)
    }
  }
  q
}

rBETINF<-function (n, mu = 0.5, sigma = 0.1, nu = 0.1, tau = 0.1) 
{
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greated than 0", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must greated than 0", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must greated than 0", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qBETINF(p, mu = mu, sigma = sigma, nu = nu, tau = tau)
  r
}

BETINF<- function (mu.link = "logit", sigma.link = "log", nu.link = "log", 
                   tau.link = "log") 
{
  mstats <- checklink("mu.link", "BETINF", substitute(mu.link), 
                c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "BETINF", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink("nu.link", "BETINF", substitute(nu.link), 
                c("inverse", "log", "identity", "own"))
  tstats <- checklink("tau.link", "BETINF", substitute(tau.link), 
                c("inverse", "log", "identity", "own"))
  structure(list(family = c("BETINF", "Beta Inflated"),
                 parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE, 
                          tau = TRUE), nopar = 4, type = "Mixed", 
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 nu.link = as.character(substitute(nu.link)),
                 tau.link = as.character(substitute(tau.link)), 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 nu.linkfun = vstats$linkfun, 
                 tau.linkfun = tstats$linkfun, 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 nu.linkinv = vstats$linkinv, 
                 tau.linkinv = tstats$linkinv, 
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 nu.dr = vstats$mu.eta, 
                 tau.dr = tstats$mu.eta, 
                 dldm = function(y, mu, sigma) {
                   p <- mu * sigma
                   q <- (1 - mu) * sigma
                   mustart <- digamma(p) - digamma(q)
                   ystart <- log(y) - log(1 - y)
                   dldm <- ifelse(((y == 0) | (y == 1)), 0, sigma *
                                    (ystart - mustart))
                   dldm
                 }, d2ldm2 = function(y, mu, sigma) {
                   p <- mu * sigma
                   q <- (1 - mu) * sigma
                   d2ldm2 <- ifelse(((y == 0) | (y == 1)), 0, 
                                    -sigma^2 * (trigamma(p) + trigamma(q)))
                   d2ldm2
                 }, dldd = function(y, mu, sigma) {
                   p <- mu * sigma
                   q <- (1 - mu) * sigma
                   mustart <- digamma(p) - digamma(q)
                   ystart <- log(y) - log(1 - y)
                   dldd <- ifelse(((y == 0)| (y == 1)), 0, mu *
                                    (ystart - mustart) + 
                                    log(1 - y) - digamma(q) + digamma(sigma))
                   dldd
                 }, d2ldd2 = function(y, mu, sigma) {
                   p <- mu * sigma
                   q <- (1 - mu) * sigma
                   d2ldd2 <- ifelse(((y == 0)| (y == 1)), 0, 
                                    -(mu^2 * trigamma(p) + (1 - mu)^2 *
                                        trigamma(q) -  trigamma(sigma)))
                   d2ldd2
                 }, dldv = function(y, nu, tau) {
                   dldv <- ifelse(y == 0, (1/nu), 0) -
                     (1/(1 + nu +  tau))
                   dldv
                 }, d2ldv2 = function(nu, tau) {
                   d2ldv2 <- -(1 + tau)/(nu * 
                                        ((1 + nu + tau)^2))
                   d2ldv2
                 }, dldt = function(y, nu, tau) {
                   dldt <- ifelse(y == 1, (1/tau), 0) - 
                     (1/(1 + nu +  tau))
                   dldt
                 }, d2ldt2 = function(nu, tau) {
                   d2ldt2 <- -(1 + nu)/(tau * ((1 + nu + tau)^2))
                   d2ldt2
                 }, d2ldmdd = function(y, mu, sigma) {
                   p <- mu * sigma
                   q <- (1 - mu) * sigma
                   d2ldmdd <- ifelse(((y == 0) | (y == 1)), 0, 
                            -sigma * (trigamma(p) * mu) - 
                              (trigamma(q) * (1 - mu)))
                   d2ldmdd
                 }, d2ldmdv = function(y) {
                   d2ldmdv <- rep(0, length(y))
                   d2ldmdv
                 }, d2ldmdt = function(y) {
                   d2ldmdt <- rep(0, length(y))
                   d2ldmdt
                 }, d2ldddv = function(y) {
                   d2ldddv <- rep(0, length(y))
                   d2ldddv
                 }, d2ldddt = function(y) {
                   d2ldddt <- rep(0, length(y))
                   d2ldddt
                 }, d2ldvdt = function(nu, tau) {
                   d2ldvdt <- 1/((1 + nu + tau)^2)
                   d2ldvdt
                 }, G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2 * 
                   dBETINF(y, mu, sigma, nu, tau, log = TRUE), 
                 rqres = expression(rqres(pfun = "pBETINF", type = "Mixed", 
                        mass.p = c(0, 1), prob.mp = cbind(nu/(1 + nu + tau), 
                        tau/(1 + nu + tau)), y = y, mu = mu,
                        sigma = sigma, nu = nu, tau = tau)), 
                 mu.initial = expression(mu <- (y + mean(y))/2), 
                 sigma.initial = expression(sigma <- rep(0.5,  length(y))),
                 nu.initial = expression(nu <- rep(0.3, length(y))), 
                 tau.initial = expression(tau <- rep(0.3, length(y))), 
                 mu.valid = function(mu) all(mu > 0 &  mu < 1),
                 sigma.valid = function(sigma) all(sigma >  0), 
                 nu.valid = function(nu) all(nu > 0), 
                 tau.valid = function(tau) all(tau > 0), 
                 y.valid = function(y) all(y >= 0 & y <= 1)), 
            class = c("gamlss.family", "family"))
}

loglog  <- function()
{
  linkfun <- function(mu) { -log(-log(mu))} 
  linkinv <- function(eta) { 
    thresh <- log(-log(.Machine$double.eps))
    eta <- pmin(thresh, pmax(eta, -thresh))
    exp(-exp(-eta))}
  mu.eta <- function(eta) pmax(exp(-eta) * exp(-exp(-eta)), 
                               .Machine$double.eps)
  valideta <- function(eta) TRUE
  link <- "loglog"
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
                 valideta = valideta, name = link), class = "link-gamlss")
}

## Example on how to use the GAMLSS with our distribution
## In mod1, y is the response variable, x is a covariate (continuous or discrete) 
## in the parametric part of the model, z is a continuous covariate which relation   
## with the response is unknown and will approximated by a P-spline. 
## In GAMLSS pb() is the function for P-spline, inter is the number of knots we want 
## to settle and lambda is the smooth parameter. 
## If the inter and lambda parameters are not set by the user, they will be estimated 
## by the function. sigma.formula is the linear predictor related to precision parameter (phi) 
## of the ZOAB distribution, nu.formula and tau.formula  are the linear predictors related 
## to probability of occurence of zero and one, respectively.
##  Inside the family distribution, we set the link functions, here I set just for mu and phi.           
# mod1=gamlss(y~1+x1+pb(z,control=pb.control(inter=50),lambda = 100),
#            sigma.formula =~1+x2,nu.formula =~1+z0,tau.formula =~1+z1, 
#             family = BETINF(mu.link="logit",sigma.link = "log"))

## If the link function is not implemented in GAMLSS, as the loglog link, we can implement
## and use as in mod2 example.            
# mod2=gamlss(y~1+x+pb(z,control=pb.control(inter=50),lambda = 100),
#            sigma.formula =~1+e,nu.formula =~1+fl,tau.formula =~1+m, 
#             family = BETINF(mu.link=loglog(),sigma.link = "log"))            





            
####################################################################################
########################### Diagnostic tools #######################################
####################################################################################
            
                    
## Residuals
RQR.ZOAB=function(mod){
  y=mod$y
  tau=mod$tau.coefficients
  rho=mod$nu.coefficients
  Z0=mod$nu.x
  Z1=mod$tau.x
  predf<-Z0%*%rho
  predm<-Z1%*%tau
  p0=exp(predf)/(1+exp(predf)+exp(predm))
  p1=exp(predm)/(1+exp(predf)+exp(predm))
  mu=mod$mu.fv
  phi=mod$sigma.fv
  n0<- length(which(y==0))
  ut0<- runif(n0,0,p0[y==0])
  rQ0<- qnorm(ut0)
  n1<- length(which(y==1))
  ut1<- runif(n1,1-p1[y==1],1)
  rQ1<- qnorm(ut1)
  ut<- pBETINF(y, mu, phi, nu = p0/(1-p0-p1), tau = p1/(1-p0-p1))
  rQ<- qnorm(ut)
  resiquan<-numeric()
  resiquan[y==0]<- rQ0
  resiquan[y==1]<- rQ1
  resiquan[which(y>0 & y<1)]<- rQ[which(y>0 & y<1)]
  return(resiquan)
}

plotres=function(mod){
  res=RQR.ZOAB(mod)
  smp <- data.frame(norm = res)
  s1<-ggplot(data = smp, mapping = aes(sample = norm))  +
    stat_qq_band(conf = 0.95) +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Quantile Residuals") +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "grey"),
      text=element_text(size=20,family="serif")
    )
  s2<-ggplot()+
    geom_point(aes(x=seq_along(res),y=res))+ 
    xlab("Index") + ylab("Quantile Residuals") +
    geom_hline(yintercept=-2,linetype="dashed")+ 
    geom_hline(yintercept=2,linetype="dashed")+ 
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "grey"),
      text=element_text(size=20,family="serif")
    )
  datresi=data.frame(res=res)
  s3<-ggplot(datresi, aes(x = res, y = after_stat(density))) + 
    geom_histogram(fill = "grey", color = "black",bins = 10) +
    xlab("Quantile Residuals") + ylab("density") +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "grey"),
      text=element_text(size=20,family="serif")
    )
  
  print(ggarrange(s1,s2,s3))
}
            



## Standard error discrete part
stderror.disc=function(mod,conf=0.95){
  tau=mod$tau
  rho=mod$rho
  Z0=mod$nu.x
  Z1=mod$tau.x
  k0=length(rho)
  k1=length(tau)
  n=mod$ntol
  
  predf<-Z0%*%rho
  predm<-Z1%*%tau
  p0=exp(predf)/(1+exp(predf)+exp(predm))
  p1=exp(predm)/(1+exp(predf)+exp(predm))
  zast<-ifelse(mod$y==0 | mod$y==1, 1,0)
  Zast<-diag(c(zast))
  
  IFrho<- t(Z0)%*%diag(as.vector(p0*(1-p0)),n)%*%Z0
  IFtau<- t(Z1)%*%diag(as.vector(p1*(1-p1)),n)%*%Z1
  IFpr<- t(Z0)%*%diag(as.vector(-p0*p1),n)%*%Z1
  IF<-matrix(rbind(cbind(IFrho, IFpr),cbind(t(IFpr), IFtau)),nrow=k0+k1,ncol=k0+k1)
  
  coefdisc=c(rho,tau)
  IF1=solve(IF)
  EPdisc=sqrt(diag(IF1))
  VARdisc=diag(IF1)
  valorpdisc = 0
  walddisc = 0
  for(i in 1:length(coefdisc)){
    walddisc[i] = wald.test(VARdisc[i], coefdisc[i], Terms = 1)$result$chi2[1]
    valorpdisc[i] = wald.test(VARdisc[i], coefdisc[i], Terms = 1)$result$chi2[3]
  }
  IC<-function(nconf,param, EP){
    lower<- c(param)-qnorm(1-nconf/2)*EP
    upper<- c(param)+qnorm(1-nconf/2)*EP
    obj.out <- list(IC=cbind(cbind(lower,upper)))
    return(obj.out)
  }
  
  interval=rbind(IC(1-conf, rho, EPdisc[1:k0])$IC,
                 IC(1-conf, tau, EPdisc[(k0+1):(k0+k1)])$IC)
  
  coefficientsd = data.frame(coefdisc, EPdisc, round(valorpdisc,digits = 4),interval)
  colnames(coefficientsd) = c("Estimate","Std.err", "Pr(>|W|)","IC-lower","IC-upper")
  return(printCoefmat(as.matrix(coefficientsd), digits = 4))
}


gmu=function(mod,alpha=NULL){
  x = mod$mu.link
  mu=mod$mu.fv
  
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
gphi=function(mod){
  x=mod$sigma.link
  phi=mod$sigma.fv
  if(x=="log"){
    dphi=(1/phi)^(-1)
    d2phi2=-1/phi^2
  }
  if(x=="1/x^2"){
    dphi=(-(2/phi^3))^(-1)
    d2phi2=6/phi^4
  }
  if(x=="sqrt"){
    dphi=((1/2)*phi^(-1/2))^(-1)
    d2phi2=(-1/4)*phi^(-3/2)
  }
  ddphi<-list("dphi"=dphi,"d2phi2"=d2phi2)
}            


## Standard error continuous part            
stderror.cont=function(mod,conf=0.95){
degree=3
order=2
n=mod$N
s=ncol(mod$sigma.x)
mu=mod$mu.fv
phi=mod$sigma.fv
zast<-ifelse(mod$y==0 | mod$y==1, 1,0)

if(is.null(mod$mu.s)==FALSE){
  p=ncol(mod$mu.x)-ncol(mod$mu.s)  
  nq=ncol(mod$mu.s)
  qj=0
  for(i in 1:nq){
    qj[i]=length(mod$mu.coefSmo[[i]][11]$knots)
  }
  X=as.matrix(mod$mu.x[,-which(colnames(mod$mu.x) %in% colnames(mod$mu.s))])
  bet=matrix(mod$mu.coefficients[-which(colnames(mod$mu.x) %in% colnames(mod$mu.s))], nrow=p)  
  gamms=list(NULL)
  for(i in 1:nq){
    gamms[[i]]=mod$mu.coefSmo[[i]][1]$coef
  }
  gamm=Reduce("rbind",gamms)
  Bs<-list(NULL)
  Ds<-list(NULL)
  for (i in 1:nq){
    ajus=pb(mod$mu.x[,which(colnames(mod$mu.x) %in% colnames(mod$mu.s)[i])],
            degree=degree,order=order,
            control=pb.control(inter=(length(mod$mu.coefSmo[[1]]$knots)-3)),
            lambda = mod$mu.lambda[i])
    Bs[[i]]=attr(ajus,"X")
    Ds[[i]]=t(attr(ajus,"D"))%*%attr(ajus,"D")*mod$mu.lambda[i]
  }
  
  B=Reduce("cbind",Bs)
  D=0.5*bdiag(Ds)
  
} else{
  p=ncol(mod$mu.x)
  bet=matrix(mod$mu.coefficients,nrow=p)
  X=mod$mu.x}

kapp=matrix(mod$sigma.coefficients,nrow=s)
N=mod$sigma.x
T1=diag(gmu(mod)$dmu,n)
T2=diag(gphi(mod)$dphi,n)

W1=(phi^2*(trigamma(mu*phi)+trigamma((1-mu)*phi)))*gmu(mod)$dmu^2
W=diag(as.numeric((1-zast)*W1),n)
Pstar1=(trigamma(mu*phi)*mu^2+trigamma((1-mu)*phi)*(1-mu)^2-trigamma(phi))*gphi(mod)$dphi^2
Pstar=diag(as.numeric((1-zast)*Pstar1),n)
W=diag(as.numeric((1-zast)*Pstar1),n)
Wstar1=(phi*trigamma(mu*phi)*mu-phi*trigamma((1-mu)*phi)*
          (1-mu))*gmu(mod)$dmu*gphi(mod)$dphi
Wstar=diag(as.numeric((1-zast)*Wstar1),n)


Kbb= t(X)%*%W%*%X
Kkk= t(N)%*%Pstar%*%N
Kgg= t(B)%*%W%*%B+D
Kbk= t(X)%*%Wstar%*%N
Kgk= t(N)%*%Wstar%*%B
Kbg= t(X)%*%W%*%B
Kggl=t(B)%*%W%*%B
Ktt=rbind(cbind(Kbb,Kbk,Kbg),cbind(t(Kbk),Kkk,Kgk),cbind(t(Kbg),t(Kgk),Kgg))
Kttinv=ginv(as.matrix(Ktt))
coef=c(bet,kapp)
EP=sqrt(diag(Kttinv)[1:length(coef)])
VAR=diag(Kttinv)[1:length(coef)]
valorp = 0
wald = 0
for(i in 1:length(coef)){
  wald[i] = wald.test(VAR[i], coef[i], Terms = 1)$result$chi2[1]
  valorp[i] = wald.test(VAR[i], coef[i], Terms = 1)$result$chi2[3]
}

IC<-function(nconf,param, EP){
  lower<- c(param)-qnorm(1-nconf/2)*EP
  upper<- c(param)+qnorm(1-nconf/2)*EP
  obj.out <- list(IC=cbind(cbind(lower,upper)))
  return(obj.out)
}

interval=rbind(IC(1-conf, bet, EP[1:p])$IC,
               IC(1-conf, kapp, EP[(p+1):(s+p)])$IC)


coefficients = data.frame(coef, EP,  round(valorp,4), interval)
colnames(coefficients) = c("Estimate","Std.err", "Pr(>|W|)","IC-lower","IC-upper")
return(printCoefmat(as.matrix(coefficients), digits = 4))
}



#Local influence discrete part
#Case-weight perturbation
# if you want to highlight a point, use the argument ref to do so
localdisc=function(mod,ref=0.5){
tau=mod$tau.coefficients
rho=mod$nu.coefficients
Z0=mod$nu.x
Z1=mod$tau.x
k0=length(rho)
k1=length(tau)

predf<-Z0%*%rho
predm<-Z1%*%tau
p0=exp(predf)/(1+exp(predf)+exp(predm))
p1=exp(predm)/(1+exp(predf)+exp(predm))
zast<-ifelse(mod$y==0 | mod$y==1, 1,0)
Zast<-diag(c(zast))

Urho<- t(Z0)%*%(zast*(1-mod$y)-p0)
Utau<- t(Z1)%*%(zast*mod$y-p1)
Uesc<-matrix(c(Urho,Utau), ncol=1)

Urr<- -t(Z0)%*%diag(as.vector(p0*(1-p0)),n)%*%Z0
Uttau<- -t(Z1)%*%diag(as.vector(p1*(1-p1)),n)%*%Z1
Urtau<- -t(Z0)%*%diag(as.vector(-p0*p1),n)%*%Z1
Udisc<-matrix(rbind(cbind(Urr,Urtau ),cbind(t(Urtau), Uttau)),nrow=k0+k1,ncol=k0+k1)
Udiscinv=ginv(as.matrix(-Udisc))

deltarhoc=t(Z0)%*%diag(as.numeric(zast*(1-mod$y)-p0),n)
deltatauc=t(Z1)%*%diag(as.numeric(zast*mod$y-p1),n)
deltacasosd=rbind(deltarhoc,deltatauc)
Zcdi=t(deltacasosd)%*%Udiscinv%*%deltacasosd
Cmaxcdi=abs(eigen(Zcdi)$vectors[,1])
s1<-qplot(seq_along(Cmaxcdi),Cmaxcdi,geom = "point",label = seq(1,length(Cmaxcdi),1))+ 
  xlab("Index") + ylab(expression(l[max])) +
  geom_text_repel(aes(label=ifelse(((Cmaxcdi)> ref),paste0(1:n),""))) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=20,family="serif")
  )
  print(s1)
}


#Local influence continuous part
#Case-weight perturbation
# if you want to highlight a point, use the argument ref to do so

localcont=function(mof,ref=0.5){
  degree=3
  order=2
  n=mod$N
  s=ncol(mod$sigma.x)
  mu=mod$mu.fv
  phi=mod$sigma.fv
  zast<-ifelse(mod$y==0 | mod$y==1, 1,0)
  
  if(is.null(mod$mu.s)==FALSE){
    p=ncol(mod$mu.x)-ncol(mod$mu.s)  
    nq=ncol(mod$mu.s)
    qj=0
    for(i in 1:nq){
      qj[i]=length(mod$mu.coefSmo[[i]][11]$knots)
    }
    X=as.matrix(mod$mu.x[,-which(colnames(mod$mu.x) %in% colnames(mod$mu.s))])
    bet=matrix(mod$mu.coefficients[-which(colnames(mod$mu.x) %in% colnames(mod$mu.s))], nrow=p)  
    gamms=list(NULL)
    for(i in 1:nq){
      gamms[[i]]=mod$mu.coefSmo[[i]][1]$coef
    }
    gamm=Reduce("rbind",gamms)
    Bs<-list(NULL)
    Ds<-list(NULL)
    for (i in 1:nq){
      ajus=pb(mod$mu.x[,which(colnames(mod$mu.x) %in% colnames(mod$mu.s)[i])],
              degree=degree,order=order,
              control=pb.control(inter=(length(mod$mu.coefSmo[[1]]$knots)-3)),
              lambda = mod$mu.lambda[i])
      Bs[[i]]=attr(ajus,"X")
      Ds[[i]]=t(attr(ajus,"D"))%*%attr(ajus,"D")*mod$mu.lambda[i]
    }
    
    B=Reduce("cbind",Bs)
    D=0.5*bdiag(Ds)
    
  } else{
    p=ncol(mod$mu.x)
    bet=matrix(mod$mu.coefficients,nrow=p)
    X=mod$mu.x}
  
  kapp=matrix(mod$sigma.coefficients,nrow=s)
  N=mod$sigma.x  
  T1=diag(gmu(mod)$dmu,n)
  T2=diag(gphi(mod)$dphi,n)
  
  
  Phi=diag(mod$sigma.fv,n)
  logy=ifelse(mod$y>0 & mod$y<1,log(mod$y),0)
  logyc=ifelse(mod$y>0 & mod$y<1,log(1-mod$y),0)
  Q= Qstar=R=P=W=Wstar=diag(0,n)
  ystar=mustar=u=0
  ystar=logy-logyc
  mustar=digamma(mu*phi)-digamma((1-mu)*phi)
  u=mu*(ystar-mustar)+logyc-digamma((1-mu)*phi)+digamma(phi)
  Q1=(phi^2*(trigamma(mu*phi)+trigamma((1-mu)*phi))+phi*(ystar-mustar)*gmu(mod)$dmu*gmu(mod)$d2mu2)*gmu(mod)$dmu^2
  Q=diag(as.numeric((1-zast)*Q1),n)
  Qstar1=(phi*(trigamma(mu*phi)*mu-trigamma((1-mu)*phi)*(1-mu))-(ystar-mustar))*gmu(mod)$dmu*gphi(mod)$dphi
  Qstar=diag(as.numeric((1-zast)*Qstar1),n)
  R1=(trigamma(mu*phi)*mu^2+trigamma((1-mu)*phi)*(1-mu)^2-trigamma(phi)+u*gphi(mod)$dphi*gphi(mod)$d2phi2)*gphi(mod)$dphi^2
  R=diag(as.numeric((1-zast)*R1))
  
  Ubb=-t(X)%*%Q%*%X
  Ukk=-t(N)%*%R%*%N
  Ugg=-t(B)%*%Q%*%B-D
  Ubk=-t(X)%*%Qstar%*%N
  Ugk=-t(B)%*%Qstar%*%N
  Ubg=-t(X)%*%Q%*%B
  Utt=rbind(cbind(Ubb,Ubg,Ubk),cbind(t(Ubg),Ugg,Ugk),cbind(t(Ubk),t(Ugk),Ukk))
  Uttinv=ginv(as.matrix(-Utt))
  
  deltabc=t(X)%*%T1%*%Phi%*%diag((ystar-mustar),n)
  deltagc=t(B)%*%T1%*%Phi%*%diag((ystar-mustar),n)
  deltakc=t(N)%*%T2%*%diag(u,n)
  deltacasos=rbind(deltabc,deltagc,deltakc)
  Zc=t(deltacasos)%*%Uttinv%*%deltacasos
  Cmaxc=abs(eigen(Zc)$vectors[,1])
  s1<-qplot(seq_along(Cmaxc),Cmaxc,geom = "point",label = seq(1,length(Cmaxc),1))+ 
    xlab("Index") + ylab(expression(l[max])) +
    geom_text_repel(aes(label=ifelse(((Cmaxc)> ref),paste0(1:n),""))) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "grey"),
      text=element_text(size=20,family="serif")
    )
  print(s1)
  
}
            

