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
  #if(any(mu <= 0)|any(mu >= 1)) stop(paste("mu must be between 0 and 1","\n",""))
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




            
##RQR
RQR=function(mod){
  y=mod$y
  tau=mod$tau.coefficients
  rho=mod$nu.coefficients
  Fl=mod$nu.x
  M=mod$tau.x
  predf<-Fl%*%rho
  predm<-M%*%tau
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




