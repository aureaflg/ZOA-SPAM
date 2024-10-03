library(gamlss)
library(purrr)
library(ggplot2)
library(Matrix)
library(matrixcalc)
library(GoFKernel)
library(pracma)
library(rootSolve)
library(ggrepel)


dSIM<-function (x, mu = 0.5, sigma = 1, log = FALSE) 
{
    #if(any(0.99999<mu & mu<1.000001)){
    #  mu[which(0.99999<mu & mu<1.000001)]=0.9999
    #}
    if (any(mu <= 0)) 
        stop(paste("mu must be between 0 and 1", "\n",""))
    if (any(sigma <= 0)) 
        stop(paste("sigma must be positive", "\n",""))
    if (any(x <= 0) || any(x >= 1)) 
        stop(paste(" must be between 0 and 1", "\n",""))
    logpdf <- -((x - mu)/(mu * (1 - mu)))^2/(2 * x * (1 - x) * 
        sigma) - (log(2 * pi * sigma) + 3 * (log(x) + log(1 - 
        x)))/2
    if (!log) 
        logpdf <- exp(logpdf)
    logpdf
}

pSIM<- function(q, mu=0.5, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(q <= 0) || any(q >= 1)) 
        stop(paste("q must be between 0 and 1", "\n", 
            ""))
    if (any(mu <= 0) || any(mu >= 1)) 
        stop(paste("mu must be between 0 and 1", "\n", 
            ""))
    if (any(sigma <= 0)) 
        stop(paste("sigma must be positive", "\n", 
            ""))
    lp <- pmax.int(length(q), length(mu), length(sigma))
    q <- rep(q, length = lp)
    sigma <- rep(sigma, length = lp)
    mu <- rep(mu, length = lp)
    zero <- rep(0, length = lp)
    pdf <- function(x, mu, sigma) 1/sqrt(2 * pi * sigma * (x * 
        (1 - x))^3) * exp(-1/2/sigma * (x - mu)^2/(x * (1 - 
        x) * mu^2 * (1 - mu)^2))
    cdfun <- function(upper, mu, sigma) {
        int <- integrate(pdf, lower = 0, upper = upper, mu, sigma)
        int$value
    }
    Vcdf <- Vectorize(cdfun)
    cdf <- Vcdf(upper = q, mu = mu, sigma = sigma)
    if (lower.tail == TRUE) 
        cdf <- cdf
    else cdf <- 1 - cdf
    if (log.p == FALSE) 
        cdf <- cdf
    else cdf <- log(cdf)
    cdf
}

qSIM<-function(p,mu=0.5,sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(sigma <= 0)) 
            stop(paste("sigma must be positive", "\n", 
                ""))
        if (any(mu <= 0) || any(mu >= 1)) 
            stop(paste("mu must be between 0 and 1", "\n", 
                ""))
        if (log.p == TRUE) 
            p <- exp(p)
        else p <- p
        if (lower.tail == TRUE) 
            p <- p
        else p <- 1 - p
        if (any(p < 0) | any(p > 1)) 
            stop(paste("p must be between 0 and 1", "\n", 
                ""))
        lp <- max(length(p), length(mu), length(sigma))
        p <- rep(p, length = lp)
        sigma <- rep(sigma, length = lp)
        mu <- rep(mu, length = lp)
        q <- rep(0, lp)
        h1 <- function(x, mu, sigma, p) pSIM(x, mu, sigma) - p
        uni <- function(mu, sigma, p) {
            val <- uniroot(h1, c(.Machine$double.eps, 1-.Machine$double.eps), mu = mu, sigma = sigma, 
                p = p)
            val$root
        }
        UNI <- Vectorize(uni)
        q <- UNI(mu = mu, sigma = sigma, p = p)
        q
}


dSIMINF<-function(x, mu = 0.5, sigma = 1, 
                  nu = 0.1, tau = 0.1, log = FALSE)
{ 
  if (any(mu <= 0) | any(mu >= 1) )  
    stop(paste("mu must be between 0 and 1", "\n", "")) 
  if (any(sigma <= 0))  
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0) )  
    stop(paste("nu must greated than 0", "\n", ""))           
  if (any(tau <= 0) )  
    stop(paste("tau must greated than 0", "\n", "")) 
  if (any(x < 0) | any(x > 1))  
    stop(paste("x must be 0<=x<=1, i.e. 0 to 1 inclusively", "\n", ""))  
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
    logfy[i] <- ifelse((x[i]>0 & x[i]<1), dSIM(x[i], mu=mu[i],
                                               sigma=sigma[i], log=TRUE), 0)
    logfy[i] <- ifelse((x[i]==0), log(nu[i]), logfy[i])          
    logfy[i] <- ifelse((x[i]==1), log(tau[i]) , logfy[i])
    logfy[i] <- logfy[i] - log(1+nu[i]+tau[i])   
  }
  
  if(log==FALSE) fy <- exp(logfy) else fy <- logfy
  fy
}

pSIMINF <- function(q, mu = 0.5, sigma = 1, nu = 0.1, tau = 0.1, 
                    lower.tail = TRUE, log.p = FALSE)
{     
  if (any(mu <= 0) | any(mu >= 1) )  
    stop(paste("mu must be between 0 and 1", "\n", "")) 
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0) )  
    stop(paste("nu must greated than 0", "\n", ""))           
  if (any(tau <= 0) )  
    stop(paste("tau must greated than 0", "\n", "")) 
  if (any(q < 0) | any(q > 1))  
    stop(paste("y must be 0<=y<=1, i.e. 0 to 1 inclusively", "\n", "")) 
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
      p[k]<-nu[k] + pSIM(q[k], mu = mu[k], 
                         sigma=sigma[k],lower.tail=TRUE,log.p=FALSE)
    }
    p[k] <- p[k]/(1+nu[k]+tau[k])
  }
  
  if(lower.tail==TRUE) p <- p else p=1-p
  if(log.p==FALSE) p <- p else p <- log(p)    
  p
}

qSIMINF <- function(p, mu = 0.5, sigma = 1, nu = 0.1, tau = 0.1, 
                    lower.tail = TRUE, log.p = FALSE)
{ if (any(mu <= 0) | any(mu >= 1) ) 
  stop(paste("mu must be between 0 and 1", "\n", "")) 
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0) )  
    stop(paste("nu must greated than 0", "\n", ""))           
  if (any(tau <= 0) )  
    stop(paste("tau must greated than 0", "\n", "")) 
  if (any(p < 0) | any(p > 1))  
    stop(paste("p must be between 0 and 1", "\n", ""))    
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
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
      q[k]<-qSIM((p[k]- (nu[k]/(1 + nu[k] + tau[k])))/
                   (1/(1 + nu[k] + tau[k])),
                 mu = mu[k],sigma = sigma[k],
                 lower.tail=TRUE, log.p=FALSE)
    }
  }
  q
}

rSIMINF <- function(n, mu = 0.5, sigma = 1, nu = 0.1, tau = 0.1)
{ if (any(mu <= 0) | any(mu >= 1) )  
  stop(paste("mu must be between 0 and 1", "\n", "")) 
  if (any(sigma <= 0)) 
    stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu <= 0) )  
    stop(paste("nu must greated than 0", "\n", ""))           
  if (any(tau <= 0) )  
    stop(paste("tau must greated than 0", "\n", "")) 
  if (any(n <= 0))  
    stop(paste("n must be a positive integer", "\n", ""))  
  n <- ceiling(n)
  p <- runif(n)
  r <- qSIMINF(p, mu=mu, sigma=sigma, nu=nu, tau=tau)
  r
}


SIMINF <- function (mu.link = "logit", sigma.link = "logit", 
                  nu.link = "log", tau.link = "log")
{
  mstats <- checklink("mu.link", "SIMINF", substitute(mu.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "SIMINF", substitute(sigma.link), 
                      c("inverse", "log", "identity", "sqrt", "own"))   
  vstats <- checklink("nu.link", "SIMINF", substitute(nu.link),    
                      c("inverse", "log", "identity", "own"))
  tstats <- checklink("tau.link", "SIMINF", substitute(tau.link),   
                      c("inverse", "log", "identity", "own")) 
  structure(
    list(family = c("SIMINF", "Simplex Inflated"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
         nopar = 4, 
         type = "Mixed",
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
         dldm = function(y,mu,sigma) {
         dldm <- ifelse(((y==0)|(y==1)),0,-((mu - y) * (mu^2 + y - 2 * y * mu))/(sigma * 
                                          y * (y - 1) * mu^3 * (mu - 1)^3))
         dldm 
         },
         d2ldm2 = function(y,mu,sigma) {
         d <- 1/(mu^3*(1-mu)^3)
         d1<-(3*sigma)/(mu*(1-mu))
         d2ldm2 <- ifelse( ((y==0)|(y==1)), 0, -((1/sigma)*(d1+d)))
         d2ldm2 
         },
         dldd = function(y,mu,sigma) {
         dldd <- ifelse( ((y==0)|(y==1)), 0 , 
                         ((y - mu)^2)/((mu^2) * (1 - mu)^2 * y * (1 - y) * 2
                                       * sigma^2) - (1/(2 * sigma)))
         dldd 
         }, 
         d2ldd2 = function(y,mu,sigma,nu,tau) { 
           d2ldd2 <- ifelse( ((y==0)|(y==1)), 0 , 
                             -(1/(2*sigma^2)))
           d2ldd2
         }, 
         dldv = function(y,nu,tau)  {dldv <- ifelse(y==0,(1/nu),0) -(1/(1+nu+tau))
         dldv
         }, 
         d2ldv2 = function(nu,tau) {d2ldv2 <- -(1+tau)/(nu*((1+nu+tau)^2))
         d2ldv2
         },
         dldt = function(y,nu,tau) { dldt <- ifelse(y==1,(1/tau),0) -(1/(1+nu+tau))
         dldt
         }, 
         
         d2ldt2 = function(nu,tau) {
           d2ldt2 <- -(1+nu)/(tau*((1+nu+tau)^2))
           d2ldt2
         },
         d2ldmdd = function(y,mu,sigma) { 
         d2ldmdd <- rep(0,length(y))
         d2ldmdd 
         },
         d2ldmdv = function(y) {
           d2ldmdv <- rep(0,length(y))
           d2ldmdv
         },
         
         d2ldmdt = function(y) {
           d2ldmdt <- rep(0,length(y))
           d2ldmdt
         },
         
         d2ldddv = function(y) {
           d2ldddv <- rep(0,length(y))
           d2ldddv
         },
         
         d2ldddt = function(y) {
           d2ldddt <- rep(0,length(y)) 
           d2ldddt 
         },
         
         d2ldvdt = function(nu,tau) {
           d2ldvdt <- 1/((1+nu+tau)^2)
           d2ldvdt  
         },
         
         G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
           -2*dSIMINF(y,mu,sigma,nu,tau,log=TRUE),                     
         rqres = expression(rqres(pfun="pSIMINF", type="Mixed",  mass.p=c(0,1),  
                                  prob.mp=cbind(nu/(1+nu+tau),tau/(1+nu+tau)), y=y, mu=mu, 
                                  sigma=sigma, nu=nu, tau=tau)),
         mu.initial = expression(mu <- (y+mean(y))/2),    #(y+mean(y))/2),# rep(mean(y),length(y)) 
         sigma.initial = expression(sigma <- rep(0.5, length(y))),
         nu.initial = expression(nu <- rep(0.3, length(y))), 
         tau.initial = expression(tau <-rep(0.3, length(y))), 
         mu.valid = function(mu) all(mu > 0 & mu < 1) , 
         sigma.valid = function(sigma)  all(sigma > 0), 
         nu.valid = function(nu)  all(nu > 0) , 
         tau.valid = function(tau) all(tau > 0), 
         y.valid = function(y)  all(y >= 0 & y <= 1)
    ),
    class = c("gamlss.family","family"))
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
  ut<- pSIMINF(y, mu, phi, nu = p0/(1-p0-p1), tau = p1/(1-p0-p1))
  rQ<- qnorm(ut)
  resiquan<-numeric()
  resiquan[y==0]<- rQ0
  resiquan[y==1]<- rQ1
  resiquan[which(y>0 & y<1)]<- rQ[which(y>0 & y<1)]
  return(resiquan)
}




##Standard Error discrete part
n=mod$N
tau=mod$tau.coefficients
rho=mod$nu.coefficients
Fl=mod$nu.x
M=mod$tau.x
k0=length(rho)
k1=length(tau)

predf<-Fl%*%rho
predm<-M%*%tau
p0=exp(predf)/(1+exp(predf)+exp(predm))
p1=exp(predm)/(1+exp(predf)+exp(predm))
zast<-ifelse(mod$y==0 | mod$y==1, 1,0)
Zast<-diag(c(zast))

IFrho<- t(Fl)%*%diag(as.vector(p0*(1-p0)),n)%*%Fl
IFtau<- t(M)%*%diag(as.vector(p1*(1-p1)),n)%*%M
IFpr<- t(Fl)%*%diag(as.vector(-p0*p1),n)%*%M
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
coefficientsd = data.frame(coefdisc, EPdisc, walddisc, valorpdisc)
colnames(coefficientsd) = c("Estimate","Std.err", "Wald", "Pr(>|W|)")
printCoefmat(as.matrix(coefficientsd), digits = options()$digits)

##Standard Error continuous part
degree=3
order=2
n=mod$N
mu=mod$mu.fv
phi=mod$sigma.fv
zast<-ifelse(mod$y==0 | mod$y==1, 1,0)
s=ncol(mod$sigma.x)
p=ncol(mod$mu.x)-ncol(mod$mu.s)

nq=ncol(mod$mu.s)
qj=0
for(i in 1:nq){
  qj[i]=length(mod$mu.coefSmo[[i]][11]$knots)
}

X=as.matrix(mod$mu.x[,-which(colnames(mod$mu.x) %in% colnames(mod$mu.s))])
bet=matrix(mod$mu.coefficients[-which(colnames(mod$mu.x) %in% colnames(mod$mu.s))], nrow=p)
kapp=matrix(mod$sigma.coefficients,nrow=s)
gamms=list(NULL)
for(i in 1:nq){
  gamms[[i]]=mod$mu.coefSmo[[i]][1]$coef
}
gamm=Reduce("rbind",gamms)

N=mod$sigma.x
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
D=bdiag(Ds)

gmu=function(mod,alpha=NULL){
  x = mod$mu.link
  mu=mod$mu.fv
  
  if(x=="logit"){
    dmu=(1/mu+1/(1-mu))^(-1)
    d2mu2=-1/mu^2+1/(1-mu)^2
  }
  if(x=="cloglog"){
    dmu=(1/(log(1-mu)*(1-mu)))^(-1)
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

T1=diag(gmu(mod)$dmu,n)
T2=diag(gphi(mod)$dphi,n)

W1<- (1/phi)*(((3*phi)/(mu*(1-mu)))+(1/(mu^3*(1-mu)^3)))*gmu(mod)$dmu^2
W<- diag(as.numeric((1-zast)*W1),n)
P1 <- as.vector((1/(2*phi^2))*gphi(mod)$dphi^2)
P<- diag(as.numeric((1-zast)*P1),n)


Kbb= t(X)%*%W%*%X
Kkk= t(N)%*%P%*%N
Kgg= t(B)%*%W%*%B+D
Kbk= matrix(0,ncol=s,nrow=p)
Kgk= matrix(0,ncol=qj,nrow=s)
Kbg= t(X)%*%W%*%B
Kggl=t(B)%*%W%*%B
Ktt=rbind(cbind(Kbb,Kbk,Kbg),cbind(t(Kbk),Kkk,Kgk),cbind(t(Kbg),t(Kgk),Kgg))
Kttinv=ginv(as.matrix(Ktt))
coef=c(mod$mu.coefficients[-5])
EP=sqrt(diag(Kttinv)[1:length(coef)])
VAR=diag(Kttinv)[1:length(coef)]
valorp = 0
wald = 0
for(i in 1:length(coef)){
  wald[i] = wald.test(VAR[i], coef[i], Terms = 1)$result$chi2[1]
  valorp[i] = wald.test(VAR[i], coef[i], Terms = 1)$result$chi2[3]
}
coefficients = data.frame(coef, EP, valorp)
colnames(coefficients) = c("Estimate","Std.err",  "Pr(>|W|)")
printCoefmat(as.matrix(coefficients), digits = options()$digits)

##Residuals
smp <- data.frame(norm = modolos$residuals)
envesim<-ggplot(data = smp, mapping = aes(sample = norm))  +
  stat_qq_band(conf = 0.95) +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Quantile Residuals ZOAS-SPAM") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )
qplot(seq_along(mod$residuals),mod$residuals)+ 
  xlab("Index") + ylab("Quantile Residuals") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

##Diagnostic
#Hessian matrix for discrete part
tau=mod$tau.coefficients
rho=mod$nu.coefficients
Fl=mod$nu.x
M=mod$tau.x
k0=length(rho)
k1=length(tau)

predf<-Fl%*%rho
predm<-M%*%tau
p0=exp(predf)/(1+exp(predf)+exp(predm))
p1=exp(predm)/(1+exp(predf)+exp(predm))
zast<-ifelse(mod$y==0 | mod$y==1, 1,0)
Zast<-diag(c(zast))

Urho<- t(Fl)%*%(zast*(1-mod$y)-p0)
Utau<- t(M)%*%(zast*mod$y-p1)
Uesc<-matrix(c(Urho,Utau), ncol=1)

Urr<- -t(Fl)%*%diag(as.vector(p0*(1-p0)),n)%*%Fl
Uttau<- -t(M)%*%diag(as.vector(p1*(1-p1)),n)%*%M
Urtau<- -t(Fl)%*%diag(as.vector(-p0*p1),n)%*%M
Udisc<-matrix(rbind(cbind(Urr, Urtau),cbind(t(Urtau), Uttau)),nrow=k0+k1,ncol=k0+k1)
Udiscinv=ginv(as.matrix(-Udisc))


#Local influence
#Case-weight perturbation
deltarhoc=t(Fl)%*%diag(as.numeric(zast*(1-mod$y)-p0),n)
deltatauc=t(M)%*%diag(as.numeric(zast*mod$y-p1),n)
deltacasosd=rbind(deltarhoc,deltatauc)
Zcdi=t(deltacasosd)%*%Udiscinv%*%deltacasosd
Cmaxcdi=abs(eigen(Zcdi)$vectors[,1])
qplot(seq_along(Cmaxcdi),Cmaxcdi,geom = "point",label = seq(1,length(Cmaxcdi),1))+ 
  xlab("Index") + ylab(expression(l[max])) +
  geom_text_repel(aes(label=ifelse(((Cmaxcdi)> 0.2),paste0(1:n),""))) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )



#Hessian matrix for continuous part
Phistar=diag(1/mod$sigma.fv,n)
phi=mod$sigma.fv
mu=mod$mu.fv
y=mod$y
Q=Qstar=R=U=S=diag(0,n)
d=ul=0
a=0
d=ifelse(mod$y==0 | mod$y==1,0,(y-mu)^2/(y*(1-y)*(mu^2)*(1-mu)^2))
a=-1/2*phi+d*(1/(2*phi^2))
U=diag((1/(mu*(1-mu)))*(d + (1/(mu^2*(1-mu)^2))),n)
ul=(2*(y-mu)*U)/(mu*(1-mu)) + (3-6*mu)/(mu^4*(1-mu)^4) + ((1-2*mu)*d)/(mu^2*(1-mu)^2)
Q=diag(as.vector((1-zast)*((1/phi)*((diag(U)-(y-mu)*ul)+diag(U)*(y-mu)*gmu(mod)$dmu*gmu(mod)$d2mu2)*gmu(mod)$dmu^2)),n)
Qstar=diag((1-zast)*(1/phi^2)*diag(U)*(y-mu)*gmu(mod)$dmu*gphi(mod)$dphi,n)
S=diag((1-zast)*((1/(2*phi^2)-1/phi^3*d+a*gphi(mod)$dphi*gphi(mod)$d2phi2)*gphi(mod)$dphi^2),n)

Ubb=-t(X)%*%Q%*%X
Ukk=-t(N)%*%S%*%N
Ugg=-t(B)%*%Q%*%B-D
Ubk= -t(X)%*%Qstar%*%N
Ugk=-t(B)%*%Qstar%*%N
Ubg=-t(X)%*%Q%*%B
Uggl=-t(B)%*%Q%*%B
Utt=rbind(cbind(Ubb,Ubg,Ubk),cbind(t(Ubg),Ugg,Ugk),cbind(t(Ubk),t(Ugk),Ukk))
Uttinv=ginv(as.matrix(-Utt))


#Influencia local
#perturbacao de casos
deltabc=t(X)%*%T1%*%U%*%Phistar%*%diag((y-mu),n)
deltagc=t(B)%*%T1%*%U%*%Phistar%*%diag((y-mu),n)
deltakc=t(N)%*%T2%*%diag(a,n)
deltacasos=rbind(deltabc,deltagc,deltakc)
Zc=t(deltacasos)%*%Uttinv%*%deltacasos
Cmaxc=Re(eigen(Zc)$vectors[,1])/sum(abs(eigen(Zc)$vectors[,1]))
qplot(seq_along(Cmaxc),Cmaxc,geom = "point",label = seq(1,length(Cmaxc),1))+ 
  xlab("Index") + ylab("dmáx") +
  geom_text_repel(aes(label=ifelse(((Cmaxc)> 0.24),paste0(1:n),""))) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )


