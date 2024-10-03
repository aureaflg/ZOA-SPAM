library(gamlss)
library(matrixcalc)
library(ggplot2)
library(stringr)
library(Matrix)
library(GoFKernel)
library('VGAM')
library(MASS)
source("C:/Users/aiflg/Documents/Códigos dissertação/pbfake.R")
source("C:/Users/aiflg/Documents/Códigos dissertação/estimação lambda.R")
source("C:/Users/aiflg/Documents/Códigos dissertação/rbetaretangular.R")


dBRINF<-function (x, mu = 0.5, alpha=0.5, sigma = 1, nu = 0.1, tau = 0.1, 
                  log = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(alpha <= 0) | any(alpha >= 1)) 
    stop(paste("alpha must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma greated than 0", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must greated than 0", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must greated than 0", "\n", ""))
  if (any(x < 0) | any(x > 1)) 
    stop(paste("x must be 0<=x<=1, i.e. 0 to 1 inclusively", 
               "\n", ""))
  lp <- pmax.int(length(x), length(mu), length(sigma),length(alpha),
                 length(nu),length(tau))
  x <- rep(x, length = lp)
  sigma <- rep(sigma, length = lp)
  alpha <- rep(alpha, length = lp)
  mu <- rep(mu, length = lp)
  nu <- rep(nu, length = lp)
  tau <- rep(tau, length = lp)
  logfy <- rep(0, length = lp)
  for(i in 1:lp){
    logfy[i] <- ifelse((x[i]>0 & x[i]<1), dBRr(x[i], mu=mu[i], alpha=alpha[i], sigma=sigma[i], log = TRUE), 0)
    logfy[i] <- ifelse((x[i]==0), log(nu[i]), logfy[i])          
    logfy[i] <- ifelse((x[i]==1), log(tau[i]) , logfy[i])
    logfy[i] <- logfy[i] - log(1+nu[i]+tau[i])   
  }
  if (log == FALSE) 
    fy <- exp(logfy)
  else fy <- logfy
  fy
}

pBRINF<-function (q, mu = 0.5, alpha= 0.5, sigma = 1, nu = 0.1, tau = 0.1, lower.tail = TRUE, 
                  log.p = FALSE) {
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(alpha <= 0) | any(alpha >= 1)) 
    stop(paste("alpha must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greated than 0", "\n", ""))
  if (any(nu <= 0)) 
    stop(paste("nu must greated than 0", "\n", ""))
  if (any(tau <= 0)) 
    stop(paste("tau must greated than 0", "\n", ""))
  if (any(q < 0) | any(q > 1)) 
    stop(paste("y must be 0<=y<=1, i.e. 0 to 1 inclusively", 
               "\n", ""))
  lp <- pmax.int(length(q), length(mu), length(sigma),length(alpha),
                 length(nu),length(tau))
  q <- rep(q, length = lp)
  sigma <- rep(sigma, length = lp)
  alpha <- rep(alpha, length = lp)
  mu <- rep(mu, length = lp)
  nu <- rep(nu, length = lp)
  tau <- rep(tau, length = lp)
  p <- rep(0, length = lp)
  for(k in 1:lp){
    if(q[k]==0){p[k]<-nu[k]}
    if(q[k]==1){p[k]<-1+nu[k]+tau[k]}
    if(q[k]>0 & q[k]<1){
      p[k]<-nu[k] + pBRr(q[k], mu = mu[k], alpha=alpha[k],
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

qBRINF<-function (p, mu = 0.5, alpha=0.5, sigma = 1, nu = 0.1, tau = 0.1, lower.tail = TRUE, 
                  log.p = FALSE) {
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(alpha <= 0) | any(alpha >= 1)) 
    stop(paste("alpha must be between 0 and 1", "\n", ""))
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
  lp <- pmax.int(length(p), length(mu), length(sigma),length(alpha),
                 length(nu),length(tau))
  p <- rep(p, length = lp)
  sigma <- rep(sigma, length = lp)
  alpha <- rep(alpha, length = lp)
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
      q[k]<-qBRr((p[k]- (nu[k]/(1 + nu[k] + tau[k])))/
                   (1/(1 + nu[k] + tau[k])),
                 mu = mu[k],alpha=alpha[k],sigma = sigma[k],
                 lower.tail=TRUE, log.p=FALSE)
    }
  }
  q
}

rBRINF<-function (n, mu = 0.5, alpha=0.5, phi = 1, p0 = 0.1, p1 = 0.1) 
{
  y<-rep(0,n)
  lp <- pmax.int(length(y), length(mu), length(sigma),length(alpha),
                 length(p0),length(p1))
  
  phi <- rep(phi, length = lp)
  alpha <- rep(alpha, length = lp)
  mu <- rep(mu, length = lp)
  p0 <- rep(p0, length = lp)
  p1 <- rep(p1, length = lp)
  
  ybetaR<-rBRr(n,mu=mu,alpha=0.7,phi=phi)
  aux_unif1<-runif(n)
  tau<-p0+p1
  eta<-p1/tau
  sim_ber<-mapply(rbinom,1,1,eta)
  y<-ifelse(aux_unif1<tau,sim_ber,ybetaR)
  y
}

pBRINF2<-function (q, mu = 0.5, alpha= 0.5, sigma = 1, 
                   p0 = 0.1, p1 = 0.1, lower.tail = TRUE, 
                   log.p = FALSE) {
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(alpha <= 0) | any(alpha >= 1)) 
    stop(paste("alpha must be between 0 and 1", "\n", ""))
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
                 length(alpha),
                 length(p0),length(p1))
  q <- rep(q, length = lp)
  sigma <- rep(sigma, length = lp)
  alpha <- rep(alpha, length = lp)
  mu <- rep(mu, length = lp)
  nu <- rep(nu, length = lp)
  tau <- rep(tau, length = lp)
  p <- rep(0, length = lp)
  for(k in 1:lp){
    if(q[k]==0 & q[k]==1){
      p[k]=(p0[k]+p1[k])*pbinom(q[k],1,p1[k]/(p0[k]+p1[k]))
    }
    if(q[k]>0 & q[k]<1){
      p[k] = (1-p0[k]-p1[k])*pBRr(q[k], mu = mu[k], alpha=alpha[k],
                                  sigma=sigma[k],lower.tail=TRUE,log.p=FALSE)
    }
  }
  if (lower.tail == TRUE) p <- p
  else p = 1 - p
  if (log.p == FALSE) p <- p
  else p <- log(p)
  p
}

pZOABR<- function(y,p0,p1,mu,sigma,alpha){
  part<-which(y>0 & y<1)
  y<-y[part]
  p0<-p0[part]
  p1<-p1[part]
  mu<-mu[part]
  sigma<-sigma[part]
  epi<- 1-sqrt(1- 4*alpha*mu*(1-mu))
  delta<- (mu - 0.5*epi)/(1-epi)
  BR<- epi*punif(y,0,1) + (1-epi)*pbeta(y,delta*sigma,
                                        (1-delta)*sigma)
  tau<- p0+p1
  eta<- p1/tau
  FD<- tau*pbinom(y,1,eta) + (1-tau)*BR
  return(FD)
}

qBRINF<-function (p, mu = 0.5, alpha=0.5, sigma = 1, nu = 0.1, tau = 0.1, lower.tail = TRUE, 
                  log.p = FALSE) {
  if (any(mu <= 0) | any(mu >= 1)) 
    stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(alpha <= 0) | any(alpha >= 1)) 
    stop(paste("alpha must be between 0 and 1", "\n", ""))
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
  lp <- pmax.int(length(p), length(mu), length(sigma),length(alpha),
                 length(nu),length(tau))
  p <- rep(p, length = lp)
  sigma <- rep(sigma, length = lp)
  alpha <- rep(alpha, length = lp)
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
      q[k]<-qBRr((p[k]- (nu[k]/(1 + nu[k] + tau[k])))/
                   (1/(1 + nu[k] + tau[k])),
                 mu = mu[k],alpha=alpha[k],sigma = sigma[k],
                 lower.tail=TRUE, log.p=FALSE)
    }
  }
  q
}








###Algorithm for discrete part
EF<- function(y,Fl,M, iter,trace) {
  
  k0<- ncol(Fl)
  k1<- ncol(M)
  n<- length(y)
  
  if(k0>1 & k1>1){ #to model the probability of zeros and ones
    
    con<-gamlss.control(trace=FALSE)
    fit=gamlss(y~ 1, nu.formula=~ -1+Fl,tau.formula=~ -1+M,
               family=BEINF,control=con)
    
    #Initializes rho, tau beta estimators with MV
    rho.ini<-as.numeric((fit$nu.coefficients))
    tau.ini<-as.numeric((fit$tau.coefficients))
    
    #Initial values for rho, tau
    th<-matrix(c(rho.ini, tau.ini),ncol=1)
    TH=th
    
    rho<- cbind(th[1:k0])
    tau<- cbind(th[(k0+1):(k0+k1)])
    
    criterio=1
    count=0
    
    while(criterio > 0.00000001){
      count=count+1
      
      # Link functions
      predf<-Fl%*%rho
      predm<-M%*%tau
      p0=exp(predf)/(1+exp(predf)+exp(predm))
      p1=exp(predm)/(1+exp(predf)+exp(predm))
      zast<-ifelse(y==0 | y==1, 1,0)
      Zast<-diag(c(zast))
      
      Urho<- t(Fl)%*%(zast*(1-y)-p0)
      Utau<- t(M)%*%(zast*y-p1)
      U<-matrix(c(Urho,Utau), ncol=1)
      
      IFrho<- t(Fl)%*%diag(as.vector(p0*(1-p0)),n)%*%Fl
      IFtau<- t(M)%*%diag(as.vector(p1*(1-p1)),n)%*%M
      IFpr<- t(Fl)%*%diag(as.vector(-p0*p1),n)%*%M
      IF<-matrix(rbind(cbind(IFrho, IFpr),cbind(t(IFpr), IFtau)),nrow=k0+k1,ncol=k0+k1)
      
      th<-th + ginv(IF)%*%U
      criterio=sqrt(t(TH-th)%*%(TH-th))
      TH=th
      
      rho<- cbind(th[1:k0])
      tau<- cbind(th[(k0+1):(k0+k1)])
      
      if(trace==T){
        print(t(rho))
        print(t(tau))
        print(count)
      }
      
      if (count==iter)
      {
        break
      }
    }
    
    predf<-Fl%*%rho
    predm<-M%*%tau
    p0=exp(predf)/(1+exp(predf)+exp(predm))
    p1=exp(predm)/(1+exp(predf)+exp(predm))
    obj.out <- list(rho=rho, tau=tau, EP=sqrt(diag(solve(IF))), 
                    p0=p0,p1=p1,iter=count)
    #class(obj.out) <- "EF.ZOABR"
    
  } else{ #if we don't model the probability of zeros and ones
    
    con<-gamlss.control(trace=FALSE)
    fit=gamlss(y~ 1,family=BEINF,control=con)
    
    #Initializes rho, tau beta estimators with MV
    nu<- as.numeric(fit$nu.fv[1])
    tau<- as.numeric(fit$tau.fv[1])
    p0.ini<- nu/(1+nu+tau)
    p1.ini<- tau/(1+nu+tau)
    
    #Initial values for rho, tau
    th<-matrix(c(p0.ini, p1.ini),ncol=1)  
    TH=th
    
    p0<- th[1]
    p1<- th[2]
    
    criterio=1
    count=0
    
    while(criterio > 0.00001)
    {
      count=count+1
      zast<-ifelse(y==0 | y==1, 1,0)
      
      Up0<-  sum((zast*(1-y))/p0 - (1-zast)/(1-p0-p1))
      Up1<- sum((zast*y)/p1 - (1-zast)/(1-p0-p1))
      
      IFp0<- n*(1/p0 + 1/(1-p0-p1))
      IFp1<- n*(1/p1 + 1/(1-p0-p1))
      IFp01<- n*(1/(1-p0-p1))
      
      U<-matrix(c(Up0,Up1), ncol=1)
      
      IF<-matrix(c(IFp0, IFp01, IFp01, IFp1),nrow=2,ncol=2)
      
      th<-th + solve(IF)%*%U
      criterio=sqrt(t(TH-th)%*%(TH-th))
      TH=th
      
      p0<- th[1]
      p1<- th[2]
      
      if(trace==T){
        print(p0)
        print(p1)
        print(count)
      }
      
      if (count==iter)
      {
        break
      }
    }
    obj.out<-(list(p0=p0, p1=p1, EP=sqrt(diag(ginv(IF))), iter=count))
    #class(obj.out) <- "EF.ZOABR"
    
    
  }
  
  return(obj.out)  
  
}

##Expectation of log-likelihood function
Q<-function(param,v,y,N,X,gamm,D,B,lamb,lig,lig.phi){
  p=ncol(X)
  s=ncol(N)
  bet=param[1:p]
  kapp=param[(1+p):(s+p)]
  alpha=param[s+p+1]
  eta1=X%*%bet+B%*%gamm
  if(lig=='logit'){
    mu<- logitlink(eta1, inverse = TRUE)
  }else if(lig=='probit'){
    mu<- probitlink(eta1, inverse = TRUE)
  } else if(lig=='cloglog'){
    mu<- clogloglink(eta1, inverse = TRUE)
  } else if(lig=='cauchit'){
    mu<- cauchitlink(eta1, inverse = TRUE)
  } else if(lig=='loglog'){
    mu<- exp(-exp(-eta1))
  } else {print('Link function not defined')
  } 
  eta2=N%*%kapp
  if(lig.phi=='log'){
    phi<-exp(eta2)
  }else if(lig.phi=='sqrt'){
    phi<-eta2^2
  } else if(lig.phi=='identity'){
    phi<-eta2
  }  else {print('Link function not defined')
  } 
  epi<-1-sqrt(1-4*alpha*mu*(1-mu))
  delta<- (mu-0.5*epi)/(1-epi)
  a=delta*phi
  b=(1-delta)*phi
  Dlamb=lamb*D
  zast<-ifelse(y==0 | y==1, 1,0)
  dens<-ifelse(y>0 & y<1, dbeta(y,a,b,log=TRUE),0)
  Q=sum((1-zast)*(v*log(epi)+(1-v)*log(1-epi)+(1-v)*dens))-
    0.5*t(gamm)%*%Dlamb%*%gamm
  return(-Q)
}

###Algoritmo EM
samzoabr<-function(formula.mu=formula,formula.phi=~1,
                formula.nu=~1,formula.tau=~1,
                data=NULL,mu.link,phi.link="log",
                iter,trace,knots=20,degreepb=3,order=2,lambda=NULL,
                method="GAIC",k=2,c=FALSE){
  
  mcall <- if(is.null(data)){ terms(formula.mu, specials = "pb") 
  }else terms(formula.mu, specials = "pb", data = data)    
  mu.X <- if(is.null(data)){ model.matrix(mcall) 
  }else model.matrix(mcall, data = data)
  label <- colnames(mu.X)#attr(mcall,"term.labels") 
  aux <- which(str_detect(label, pattern = "pb")==T)
  aux1 <- which(colnames(mu.X)==label[aux])
  X <- matrix(mu.X[,-aux1], ncol=(ncol(mu.X)-length(aux1)))
  Z <- as.matrix(mu.X[,aux1])
  y <- if(is.null(data)){ model.frame(mcall)[,1] 
  }else model.frame(mcall, data = data)[,1]
  N <- if(is.null(data)){ model.matrix(formula.phi)
  }else model.matrix(formula.phi, data = data)
  Fl <- if(is.null(data)){
    if(length(attr(model.matrix(formula.nu),"assign"))==1){
      as.matrix(rep(1,nrow(mu.X)))
    }else{model.matrix(formula.nu)}
  }else{ 
    if(length(attr(model.matrix(formula.nu, data = data),"assign"))==1){
      as.matrix(rep(1,nrow(mu.X)))
    }else{model.matrix(formula.nu, data = data)}}
  M <- if(is.null(data)){
    if(length(attr(model.matrix(formula.tau),"assign"))==1){
      as.matrix(rep(1,nrow(mu.X)))
    }else{model.matrix(formula.tau)}
  }else{ 
    if(length(attr(model.matrix(formula.tau, data = data),"assign"))==1){
      as.matrix(rep(1,nrow(mu.X)))
    }else{model.matrix(formula.tau, data = data)}}
  
  
  if(is.null(data)){
    data=cbind(model.frame(mcall),model.frame(formula.phi))
  } else data=data
  
  p<-ncol(X)
  s<-ncol(N)
  nq<-ncol(Z)
  s0<-ncol(Fl)
  s1<-ncol(M)
  ntol<-nrow(X)
  
  ##Defining the basis matrix and penalization matrix
  Bs<-list(NULL)
  Ds<-list(NULL)
  qj=rep(0,nq)
  for (au in 1:nq){
    ajus<-pbfake(as.vector(Z[,au]),inter = knots, degree = degreepb,
                 order = order)
    Bs[[au]]=attr(ajus,"X")
    Ds[[au]]=t(attr(ajus,"D"))%*%attr(ajus,"D")
    bB=t(Bs[[au]])%*%rep(1,nrow(Bs[[au]]))
    bBl=c(sqrt(sum(bB^2)),rep(0,(length(bB)-1)))
    u=bB-bBl
    H=diag(1,ncol(Bs[[au]]))-(2/sqrt(sum(u^2))^2)*u%*%t(u)
    B1=Bs[[au]]%*%H
    Bs[[au]]=B1[,-1]
    D1=H%*%Ds[[au]]
    Ds[[au]]=D1[-1,]%*%H[,-1]
    qj[au]=ncol(Bs[[au]])
  }
  
  B=Reduce("cbind",Bs)
  D=as.matrix(bdiag(Ds))
  
  if(!is.null(lambda)){
    lamb=rep(0,sum(qj))
    if(length(lambda)==nq){
      for(au2 in 1:nq){
        if(au2==1){
          lamb[1:qj[au2]]<-rep(lambda[au2],qj[au2])
          next
        }
        lamb[(sum(qj[1:(au2-1)])+1):(sum(qj[1:(au2-1)])+qj[au2])]<-rep(lambda[au2],qj[au2])
      }
    }else lamb<-rep(lambda,sum(qj))
  }else{lamb<-rep(10,sum(qj))}
  
  
  #Start values
  bet<-solve(t(X)%*%X)%*%t(X)%*%y
  gamm<-ginv(t(B)%*%B+lamb*D)%*%t(B)%*%(y-X%*%bet)
  kapp<-solve(t(N)%*%N)%*%t(N)%*%y
  alpha<-(2*sum(y^2))/(sum(y^2)+ntol)
  par0 <- c(bet,kapp,alpha)
  zast<-ifelse(y==0 | y==1, 1,0)
  
  #Estimating discrete part
  
  fitdisc<- EF(y,Fl,M, iter=30,trace=trace)
  
  
  #Estimating continuous part
  i=0
  dif=1
  while(dif>0.001){
    
    i=i+1
    
    eta1<-X%*%bet+B%*%gamm
    eta2<-N%*%kapp
    if(phi.link=='log'){
      phi<-exp(eta2)
    }else if(phi.link=='sqrt'){
      phi<-eta2^2
    } else if(phi.link=='identity'){
      phi<-eta2
    }  else {print('Link function not defined')
    } 
    
    if(mu.link=='logit'){
      mu<- logitlink(eta1, inverse = TRUE)
    }else if(mu.link=='probit'){
      mu<- probitlink(eta1, inverse = TRUE)
    } else if(mu.link=='cloglog'){
      mu<- clogloglink(eta1, inverse = TRUE)
    } else if(mu.link=='cauchit'){
      mu<- cauchitlink(eta1, inverse = TRUE)
    } else if(mu.link=='loglog'){
      mu<- exp(-exp(-eta1))
    } else {print('Link function not defined')
    } 
    epi<-1-sqrt(1-4*alpha*mu*(1-mu))
    delta<-(mu-0.5+0.5*(1-epi))/(1-epi)
    a<-delta*phi
    b<-(1-delta)*phi
    
    ###E step
    v<-epi^(1-zast)/((epi+(1-epi)*dbeta(y,a,b))^(1-zast))
    
    ##M step
    gammnew=nlminb(gamm, Q, v = v, y = y, 
                    N = N, X = X, param = par0, D = D, B = B,lamb = lamb,
                    lig=mu.link,lig.phi=phi.link,lower=c(rep(-Inf,sum(qj))), 
                    upper=c(rep(Inf,sum(qj))))$par
    gamm=gammnew
    resu=nlminb(par0, Q, v = v, y = y, 
                N = N, X = X, gamm = gamm, D = D, B = B,lamb = lamb,
                lig=mu.link,lig.phi=phi.link,lower=c(rep(-Inf,s+p),0.001), 
                upper=c(rep(Inf,s+p),0.999))$par
    
    bet=resu[1:p]
    kapp=resu[(1+p):(s+p)]
    alpha=resu[s+p+1]
    dif=sqrt(t(resu-par0)%*%(resu-par0))
    par0=c(bet,kapp,alpha)
    
    if(is.null(lambda)){
      statuslam=1
    }else{statuslam=0}
    quant=estlambdainf(lambda=lamb,y=y,zast=zast,X=X,N=N,Ds=Ds,D=D,
                       Bs=Bs,B=B,qj=qj,bet=bet,gamm=gamm,kapp=kapp,
                       alpha=alpha,method=method,k=k,c=c,
                       statuslam=statuslam,lig=mu.link,lig.phi=phi.link,
                       s0=s0,s1=s1)
    
    
    lamb0=quant$lambda
    lamb=lamb0 
    #Stop if exceeds requested number of iterations
    
    if (i==iter)
    {
      break
    }
    if (trace==T){
      print(bet)
      print(kapp)
      print(alpha)
      print(i)
    }
    
  }
  quant1=estlambdainf(lambda=lamb,y=y,zast=zast,X=X,N=N,Ds=Ds,D=D,Bs=Bs,B=B,
                      qj=qj,bet=bet,gamm=gamm,kapp=kapp,alpha=alpha,
                      method=method,k=k,c=c,statuslam=0,
                      lig=mu.link,lig.phi=phi.link,s0=s0,s1=s1)
  mu.s=matrix(0,ncol=nq,nrow=ntol)
  mu.coefSmo = list()
  for(au3 in 1:nq){
    mu.coefSmo[[au3]] = list()
    if(au3==1){
      mu.s[,1]<-B[,(1:qj[au3])]%*%gamm[(1:qj[au3])]
      mu.coefSmo[[1]]$edf=quant$edf
      mu.coefSmo[[1]]$knots=knots
      mu.coefSmo[[1]]$lambda=lamb[1]
      mu.coefSmo[[1]]$coef=gamm[(1:qj[au3])]
      mu.coefSmo[[1]]$Bk=B[,(1:qj[au3])]
      mu.coefSmo[[1]]$Dk= Ds[[1]]
      mu.coefSmo[[1]]$B=B
      mu.coefSmo[[1]]$D=D
      mu.coefSmo[[1]]$gamm=gamm
      next
    }
    mu.s[,au3]<-B[,(sum(qj[1:(au3-1)])+1):(sum(qj[1:(au3-1)])+qj[au3])]%*%
      gamm[((sum(qj[1:(au3-1)])+1):(sum(qj[1:(au3-1)])+qj[au3]))]
    mu.coefSmo[[au3]]$edf=quant$edf
    mu.coefSmo[[au3]]$knots=knots 
    mu.coefSmo[[au3]]$lambda=lamb[(sum(qj[1:(au3-1)])+1)]
    mu.coefSmo[[au3]]$coef=gamm[((sum(qj[1:(au3-1)])+1):(sum(qj[1:(au3-1)])+qj[au3]))]
    mu.coefSmo[[au3]]$Bk=B[,(sum(qj[1:(au3-1)])+1):(sum(qj[1:(au3-1)])+qj[au3])]
    mu.coefSmo[[au3]]$Dk= Ds[[au3]]
    mu.coefSmo[[au3]]$B=B
    mu.coefSmo[[au3]]$D=D
    mu.coefSmo[[au3]]$gamm=gamm
  }
  
  phi.lp<-N%*%kapp
  if(phi.link=='log'){
    phi.fv<-exp(phi.lp)
  }else if(phi.link=='sqrt'){
    phi.fv<-phi.lp^2
  } else if(phi.link=='identity'){
    phi.fv<-phi.lp
  }  else {print('Link function not defined')
  } 
 
  mu.lp<-X%*%bet+B%*%gamm
  if(mu.link=='logit'){
    mu.fv<- logitlink(mu.lp, inverse = TRUE)
  }else if(mu.link=='probit'){
    mu.fv<- probitlink(mu.lp, inverse = TRUE)
  } else if(mu.link=='cloglog'){
    mu.fv<- clogloglink(mu.lp, inverse = TRUE)
  } else if(mu.link=='cauchit'){
    mu.fv<- cauchitlink(mu.lp, inverse = TRUE)
  } else if(mu.link=='loglog'){
    mu.fv<- exp(-exp(-mu.lp))
  } else {print('Link function not defined')
  } 
  
  if(ncol(Fl)<2 & ncol(M)<2){
    p0<-fitdisc$p0
    p1<-fitdisc$p1
    n0<- length(which(y==0))
    ut0<- runif(n0,0,p0)
    rQ0<- qnorm(ut0)
    n1<- length(which(y==1))
    ut1<- runif(n1,1-p1,1)
    rQ1<- qnorm(ut1)
    ut<- pBRINF(y, mu.fv, alpha, phi.fv, nu=p0/(1-p0-p1), 
                tau=p1/(1-p0-p1))
    rQ<- qnorm(ut)
  }
  if(ncol(Fl)>1 & ncol(M)>1){
    p0<-fitdisc$p0
    p1<-fitdisc$p1
    n0<- length(which(y==0))
    ut0<- runif(n0,0,p0[y==0])
    rQ0<- qnorm(ut0)
    n1<- length(which(y==1))
    ut1<- runif(n1,1-p1[y==1],1)
    rQ1<- qnorm(ut1)
    ut<- pBRINF(y, mu.fv, alpha, phi.fv, nu=p0/(1-p0-p1), 
                tau=p1/(1-p0-p1))
    rQ<- qnorm(ut)
  }
  resiquan<-numeric()
  resiquan[y==0]<- rQ0
  resiquan[y==1]<- rQ1
  resiquan[which(y>0 & y<1)]<- rQ[which(y>0 & y<1)]
  
  return(list(bet=bet,alpha=alpha,kapp=kapp,
              rho=fitdisc$rho,tau=fitdisc$tau,
              EPdisc=fitdisc$EP,
              niter=i,mu.fv=mu.fv,mu.s=mu.s,
              phi.fv=phi.fv,p0.fv=fitdisc$p0,
              p1.fv=fitdisc$p1,
              mu.link=mu.link,phi.link=phi.link,
              ntol=ntol,mu.formula=formula.mu,
              phi.formula=formula.phi, 
              nu.formula=formula.nu,
              tau.formula=formula.tau,
              mu.x=X,sigma.x=N,nu.x=Fl,
              tau.x=M,data=data,
              y=y,mu.coefSmo=mu.coefSmo,
              mu.lp=mu.lp,phi.lp=phi.lp,
              aic=quant1$aic,
              bic=quant1$bic,aicc=quant1$aicc,
              hqic=quant1$hqic,sabic=quant1$sabic))
}



RQR=function(mod){
  p0<-mod$p0.fv
  p1<-mod$p1.fv
  y<-mod$y
  mu<-mod$mu.fv
  phi<-mod$phi.fv
  alpha<-mod$alpha
  n0<-length(which(y==0))
  ut0<-runif(n0,0,p0[y==0])
  rQ0<-qnorm(ut0)
  n1<-length(which(y==1))
  ut1<-runif(n1,1-p1[y==1],1)
  rQ1<-qnorm(ut1)
  ut<-pBRINF(y, mu, alpha, phi, nu=p0/(1-p0-p1), 
              tau=p1/(1-p0-p1))
  rQ<- qnorm(ut)
  resiquan<-numeric()
  resiquan[y==0]<- rQ0
  resiquan[y==1]<- rQ1
  resiquan[which(y>0 & y<1)]<- rQ[which(y>0 & y<1)]
  return(resiquan)
}




sim = 100 # simulated envelope Monte Carlo (MC) replica number
e<- matrix(0,588,sim)
for(i in 1:sim) #loop of the MC
{
  e[,i]<-sort(rnorm(588,0,1))
}

e1<- numeric(n)
e2<- numeric(n)
for(i in 1:n){
  eo <- sort(e[i,])
  e1[i] <- (eo[2]+eo[3])/2
  e2[i] <- (eo[97]+eo[98])/2
}
#
med <- apply(e,1,mean)

faixa<-range(residuoQ,e1,e2)
qqnorm(residuoQ, pch = 19,ylim=faixa,main="(d)",xlab="Theoretical Quantiles",ylab="Randomized Quantile Residuals")
par(new=TRUE)
qqnorm(e1,axes=F,type="l",ylim=faixa,main=" ",xlab="",ylab="")
par(new=TRUE)
qqnorm(e2,axes=F,type="l",ylim=faixa,main=" ",xlab="",ylab="")
par(new=TRUE)
qqnorm(med,axes=F,type="l",lty=2,ylim=faixa,main=" ",xlab="",ylab="")





##Ajustando
dado=read.table("dados inf 100.txt")
x=dado$x;e=dado$e;z=dado$z;m=dado$m;fl=dado$w
phi=dado$phi;p0=dado$p0;p1=dado$p1;eta1=dado$eta1
tau=dado$tau;nu=dado$nu
alph=0.7
#eta2=1+2*simu2;phi=exp(eta2)
#eta1=-simu+d
n=100
alpha=rep(0,300)
betam=matrix(0,nrow=1,ncol=300)
kapp=rho=taum=matrix(0,nrow=2,ncol=300)
gama=matrix(0,nrow=n,ncol=300)
#logit
mu=exp(eta1)/(1+exp(eta1))
#probit
mu=pnorm(eta1)
#cloglog
mu=1-exp(-exp(eta1))
#loglog
mu=exp(-exp(-eta1))
#cauchit
mu=(1/pi)*(atan(eta1)+0.5*pi)

j=0
while(j<151){
  j=j+1
  y<-rBRINF(n, mu = mu, alpha=alph, phi = phi, p0 = p0, p1 = p1) 
  adjust<-try(samzoabr(y~-1+x+pb(z),formula.phi=~1+e,
                       formula.nu =~1+fl,formula.tau =~1+m,
                       mu.link="logit",iter=30,trace=T,
                       knots=50,lambda=100))
  if(inherits(adjust, "try-error"))
  {
    j = j-1
    next
  }
  betam[j]=adjust$bet
  kapp[,j]=adjust$kapp
  rho[,j]=adjust$rho
  taum[,j]=adjust$tau
  alpha[j]=adjust$alpha
  gama[,j]=adjust$mu.s
  cat("iteracao"," ",j, "\n")
}
dat=data.frame(beta1=betam[1,],kapp0=kapp[1,],
               kapp1=kapp[2,],rho0=rho[1,],
               rho1=rho[2,],taum0=taum[1,],
               taum1=taum[2,],alpha=alpha,gama=t(gama))
  
write.table(dat,file = "Simulacao zoabr logit 100.txt")
  

dataa=data.frame(mu=adjust1$mu.s,
                 g=factor(rep(1,500)),dd=z,d=cos(z))

dataa=data.frame(mu=vec(gama[,1:150]),g=factor(rep(c(1:150), each=100)),dd=z,d=cos(z))
ggplot(dataa,aes(x=dd,y=mu,group=g))+
  geom_line()+
  geom_line(aes(x = dd, y = d, group = "none", colour = "red"),size=1.2)

# while(any(y>0.9999999)){
#   y[which(y>0.9999999)]=rbeta(length(which(y>0.9999999)),
#                               delta[which(y>0.9999999)]*phi[which(y>0.9999999)],
#                               (1-delta[which(y>0.9999999)])*phi[which(y>0.9999999)])
# }
# while(any(y<1e-200 )){
#   y[which(y<1e-200)]=rbeta(length(which(y<1e-200)),
#                               delta[which(y<1e-200)]*phi[which(y<1e-200)],
#                               (1-delta[which(y<1e-200)])*phi[which(y<1e-200)])
# }