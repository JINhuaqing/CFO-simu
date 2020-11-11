#setwd("C:/Users/Phase I/")

rm(list = ls(all = TRUE))
library(MCMCpack)
library(cubature)
library(e1071)

#true toxicity/efficacy probability
tprob1=c(0.02, 0.03, 0.05, 0.06, 0.07, 0.08)
eprob1=c(0.40, 0.25, 0.20, 0.15, 0.10, 0.05)

eprob2=c(0.25, 0.46, 0.30, 0.25, 0.16, 0.10)
tprob2=c(0.03, 0.10, 0.15, 0.20, 0.28, 0.35)

eprob3=c(0.10, 0.15, 0.35, 0.18, 0.12, 0.07)
tprob3=c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)

eprob4=c(0.05, 0.15, 0.30, 0.45, 0.35, 0.30)
tprob4=c(0.02, 0.04, 0.06, 0.10, 0.18, 0.40)

eprob5=c(0.15, 0.25, 0.33, 0.47, 0.60, 0.40)
tprob5=c(0.01, 0.02, 0.03, 0.05, 0.06, 0.07)


eprob6=c(0.05, 0.15, 0.30, 0.35, 0.40, 0.50)
tprob6=c(0.02, 0.05, 0.15, 0.20, 0.25, 0.40)


#default skeleton prob
prob6=c(0.05, 0.12, 0.19, 0.26, 0.33, 0.40)
prob5=c(0.12, 0.19, 0.26, 0.33, 0.40, 0.33)
prob4=c(0.19, 0.26, 0.33, 0.40, 0.33, 0.26)
prob3=c(0.26, 0.33, 0.40, 0.33, 0.26, 0.19)
prob2=c(0.33, 0.40, 0.33, 0.26, 0.19, 0.12)
prob1=c(0.40, 0.33, 0.26, 0.19, 0.12, 0.05)


response.tox=tprob1
response.eff=eprob1


mu_para_1=0;
mu_para_2=10;
tau_para_1=0.1;
tau_para_2=0.1;


csize=3; #cohort size
ncycle=10; #simulation times
total.s1=24; # sample size in stage 1
total.ss=60; # sample size in both stage 2
SKL=rep(0,length(response.eff));
Dose.tox=NULL;
Dose.eff=NULL;
Dose.npat=NULL; # num of sub for the all trial 
Dose.npat1=NULL;# num of sub for the trial of stage 1

skeleton=matrix(c(prob1,prob2,prob3,prob4,prob5,prob6),length(prob1),length(prob1))
OBD=rep(0, length(response.eff))
MTD=rep(0, length(response.tox))


Deci=rep(0,2)#record the decision choice: Deci[1]: parametric; Deci[2]: nonparametric
mtd.Deci=rep(0,2)
theta_U=0.7;
theta_L=0.4;
theta_T=0.95;
phi_T=0.3;

## pava is the pool-adjacent-violator algorithm to perform isotonic transformation for the posterior means
pava <- function (x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1) 
        return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
        stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol))) 
          break
        i <- min((1:(n - 1))[viol])
        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i + 1]
        ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
        x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
        lvlsets[ilvl] <- lvl1
    }
    x
}


## betavar computes variances of beta distributions 
betavar<-function(a,b){
    variance <- a*b/((a+b)^2*(a+b+1))
    return(variance)
}

## posttoxf1 computes the inner part for integration of posterior mean based on probit model
posttoxf1<- function(mu,tau,dose.tox,dose.npat,dose.level,j){ 
	    return(pnorm(dose.level[j], mu, sqrt(tau))*posterior1(mu,tau,dose.tox,dose.npat,dose.level))
}


## posterior1 computes posterior function without normalizing based on probit model
posterior1<- function(mu,tau,dose.tox,dose.npat,dose.level){ 
    lik=1;
	  for(i in 1:length(dose.level))
	  {
	  	  pi = pnorm(dose.level[i], mu, sqrt(tau))
	  	  lik = lik*pi^dose.tox[i]*(1-pi)^(dose.npat[i]-dose.tox[i]);
	  }
	  return(lik*dnorm(mu,mu_para_1,sqrt(mu_para_2))*dgamma(tau,tau_para_1,tau_para_2));
}

posterior <- function(alpha, p, y, d){
    sigma2 = 2;
    lik=1;
    for(i in 1:length(y))
    {
        pi = p[d[i]]^(exp(alpha));
        lik = lik*pi^y[i]*(1-pi)^(1-y[i]);
    }
    return(lik*exp(-0.5*alpha*alpha/sigma2));
}


## loglikepmd computes the posterior mean deviance
loglikepmd<- function(mu,tau,dose.tox,dose.npat,dose.level,marg){ 
    lik=1;
	  for(i in 1:length(dose.level))
	  {
	  	  pi = pnorm(dose.level[i], mu, sqrt(tau));
	  	  pi=max(0.01,min(0.99,pi));
	  	  lik = lik*pi^dose.tox[i]*(1-pi)^(dose.npat[i]-dose.tox[i]);
	  }
       
	  return(-2*log(lik)*lik*dnorm(mu,mu_para_1,sqrt(mu_para_2))*dgamma(tau,tau_para_1,tau_para_2)/marg);
}


## loglikedpm computes the deviance of posterior mean
loglikedpm<- function(postmean.mu,postmean.tau,dose.tox,dose.npat,dose.level){ 
    lik=1;
    for(i in 1:length(dose.level))
    {
    	  pi = pnorm(dose.level[i], postmean.mu, sqrt(postmean.tau));
    	  lik = lik*pi^dose.tox[i]*(1-pi)^(dose.npat[i]-dose.tox[i]);
    }
    return(-2*log(lik));
}


## postmeanmu computes the inner part for intergration of the posterior mean of mu
postmeanmu<- function(mu,tau,dose.tox,dose.npat,dose.level){ 
    lik=1;
	  for(i in 1:length(dose.level))
	  {
	  	  pi = pnorm(dose.level[i], mu, sqrt(tau))
	  	  lik = lik*pi^dose.tox[i]*(1-pi)^(dose.npat[i]-dose.tox[i]);
	  }
	  return(mu*lik*dnorm(mu,mu_para_1,sqrt(mu_para_2))*dgamma(tau,tau_para_1,tau_para_2));
}


## postmeantau computes the inner part for intergration of the posterior mean of tau
postmeantau<- function(mu,tau,dose.tox,dose.npat,dose.level){ 
     lik=1;
	for(i in 1:length(dose.level))
	{
		pi = pnorm(dose.level[i], mu, sqrt(tau))
		lik = lik*pi^dose.tox[i]*(1-pi)^(dose.npat[i]-dose.tox[i]);
	}
	return(tau*lik*dnorm(mu,mu_para_1,sqrt(mu_para_2))*dgamma(tau,tau_para_1,tau_para_2));
}


posttoxf <- function(alpha, p, y, d, j) { p[j]^(exp(alpha))*posterior(alpha, p, y, d); }


logbeta<-function(pi,alpha0, beta0, dose.tox, dose.npat,j){
  lik=pi^dose.tox[j]*(1-pi)^(dose.npat[j]-dose.tox[j]);
  return(-2*log(lik)*dbeta(pi,alpha0[j]+dose.tox[j],dose.npat[j]-dose.tox[j]+beta0[j]));

}


prior1<- function(mu,tau,startdose){ 
	return(pnorm(startdose, mu, sqrt(tau))*dnorm(mu,mu_para_1,sqrt(mu_para_2))*dgamma(tau,tau_para_1,tau_para_2))
}


priormean <- function(x){
marginal=integrate(function(tau)
              {
                sapply(tau,function(tau)
                       {
                         integrate(function(mu) prior1(mu,tau,x), -Inf, Inf)$value
                       })
              }, 0, Inf)$value;
return(marginal)
}


post <- function(alpha,pp,dose.eff,dose.npat) {
  like=1;
  for (ii in 1: length(dose.npat)) {
    like1=(pp[ii]^exp(alpha))^dose.eff[ii]*(1-pp[ii]^exp(alpha))^(dose.npat[ii]-dose.eff[ii]);
    like=like*like1;
  }
  return(like*dnorm(alpha,0,sqrt(2)))
}


postmean <- function(alpha,pp,dose.eff,dose.npat,j,normcons1){
  like=1;
  for (ii in 1: length(dose.npat)) {
    like1=(pp[ii]^exp(alpha))^dose.eff[ii]*(1-pp[ii]^exp(alpha))^(dose.npat[ii]-dose.eff[ii]);
    like=like*like1; 
  }
  return(pp[j]^exp(alpha)*like*dnorm(alpha,0,sqrt(2))/normcons1)
}


#===================pre-processing===================
### Compute the effective sample size for beta prior. (d_j and prior variance for probit model)

prior_mean<-seq(0.1,0.1*length(response.tox),0.1);

prior2<- function(mu,tau,dose.level,prior_mean,i)
{ 
	return((pnorm(dose.level[i], mu, sqrt(tau))-prior_mean[i])^2*dnorm(mu,mu_para_1,sqrt(mu_para_2))*dgamma(tau,tau_para_1,tau_para_2))
}


#compute d_j
dose.level=0*prior_mean;
for (i in 1:length(response.tox)){
dose.level[i]=uniroot(function(x) {priormean(x)-prior_mean[i]},c(-20,20),tol=0.01)$root
}


priorvar=0*prior_mean;
for (i in 1:length(response.tox)){
priorvar[i] <-integrate(function(tau)
              {
                sapply(tau,function(tau)
                       {
                         integrate(function(mu) prior2(mu,tau,dose.level, prior_mean,i), -Inf, Inf)$value
                       })
              }, 0, Inf)$value;
}


var_p<-priorvar;
mean_p<-prior_mean;
alpha0=mean_p^2*(1-mean_p)/(var_p)-mean_p;
beta0=alpha0*(1-mean_p)/mean_p;
alpha0
beta0


#===============Main program===============
STOP=0;

for (simu in 1:ncycle) {
    print(simu)
    
    ##############initial values################
    dose.up=length(dose.level);
    dose.curr=1;##start from the first dose level
    dose=NULL;
    y=NULL;
    dose.tox=rep(0,length(dose.level));  #number of toxicities at each dose level
    dose.npat=rep(0,length(dose.level)); #number of patients at each dose level
    dose.eff=rep(0,length(dose.level));  #number of patients at each dose level
    
    stop1=0;
    
    
    #==============Stage I==============
    while ((sum(dose.npat)<total.s1)&&(stop1==0)){
    
      #dose=c(dose,rep(dose.level[dose.curr], csize));
      tox.curr=rbinom(csize, 1, response.tox[dose.curr]);
      dose.tox[dose.curr]=dose.tox[dose.curr]+sum(tox.curr);
      dose.npat[dose.curr]=dose.npat[dose.curr]+csize;
      y= c(y, tox.curr)
      dose = c(dose, rep(dose.curr, csize))
      
      ##use nonparametric model for decision###
      
      pp.np=1-pbeta(phi_T,alpha0[dose.curr]+dose.tox[dose.curr],beta0[dose.curr]+dose.npat[dose.curr]-dose.tox[dose.curr])
      
      if (pp.np>theta_U) {decision.np=-1};           
      if (pp.np<theta_L) {decision.np=1}; 
      if ((pp.np<=theta_U) && (pp.np>=theta_L)) {decision.np=0};
      
      
      normcons=integrate(function(tau)
                    {
                      sapply(tau,function(tau)
                             {
                               integrate(function(mu) posterior1(mu,tau,dose.tox,dose.npat,dose.level), -Inf, Inf)$value
                             })
                    }, 0, Inf)$value;
      
      j=dose.curr;
      pp.p=integrate(function(tau)
       	             {
          	             sapply(tau,function(tau)
            	                 {
            	                   integrate(function(mu) posterior1(mu,tau,dose.tox,dose.npat,dose.level), -Inf, dose.level[j]-sqrt(tau)*qnorm(phi_T))$value
            	                 })
            	        }, 0, Inf)$value/normcons;
    
      ##stopping rules
      
      if((pp.p>=theta_T)||(pp.np>=theta_T)){dose.up=dose.curr-1;}
      if(dose.up<1){stop1=1;break}
      
      if (pp.p>theta_U) {decision.p=-1};           
      if (pp.p<theta_L) {decision.p=1}; 
      if ((pp.p<=theta_U) && (pp.p>=theta_L)) {decision.p=0};
      
      
      ##make final decision for toxicity using the posterior model probability##
      ##posterior model probability for parametric model##
    
      normcons.p=integrate(function(tau)
                    {
                      sapply(tau,function(tau)
                             {
                               integrate(function(mu) posterior1(mu,tau,dose.tox,dose.npat,dose.level), -Inf, Inf)$value
                             })
                    }, 0, Inf)$value;
      
      ##posterior model probability for nonparametric model##
      
      normcons.np=1;
      for (j in 1:length(dose.level)){
      priorcons=gamma(alpha0[j]+beta0[j])/(gamma(alpha0[j])*gamma(beta0[j]));
      normcons.np=normcons.np*gamma(alpha0[j]+dose.tox[j])*gamma(beta0[j]+dose.npat[j]-dose.tox[j])/gamma(alpha0[j]+beta0[j]+dose.npat[j])*priorcons;
      #choose(dose.npat[j],dose.tox[j])
      }
      #print(c(normcons.np,normcons.p))
      
      if (normcons.p>normcons.np){decision<-decision.p; Deci[1]=Deci[1]+1}
      else {decision<-decision.np;Deci[2]=Deci[2]+1}
      #print(c(normcons.np,normcons.p))
      dose.curr=max(min(dose.curr+decision,dose.up),1);
      
    }#while
  
  
    if(stop1==0){
    
    normcons.p=integrate(function(tau)
                  {
                    sapply(tau,function(tau)
                           {
                             integrate(function(mu) posterior1(mu,tau,dose.tox,dose.npat,dose.level), -Inf, Inf)$value
                           })
                  }, 0, Inf)$value;
    
    ##posterior model probability for nonparametric model
    
    normcons.np=1;
    for (j in 1:length(dose.level)){
    priorcons=gamma(alpha0[j]+beta0[j])/(gamma(alpha0[j])*gamma(beta0[j]));
    normcons.np=normcons.np*gamma(alpha0[j]+dose.tox[j])*gamma(beta0[j]+dose.npat[j]-dose.tox[j])/gamma(alpha0[j]+beta0[j]+dose.npat[j])*priorcons;
    #choose(dose.npat[j],dose.tox[j])
    }
    #print(c(normcons.np,normcons.p))
    
    
    if (normcons.p>=normcons.np){
    pi_hat=rep(0,length(response.tox));
    
    for(j in 1:length(response.tox)) { 
       pi_hat[j]=integrate(function(tau)
     	             {
        	             sapply(tau,function(tau)
          	                 {
          	                   integrate(function(mu) posttoxf1(mu,tau,dose.tox,dose.npat,dose.level,j), -Inf, Inf)$value
          	                 })
          	        }, 0, Inf)$value/normcons.p;
                     }
    
    dose.star=max(which(abs(pi_hat-phi_T)==min(abs(pi_hat-phi_T))));
    
    mtd.Deci[1]=mtd.Deci[1]+1;
    MTD[dose.star]=MTD[dose.star]+1;}
    else {
    dose.mean=(dose.tox+alpha0)/(dose.npat+alpha0+beta0);
    dose.pava=pava(dose.mean);
    dose.star=max(which.min(abs(dose.pava-phi_T)));
    MTD[dose.star]=MTD[dose.star]+1;
    mtd.Deci[2]=mtd.Deci[2]+1;
    }#else
    
    }
    #Dose.npat=Dose.npat+dose.npat;
    STOP=STOP+stop1;
    Dose.npat1=rbind(Dose.npat1,dose.npat);
    
    dose.curr = dose.star

    #==============Stage II==============
    while ((sum(dose.npat)<total.ss) & stop1==0){
      eff.curr=rbinom(csize, 1, response.eff[dose.curr]);
      tox.curr=rbinom(csize, 1, response.tox[dose.curr]);
      dose.eff[dose.curr]=dose.eff[dose.curr]+sum(eff.curr);
      
      dose.tox[dose.curr]=dose.tox[dose.curr]+sum(tox.curr);
      dose.npat[dose.curr]=dose.npat[dose.curr]+csize;
      y = c(y, tox.curr)
      dose = c(dose, rep(dose.curr, csize))
     
      
      ##monitor toxicity
      pi.hat=rep(0, length(dose.eff))
      
      normcons=integrate(function(tau)
      {
        sapply(tau,function(tau)
        {
          integrate(function(mu) posterior1(mu,tau,dose.tox,dose.npat,dose.level), -Inf, Inf)$value
        })
      }, 0, Inf)$value;
      
      for(j in 1:length(dose.eff)) { 
        pi.hat[j]=integrate(function(tau)
        {
          sapply(tau,function(tau)
          {
            integrate(function(mu) posterior1(mu,tau,dose.tox,dose.npat,dose.level), -Inf, dose.level[j]-sqrt(tau)*qnorm(phi_T))$value
          })
        }, 0, Inf)$value/normcons;
    
      }
      
      
      diff = abs(pi.hat-phi_T);
      
      adm.dose = min(which(diff==min(diff)))
      
      
      ##use CRM for six skeletons
      Mk=rep(1/adm.dose,adm.dose);
      post.prob=rep(0,adm.dose);
      post.mean=rep(0,adm.dose);
      for (skel in 1:adm.dose){
        pp=skeleton[skel,];
        #like=prod((pp^exp(alpha))^dose.eff)*prod((1-pp^exp(alpha))^(dose.npat-dose.eff));
        post.prob[skel]=integrate(post,lower=-Inf,upper=Inf,pp,dose.eff,dose.npat)$value;
      }
      
      Mk.post=post.prob*Mk/sum(post.prob*Mk)
      #skle.sele=max(which.max(Mk.post));
      skle.sele=rdiscrete(1, Mk.post, values = 1:length(Mk.post))
      SKL[skle.sele]=SKL[skle.sele]+1;
      
      dose.curr=skle.sele;
    }#while
    
    if (stop1==0){
      dose.obd=which.max(Mk.post)
      OBD[dose.obd]=OBD[dose.obd]+1;
      
    }
    Dose.tox=rbind(Dose.tox, dose.tox);
    Dose.eff=rbind(Dose.eff,dose.eff);
    Dose.npat=rbind(Dose.npat,dose.npat);

}#for (simu in 1:ncycle)


#=============output=============

Selection <- 100*OBD/ncycle
Allocation <- colMeans(Dose.npat)
effs.nums <- colSums(Dose.eff)/ncycle
toxs.nums <- colSums(Dose.tox)/ncycle
tol.effs <- sum(effs.nums)
tol.toxs <- sum(toxs.nums)
tol.Subjs <- sum(Allocation)
errStop <- 1 - sum(OBD)/ncycle

Results <- list(
  Selection=Selection, 
  Allocation=Allocation,
  effs.nums=effs.nums, 
  toxs.nums=toxs.nums, 
  tol.effs=tol.effs,
  tol.toxs=tol.toxs,
  tol.Subjs=tol.Subjs,
  errStop=errStop)

i <- 1
fname <- paste0("AdaModel", i, ".RData")
save(Results, file=fname)
load("test.RData")
