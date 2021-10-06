# this file is to compare CRM with prior variance 2 and 1.34, I use 1.34 
library(dfcrm)


posterior2 <- function(alpha, p, y, d) {
    sigma2 = 2 #use 2 before, but in most of papers, they adopt variance = 1.34
    lik=1;
    for(i in 1:length(y))
    {
        pi = p[d[i]]^(exp(alpha));
        lik = lik*pi^y[i]*(1-pi)^(1-y[i]);
    }
    return(lik*exp(-0.5*alpha*alpha/sigma2));
}

# used to calculate the posterior mean of pi
posttoxf2 <- function(alpha, p, y, d, j) { 
    p[j]^(exp(alpha))*posterior2(alpha, p, y, d); 
}

crm2.simu.fn <-function(target = 0.30, ## Target toxicity pr
              p.true, p.prior, init.level=1, cohortsize=1, ncohort=12, p.eli = 0.95){
    
    
    
    ndose <- length(p.true)
    if (missing(p.prior)){
        p.prior <- getprior(0.05, target, ceiling(ndose/2), ndose)
        # p.prior <- getprior(0.12, target, ceiling(ndose/2), ndose)
    }
    pts=rep(0,ndose);
    dlt=rep(0,ndose);
    pi.hat = numeric(ndose); # estimate of toxicity prob
    nstop = 0; # number of trial stopped due to high toxicity 
    t.start=Sys.time();
        
    y=NULL;  #binary outcome
    d=NULL;  #dose level
    dose.curr = init.level
    stop=0; #indicate if trial stops early
    for(i in 1:ncohort)
    {
        
        # generate data for the new patient
        y = c(y, rbinom(cohortsize, 1, p.true[dose.curr]));
        d = c(d, rep(dose.curr, cohortsize));
        
        # calculate posterior mean of toxicity probability at each dose leavel
        marginal=integrate(posterior2,lower=-Inf,upper=Inf,p.prior,y,d)$value
        for(j in 1:ndose) { pi.hat[j] = integrate(posttoxf2,lower=-Inf,upper=Inf,p.prior,y,d,j)$value/marginal;}
        
        # calculate pr(pi_1>target)
        p.overtox = integrate(posterior2,lower=-Inf,upper=log(log(target)/log(p.prior[1])),p.prior,y,d)$value/marginal;	
        if(p.overtox>p.eli) { stop=1; break;}
        
        diff = abs(pi.hat-target);
        dose.best =which.min(diff);
        if(dose.best>dose.curr && dose.curr != ndose) dose.curr = dose.curr+1;
        if(dose.best<dose.curr && dose.curr != 1) dose.curr = dose.curr-1;
    }
    if(stop==1) 
    {
        MTD <- 99 
    }
    else 
    { 
        MTD <- dose.best
    }

    for (j in 1:ndose){
        pts[j] <- sum(d==j)
        dlt[j] <- sum(y[d==j])
    }	
        
        
    dose.ns <- pts
    dlt.ns <- dlt
    list(MTD=MTD, dose.ns=dose.ns, DLT.ns=dlt.ns, p.true=p.true, target=target)
}
