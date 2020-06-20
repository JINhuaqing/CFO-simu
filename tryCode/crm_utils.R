library(dfcrm)


posterior <- function(alpha, p, y, d) {
    sigma2 = 2;
    lik=1;
    for(i in 1:length(y))
    {
        pi = p[d[i]]^(exp(alpha));
        lik = lik*pi^y[i]*(1-pi)^(1-y[i]);
    }
    return(lik*exp(-0.5*alpha*alpha/sigma2));
}

# used to calculate the posterior mean of pi
posttoxf <- function(alpha, p, y, d, j) { p[j]^(exp(alpha))*posterior(alpha, p, y, d); }

crm.simu.fn <-function(target = 0.30, ## Target toxicity pr
              p.true, cohortsize=1, ncohort=12, p.eli = 0.95, ntrial=1){
    
    
    
    ndose <- length(p.true)
    p.prior <- getprior(0.05, target, ceiling(ndose/2), ndose)
    ndose = length(p.prior);
    sel=rep(0,ndose);
    pts=rep(0,ndose);
    dlt=rep(0,ndose);
    tox=0;
    pi.hat = numeric(ndose); # estimate of toxicity prob
    dose.select=numeric(ndose);  # dose selection
    ntox=rep(0, ndose);  # number of toxicity at each dose level
    ntrted=rep(0, ndose); # number of patient at each dose level
    nstop = 0; # number of trial stopped due to high toxicity 
    t.start=Sys.time();
    d.mtd=which.min(abs(p.true-target));  
    for(trial in 1:ntrial)
    {
        
        y=NULL;  #binary outcome
        d=NULL;  #dose level
        dose.curr = ceiling(ndose/2);  # current dose level
        stop=0; #indicate if trial stops early
        for(i in 1:ncohort)
        {
            
            # generate data for the new patient
            y = c(y, rbinom(cohortsize, 1, p.true[dose.curr]));
            d = c(d, rep(dose.curr, cohortsize));
            
            # calculate posterior mean of toxicity probability at each dose leavel
            marginal=integrate(posterior,lower=-Inf,upper=Inf,p.prior,y,d)$value
            for(j in 1:ndose) { pi.hat[j] = integrate(posttoxf,lower=-Inf,upper=Inf,p.prior,y,d,j)$value/marginal;}
            
            # calculate pr(pi_1>target)
            p.overtox = integrate(posterior,lower=-Inf,upper=log(log(target)/log(p.prior[1])),p.prior,y,d)$value/marginal;	
            if(p.overtox>p.eli) { stop=1; break;}
            
            diff = abs(pi.hat-target);
            dose.best = min(which(diff==min(diff)));
            if(dose.best>dose.curr && dose.curr != ndose) dose.curr = dose.curr+1;
            if(dose.best<dose.curr && dose.curr != 1) dose.curr = dose.curr-1;
        }
        if(stop==1) { nstop=nstop+1; }
        else { dose.select[dose.best] = dose.select[dose.best]+1; 
        sel[dose.best]=sel[dose.best]+1/ntrial*100}
        if(sum(y)/length(d)>target){tox=tox+1/ntrial*100}
        for (j in 1:ndose){
            pts[j]<-pts[j]+sum(d==j)/ntrial
            dlt[j]<-dlt[j]+sum(y[d==j])/ntrial
        }	
        
        
    }
    pts<-round(pts,1)
    dlt<-round(dlt,1)
    MTD <- which.max(sel)
    dose.ns <- pts
    dlt.ns <- dlt
    list(MTD=MTD, dose.ns=dose.ns, DLT.ns=dlt.ns, p.true=p.true, target=target)
}
