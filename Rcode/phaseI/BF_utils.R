library(TruncatedDistributions)
library(BOIN)
source("utilities.R")

# truncated Beta sampler
tbeta.sampler.low <- function(pcs, alp, bet, pU=0.9){
    n <- length(pcs)
    rtbeta(n, alp, bet, pcs, pU)
}

tbeta.sampler.up <- function(pcs, alp, bet, pL=0){
    n <- length(pcs)
    rtbeta(n, alp, bet, pL, pcs)
}
# Function to make a decision based on the 
# Pr(\theta_R>\phi|D)/Pr(\theta_C>\phi|D) and 
# Pr(\theta_L<\phi|D)/Pr(\theta_C<\phi|D).
make.decision.BF.fn <- function(phi, ys, ns, alp.prior, bet.prior, over.doses, diag=FALSE){
   if (over.doses[2]==1){
       if (diag){
           rev <- list(final.action=1,
                       p1.sps=NULL,
                       p2.sps=NULL,
                       p3.sps=NULL,
                       BFs=c(NULL, NULL))
           return(rev)
       }else{
           return(1)
       }
   }else{
       alps <- ys + alp.prior
       bets <- ns - ys + bet.prior
       p2.under <- pbeta(phi, alps[2], bets[2])
       p2.over <- 1 - p2.under

       p2.sps <- rbeta(10000, alps[2], bets[2])
       p1.sps <- NULL
       p3.sps <- NULL
       if (!is.na(ys[1])){
           p1.sps <- tbeta.sampler.up(p2.sps, alps[1], bets[1])
           p1.under <- mean(p1.sps<=phi)
           p.under.BF <- p1.under/p2.under
       }else{
           p.under.BF <- 0.1
       }
       if (!(is.na(ys[3]) | over.doses[3]==1) ){
           p3.sps <- tbeta.sampler.low(p2.sps, alps[3], bets[3])
           p3.over <- mean(p3.sps>=phi)
           p.over.BF <- p3.over/p2.over
       }else{
           p.over.BF <- 0.1
       }
       CV1 <- 0.5
       CV2 <- 0.5
       dec.over <- (log10(p.over.BF)>CV1) 
       dec.under <-  (log10(p.under.BF)>CV2)
       if (dec.over & (!dec.under)){
           final.action <- 3
       }else if (dec.under & (!dec.over)){
           final.action <- 1
       }else{
           final.action <- 2
       }
       if (diag){
           rev <- list(final.action=final.action, 
                       p1.sps=p1.sps, 
                       p2.sps=p2.sps,
                       p3.sps=p3.sps,
                       BFs=c(p.under.BF, p.over.BF))
           return(rev)
       }else{
           return(final.action)
       }
   }
}

overdose.fn <- function(phi, add.args=list()){
    args <- c(list(phi=phi), add.args)
    y <- add.args$y
    n <- add.args$n
    alp.prior <- add.args$alp.prior
    bet.prior <- add.args$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    if ((pp >= 0.95) & (add.args$n>=3)){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

# Simulation function for BF method
BF.simu.fn <- function(phi, p.true, ncohort=12,
                              cohortsize=1, add.args=list()){
    # phi: Target DIL rate
    # p.true: True DIL rates under the different dose levels
    # ncohort: The number of cohorts
    # cohortsize: The sample size in each cohort
    # alp.prior, bet.prior: prior parameters
    # add.args, list of argments. 
    #    list(alp.prior=, bet.prior)
    earlystop <- 0
    ndose <- length(p.true)
    cidx <- ceiling(ndose/2)
    
    tys <- rep(0, ndose) # number of responses for different doses.
    tns <- rep(0, ndose) # number of subject for different doses.
    tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
    
    
    
    
    for (i in 1:ncohort){
        pc <- p.true[cidx] 
        
        # sample from current dose
        cres <- rbinom(cohortsize, 1, pc)
        
        # update results
        tys[cidx] <- tys[cidx] + sum(cres)
        tns[cidx] <- tns[cidx] + cohortsize
        
        
        
        cy <- tys[cidx]
        cn <- tns[cidx]
        
        add.args.od <- c(list(y=cy, n=cn, tys=tys, tns=tns, cidx=cidx), add.args)
        
        
        if (overdose.fn(phi, add.args.od)){
            tover.doses[cidx:ndose] <- 1
        }
        
        if (tover.doses[1] == 1){
            earlystop <- 1
            break()
        }
        
        
        # the results for current 3 dose levels
        if (cidx!=1){
            cys <- tys[(cidx-1):(cidx+1)]
            cns <- tns[(cidx-1):(cidx+1)]
            cover.doses <- tover.doses[(cidx-1):(cidx+1)]
        }else{
            cys <- c(NA, tys[1:(cidx+1)])
            cns <- c(NA, tns[1:(cidx+1)])
            cover.doses <- c(NA, tover.doses[1:(cidx+1)])
        }
        
        alp.prior <- add.args$alp.prior
        bet.prior <- add.args$bet.prior
        idx.chg <- make.decision.BF.fn(phi, cys, cns, alp.prior, bet.prior, cover.doses) - 2
        
        
        cidx <- idx.chg + cidx
        
    }
    
    
    if (earlystop==0){
        MTD <- select.mtd(phi, tns, tys)$MTD
    }else{
        MTD <- 99
    }
    list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi)
}

#BF.simu.fn(0.3, c(0.1, 0.3, 0.4), add.args = list(alp.prior=1, bet.prior=1))
