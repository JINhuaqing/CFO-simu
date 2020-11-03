rm(list=ls())
# setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
source("ORM_utils.R")
source("utilities.R")
library(magrittr)

make.move.fn <- function(ps, m=10){
    # ps: Output from move.dose.probs.fn 
    # m: number of samples to draw for the majority vote
    # output: 
    #    action: D--1, S--2, E--3
    nlevel <- length(ps)
    cps <- cumsum(ps)
    rvs <- runif(m)
    locs <- sapply(rvs, function(rv)rv<= cps)
    actions <- apply(locs, 2, which.max)
    if (m==1){
        final.action <- actions
    }else{
        res <- sapply(1:nlevel, function(i)actions==i)
        res <- colSums(res)
        final.action <- which.max(res)
    }
    final.action
}


under.eff.fn <- function(phiE, add.args=list())
    {
    args <- c(list(phi=phiE), add.args)
    x <- add.args$x
    n <- add.args$n
    alp.prior <- add.args$alp.prior.eff
    bet.prior <- add.args$bet.prior.eff
    ppE <- 1 - post.prob.fn(phiE, x, n, alp.prior, bet.prior)
    if ((ppE >= 0.95) & (n>3)){
        return(TRUE)
    }else{
        return(FALSE)
    }
}


# The fn to return the probs of the movement based on the efficacy in the admissible set
# In fact, we don't need so complicated form. 
# We can simply select the next dose via comparing x/n for different doses
ORM.Eff.move.probs <- function(ad.xs, ad.ns, alp.prior, bet.prior){
    alps <- ad.xs + alp.prior
    bets <- ad.ns - ad.xs + bet.prior
    nd <- length(ad.xs)
    
    Nsps <- 100000
    sps.list <- list() 
    for (i in 1:nd){
        sps.list[[i]] <- rbeta(Nsps, alps[i], bets[i])
    }
    
    spss <- do.call(rbind, sps.list)
    argMaxs <- apply(spss, 2, which.max)
    probs <- as.vector(table(argMaxs))/Nsps
    
    probs
}

# Simulation function for ORM-phase I/II
ORM.Eff.simu.fn <- function(phi, phiE, p.true, pE.true, ncohort=10, init.level=1,  cohortsize=3, add.args=list()){
    # phi: Target DIL rate
    # phiE: The minimal efficacy rate, only used for early stop
    # p.true: True DIL rates under the different dose levels
    # pE.true: True efficacy probs under the different dose levels
    # ncohort: The number of cohorts
    # cohortsize: The sample size in each cohort
    # alp.prior, bet.prior: prior parameters
    
    earlystop <- 0
    ndose <- length(p.true)
    cidx <- init.level
    
    tys <- rep(0, ndose) # number of DLT responses for different doses.
    txs <- rep(0, ndose) # number of efficacy responses for different doses.
    tns <- rep(0, ndose) # number of subject for different doses.
    tover.doses <- rep(0, ndose) # Whether each dose is too toxic or not, 1 yes.
    tunder.effs <- rep(0, ndose) # Whether the dose is not efficacious or not, 1 yes
    # if a dose is not efficacious enough or it is too toxic, it is would be eliminated from the admissible set.
    
    for (i in 1:ncohort){
        pc <- p.true[cidx] 
        pEc <- pE.true[cidx] 
        
        # sample from current dose
        cres <- rbinom(cohortsize, 1, pc)
        cEres <- rbinom(cohortsize, 1, pEc)
        
        # update results
        tys[cidx] <- tys[cidx] + sum(cres)
        txs[cidx] <- txs[cidx] + sum(cEres)
        tns[cidx] <- tns[cidx] + cohortsize
        
        #if (i == ncohort){
        #   break()
        #}
        
        cy <- tys[cidx]
        cx <- txs[cidx]
        cn <- tns[cidx]
        
        add.args <- c(list(y=cy, n=cn, x=cx, tys=tys, txs=txs, tns=tns, cidx=cidx), add.args)
        
        if (overdose.fn(phi, add.args)){
            tover.doses[cidx:ndose] <- 1
        }
        
        if (under.eff.fn(phiE, add.args)){
            tunder.effs[cidx] <- 1
        }
        
        
        if ((tover.doses[1] == 1) | (sum(tunder.effs)==ndose)){
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
        
        # The up.idx of the admissible set
        up.idx <- make.decision.ORM.fn(phi, cys, cns, add.args$alp.prior, add.args$bet.prior, cover.doses) - 2 + cidx
        if (up.idx == 1){
            cidx <- 1
        }else{
            if (cidx == 1){
                low.idx <- 1
            }else{
                low.idx <- cidx - 1
            }
            ad.xs <- txs[low.idx:up.idx]
            ad.ns <- tns[low.idx:up.idx]
            probs <- ORM.Eff.move.probs(ad.xs, ad.ns, add.args$alp.prior.eff, add.args$bet.prior.eff)
            #cidx <- which.max(probs)
            cidx <- make.move.fn(probs, m=10)
            
        }
        
    }
    
    
    if (earlystop==0){
        MTD <- select.mtd(phi, tns, tys)$MTD
        OBD.probs <- ORM.Eff.move.probs(txs[1:MTD], tns[1:MTD], add.args$alp.prior.eff, add.args$bet.prior.eff)
        OBD <- which.max(OBD.probs)
        
    }else{
        OBD <- 99
    }

    list(OBD=OBD, dose.ns=tns, eff.ns=txs, DLT.ns=tys, pE.true=pE.true, min.eff=phiE, 
         p.true=p.true, target=phi, over.doses=tover.doses, under.eff=tunder.effs)
    
}
    
    
    
phi <- 0.3
phiE <- 0.4
#p.true <- c(0.25, 0.35, 0.4, 0.45, 0.5)
#pE.true <- c(0.3, 0.35, 0.45, 0.65, 0.75)
#p.true <- c(0.1, 0.12, 0.15, 0.2, 0.25)
#pE.true <- c(0.25, 0.35, 0.6, 0.6, 0.6)
#p.true <- c(0.05, 0.1, 0.25, 0.3, 0.5)
#pE.true <- c(0.2, 0.4, 0.6, 0.70, 0.55)
p.true <- c(0.05, 0.07, 0.1, 0.15, 0.35)
pE.true <- c(0.1, 0.2, 0.35, 0.5, 0.55)

add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)


ress <- list()
for (i in 1:100){
    print(i)
    res <- ORM.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, ncohort=12);res
    ress[[i]] <- res
}


phase12.post.fn(ress)
