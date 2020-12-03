library(Iso)
ns <- c(3, 5, 24, 25)
ys <- c(0, 0, 4, 4)
xs <- c(0, 0, 2, 0)

pE <- (ys+0.1)/(ns+0.2);pE

pT <- (xs+0.1)/(ns+0.2);pT
pT <- pava(pT, ns+0.2);pT

source("ORM_utils.R")
source("./phaseI_II/ORM_eff_utils.R")


# Simulation function for ORM-phase I/II
ORM.Eff.simu.fn.realdata <- function(phi, phiE, p.true, pE.true, ncohort=10, init.level=1,  cohortsize=3, add.args=list(),
                            ph1=0){
    # phi: Target DIL rate
    # phiE: The minimal efficacy rate, only used for early stop
    # p.true: True DIL rates under the different dose levels
    # pE.true: True efficacy probs under the different dose levels
    # ncohort: The number of cohorts
    # cohortsize: The sample size in each cohort
    # alp.prior, bet.prior: prior parameters
    
    dose.seq <- c()
    upidx.seq <- c()
    tox.seq <- NULL
    eff.seq <- NULL
    ttox.seq <- NULL
    teff.seq <- NULL
    tover.seq <- NULL
    tunder.seq <- NULL
    
    tadd.args <- add.args
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
        dose.seq <- c(dose.seq, cidx)
        
        pc <- p.true[cidx] 
        pEc <- pE.true[cidx] 
        
        # sample from current dose
        cres <- rbinom(cohortsize, 1, pc)
        cEres <- rbinom(cohortsize, 1, pEc)
        tox.seq <- rbind(tox.seq, cres)
        eff.seq <- rbind(eff.seq, cEres)
        
        # update results
        tys[cidx] <- tys[cidx] + sum(cres)
        txs[cidx] <- txs[cidx] + sum(cEres)
        tns[cidx] <- tns[cidx] + cohortsize
        
        ttox.seq <- rbind(ttox.seq, tys)
        teff.seq <- rbind(teff.seq, txs)
        #if (i == ncohort){
        #   break()
        #}
        
        cy <- tys[cidx]
        cx <- txs[cidx]
        cn <- tns[cidx]
        
        add.args <- c(list(y=cy, n=cn, x=cx, tys=tys, txs=txs, tns=tns, cidx=cidx), tadd.args)
        
        if (overdose.fn(phi, add.args)){
            tover.doses[cidx:ndose] <- 1
        }
        
        if (under.eff.fn(phiE, add.args)){
            tunder.effs[cidx] <- 1
        }else{
            tunder.effs[cidx] <- 0
        }
        
        tover.seq <- rbind(tover.seq, tover.doses)
        tunder.seq <- rbind(tunder.seq, tunder.effs)
        
        if ((tover.doses[1] == 1) | (sum(tunder.effs[tover.doses==0])==sum(tover.doses==0))){
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
        upidx.seq <- c(upidx.seq, up.idx)
        if (i<=ph1){
            cidx <- up.idx
        }else{
            
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
                if (length(ad.xs)==1){
                    cidx <- low.idx 
                }else{
                    probs <- ORM.Eff.move.probs(ad.xs, ad.ns, add.args$alp.prior.eff, add.args$bet.prior.eff)
                    cidx <- which.max(probs) + low.idx -1
                    #cidx <- make.move.fn(probs, m=1) + low.idx - 1
                }
                
            }
        }
        
    }
    
    
    if (earlystop==0){
        MTD <- select.mtd(phi, tns, tys)$MTD
        if ( (MTD==99) | (sum(tunder.effs[1:MTD]) == MTD)){
            OBD <- 99
        }else{
            OBD.probs <- ORM.Eff.move.probs(txs[1:MTD], tns[1:MTD], add.args$alp.prior.eff, add.args$bet.prior.eff)
            OBD <- which.max(OBD.probs)
        }
        
    }else{
        OBD <- 99
    }
    
    list(OBD=OBD, dose.ns=tns, eff.ns=txs, DLT.ns=tys, pE.true=pE.true, min.eff=phiE, 
         p.true=p.true, target=phi, over.doses=tover.doses, under.eff=tunder.effs, 
         dose.seq = dose.seq ,
         tox.seq = tox.seq ,
         eff.seq = eff.seq ,
         ttox.seq = ttox.seq ,
         teff.seq = teff.seq ,
         tover.seq = tover.seq ,
         tunder.seq = tunder.seq,
         upidx.seq = upidx.seq
    )
    
}

phi <- 0.2
phiE <- 0.20
add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)
res <- ORM.Eff.simu.fn.realdata(phi=phi, phiE=0.2, p.true=pT, pE.true=pE, ncohort=19, init.level=1,  
                         cohortsize=3, add.args=add.args, ph1=0)
    

res$OBD
res$tunder.seq
res$tover.seq
res$dose.seq
res$tox.seq
res$eff.seq
res$ttox.seq
