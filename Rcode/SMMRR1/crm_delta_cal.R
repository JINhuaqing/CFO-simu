#setwd("C:/Users/JINHU/Documents/ProjectCode/CFO")
setwd("/home/r6user2/MyResearch/CFO/")
library(dfcrm)
library(parallel)
source("Rcode/phaseI/crm_utils.R")
source("Rcode/utilities.R")

gen.cal.Set.fn <- function(phi, psi, ndose){
    #args:
    #   phi: the target DLT rate
    #   psi: parameter to control
    #   ndose: number of dose levels
    pU <- phi*psi/(1-phi*(1-psi))
    pL <- phi/(psi+phi*(1-psi))
    
    sets <- list()
    for (k in 1:ndose){
        sc <- rep(pL, ndose)
        sc[k] <- phi
        if (k < ndose){
            sc[(k+1):ndose] <- pU
        }
        sets[[k]] <- sc
        
    }
    
    return(sets)
}


PCS.fn <- function(ress){
    tMTD <- MTD.level(ress[[1]]$target, ress[[1]]$p.true)
    pcs.rate <- mean(sapply(ress, function(res) res$MTD==tMTD))
    pcs.rate
}

crm.cal.single <- function(delta, phi, psi, ndose, nu, ncohort, cohortsize=3){
    #args:
    #   delta: indiff interval
    #   phi: the target DLT rate
    #   psi: parameter to control
    #   ndose: number of dose levels
    #   nu: prior MTD
    #   ncohort: number of cohorts
    #   cohortsize: cohort size
   p.prior <- getprior(delta, phi, nu, ndose) 
   cal.Set <- gen.cal.Set.fn(phi, psi, ndose)
   pcss <- c()
   nsim <- 1000
   for (k in 1:ndose){
       # for linux 
       run.fn <- function(i){
           res <- crm.simu.fn(target=phi, p.true=cal.Set[[k]], p.prior=p.prior, cohortsize=cohortsize, ncohort=ncohort)
           res
       }
       ress <- mclapply(1:nsim, run.fn, mc.cores=40)
       # for windows
       # ress <- list()
       # for (i in 1:nsim){
       #     res <- crm.simu.fn(target=phi, p.true=cal.Set[[k]], p.prior=p.prior, 
       #                    cohortsize=cohortsize, ncohort=ncohort)
       #     ress[[i]] <- res
       # }
       pcss <- c(pcss, PCS.fn(ress))
       
   }
   
   return(pcss)
}

crm.cal <- function(deltas, phi, psi, ndose, nu, ncohort, cohortsize=3){
    #args:
    #   deltas: indiff intervals
    #   phi: the target DLT rate
    #   psi: parameter to control
    #   ndose: number of dose levels
    #   nu: prior MTD
    #   ncohort: number of cohorts
    #   cohortsize: cohort size
    
    a.pcss <- c()
    for (delta in deltas){
           pcss <-  crm.cal.single(delta, phi, psi, ndose, nu, ncohort, cohortsize)
           a.pcs <- mean(pcss)
           a.pcss <- c(a.pcss, a.pcs)
        
    }
    
    a.pcss
}

phi <- 0.3
deltas <- c(0.01, phi-0.05, 0.01)
psi <- 2
ndose <- 5
ncohort <- 10

for (ndose in c(3, 5, 7, 9)){
nu <- ceiling(ndose/2)
res <- crm.cal(deltas, phi, psi, ndose, nu, ncohort)
f.Name <- paste0("./results/SMMR-R1/CRM_cal_phi", phi*100, "_ndose_", ndose, "_ncohort_", ncohort, ".RData")
save(res, file=f.Name)
}
