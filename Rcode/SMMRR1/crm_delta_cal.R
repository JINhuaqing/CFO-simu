setwd("C:/Users/JINHU/Documents/ProjectCode/CFO")
library(dfcrm)
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
       ress <- list()
       for (i in 1:nsim){
           res <- crm.simu.fn(target=phi, p.true=cal.Set[[1]], p.prior=p.prior, 
                          cohortsize=cohortsize, ncohort=ncohort)
           ress[[i]] <- res
       }
       pcss <- c(PCS.fn(ress), pcss)
       
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
           a.pcss <- c(a.pcs, a.pcss)
        
    }
    
    a.pcss


}