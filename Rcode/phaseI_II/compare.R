#rm(list=ls())
#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
#set.seed(0);
source("ORM_eff_utils.R")
source("STEIN-utils.R")
source("../utilities.R")
library(parallel)

phi <- 0.30
cohortsize = 3 # cohort size
ncohort = 10 # number of cohorts
nsimu <- 1000

#p.true <- c(0.1, 0.12, 0.15, 0.2, 0.25) #1
#pE.true <- c(0.25, 0.60, 0.4, 0.3, 0.2)

p.true <- c(0.25, 0.40, 0.45, 0.5, 0.55) #2
pE.true <- c(0.3, 0.35, 0.45, 0.65, 0.75)

#p.true <- c(0.05, 0.07, 0.1, 0.30, 0.40) # 3
#pE.true <- c(0.1, 0.2, 0.40, 0.5, 0.70)

p.true <- c(0.02, 0.05, 0.07, 0.1, 0.15) # 4
pE.true <- c(0.05, 0.08, 0.15, 0.30, 0.45)

p.true <- c(0.10, 0.13, 0.15, 0.3, 0.45, 0.5) # 5
pE.true <- c(0.30, 0.40, 0.50, 0.50, 0.50, 0.50)

#p.true <- c(0.10, 0.13, 0.25, 0.3, 0.35, 0.5) # 6
#pE.true <- c(0.30, 0.40, 0.45, 0.50, 0.60, 0.70)

psi1<-0.35 # lowest acceptable efficacy rate
psi2<-0.65 # very promising efficacy rate that leads to dose retainment

phiE <- 0.5
add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)

run.fn1 <- function(k){
    print(k)
    res <- ORM.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, cohortsize=cohortsize, ncohort=ncohort, ph1=0)
    res
}

run.fn2 <- function(k){
    print(k)
    res <- ORM.Eff.simu.fn.alter(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, cohortsize=cohortsize, ncohort=ncohort, ph1=0)
    res
}

ress <- mclapply(1:nsimu, run.fn1, mc.cores=20)
ress2 <- mclapply(1:nsimu, run.fn2, mc.cores=20)

#ress <- list()
#for (i in 1:100){
#    print(i)
#    res <- ORM.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, ncohort=12)
#    ress[[i]] <- res
#}
#
#ress2 <- list()
#for (i in 1:100){
#    print(i)
#    res <- ORM.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, ncohort=12, ph1=3)
#    ress2[[i]] <- res
#}

sum.res.stein <- get.oc(phi, pE.true, p.true, psi1, psi2, 
              ncohort, cohortsize, startdose=1, cutoff.eli=0.95, ntrial=nsimu)

sum.res.orm <- phase12.post.fn(ress)
sum.res.orm2 <- phase12.post.fn(ress2)

sum.all <- list(
                stein=sum.res.stein,
                orm=sum.res.orm,
                orm2=sum.res.orm2
                )

#save(sum.all, file="./results4.RData")

sum.fn <- function(ress){
    nams <- names(ress)
    ndose <- length(ress[[1]]$Allocation)
    res.df <- data.frame(levels=1:ndose)
    for (nam in nams){
        res.df[paste0(nam, ".Sel")] <- ress[[nam]]$Selection
    }
    for (nam in nams){
        res.df[paste0(nam, ".Allo")] <- ress[[nam]]$Allocation
    }
    res.df
}

sum.fn(sum.all)

