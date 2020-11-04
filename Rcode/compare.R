#rm(list=ls())
setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
source("ORM_eff_utils.R")
source("STEIN-utils.R")
source("utilities.R")

phi <- 0.30
cohortsize = 3 # cohort size
ncohort = 12 # number of cohorts

#p.true <- c(0.25, 0.35, 0.4, 0.45, 0.5)
#pE.true <- c(0.3, 0.35, 0.45, 0.65, 0.75)

p.true <- c(0.1, 0.12, 0.15, 0.2, 0.25)
pE.true <- c(0.25, 0.35, 0.6, 0.6, 0.6)

#p.true <- c(0.05, 0.07, 0.1, 0.15, 0.35)
#pE.true <- c(0.1, 0.2, 0.35, 0.5, 0.55)

#p.true <- c(0.02, 0.05, 0.07, 0.1, 0.15)
#pE.true <- c(0.05, 0.08, 0.15, 0.30, 0.45)

psi1<-0.35 # lowest acceptable efficacy rate
psi2<-0.65 # very promising efficacy rate that leads to dose retainment

phiE <- 0.5
add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)

ress <- list()
for (i in 1:100){
    print(i)
    res <- ORM.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, ncohort=12)
    ress[[i]] <- res
}

ress2 <- list()
for (i in 1:100){
    print(i)
    res <- ORM.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, ncohort=12, ph1=3)
    ress2[[i]] <- res
}

sum.res.stein <- get.oc(phi, pE.true, p.true, psi1, psi2, 
              ncohort, cohortsize, startdose=1, cutoff.eli=0.95, ntrial=1000)
sum.res.stein.raw <- get.oc.raw(phi, pE.true, p.true, psi1, psi2, 
              ncohort, cohortsize, startdose=1, cutoff.eli=0.95, ntrial=1000)

sum.res.orm <- phase12.post.fn(ress)
sum.res.orm2 <- phase12.post.fn(ress2)

sum.res.stein$pts
sum.res.stein.raw$pts
sum.res.orm$Allocation
sum.res.orm2$Allocation

sum.res.stein$sel
sum.res.stein.raw$sel
sum.res.orm$Selection
sum.res.orm2$Selection


