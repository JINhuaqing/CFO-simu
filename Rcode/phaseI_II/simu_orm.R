source("phaseI_II/ORM_eff_utils.R")
source("phaseI_II/STEIN-utils.R")
source("utilities.R")
library(parallel)



add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)

psi1<-0.3 # lowest acceptable efficacy rate
psi2<-0.65 # very promising efficacy rate that leads to dose retainment

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

ress1 <- mclapply(1:nsimu, run.fn1, mc.cores=ncore)
ress2 <- mclapply(1:nsimu, run.fn2, mc.cores=ncore)


sum.res.stein <- get.oc(phi, pE.true, p.true, psi1, psi2, 
              ncohort, cohortsize, startdose=1, cutoff.eli=0.95, ntrial=nsimu)

sum.res.orm <- phase12.post.fn(ress1)
sum.res.orm$p.true <- p.true
sum.res.orm$pE.true <- pE.true
sum.res.orm.alter <- phase12.post.fn(ress2)
sum.res.orm.alter$p.true <- p.true
sum.res.orm.alter$pE.true <- pE.true



