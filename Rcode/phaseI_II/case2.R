#rm(list=ls())
#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
#set.seed(0);
setwd("../")
source("phaseI_II/ORM_eff_utils.R")
source("phaseI_II/STEIN-utils.R")
source("utilities.R")
library(parallel)

#fName <- paste0("res_plateau_", 4, ".RData")
#load(fName)
#
#phase.I.II.pretty.tb(sum.all)

phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsimu <- 1000

# umbrella-shape
umbrellas <- list()
umbrellas[[1]] <- list(p.true=c(0.05, 0.1, 0.2, 0.3, 0.35),
                      pE.true=c(0.2, 0.3, 0.5, 0.35, 0.2))
umbrellas[[2]] <- list(p.true=c(0.05, 0.1, 0.15, 0.2, 0.3),
                      pE.true=c(0.2, 0.3, 0.35, 0.5, 0.2))
umbrellas[[3]] <- list(p.true=c(0.05, 0.1, 0.15, 0.2, 0.3),
                      pE.true=c(0.2, 0.3, 0.5, 0.35, 0.2))
umbrellas[[4]] <- list(p.true=c(0.05, 0.1, 0.15, 0.2, 0.3),
                      pE.true=c(0.3, 0.5, 0.45, 0.35, 0.2))
umbrellas[[5]] <- list(p.true=c(0.05, 0.2, 0.25, 0.3, 0.4),
                      pE.true=c(0.2, 0.6, 0.55, 0.35, 0.2))
umbrellas[[6]] <- list(p.true=c(0.05, 0.22, 0.25, 0.3, 0.4),
                      pE.true=c(0.2, 0.6, 0.55, 0.35, 0.2))
umbrellas[[7]] <- list(p.true=c(0.10, 0.22, 0.25, 0.3, 0.4),
                      pE.true=c(0.3, 0.6, 0.55, 0.35, 0.2))
umbrellas[[8]] <- list(p.true=c(0.05, 0.1, 0.15, 0.2, 0.3),
                      pE.true=c(0.2, 0.6, 0.7, 0.55, 0.2))
umbrellas[[9]] <- list(p.true=c(0.05, 0.1, 0.15, 0.2, 0.3),
                      pE.true=c(0.2, 0.55, 0.7, 0.55, 0.2))
# 7, 8, 9 is good
idx <- 9
p.true <- umbrellas[[idx]]$p.true
pE.true <- umbrellas[[idx]]$pE.true


psi1<-0.3 # lowest acceptable efficacy rate
psi2<-0.65 # very promising efficacy rate that leads to dose retainment

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

ress1 <- mclapply(1:nsimu, run.fn1, mc.cores=8)
ress2 <- mclapply(1:nsimu, run.fn2, mc.cores=8)

sum.res.stein <- get.oc(phi, pE.true, p.true, psi1, psi2, 
              ncohort, cohortsize, startdose=1, cutoff.eli=0.95, ntrial=nsimu)

sum.res.orm <- phase12.post.fn(ress1)
sum.res.orm$p.true <- p.true
sum.res.orm$pE.true <- pE.true
sum.res.orm.alter <- phase12.post.fn(ress2)
sum.res.orm.alter$p.true <- p.true
sum.res.orm.alter$pE.true <- pE.true
source("phaseI_II/simu_Ada.R")

sum.all <- list(
                stein=sum.res.stein,
                orm=sum.res.orm,
                orm.alter=sum.res.orm.alter,
                ada=sum.res.ada
                )

fName <- paste0("res_umbrella_", idx, ".RData")
save(sum.all, file=fName)

print(OBD.level(phi, phiE, p.true, pE.true))
phase.I.II.pretty.tb(sum.all)

