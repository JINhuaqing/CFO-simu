#rm(list=ls())
#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
#set.seed(0);
setwd("../")
source("phaseI_II/ORM_eff_utils.R")
source("phaseI_II/STEIN-utils.R")
source("utilities.R")
library(parallel)

phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsimu <- 1000

# plateau
plateaus <- list()
plateaus[[1]] <- list(p.true=c(0.05, 0.1, 0.2, 0.3, 0.35),
                      pE.true=c(0.2, 0.3, 0.5, 0.5, 0.5))
plateaus[[2]] <- list(p.true=c(0.05, 0.1, 0.3, 0.5, 0.6),
                      pE.true=c(0.2, 0.3, 0.5, 0.5, 0.5))
plateaus[[3]] <- list(p.true=c(0.05, 0.1, 0.3, 0.5, 0.6),
                      pE.true=c(0.1, 0.2, 0.6, 0.6, 0.6))
idx <- 1
p.true <- plateaus[[idx]]$p.true
pE.true <- plateaus[[idx]]$pE.true


psi1<-0.3 # lowest acceptable efficacy rate
psi2<-0.8 # very promising efficacy rate that leads to dose retainment

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

ress1 <- mclapply(1:nsimu, run.fn1, mc.cores=20)
ress2 <- mclapply(1:nsimu, run.fn2, mc.cores=20)

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

fName <- paste0("res_plateau_", idx, ".RData")
save(sum.all, file=fName)

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

