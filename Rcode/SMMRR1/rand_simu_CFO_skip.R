# compare CFO and CFO allowing dose skipping
setwd("/home/r6user2/MyResearch/CFO/Rcode/")
source("phaseI_II/ORM_eff_utils.R")
source("SMMRR1/ORMConv_eff_utils.R")
source("utilities.R")
library(parallel)

phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsim <- 200
ncore <- 40
ndose <- 5
mu1 <- 0.52
mu2 <- mu1


add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)

nsim1 <- ceiling(nsim/2)
f.Name.um <- paste0("../results/SMMR-R1/phase12_rc_umbrella_nsim", nsim1, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")
f.Name.pl <- paste0("../results/SMMR-R1/phase12_rc_plateau_nsim", nsim1, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")
load(f.Name.um)
load(f.Name.pl)



run.fn1 <- function(k){
    print(c(1, k))
    if (k <= nsim1){
        sc <- scs.um[[k]]
    }else{
        sc <- scs.pl[[k-nsim1]]
    }
    p.true <- sc$ps
    pE.true <- sc$qs
    
    res <- ORM.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, cohortsize=cohortsize, ncohort=ncohort, ph1=0)
    res$k.OBD <- sc$k.OBD
    res$k.MTD <- sc$k.MTD
    res
}

run.fn2 <- function(k){
    print(c(2, k))
    if (k <= nsim1){
        sc <- scs.um[[k]]
    }else{
        sc <- scs.pl[[k-nsim1]]
    }
    p.true <- sc$ps
    pE.true <- sc$qs
    res <- ORM2.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, cohortsize=cohortsize, ncohort=ncohort, ph1=0)
    res$k.OBD <- sc$k.OBD
    res$k.MTD <- sc$k.MTD
    res
}

ress1 <- mclapply(1:nsim, run.fn1, mc.cores=ncore)
ress2 <- mclapply(1:nsim, run.fn2, mc.cores=ncore)


save.f.Name1 <- paste0("../results/SMMR-R1/phase12_CFO_nsim", nsim, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")
save.f.Name2 <- paste0("../results/SMMR-R1/phase12_CFO2_nsim", nsim, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")

#sum.res.orm <- phase12.post.fn(ress1)
#sum.res.orm.alter <- phase12.post.fn(ress2)



