setwd("/home/r6user2/MyResearch/CFO/Rcode/")
source("SMMRR1/WT_utils.R")
#source("SMMRR1/MADA_utils.R")
source("utilities.R")
library(parallel)

phi <- 0.30
psi1 = phiE = 0.30
psi2 <- 0.65
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsim <- 200
ncore <- 40
ndose <- 5
mu1 <- 0.52
mu2 <- mu1
seeds <- 1:nsim


add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)

nsim1 <- ceiling(nsim/2)
f.Name.um <- paste0("../results/SMMR-R1/phase12_rc_umbrella_nsim", nsim1, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")
f.Name.pl <- paste0("../results/SMMR-R1/phase12_rc_plateau_nsim", nsim1, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")
load(f.Name.um)
load(f.Name.pl)



run.fn1 <- function(k){
    set.seed(seeds[k])
    print(c(1, k))
    if (k <= nsim1){
        sc <- scs.um[[k]]
    }else{
        sc <- scs.pl[[k-nsim1]]
    }
    p.true <- sc$ps
    pE.true <- sc$qs
    
    #res <- MADA.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, cohortsize=cohortsize, ncohort=ncohort, init.level=1)
    res <- WT.simu.fn(tul=phi, ell=phiE, p0=p.true, q0=pE.true, cohortsize=cohortsize, ncohort=ncohort, start.comb=1)
    res$k.OBD <- sc$k.OBD
    res$k.MTD <- sc$k.MTD
    res
}


ress1 <- mclapply(1:nsim, run.fn1, mc.cores=ncore)



save.f.Name1 <- paste0("../results/SMMR-R1/phase12_STEIN_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")

save(ress1, file=save.f.Name1)

phase12.post.process.random.all(ress1)
