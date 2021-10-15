setwd("/home/r6user2/MyResearch/CFO/Rcode/")
source("SMMRR1/WT_utils.R")
source("SMMRR1/MADA_utils.R")
source("SMMRR1/STEIN_single_utils.R")
source("phaseI_II/ORM_eff_utils.R")
source("SMMRR1/CFO_effTox_utils.R")
source("utilities.R")
library(parallel)

phi <- 0.30
psi1 = phiE = 0.30
psi2 <- 0.65
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsim <- 5000
ncore <- 40
ndose <- 5
seeds <- 1:nsim
add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)

# Target = 0.3
# dose 5, mu1=mu2=0.23, 0.05
# dose 5, mu1=mu2=0.38, 0.07 
# dose 5, mu1=mu2=0.53, 0.1
# dose 5, mu1=mu2=0.71, 0.15
mu1 <- 0.23
typ <- "um"

for (mu1 in c(0.23, 0.38, 0.53, 0.71)){
for (typ in c("um", "pl")){
mu2 <- mu1



f.Name.um <- paste0("../results/SMMR-R1/phase12_rc_umbrella_nsim", nsim, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")
f.Name.pl <- paste0("../results/SMMR-R1/phase12_rc_plateau_nsim", nsim, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")
load(f.Name.um)
load(f.Name.pl)

if (typ=="um"){
    scs <- scs.um
}else {
    scs <- scs.pl
}


run.fn.STEIN <- function(k){
    set.seed(seeds[k])
    print(c(1, k))
    sc <- scs[[k]]
    p.true <- sc$ps
    pE.true <- sc$qs
    
    res <- STEIN.simu.fn(phi, pT.true=p.true, pE.true=pE.true, cohortsize=cohortsize, ncohort=ncohort, startdose=1, psi1=psi1, psi2=psi2)
    res$k.OBD <- sc$k.OBD
    res$k.MTD <- sc$k.MTD
    res
}

run.fn.WT <- function(k){
    set.seed(seeds[k])
    print(c(2, k))
    sc <- scs[[k]]
    p.true <- sc$ps
    pE.true <- sc$qs
    
    res <- WT.simu.fn(tul=phi, ell=phiE, p0=p.true, q0=pE.true, cohortsize=cohortsize, ncohort=ncohort, start.comb=1)
    res$k.OBD <- sc$k.OBD
    res$k.MTD <- sc$k.MTD
    res
}

run.fn.MADA <- function(k){
    set.seed(seeds[k])
    print(c(3, k))
    sc <- scs[[k]]
    p.true <- sc$ps
    pE.true <- sc$qs
    
    res <- MADA.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, cohortsize=cohortsize, ncohort=ncohort, init.level=1)
    res$k.OBD <- sc$k.OBD
    res$k.MTD <- sc$k.MTD
    res
}

run.fn.CFO <- function(k){
    set.seed(seeds[k])
    print(c(4, k))
    sc <- scs[[k]]
    p.true <- sc$ps
    pE.true <- sc$qs
    
    res <- ORM.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, cohortsize=cohortsize, ncohort=ncohort, ph1=0)
    res$k.OBD <- sc$k.OBD
    res$k.MTD <- sc$k.MTD
    res
}

run.fn.CFO2 <- function(k){
    set.seed(seeds[k])
    print(c(5, k))
    sc <- scs[[k]]
    p.true <- sc$ps
    pE.true <- sc$qs
    res <- ORM2.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, cohortsize=cohortsize, ncohort=ncohort, ph1=0)
    res$k.OBD <- sc$k.OBD
    res$k.MTD <- sc$k.MTD
    res
}

prefix <- paste0("../results/SMMR-R1/phase12_", typ)
save.f.Name.STEIN <- paste0(prefix, "_STEIN2_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.CFO <- paste0(prefix, "_CFO_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.CFO2 <- paste0(prefix, "_CFO2_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.WT <- paste0(prefix, "_WT_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.MADA <- paste0(prefix, "_MADA_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")

ress.STEIN <- mclapply(1:nsim, run.fn.STEIN, mc.cores=ncore)
save(ress.STEIN, file=save.f.Name.STEIN)
# ress.WT <-    mclapply(1:nsim, run.fn.WT, mc.cores=ncore)
# save(ress.WT, file=save.f.Name.WT)
# ress.MADA <-  mclapply(1:nsim, run.fn.MADA, mc.cores=ncore)
# save(ress.MADA, file=save.f.Name.MADA)
# ress.CFO <-   mclapply(1:nsim, run.fn.CFO, mc.cores=ncore)
# save(ress.CFO, file=save.f.Name.CFO)
# ress.CFO2 <-  mclapply(1:nsim, run.fn.CFO2, mc.cores=ncore)
# save(ress.CFO2, file=save.f.Name.CFO2)

#phase12.post.process.random.all(ress.WT, ress.STEIN, ress.MADA, 
#                                ress.CFO, ress.CFO2, 
#                                names=c("WT", "STEIN", "MADA", "CFO", "CFO2"))

}
}
