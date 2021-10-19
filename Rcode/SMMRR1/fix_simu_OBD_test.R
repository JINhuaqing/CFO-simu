# For OBD simulation under fixed scs
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

scs <- list()
scs[[1]] <- list(p.true=c(0.05, 0.1, 0.3, 0.5, 0.6),
                      pE.true=c(0.2, 0.3, 0.5, 0.5, 0.5))
scs[[2]] <- list(p.true=c(0.15, 0.25, 0.3, 0.35, 0.4),
                      pE.true=c(0.2, 0.5, 0.5, 0.5, 0.5))
scs[[3]] <- list(p.true=c(0.10, 0.22, 0.25, 0.3, 0.4),
                      pE.true=c(0.3, 0.6, 0.55, 0.35, 0.2))
scs[[4]] <- list(p.true=c(0.05, 0.15, 0.25, 0.40, 0.45),
                      pE.true=c(0.08, 0.17, 0.45, 0.30, 0.25))
scs[[5]] <- list(p.true=c(0.05, 0.07, 0.1, 0.12, 0.16),
                      pE.true=c(0.35, 0.45, 0.5, 0.55, 0.75))
scs[[6]] <- list(p.true=c(0.25, 0.40, 0.5, 0.55, 0.6),
                pE.true=c(0.15, 0.25, 0.5, 0.5, 0.5))
scs[[7]] <- list(p.true=c(0.40, 0.5, 0.55, 0.6, 0.7),
                 pE.true=c(0.15, 0.25, 0.5, 0.5, 0.5))


idx <- 3

p.true <- scs[[idx]]$p.true
pE.true <- scs[[idx]]$pE.true

# test whether OBD = max U in STEIN
# u<-pE.true-0.33*p.true-1.09*(p.true>phi) 	
# print(u)
# print(which.max(u))
# print(OBD.level(phi, phiE, p.true, pE.true))

for (idx in 3){

p.true <- scs[[idx]]$p.true
pE.true <- scs[[idx]]$pE.true


run.fn.STEIN <- function(k){
    set.seed(seeds[k])
    print(c(1, k))
    res <- STEIN.simu.fn(phi, pT.true=p.true, pE.true=pE.true, cohortsize=cohortsize, ncohort=ncohort, startdose=1, psi1=psi1, psi2=psi2)
    res
}

run.fn.WT <- function(k){
    set.seed(seeds[k])
    print(c(2, k))
    
    res <- WT.simu.fn(tul=phi, ell=phiE, p0=p.true, q0=pE.true, cohortsize=cohortsize, ncohort=ncohort, start.comb=1)
    res
}

run.fn.MADA <- function(k){
    set.seed(seeds[k])
    print(c(3, k))
    
    res <- MADA.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, cohortsize=cohortsize, ncohort=ncohort, init.level=1)
    res
}

run.fn.CFO <- function(k){
    set.seed(seeds[k])
    print(c(4, k))
    
    res <- ORM.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, cohortsize=cohortsize, ncohort=ncohort, ph1=0)
    res
}

run.fn.CFO2 <- function(k){
    set.seed(seeds[k])
    print(c(5, k))

    res <- ORM2.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, cohortsize=cohortsize, ncohort=ncohort, ph1=0)
    res
}

prefix <- paste0("../results/SMMR-R1/phase12_Fixed_", idx)
save.f.Name.STEIN <- paste0(prefix, "_STEIN_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.CFO <- paste0(prefix, "_CFO_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.CFO2 <- paste0(prefix, "_CFO2_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.WT <- paste0(prefix, "_WT_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.MADA <- paste0(prefix, "_MADA_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")

ress.STEIN <- mclapply(1:nsim, run.fn.STEIN, mc.cores=ncore)
#save(ress.STEIN, file=save.f.Name.STEIN)
#ress.WT <-    mclapply(1:nsim, run.fn.WT, mc.cores=ncore)
#save(ress.WT, file=save.f.Name.WT)
#ress.MADA <-  mclapply(1:nsim, run.fn.MADA, mc.cores=ncore)
#save(ress.MADA, file=save.f.Name.MADA)
# ress.CFO <-   mclapply(1:nsim, run.fn.CFO, mc.cores=ncore)

#sum.res.WT <- phase12.post.fn(ress.WT)
sum.res.STEIN <- phase12.post.fn(ress.STEIN)

sum.all <- list(
#CFO=sum.res.CFO ,
#CFO2=sum.res.CFO2,
#MADA=sum.res.MADA,
STEIN=sum.res.STEIN)
#WT=sum.res.WT)

print(OBD.level(phi, phiE, p.true, pE.true))
print(phase.I.II.pretty.tb(sum.all))
}
