setwd("../")
source("utilities.R")
source("phaseI_II/ORM_eff_utils.R")
library(parallel)


phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsimu <- 1000
ncore <- 10

# No OBD
ps <- list()
ps[[1]] <- list(p.true=c(0.25, 0.40, 0.5, 0.55, 0.6),
                pE.true=c(0.15, 0.25, 0.3, 0.45, 0.55))
ps[[2]] <- list(p.true=c(0.25, 0.40, 0.5, 0.55, 0.6),
                pE.true=c(0.15, 0.25, 0.5, 0.5, 0.5))
ps[[3]] <- list(p.true=c(0.25, 0.40, 0.5, 0.55, 0.6),
                pE.true=c(0.15, 0.25, 0.5, 0.3, 0.1))
ps[[4]] <- list(p.true=c(0.1, 0.15, 0.25, 0.30, 0.45),
                pE.true=c(0.05, 0.10, 0.17, 0.25, 0.35))
ps[[5]] <- list(p.true=c(0.22, 0.45, 0.55, 0.65, 0.70),
                pE.true=c(0.03, 0.10, 0.20, 0.35, 0.40))
 
idx <- 1
p.true <- ps[[idx]]$p.true
pE.true <- ps[[idx]]$pE.true

add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)

for (i in 1:5000){
    set.seed(i)
    res <- ORM.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, cohortsize=cohortsize, ncohort=ncohort, ph1=0)
    print(i)
    print(res)
}
#source("phaseI_II/simu_efftox.R")
#fName <- paste0("../results/efftox_res_noOBD_", nsimu, "_", idx, ".RData")
#save(sum.res.efftox, file=fName)
