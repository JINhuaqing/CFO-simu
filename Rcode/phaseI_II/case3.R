setwd("../")
source("utilities.R")
library(parallel)


phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsimu <- 5000
ncore <- 20

# increasing
increas <- list()
increas[[1]] <- list(p.true=c(0.05, 0.07, 0.1, 0.12, 0.16),
                      pE.true=c(0.35, 0.45, 0.5, 0.55, 0.75))
increas[[2]] <- list(p.true=c(0.03, 0.16, 0.27, 0.45, 0.55),
                      pE.true=c(0.15, 0.38, 0.45, 0.6, 0.7))
increas[[3]] <- list(p.true=c(0.03, 0.22, 0.30, 0.45, 0.55),
                      pE.true=c(0.15, 0.30, 0.45, 0.6, 0.7))
increas[[4]] <- list(p.true=c(0.10, 0.30, 0.35, 0.45, 0.55),
                      pE.true=c(0.25, 0.40, 0.45, 0.6, 0.7))
# 
idx <- 1
p.true <- increas[[idx]]$p.true
pE.true <- increas[[idx]]$pE.true

#source("phaseI_II/simu_efftox.R")
#fName <- paste0("../results/efftox_res_increa_", nsimu, "_", idx, ".RData")
#save(sum.res.efftox, file=fName)

source("phaseI_II/simu_orm.R")
source("phaseI_II/simu_Ada.R")
sum.all <- list(
                stein=sum.res.stein,
                orm=sum.res.orm,
#                efftox=sum.res.efftox,
                ada=sum.res.ada
                )

fName <- paste0("../results/res_increa_", nsimu, "_", idx, ".RData")
save(sum.all, file=fName)
print(OBD.level(phi, phiE, p.true, pE.true))
phase.I.II.pretty.tb(sum.all)



