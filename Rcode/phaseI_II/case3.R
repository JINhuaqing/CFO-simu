setwd("../")
source("utilities.R")
library(parallel)


phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsimu <- 1000
ncore <- 20

# increasing
increas <- list()
increas[[1]] <- list(p.true=c(0.05, 0.07, 0.1, 0.12, 0.16),
                      pE.true=c(0.35, 0.45, 0.5, 0.55, 0.75))
idx <- 1
p.true <- increas[[idx]]$p.true
pE.true <- increas[[idx]]$pE.true

#source("phaseI_II/simu_efftox.R")
#fName <- paste0("../results/efftox_res_increa_", nsimu, "_", idx, ".RData")
#save(sum.res.efftox, file=fName)

#source("phaseI_II/simu_orm.R")
source("phaseI_II/simu_Ada.R")
sum.all <- list(
#                stein=sum.res.stein,
#                orm=sum.res.orm,
#                efftox=sum.res.efftox,
                ada=sum.res.ada
                )

fName <- paste0("../results/ada_res_increa_", nsimu, "_", idx, ".RData")
#fName <- paste0("../results/res_increa_", nsimu, "_", idx, ".RData")
save(sum.all, file=fName)
print(OBD.level(phi, phiE, p.true, pE.true))
phase.I.II.pretty.tb(sum.all)



