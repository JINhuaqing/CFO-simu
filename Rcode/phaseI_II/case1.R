setwd("../")
source("utilities.R")
library(parallel)


phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsimu <- 1000
ncore <- 10

# plateau
plateaus <- list()
plateaus[[1]] <- list(p.true=c(0.05, 0.1, 0.3, 0.5, 0.6),
                      pE.true=c(0.2, 0.3, 0.5, 0.5, 0.5))
plateaus[[2]] <- list(p.true=c(0.15, 0.25, 0.3, 0.35, 0.4),
                      pE.true=c(0.2, 0.5, 0.5, 0.5, 0.5))
idx <- 1
p.true <- plateaus[[idx]]$p.true
pE.true <- plateaus[[idx]]$pE.true


#source("phaseI_II/simu_efftox.R")
#fName <- paste0("../results/efftox_res_plateau_", nsimu, "_", idx, ".RData")
#save(sum.res.efftox, file=fName)


#source("phaseI_II/simu_orm.R")
source("phaseI_II/simu_Ada.R")
sum.all <- list(
#                stein=sum.res.stein,
#                orm=sum.res.orm,
#                orm.alter=sum.res.orm.alter,
#                efftox=sum.res.efftox,
                ada=sum.res.ada
                )

fName <- paste0("../results/ada_res_plateau_", nsimu, "_", idx, ".RData")
#fName <- paste0("../results/res_plateau_", nsimu, "_", idx, ".RData")
save(sum.all, file=fName)
print(OBD.level(phi, phiE, p.true, pE.true))
phase.I.II.pretty.tb(sum.all)


