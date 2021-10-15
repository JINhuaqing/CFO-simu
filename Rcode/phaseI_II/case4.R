setwd("../")
source("utilities.R")
library(parallel)


phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsimu <- 1000
ncore <- 20

# No OBD
ps <- list()
ps[[1]] <- list(p.true=c(0.25, 0.40, 0.5, 0.55, 0.6),
                pE.true=c(0.15, 0.25, 0.5, 0.5, 0.5))
 
idx <- 1
p.true <- ps[[idx]]$p.true
pE.true <- ps[[idx]]$pE.true

#source("phaseI_II/simu_efftox.R")
#fName <- paste0("../results/efftox_res_noOBD_", nsimu, "_", idx, ".RData")
#save(sum.res.efftox, file=fName)

#source("phaseI_II/simu_orm.R")
source("phaseI_II/simu_Ada.R")
sum.all <- list(
#                stein=sum.res.stein,
#                orm=sum.res.orm,
#                efftox=sum.res.efftox,
                ada=sum.res.ada
                )

fName <- paste0("../results/ada_res_noOBD_", nsimu, "_", idx, ".RData")
#fName <- paste0("../results/res_noOBD_", nsimu, "_", idx, ".RData")
save(sum.all, file=fName)
print(OBD.level(phi, phiE, p.true, pE.true))
phase.I.II.pretty.tb(sum.all)



