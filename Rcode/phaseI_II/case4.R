setwd("../")
source("utilities.R")
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
 
idx <- 3
p.true <- ps[[idx]]$p.true
pE.true <- ps[[idx]]$pE.true

#source("phaseI_II/simu_efftox.R")
#fName <- paste0("../results/efftox_res_noOBD_", nsimu, "_", idx, ".RData")
#save(sum.res.efftox, file=fName)
source("phaseI_II/simu_orm.R")
source("phaseI_II/simu_Ada.R")
sum.all <- list(
                stein=sum.res.stein,
                orm=sum.res.orm,
#                orm.alter=sum.res.orm.alter,
#                efftox=sum.res.efftox,
                ada=sum.res.ada
                )

fName <- paste0("../results/res_noOBD_", nsimu, "_", idx, ".RData")
save(sum.all, file=fName)
print(OBD.level(phi, phiE, p.true, pE.true))
phase.I.II.pretty.tb(sum.all)



