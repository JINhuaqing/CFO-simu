setwd("../")
source("utilities.R")
library(parallel)


phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsimu <- 1000
ncore <- 10

# umbrella-shape
umbrellas <- list()
umbrellas[[1]] <- list(p.true=c(0.05, 0.1, 0.2, 0.3, 0.35),
                      pE.true=c(0.2, 0.3, 0.5, 0.35, 0.2))
umbrellas[[2]] <- list(p.true=c(0.05, 0.1, 0.15, 0.2, 0.3),
                      pE.true=c(0.2, 0.3, 0.35, 0.5, 0.2))
umbrellas[[3]] <- list(p.true=c(0.05, 0.1, 0.15, 0.2, 0.3),
                      pE.true=c(0.2, 0.3, 0.5, 0.35, 0.2))
umbrellas[[4]] <- list(p.true=c(0.05, 0.1, 0.15, 0.2, 0.3),
                      pE.true=c(0.3, 0.5, 0.45, 0.35, 0.2))
umbrellas[[5]] <- list(p.true=c(0.05, 0.2, 0.25, 0.3, 0.4),
                      pE.true=c(0.2, 0.6, 0.55, 0.35, 0.2))
umbrellas[[6]] <- list(p.true=c(0.05, 0.22, 0.25, 0.3, 0.4),
                      pE.true=c(0.2, 0.6, 0.55, 0.35, 0.2))
umbrellas[[7]] <- list(p.true=c(0.10, 0.22, 0.25, 0.3, 0.4),
                      pE.true=c(0.3, 0.6, 0.55, 0.35, 0.2))
umbrellas[[8]] <- list(p.true=c(0.05, 0.1, 0.15, 0.2, 0.3),
                      pE.true=c(0.2, 0.6, 0.7, 0.55, 0.2))
umbrellas[[9]] <- list(p.true=c(0.05, 0.1, 0.15, 0.2, 0.3),
                      pE.true=c(0.2, 0.55, 0.7, 0.55, 0.2))
umbrellas[[10]] <- list(p.true=c(0.05, 0.15, 0.30, 0.45, 0.60),
                      pE.true=c(0.08, 0.46, 0.25, 0.20, 0.1))
umbrellas[[11]] <- list(p.true=c(0.15, 0.25, 0.40, 0.45, 0.50),
                      pE.true=c(0.15, 0.45, 0.30, 0.25, 0.2))
umbrellas[[12]] <- list(p.true=c(0.05, 0.15, 0.25, 0.40, 0.45),
                      pE.true=c(0.08, 0.17, 0.45, 0.30, 0.25))
# 
idx <- 12
p.true <- umbrellas[[idx]]$p.true
pE.true <- umbrellas[[idx]]$pE.true

#source("phaseI_II/simu_efftox.R")
#fName <- paste0("../results/efftox_res_umbrella_", nsimu, "_", idx, ".RData")
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

fName <- paste0("../results/res_umbrella_", nsimu, "_", idx, ".RData")
save(sum.all, file=fName)
print(OBD.level(phi, phiE, p.true, pE.true))
phase.I.II.pretty.tb(sum.all)



