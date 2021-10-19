setwd("../")
source("utilities.R")


phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
ndose <- 5
nsim <- 5000

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
idx <- 7 # 5, 7
p.true <- scs[[idx]]$p.true
pE.true <- scs[[idx]]$pE.true
prefix <- paste0("../results/SMMR-R1/phase12_Fixed_", idx)
save.f.Name.STEIN <- paste0(prefix, "_STEIN2_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.CFO <- paste0(prefix, "_CFO_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.CFO2 <- paste0(prefix, "_CFO2_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.WT <- paste0(prefix, "_WT_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.MADA <- paste0(prefix, "_MADA_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")

load(save.f.Name.CFO)
load(save.f.Name.CFO2)
load(save.f.Name.MADA)
load(save.f.Name.STEIN)
load(save.f.Name.WT)

sum.res.CFO <- phase12.post.fn(ress.CFO)
sum.res.CFO2 <- phase12.post.fn(ress.CFO2)
sum.res.MADA <- phase12.post.fn(ress.MADA)
sum.res.WT <- phase12.post.fn(ress.WT)
sum.res.STEIN <- phase12.post.fn(ress.STEIN)

sum.all <- list(
#CFO=sum.res.CFO ,
CFO2=sum.res.CFO2,
MADA=sum.res.MADA,
STEIN=sum.res.STEIN,
WT=sum.res.WT
)

print(OBD.level(phi, phiE, p.true, pE.true))
phase.I.II.pretty.tb(sum.all)[, -8]


