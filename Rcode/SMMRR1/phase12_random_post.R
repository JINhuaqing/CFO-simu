rm(list=ls())
#setwd("C:/Users/JINHU/Documents/ProjectCode/CFO/Rcode")
setwd("/home/r6user2/Documents/TQ/CFO/Rcode")
source("utilities.R")


phi <- 0.30
psi1 = phiE = 0.30
psi2 <- 0.65
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsim <- 5000
ncore <- 10
ndose <- 5

# Target = 0.3
# dose 5, mu1=mu2=0.23, 0.05
# dose 5, mu1=mu2=0.38, 0.07 
# dose 5, mu1=mu2=0.53, 0.1
# dose 5, mu1=mu2=0.71, 0.15
mu1 <- 0.38
typ <- "pl"
mu2 <- mu1

# random scenario results
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
load(save.f.Name.STEIN)
load(save.f.Name.WT)
load(save.f.Name.CFO)
load(save.f.Name.MADA)
load(save.f.Name.CFO2)
phase12.post.process.random.all(ress.WT, ress.STEIN, ress.MADA, 
                                ress.CFO, ress.CFO2, names=c("WT", "STEIN", "MADA", "CFO", "CFO2"))

# 
fadf
# CFO vs CFOconv
nsim <- 10000
save.f.Name1 <- paste0("../results/SMMR-R1/phase12_CFO_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name2 <- paste0("../results/SMMR-R1/phase12_CFOcov_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")


load(save.f.Name1)
load(save.f.Name2)

phase12.post.process.random.all(ress1, ress2, names=c("CFO", "CFO-conv"))



