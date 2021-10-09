rm(list=ls())
setwd("C:/Users/JINHU/Documents/ProjectCode/CFO/Rcode")
source("utilities.R")


phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsim <- 200
ncore <- 40
ndose <- 5
mu1 <- 0.52
mu2 <- mu1


save.f.Name1 <- paste0("../results/SMMR-R1/phase12_CFO_nsim", nsim, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")
save.f.Name2 <- paste0("../results/SMMR-R1/phase12_CFO2_nsim", nsim, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")

load(save.f.Name1)
load(save.f.Name2)

phase12.post.process.random.all(ress1, ress2, names=c("CFO", "CFO-conv"))



