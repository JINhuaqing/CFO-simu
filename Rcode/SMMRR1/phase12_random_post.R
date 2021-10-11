rm(list=ls())
#setwd("C:/Users/JINHU/Documents/ProjectCode/CFO/Rcode")
setwd("/home/r6user2/Documents/TQ/CFO/Rcode")
source("utilities.R")


phi <- 0.30
phiE <- 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsim <- 5000
ncore <- 40
ndose <- 5
ncohort <- 20
cohortsize <- 3
mu1 <- 0.52
mu2 <- mu1

save.f.Name1 <- paste0("../results/SMMR-R1/phase12_CFO_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name2 <- paste0("../results/SMMR-R1/phase12_CFOcov_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")


load(save.f.Name1)
load(save.f.Name2)

phase12.post.process.random.all(ress1, ress2, names=c("CFO", "CFO-conv"))



