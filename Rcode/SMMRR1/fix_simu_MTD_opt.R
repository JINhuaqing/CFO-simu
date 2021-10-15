# simulation for MTD under fixed scinarios
# for optimal design
#setwd("../")
setwd("../../Rcode/")
#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
library(magrittr)
library(parallel)

source("utilities.R")
source("SMMRR1/optDes_utils.R")


target <- 0.33
ncohort <- 10
cohortsize <- 3
init.level <- 1
nsimu <- 5000

add.args <- list(alp.prior=target, bet.prior=1-target)
p.trues <- list()
p.trues[[1]] <- c(0.33, 0.45, 0.58, 0.70, 0.80) 
p.trues[[2]] <- c(0.18, 0.33, 0.52, 0.60, 0.70) 
p.trues[[3]] <- c(0.12, 0.20, 0.33, 0.40, 0.50) 
p.trues[[4]] <- c(0.01, 0.02, 0.03, 0.33, 0.50) 
p.trues[[5]] <- c(0.00, 0.00, 0.05, 0.10, 0.33) 
p.trues[[6]] <- c(0.45, 0.55, 0.65, 0.75, 0.85)

idx <- 6
p.true <- p.trues[[idx]]
N <- ncohort * cohortsize
set.seed(2021)
phaseI.opt.simu.fn(nsimu, N, target, p.true)



