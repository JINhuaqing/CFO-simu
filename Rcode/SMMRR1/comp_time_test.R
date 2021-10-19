# To calculate the computation time 
setwd("../../Rcode/")
library(magrittr)
library(rbenchmark)
source("utilities.R")
source("ORM_utils.R")
source("SMMRR1/CFO_effTox_utils.R")


# time for phase I trial
target <- 0.33
ncohort <- 10
cohortsize <- 3
init.level <- 1

add.args <- list(alp.prior=target, bet.prior=1-target)
p.true <- c(0.00, 0.00, 0.05, 0.10, 0.33) 


benchmark(
CFOph1 = ORM.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize,
                                init.level=init.level,  add.args=add.args)
)


# time for phase I/II trial
phi <- 0.30
psi1 = phiE = 0.30
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)

p.true <- c(0.05, 0.07, 0.1, 0.12, 0.16)
pE.true <- c(0.35, 0.45, 0.5, 0.55, 0.75)
benchmark(
CFOph12 = ORM2.Eff.simu.fn(phi, phiE, p.true=p.true, pE.true=pE.true, add.args=add.args, cohortsize=cohortsize, ncohort=ncohort, ph1=0)
)
