# simulation for OBD under fixed scinarios
# for optimal design
#setwd("../")
setwd("../../Rcode/")
#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
library(magrittr)
library(parallel)

source("utilities.R")
source("SMMRR1/optDes_utils.R")


phi <- 0.30
psi1 = phiE = 0.30
psi2 <- 0.65
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsim <- 5000
ncore <- 40
ndose <- 5
seeds <- 1:nsim
add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)

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


#sapply(scs, function(x) {100*c(mean(x$p.true), mean(x$pE.true))})

idx <- 7

p.true <- scs[[idx]]$p.true
pE.true <- scs[[idx]]$pE.true

N <- ncohort * cohortsize
set.seed(2021)
round(phase12.opt.simu.fn(nsim, N, phi, phiE, p.true, pE.true), 1)



