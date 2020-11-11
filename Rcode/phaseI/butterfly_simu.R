rm(list=ls())
library(parallel)
#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/tryCode")
library(magrittr)
source("utilities.R")
source("butterfly_utils.R")
target <- 0.3

p.true1 <- c(0.1, 0.2, 0.3)
p.true2 <- c(0.05, 0.22, 0.38)

p.true3 <- c(0.2, 0.3, 0.4)
p.true4 <- c(0.18, 0.3, 0.45)

p.true5 <- c(0.07, 0.13, 0.21)
p.true6 <- c(0.04, 0.1, 0.2)

ncohort <- 12
cohortsize <- 1
m <- 2 
p.true <- p.true3
prefix <- paste0("m=", m)
tidx <- which.min(abs(p.true-target))
type <- "CRM"
if (m==50){
    suffix <- "[5pt]"
}else{
    suffix <- NULL
}
#add.args <- list(p.prior=c(0.1, 0.2, 0.3))
#add.args <- list(alp.prior=0.1, bet.prior=0.1)
add.args <- list(alp.prior=0.1, bet.prior=0.1, p.prior=c(0.1, 0.2, 0.3))

#butterfly.simu.fn(target, p.true1, type="BB+CRM", add.args=add.args)
#butterfly.simu.fn(target, p.true2, type="CRM", add.args=list(p.prior=c(0.1, 0.2, 0.3)))

#res <- nsimu.fn(target, p.true2, ncohort=ncohort, cohortsize=cohortsize, nsimu=1000, m=10)
#res <- post.process(res)
#latex.out.fn(res, prefix, tidx)


# Simulation code for linux 
run.fn <- function(k){
    print(k)
    res <- butterfly.simu.fn(target, p.true, type=type, ncohort=ncohort, cohortsize=cohortsize, m=m, add.args=add.args)
    res
}
results <- mclapply(1:1000, run.fn, mc.cores=20)
res <- post.process.raw(results)
latex.out.fn(res, prefix, tidx, suffix)

