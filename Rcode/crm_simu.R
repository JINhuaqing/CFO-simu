library(magrittr)
library(parallel)
source("utilities.R")
source("crm_utils.R")

target <- 0.2
prior <- c(0.1, 0.2, 0.3)
rate <- 0.1
cycle <- 1

p.true1 <- c(0.1, 0.2, 0.3)
p.true2 <- c(0.05, 0.22, 0.38)

p.true3 <- c(0.2, 0.3, 0.4)
p.true4 <- c(0.18, 0.3, 0.45)

p.true5 <- c(0.07, 0.13, 0.21)
p.true6 <- c(0.04, 0.1, 0.2)

ncohort <- 10
cohortsize <- 3

crm.simu.fn(target, p.true1, ncohort=ncohort, cohortsize=cohortsize)


p.prior <- c(0.1, 0.2, 0.3)
target <- 0.3
ncohort <- 10
cohortsize <- 3
mu <- 0.3
run.fn <- function(k){
    print(k)
    p.true.all <- gen.rand.doses(3, target, mu1=mu, mu2=mu)
    p.true <- p.true.all$p.true
    tmtd <- p.true.all$mtd.level

    crm.res <- crm.simu.fn(target=target, p.true=p.true, p.prior=p.prior, cohortsize=cohortsize, ncohort=ncohort)

    ress <- list(
                 crm = crm.res, 
                 crm = crm.res, 
                 paras=list(p.true=p.true, 
                         mtd=tmtd, 
                         p.prior=p.prior,
                         target=target,
                         ncohort=ncohort,
                         cohortsize=cohortsize
                         )
                 )
    ress
}

nsimu <- 10000
results <- mclapply(1:nsimu, run.fn, mc.cores=20)
post.process.random(results)

