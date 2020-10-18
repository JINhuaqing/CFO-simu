library(magrittr)
library(parallel)
source("utilities.R")
source("crm_utils.R")
source("boin_utils.R")
source("ORM_utils.R")
source("BMS_utils.R")


set.seed(1)
target <- 0.25
ncohort <- 10
cohortsize <- 3
init.level <- 1

p.prior <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
add.args <- list(alp.prior=target, bet.prior=1-target, p.prior=p.prior, m=50)
add.args1 <- list(alp.prior=target, bet.prior=1-target, p.prior=p.prior, m=1)
add.args2 <- list(alp.prior=target, bet.prior=1-target, p.prior=p.prior, m=10)
add.args3 <- list(alp.prior=target, bet.prior=1-target, p.prior=p.prior, m=50)
add.args4 <- list(alp.prior=target, bet.prior=1-target, p.prior=p.prior, m=Inf)


run.fn <- function(k){
    print(k)
    p.true <- c(0.25, 0.35, 0.5, 0.6, 0.7, 0.8)
    tmtd <- 1

    bms1.res <- BMS.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, 
                            init.level=init.level, add.args=add.args1)
    orm.res <- ORM.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize,
                            init.level=init.level,  add.args=add.args)
    bms2.res <- BMS.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize,
                           init.level=init.level,  add.args=add.args2)
    bms3.res <- BMS.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize,
                           init.level=init.level, add.args=add.args3)
    bms4.res <- BMS.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize,
                           init.level=init.level, add.args=add.args4)
    crm.res <- crm.simu.fn(target=target, p.true=p.true, p.prior=p.prior, 
                          init.level=init.level, cohortsize=cohortsize, ncohort=ncohort)
    boin.res <- boin.simu.fn(target=target, p.true=p.true, ncohort=ncohort, 
                            init.level=init.level, cohortsize=corhortsize)

    ress <- list(
                 orm=orm.res,
                 bms1 = bms1.res, 
                 bms2 = bms2.res, 
                 bms3 = bms3.res, 
                 bms4 = bms4.res, 
                 crm = crm.res, 
                 boin = boin.res, 
                 paras=list(p.true=p.true, 
                         mtd=tmtd, 
                         p.prior=p.prior,
                         add.args=add.args,
                         target=target,
                         ncohort=ncohort,
                         cohortsize=cohortsize,
                         m=m)
                 )
    ress
}


nsimu <- 10000
file.name <- paste0("../results/", "Fix_simu", nsimu, "Level_6_fixed1_corhortsize3",  ".RData")
results <- mclapply(1:nsimu, run.fn, mc.cores=15)
post.process.random(results)
save(results, file=file.name)
