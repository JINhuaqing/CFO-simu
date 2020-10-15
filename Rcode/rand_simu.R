library(magrittr)
library(parallel)
source("utilities.R")
source("crm_utils.R")
source("boin_utils.R")
source("BF_utils.R")
source("ORM_utils.R")
source("BMS_utils.R")


#set.seed(10)
target <- 0.3
ncohort <- 10
cohortsize <- 3
m <- 50

p.prior <- c(0.1, 0.2, 0.3, 0.4, 0.5)#, 0.6, 0.7)
add.args <- list(alp.prior=target, bet.prior=1-target, p.prior=p.prior, m=m)
add.args1 <- list(alp.prior=target, bet.prior=1-target, p.prior=p.prior, m=1)
add.args2 <- list(alp.prior=target, bet.prior=1-target, p.prior=p.prior, m=10)
add.args3 <- list(alp.prior=target, bet.prior=1-target, p.prior=p.prior, m=50)
add.args4 <- list(alp.prior=target, bet.prior=1-target, p.prior=p.prior, m=Inf)

## Target = 0.3
# dose 3, mu1=0.55, mu2=0.40, 0.1
# dose 3, mu1=mu2=0.30, 0.07
# dose 3, mu1=mu2=0.46, 0.1
# dose 3, mu1=mu2=0.64, 0.15

# dose 5, mu1=0.60, mu2=0.50, 0.1
# dose 5, mu1=mu2=0.38, 0.07
# dose 5, mu1=mu2=0.53, 0.1
# dose 5, mu1=mu2=0.71, 0.15

# dose 7, mu1=0.70, mu2=0.50, 0.1
# dose 7, mu1=mu2=0.42, 0.07
# dose 7, mu1=mu2=0.56, 0.1
# dose 7, mu1=mu2=0.74, 0.15

mu <- 0.53
run.fn <- function(k){
    print(k)
    p.true.all <- gen.rand.doses(5, target, mu1=mu, mu2=mu)
    p.true <- p.true.all$p.true
    tmtd <- p.true.all$mtd.level

    bms1.res <- BMS.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, add.args=add.args1)
    orm.res <- ORM.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
    bms2.res <- BMS.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, add.args=add.args2)
    bms3.res <- BMS.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, add.args=add.args3)
    bms4.res <- BMS.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, add.args=add.args4)
    crm.res <- crm.simu.fn(target=target, p.true=p.true, p.prior=p.prior, cohortsize=cohortsize, ncohort=ncohort)
    boin.res <- boin.simu.fn(target=target, p.true=p.true, ncohort=ncohort, cohortsize)

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
file.name <- paste0("../results/", "Simu", nsimu, "Level_5_10_ms_corhortsize3",  ".RData")
results <- mclapply(1:nsimu, run.fn, mc.cores=15)
post.process.random(results)
save(results, file=file.name)

