library(magrittr)
library(parallel)
source("utilities.R")
source("crm_utils.R")
source("boin_utils.R")
source("butterfly_utils.R")


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

    butterfly.bf.res <- butterfly.simu.fn(target, p.true, type="BF", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
    butterfly.bms1.res <- butterfly.simu.fn(target, p.true, type="BMS", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args1)
    butterfly.bms2.res <- butterfly.simu.fn(target, p.true, type="BMS", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args2)
    butterfly.bms3.res <- butterfly.simu.fn(target, p.true, type="BMS", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args3)
    butterfly.bms4.res <- butterfly.simu.fn(target, p.true, type="BMS", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args4)
#    butterfly.odds.res <- butterfly.simu.fn(target, p.true, type="Odds", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
#    butterfly.bb.res <- butterfly.simu.fn(target, p.true, type="BB", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
#    butterfly.crm.res <- butterfly.simu.fn(target, p.true, type="CRM", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
#    butterfly.bbcrm.res <- butterfly.simu.fn(target, p.true, type="BB+CRM", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
    crm.res <- crm.simu.fn(target=target, p.true=p.true, p.prior=p.prior, cohortsize=cohortsize, ncohort=ncohort)
    boin.res <- boin.simu.fn(target=target, p.true=p.true, ncohort=ncohort, cohortsize)

    ress <- list(
                 butterfly.bf = butterfly.bf.res, 
                 butterfly.bms1 = butterfly.bms1.res, 
                 butterfly.bms2 = butterfly.bms2.res, 
                 butterfly.bms3 = butterfly.bms3.res, 
                 butterfly.bms4 = butterfly.bms4.res, 
#                 butterfly.od1ds = butterfly.odds.res, 
#                 butterfly.bb = butterfly.bb.res, 
#                 butterfly.crm = butterfly.crm.res, 
#                 butterfly.bbcrm = butterfly.bbcrm.res, 
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

#for (i in 1:20){
    #p.true.all <- gen.rand.doses(5, target, mu1=0.60, mu2=0.50)
    #p.true <- p.true.all$p.true
    #tmtd <- p.true.all$mtd.level
    #butterfly.bf.res <- butterfly.simu.fn(target, p.true, type="BF", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
    #print(butterfly.bf.res)
    #}

nsimu <- 10000
file.name <- paste0("../results/", "Simu", nsimu, "Level_5_10_ms",  ".RData")
results <- mclapply(1:nsimu, run.fn, mc.cores=20)
save(results, file=file.name)
post.process.random(results)

