library(magrittr)
library(parallel)
source("utilities.R")
source("crm_utils.R")
source("boin_utils.R")
source("butterfly_utils.R")


target <- 0.3
ncohort <- 10
cohortsize <- 3
m <- 10

add.args <- list(alp.prior=0.1, bet.prior=0.1, p.prior=c(0.1, 0.2, 0.3, 0.4, 0.5), m=m)
p.prior <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# dose 3, mu1=0.55, mu2=0.40
# dose 5, mu1=0.60, mu2=0.50
run.fn <- function(k){
    print(k)
    p.true.all <- gen.rand.doses(5, target, mu1=0.60, mu2=0.50)
    p.true <- p.true.all$p.true
    tmtd <- p.true.all$mtd.level
    butterfly.odds.res <- butterfly.simu.fn(target, p.true, type="Odds", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
    butterfly.bb.res <- butterfly.simu.fn(target, p.true, type="BB", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
    butterfly.crm.res <- butterfly.simu.fn(target, p.true, type="CRM", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
    butterfly.bbcrm.res <- butterfly.simu.fn(target, p.true, type="BB+CRM", ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
    crm.res <- crm.simu.fn(target=target, p.true=p.true, p.prior=p.prior, cohortsize=cohortsize, ncohort=ncohort)
    boin.res <- boin.simu.fn(target=target, p.true=p.true, ncohort=ncohort, cohortsize)

    ress <- list(
                 butterfly.odds = butterfly.odds.res, 
                 butterfly.bb = butterfly.bb.res, 
                 butterfly.crm = butterfly.crm.res, 
                 butterfly.bbcrm = butterfly.bbcrm.res, 
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
file.name <- paste0("../results/", "Allmethods", nsimu, "TSS_30_Level_5",  ".RData")
#file.name <- paste0("../results/", "Odds_rand", nsimu, "_m_", m, ".RData")
#file.name <- paste0("../results/", "rand", nsimu, "_m_", m, ".RData")
results <- mclapply(1:nsimu, run.fn, mc.cores=20)
save(results, file=file.name)
post.process.random(results)

