setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")

source("utilities.R")
source("./phaseI/crm_utils.R")
source("./phaseI/boin_utils.R")
source("ORM_utils.R")


target <- 0.3
ncohort <- 12
cohortsize <- 3
init.level <- 1
nsimu <- 1000

p.prior <- c(0.1, 0.2, 0.3, 0.4, 0.5)
add.args <- list(alp.prior=target, bet.prior=1-target, p.prior=p.prior)
p.true <- c(0.05, 0.10, 0.3, 0.5, 0.6)
tmtd <- 2


results <- list()
for (i in 1:nsimu){
    print(i)
    #orm.res <- ORM.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize,
    #                            init.level=init.level,  add.args=add.args)
    crm.res <- crm.simu.fn(target=target, p.true=p.true, p.prior=p.prior, 
                              init.level=init.level, cohortsize=cohortsize, ncohort=ncohort)
    boin.res <- boin.simu.fn(target=target, p.true=p.true, ncohort=ncohort, 
                                init.level=init.level, cohortsize=cohortsize)
    ress <- list(
                 #orm=orm.res,
                 crm = crm.res, 
                 boin = boin.res, 
                 paras=list(p.true=p.true, 
                             mtd=tmtd, 
                             p.prior=p.prior,
                             add.args=add.args,
                             target=target,
                             ncohort=ncohort,
                             cohortsize=cohortsize)
        )
    results[[i]] <- ress
    
}


post.process.random(results)
