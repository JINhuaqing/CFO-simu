setwd("C:/Users/JINHU/Documents/ProjectCode/CFO/Rcode")
source("utilities.R")
source("ORM_utils.R")


gen.rand.overdoses <- function(ndose, target, delta=0.1, mu2=0.55, inv.fn=qnorm, fn=pnorm){
    sigma0 <- 0.05
    sigma2 <- 0.35
    raw.p.trues <- rep(0, ndose)
    
    raw.p.trues[1] <- rnorm(1, inv.fn(target+delta), sigma0)
    
    for (j in 2:ndose){
        eps.plus <- rnorm(1, mu2, sigma2)
        raw.p.trues[j] <- raw.p.trues[j-1] + eps.plus**2
    }
    p.trues <- fn(raw.p.trues)
    p.trues
}


gen.rand.overdoses.phase12 <- function(ndose, phi, psi, delta=0.1, psi.U=0.8, mu2=0.3){
    # args:
    # ndose: number of dose levels
    # phi: DLT target rate
    # psi: the minimal acceptable efficacy rate 
    # psi.U: the upper bound of the efficacy rate
    # mu1, mu2: parameters for phase I trial
    
    k.OBD <- sample.int(ndose, 1)
    ps <- gen.rand.overdoses(ndose, phi, delta=delta, mu2=mu2)
    
    rv <- runif(1)
    q.OBD <- runif(1, psi, psi.U)
    if (rv <= 0.5){
        if (k.OBD == 1){
            qs.u <- sort(runif(ndose-1, 0, q.OBD), decreasing=TRUE)
            qs <- c(q.OBD, qs.u)
        }else if (k.OBD == ndose){
            qs.l <- sort(runif(ndose-1, 0, q.OBD))
            qs <- c(qs.l, q.OBD)
        }else {
            qs.l <- sort(runif(k.OBD-1, 0, q.OBD))
            qs.u <- sort(runif(ndose-k.OBD, 0, q.OBD), decreasing=TRUE)
            qs <- c(qs.l, q.OBD, qs.u)
        }
        
    }else{
        if (k.OBD == 1){
            qs <- rep(q.OBD, ndose)
        }else if (k.OBD == ndose){
            qs.l <- sort(runif(ndose-1, 0, q.OBD))
            qs <- c(qs.l, q.OBD)
        }else {
            qs.l <- sort(runif(k.OBD-1, 0, q.OBD))
            qs.u <- rep(q.OBD, ndose-k.OBD)
            qs <- c(qs.l, q.OBD, qs.u)
        } 
    }
    
    
    res <- list(
        qs=qs,
        ps=ps,
        k.OBD=k.OBD
    )
    res
}


target <- 0.33
ndose <- 5

ncohort <- 10
cohortsize <- 3

add.args <- list(alp.prior=target, bet.prior=1-target)

mu <- 0.69
run.fn <- function(k){
    set.seed(seeds[k])
    print(k)
    
    p.true <- gen.rand.overdoses(ndose, target, mu2=0.3)
    
    orm.res <- ORM.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
    
    ress <- list(
        orm=orm.res,
        paras=list(p.true=p.true, 
                   mtd=99, 
                   add.args=add.args,
                   target=target,
                   ncohort=ncohort,
                   cohortsize=cohortsize)
    )
    ress
}


nsimu <- 10000
seeds <- 1:nsimu
file.name <- paste0("../results/SMMR-R2/", "MTDOverDose", nsimu, ".RData")
results <- mclapply(1:nsimu, run.fn, mc.cores=15)
post.process.random(results)
save(results, file=file.name)
