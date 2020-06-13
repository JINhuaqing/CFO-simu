rm(list=ls())
library(TruncatedDistributions)
library(BOIN)
# Compute BOIN interval based on target DLT rate
BOIN.int <- function(phi, phiL, phiU){
    if (missing(phiL)) {phiL <- 0.6 * phi}
    if (missing(phiU)) {phiU <- 1.4 * phi}
    
    deltaL <- log((1-phiL)/(1-phi))/log(phi*(1-phiL)/phiL/(1-phi))
    deltaU <- log((1-phi)/(1-phiU))/log(phiU*(1-phi)/phi/(1-phiU))
    res <- c(deltaL, deltaU)
    names(res) <-  c("BOIN lb", "BOIN ub")
    res
}


# truncated Beta sampler
tbeta.sampler.low <- function(pcs, phi, alp, bet, pU=0.9){
    n <- length(pcs)
    rtbeta(n, alp, bet, pcs, pU)
}

tbeta.sampler.up <- function(pcs, phi, alp, bet, pL=0){
    n <- length(pcs)
    rtbeta(n, alp, bet, pL, pcs)
}


# probabilities of De-escalation, Stay or Escalation
move.dose.probs.fn <- function(ys, ns, alp.prior, bet.prior, BOINs, phi, over.doses){
   # ys, ns: trial results for current 3 dose levels.
   # phi: target DLT rate
   # over.doses: vector of size 3. Whether current 3 dose levels are overdosed or not. 
   #             First element is always 0. 
   if (over.doses[2]==1){
       ps <- c(1, 0, 0)
   }else{
       alps <- ys + alp.prior
       bets <- ns - ys + bet.prior
       p2 <- pbeta(BOINs[2], alps[2], bets[2]) - pbeta(BOINs[1], alps[2], bets[2])
       
       p2.sps <- rbeta(10000, alps[2], bets[2])
       if (!is.na(ys[1])){
           p1.sps <- tbeta.sampler.up(p2.sps, phi, alps[1], bets[1])
           p1 <- mean(p1.sps>=BOINs[1] & p1.sps<=BOINs[2])
       }else{
           p1 <- 0
       }
       if (!(is.na(ys[3]) | over.doses[3]==1) ){
           p3.sps <- tbeta.sampler.low(p2.sps, phi, alps[3], bets[3])
           p3 <- mean(p3.sps>=BOINs[1] & p3.sps<=BOINs[2])
       }else{
           p3 <- 0
       }
       
       ups <- c(p1, p2, p3)
       ps <- ups/sum(ups)
   }
   names(ps) <- c("D prob", "S prob", "E prob")
   ps
}


# Make a decison among De-escalation, Stay and Escalation
make.move.fn <- function(ps){
    # ps: Output from move.dose.probs.fn 
    # output: 
    #    action: D--1, S--2, E--3
    cps <- cumsum(ps)
    rv <- runif(1)
    locs <- rv <= cps
    action <- which.max(locs)
    action
}



# posterior probability of pj >= phi given data
post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)

#phi <- 0.2
#p.true <- c(0.2, 0.3, 0.4)
#ncohort <- 12
#cohortsize <- 1
#alp.prior <- 0.1
#bet.prior <- 0.1

# Simulation function
butterfly.simu.fn <- function(phi, p.true, ncohort=12, 
                                   cohortsize=1, alp.prior=0.1, bet.prior=0.1){
    # phi: Target DIL rate
    # p.true: True DIL rates under the different dose levels
    # ncohort: The number of cohorts
    # cohortsize: The sample size in each cohort
    # alp.prior, bet.prior: prior parameters
    earlystop <- 0
    ndose <- length(p.true)
    cidx <- ceiling(ndose/2)
     
    tys <- rep(0, ndose) # number of responses for different doses.
    tns <- rep(0, ndose) # number of subject for different doses.
    tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
    
    
    BOINs <- BOIN.int(phi)
    
    
    for (i in 1:ncohort){
        pc <- p.true[cidx] 
        cres <- rbinom(cohortsize, 1, pc)
        
        # update results
        tys[cidx] <- tys[cidx] + sum(cres)
        tns[cidx] <- tns[cidx] + cohortsize
        
        
        
        cy <- tys[cidx]
        cn <- tns[cidx]
        
        pp <- post.prob.fn(phi, cy, cn, alp.prior, bet.prior)
        
        if ((pp >= 0.95) & (cn>=3)){
            tover.doses[cidx:ndose] <- 1
        }
        
        if (tover.doses[1] == 1){
            earlystop <- 1
            break()
        }
        
        
        # the results for current 3 dose levels
        if (cidx!=1){
            cys <- tys[(cidx-1):(cidx+1)]
            cns <- tns[(cidx-1):(cidx+1)]
            cover.doses <- tover.doses[(cidx-1):(cidx+1)]
        }else{
            cys <- c(NA, tys[1:(cidx+1)])
            cns <- c(NA, tns[1:(cidx+1)])
            cover.doses <- c(NA, tover.doses[1:(cidx+1)])
        }
        
        ps <- move.dose.probs.fn(cys, cns, alp.prior, bet.prior, BOINs, phi, cover.doses) 
        idx.chg <- make.move.fn(ps) - 2
        cidx <- idx.chg + cidx
        #print(c(ps, idx.chg))
        #print(c(tys, tns, tover.doses))
        #print(c(cys, cns, cover.doses))
        
    }
    
        
    if (earlystop==0){
        MTD <- select.mtd(phi, tns, tys)$MTD
    }else{
        MTD <- 99
    }
    list(MTD=MTD, dose.ns=tns, DLT.ns=tys)
}

nsimu.fn <- function(phi, p.true, ncohort=12, cohortsize=1, nsimu=1000){
    MTDs <- list()
    dose.nss <- list()
    DLT.nss <- list()
    for (k in 1:nsimu){
        print(k)
        res <- butterfly.simu.fn(phi, p.true, ncohort=ncohort, cohortsize=cohortsize)
        dose.nss[[k]] <- res$dose.ns
        
        ncls <- length(res$dose.ns)
        
        MTDvec <- rep(0, ncls)
        if (res$MTD !=99){
            MTDvec[res$MTD] <- 1
            MTDs[[k]] <- MTDvec
        }
        DLT.nss[[k]] <- res$DLT.ns
    }
    
    list(MTDs=MTDs, dose=dose.nss, DLTs=DLT.nss)
}


post.process <- function(res){
    MTDs <- res$MTDs
    MTDs <- do.call(rbind, MTDs)
    MTDs.percent <- colMeans(MTDs)
    
    dose <- res$dose
    dose <- do.call(rbind, dose)
    dose.mean <- colMeans(dose)
    
    DLTs <- res$DLTs
    DLTs <- do.call(rbind, DLTs)
    DLTs.mean <- colMeans(DLTs)
    
    list(MTDs.percent=MTDs.percent, av.dose=dose.mean, av.DLT=DLTs.mean, t.dose=sum(dose.mean), t.DLT=sum(DLTs.mean))
}
    


target <- 0.2

p.true1 <- c(0.1, 0.2, 0.3)
p.true2 <- c(0.05, 0.22, 0.38)

p.true3 <- c(0.2, 0.3, 0.4)
p.true4 <- c(0.18, 0.3, 0.45)

p.true5 <- c(0.07, 0.13, 0.21)
p.true6 <- c(0.04, 0.1, 0.2)

ncohort <- 12
cohortsize <- 1

res <- nsimu.fn(target, p.true3, ncohort=ncohort, cohortsize=cohortsize, nsimu=1000)
post.process(res)
