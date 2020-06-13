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
tbeta.sampler.low <- function(pcs, phi, alp, bet, pU=0.8){
    n <- length(pcs)
    rtbeta(n, alp, bet, pcs, pU)
}

tbeta.sampler.up <- function(pcs, phi, alp, bet, pL=0){
    n <- length(pcs)
    rtbeta(n, alp, bet, pL, pcs)
}


# desicion make  prob
move.dose.probs <- function(ys, ns, alp.prior, bet.prior, BOINs, phi){
   # ys, ns: trial results for 3 dose level.
   # phi: target DLT rate
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
   if (!is.na(ys[3])){
       p3.sps <- tbeta.sampler.low(p2.sps, phi, alps[3], bets[3])
       p3 <- mean(p3.sps>=BOINs[1] & p3.sps<=BOINs[2])
   }else{
       p3 <- 0
   }
   
   
   ups <- c(p1, p2, p3)
   ps <- ups/sum(ups)
   names(ps) <- c("Dose Low", "Current Dose", "Dose Upper")
   ps
}

make.move <- function(ps){
    cps <- cumsum(ps)
    rv <- runif(1)
    locs <- rv <= cps
    which.max(locs)
}


    

# Simulation under 3 dose levels
butterfly.simu.3levels <- function(phi, p.true, ncohort=12, 
                                   cohortsize=1, alp.prior=0.1, bet.prior=0.1){
   BOINs <- BOIN.int(phi)
   pc <- p.true[2] 
   
   res.init <- rbinom(cohortsize, 1, pc)
   tys <- c(0, sum(res.init), 0)
   tns <- c(0, cohortsize, 0)
   
   ys <- c(0, sum(res.init), 0)
   ns <- c(0, cohortsize, 0)
   
   last.idx <- 2
   for (i in 2:ncohort){
      ps <- move.dose.probs(ys, ns, alp.prior, bet.prior, BOINs, phi) 
      dose.chg <- make.move(ps) - 2
      cidx <- dose.chg + last.idx
      pc <- p.true[cidx]
      
      res <- rbinom(cohortsize, 1, pc)
      tys[cidx] <- tys[cidx] + sum(res)
      tns[cidx] <- tns[cidx] + cohortsize
      
      if (cidx==1){
          ys <- c(NA, tys[1:2])
          ns <- c(NA, tns[1:2])
      }
      if (cidx==2){
          ys <- tys
          ns <- tns
      }
      if (cidx==3){
          ys <- c(tys[2:3], NA)
          ns <- c(tns[2:3], NA)
      }
      last.idx <- cidx
   }
   
       
   final.dose <- select.mtd(phi, tns, tys)$MTD
   list(MTD.dose=final.dose, dose.ns=tns, DLT.ns=tys)
}


nsimu.f <- function(phi, p.true, ncohort=12, cohortsize=1, nsimu=1000){
    MTDs <- list()
    dose.nss <- list()
    DLT.nss <- list()
    for (k in 1:nsimu){
        res <- butterfly.simu.3levels(phi, p.true, ncohort=ncohort, cohortsize=cohortsize)
        dose.nss[[k]] <- res$dose.ns
        
        ncls <- length(res$dose.ns)
        
        MTDvec <- rep(0, ncls)
        MTDvec[res$MTD] <- 1
        MTDs[[k]] <- MTDvec
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
    
    list(MTDs.percent=MTDs.percent, av.dose=dose.mean, av.DLT=DLTs.mean)
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

res <- nsimu.f(target, p.true1, ncohort=ncohort, cohortsize=cohortsize, nsimu=1000)
post.process(res)