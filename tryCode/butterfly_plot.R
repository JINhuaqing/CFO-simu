rm(list=ls())
library(magrittr)


# probabilities of De-escalation, Stay or Escalation
move.dose.post.prob.fn <- function(ys, ns, alp.prior, bet.prior, BOINs, phi){
    # ys, ns: trial results for current 3 dose levels.
    # phi: target DLT rate
    post.spss <- list()
        alps <- ys + alp.prior
        bets <- ns - ys + bet.prior
        p2 <- pbeta(BOINs[2], alps[2], bets[2]) - pbeta(BOINs[1], alps[2], bets[2])
        
        p2.sps <- rbeta(10000, alps[2], bets[2])
        post.spss$C <- p2.sps
        if (!is.na(ys[1])){
            p1.sps <- tbeta.sampler.up(p2.sps, phi, alps[1], bets[1])
            post.spss$L <- p1.sps
        }else{
            p1 <- 0
            post.spss$L <- NA
        }
        if (!is.na(ys[3])){
            p3.sps <- tbeta.sampler.low(p2.sps, phi, alps[3], bets[3])
            post.spss$R <- p3.sps
        }else{
            post.spss$R <- NA
        }
        
    post.spss
}


butterfly.simu.plot.fn <- function(phi, p.true, ncohort=12, 
                                  cohortsize=1, alp.prior=0.1, bet.prior=0.1, m=10){
        # phi: Target DIL rate
        # p.true: True DIL rates under the different dose levels
        # ncohort: The number of cohorts
        # cohortsize: The sample size in each cohort
        # alp.prior, bet.prior: prior parameters
        # m: number of samples to draw for majority vote
        earlystop <- 0
        ndose <- length(p.true)
        cidx <- ceiling(ndose/2)
        
        tys <- rep(0, ndose) # number of responses for different doses.
        tns <- rep(0, ndose) # number of subject for different doses.
        tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
        
        
        BOINs <- BOIN.int(phi)
        post.dists <- list()
        
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
            post.dists[[i]] <- move.dose.post.prob.fn(cys, cns, alp.prior, bet.prior, BOINs, phi)
            idx.chg <- make.move.fn(ps, m=m) - 2
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
        list(MTD=MTD, dose.ns=tns, DLT.ns=tys, post.dists=post.dists)
}


phi <- 0.2
p.true1 <- c(0.1, 0.2, 0.3)
p.true2 <- c(0.05, 0.22, 0.38)

p.true3 <- c(0.2, 0.3, 0.4)
p.true4 <- c(0.18, 0.3, 0.45)

p.true5 <- c(0.07, 0.13, 0.21)
p.true6 <- c(0.04, 0.1, 0.2)

ncohort <- 12
cohortsize <- 1
m <- 10
res <- butterfly.simu.plot.fn(phi, p.true3, ncohort=ncohort, cohortsize=cohortsize, m=m)

post.dists <- res$post.dists
length(post.dists)
BOINint <- BOIN.int(0.2)

setwd("C:/Users/Dell/Downloads")
png("s3.png", height=1200,width=900, res=100)
par(mfrow=c(4, 3))
for (i in 1:length(post.dists)){
    post.dist <- post.dists[[i]]
    L <- post.dist$L
    R <- post.dist$R
    C <- post.dist$C
    plot(density(C), col=2, lty=1, type="l", main="Posterior density")
    if (!sum(is.na(L))){
        lines(density(L), col=3, lty=2, type="l")
    }
    if (!sum(is.na(R))){
        lines(density(R), col=4, lty=3, type="l")
    }
    abline(v=BOINint, col=1, lty=4)
    legend("topright", c('Center', 'Left', "Right"), col=c(2, 3, 4), lty=c(1, 2, 3), cex=0.8)
}
par(mfrow=c(1, 1))
dev.off()
