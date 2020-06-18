#rm(list=ls())
library(magrittr)
library(TruncatedDistributions)


# probabilities of De-escalation, Stay or Escalation
move.dose.post.prob.fn <- function(ys, ns, alp.prior, bet.prior, BOINs){
    # ys, ns: trial results for current 3 dose levels.
    post.spss <- list()
        alps <- ys + alp.prior
        bets <- ns - ys + bet.prior
        p2 <- pbeta(BOINs[2], alps[2], bets[2]) - pbeta(BOINs[1], alps[2], bets[2])
        
        p2.sps <- rbeta(10000, alps[2], bets[2])
        post.spss$C <- p2.sps
        if (!is.na(ys[1])){
            p1.sps <- tbeta.sampler.up(p2.sps, alps[1], bets[1])
            post.spss$L <- p1.sps
        }else{
            p1 <- 0
            post.spss$L <- NA
        }
        if (!is.na(ys[3])){
            p3.sps <- tbeta.sampler.low(p2.sps, alps[3], bets[3])
            post.spss$R <- p3.sps
        }else{
            post.spss$R <- NA
        }
        
    post.spss
}

butterfly.simu.plot.fn <- function(phi, p.true, ncohort=12,  m=10,
                              cohortsize=1, type="BB", add.args=list()){
    # phi: Target DIL rate
    # p.true: True DIL rates under the different dose levels
    # ncohort: The number of cohorts
    # cohortsize: The sample size in each cohort
    # alp.prior, bet.prior: prior parameters
    # m: number of samples to draw for majority vote
    # add.args, list of argments. 
    #    BB: list(alp.prior=, bet.prior)
    #    CRM: list(p.prior)
    #    BB+CRM: list(alp.prior=, bet.prior=, p.prior)
    earlystop <- 0
    ndose <- length(p.true)
    cidx <- ceiling(ndose/2)
    
    tys <- rep(0, ndose) # number of responses for different doses.
    tns <- rep(0, ndose) # number of subject for different doses.
    tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
    
    
    BOINs <- BOIN.int(phi)
    if (type=="BB"){
        post.dists <- list()
    }
    
    
    for (i in 1:ncohort){
        pc <- p.true[cidx] 
        cres <- rbinom(cohortsize, 1, pc)
        
        # update results
        tys[cidx] <- tys[cidx] + sum(cres)
        tns[cidx] <- tns[cidx] + cohortsize
        
        
        
        cy <- tys[cidx]
        cn <- tns[cidx]
        
        add.args.od <- c(list(y=cy, n=cn, tys=tys, tns=tns, cidx=cidx), add.args)
        
        if (overdose.fn(phi, type, add.args.od)){
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
        
        if (type=="BB"){
            alp.prior <- add.args$alp.prior
            bet.prior <- add.args$bet.prior
            ps <- move.dose.probs.fn(cys, cns, alp.prior, bet.prior, BOINs, cover.doses) 
            post.dists[[i]] <- move.dose.post.prob.fn(cys, cns, alp.prior, bet.prior, BOINs)
            idx.chg <- make.move.fn(ps, m=m) - 2
        }else if (type=="CRM"){
            p.prior <- add.args$p.prior
            ps <- move.dose.probs.crm.fn(tys, tns, p.prior, BOINs, cidx, cover.doses)
            idx.chg <- make.move.fn(ps, m=m) - 2
        }else{
            alp.prior <- add.args$alp.prior
            bet.prior <- add.args$bet.prior
            p.prior <- add.args$p.prior
            ps1 <- move.dose.probs.fn(cys, cns, alp.prior, bet.prior, BOINs, cover.doses) 
            ps2 <- move.dose.probs.crm.fn(tys, tns, p.prior, BOINs, cidx, cover.doses)
            idx.chg <- mul.make.move.fn(ps1, ps2, ms=c(m/2, m/2)) - 2
        }
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
    if (type=="BB"){
        return(list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi, post.dists=post.dists))
    }else{
        return(list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi))
    }
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
add.args <- list(alp.prior=0.5, bet.prior=0.5, p.prior=c(0.1, 0.2, 0.3))# Choose this prior only for plotting
res <- butterfly.simu.plot.fn(phi, p.true3, ncohort=ncohort, cohortsize=cohortsize, m=m, type="BB", add.args=add.args)

post.dists <- res$post.dists
length(post.dists)
BOINint <- BOIN.int(phi)

#setwd("/Users/jinhuaqing/Downloads")
setwd("C:/Users/Dell/Downloads")
png("s3_m10.png", height=1200,width=900, res=100)
par(mfrow=c(4, 3))
for (i in 1:length(post.dists)){
    post.dist <- post.dists[[i]]
    L <- post.dist$L
    R <- post.dist$R
    C <- post.dist$C
    fitC <- density(C)
    maxv <- max(fitC$y)
    if (!sum(is.na(L))){
        fitL <- density(L)
        maxL <- max(fitL$y)
        if (maxv <= maxL){
            maxv <- maxL
        }
    }
    if (!sum(is.na(R))){
        fitR <- density(R)
        maxR <- max(fitR$y)
        if (maxv <= maxR){
            maxv <- maxR
        }
    }
    
    plot(density(C), col=2, lty=1, type="l", main="Posterior density", ylim=c(0, 10), xlim=c(0, 1))
    if (!sum(is.na(L))){
        lines(density(L, adjust=1), col=3, lty=2, type="l")
    }
    if (!sum(is.na(R))){
        lines(density(R, adjust=1), col=4, lty=3, type="l")
    }
    abline(v=BOINint, col=1, lty=4)
    legend("topright", c('Center', 'Left', "Right"), col=c(2, 3, 4), lty=c(1, 2, 3), cex=0.8)
}
par(mfrow=c(1, 1))
dev.off()

