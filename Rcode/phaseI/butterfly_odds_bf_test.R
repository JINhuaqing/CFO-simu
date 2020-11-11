setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
source("butterfly_utils.R")


p.true <- c(0.1, 0.2, 0.3, 0.4, 0.5)
phi <- 0.2
ncohort <- 10
cohortsize <- 3
add.args <- list(alp.prior=0.2, bet.prior=0.8, p.prior=c(0.1, 0.2, 0.3, 0.4, 0.5), m=10)


earlystop <- 0
ndose <- length(p.true)
cidx <- ceiling(ndose/2)
 
tys <- rep(0, ndose) # number of responses for different doses.
tns <- rep(0, ndose) # number of subject for different doses.
tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes




for (i in 1:ncohort){
    pc <- p.true[cidx] 
    
    # sample from current dose
    cres <- rbinom(cohortsize, 1, pc)
    
    # update results
    tys[cidx] <- tys[cidx] + sum(cres)
    tns[cidx] <- tns[cidx] + cohortsize
    
    
    
    cy <- tys[cidx]
    cn <- tns[cidx]
    
    add.args.od <- c(list(y=cy, n=cn, tys=tys, tns=tns, cidx=cidx), add.args)
    typeod <- "BB+"

    
    if (overdose.fn(phi, typeod, add.args.od)){
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
    
        alp.prior <- add.args$alp.prior
        bet.prior <- add.args$bet.prior
        dec.mk <- make.decision.odds.fn(phi, cys, cns, alp.prior, bet.prior, cover.doses, diag=TRUE)
        #dec.mk <- make.decision.BF.fn(phi, cys, cns, alp.prior, bet.prior, cover.doses, diag=TRUE)
        idx.chg <- dec.mk$final.action - 2


        print(cidx)
        cidx <- idx.chg + cidx
        print(c(cidx, earlystop))
        if (!is.null(dec.mk$oddss))
        print(log(dec.mk$oddss))
        #print(log10(dec.mk$BFs))
        if (!is.null(dec.mk$p2.sps))
        myplot(dec.mk)
}

    
if (earlystop==0){
    MTD <- select.mtd(phi, tns, tys)$MTD
}else{
    MTD <- 99
}
list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi)


myplot <- function(dec.mk){
    plot(density(dec.mk$p2.sps), lwd=2, lty=1, col=1, xlim=c(0, 1), ylim=c(0, 8))
    if (!is.null(dec.mk$p1.sps))
    lines(density(dec.mk$p1.sps), lwd=2, lty=2, col=2)
    if (!is.null(dec.mk$p3.sps))
    lines(density(dec.mk$p3.sps), lwd=2, lty=3, col=3)
    abline(v=phi, lwd=2, col=4, lty=4)
}
