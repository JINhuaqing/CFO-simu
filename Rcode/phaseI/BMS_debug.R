library(magrittr)
setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
source("butterfly_utils.R")

sample.ML <- function(phi, up.lim, N=10000){
    # phi: Target toxicity
    # N : the number of samples
    if (missing(up.lim)){
        up.lim <- 0.6
    }
    sps.all <- runif(2*N, min=phi, max=up.lim)
    sps.all <- matrix(sps.all, ncol=2)
    reorder.sps <- t(apply(sps.all, 1, sort))
    p2.sps <- reorder.sps[, 1]
    p3.sps <- reorder.sps[, 2]
    list(p2.sps=p2.sps, p3.sps=p3.sps)
}

sample.MC <- function(phi, up.lim, N=10000){
    # phi: Target toxicity
    # N : the number of samples
    if (missing(up.lim)){
        up.lim <- 0.6
    }
    
    p1.sps <- runif(N, min=0, max=phi)
    p3.sps <- runif(N, min=phi, max=up.lim)
    list(p1.sps=p1.sps, p3.sps=p3.sps)
}


sample.MR <- function(phi, N=10000){
    # phi: Target toxicity
    # N : the number of samples
    sps.all <- runif(2*N, min=0, max=phi)
    sps.all <- matrix(sps.all, ncol=2)
    reorder.sps <- t(apply(sps.all, 1, sort))
    p1.sps <- reorder.sps[, 1]
    p2.sps <- reorder.sps[, 2]
    list(p1.sps=p1.sps, p2.sps=p2.sps)
}


BMS.prob3 <- function(phi, cys, cns, model="L"){
    yL <- cys[1]
    yC <- cys[2]
    yR <- cys[3]
    
    nL <- cns[1]
    nC <- cns[2]
    nR <- cns[3]
    
    if (model=="L"){
        sps <- sample.ML(phi) 
        p2.sps <- sps$p2.sps
        p3.sps <- sps$p3.sps
        itm1 <- phi**yL * (1-phi)**(nL-yL)
        itm2.all <- p2.sps**yC*(1-p2.sps)**(nC-yC) * p3.sps**yR*(1-p3.sps)**(nR-yR)
        itm2 <- mean(itm2.all)
        return(itm1*itm2)
    }else if (model=="C"){
        sps <- sample.MC(phi) 
        p1.sps <- sps$p1.sps
        p3.sps <- sps$p3.sps
        itm1 <- phi**yC * (1-phi)**(nC-yC)
        itm2.all <- p1.sps**yL*(1-p1.sps)**(nL-yL) * p3.sps**yR*(1-p3.sps)**(nR-yR)
        itm2 <- mean(itm2.all)
        return(itm1*itm2)
    }else if (model=="R"){
        sps <- sample.MR(phi) 
        p1.sps <- sps$p1.sps
        p2.sps <- sps$p2.sps
        itm1 <- phi**yR * (1-phi)**(nR-yR)
        itm2.all <- p1.sps**yL*(1-p1.sps)**(nL-yL) * p2.sps**yC*(1-p2.sps)**(nC-yC)
        itm2 <- mean(itm2.all)
        return(itm1*itm2)
    }
    
}

BMS.prob2L <- function(phi, cys, cns, model="L"){
    yL <- cys[1]
    yC <- cys[2]
    
    nL <- cns[1]
    nC <- cns[2]
    
    if (model=="L"){
        p2.sps <- runif(10000, min=phi, max=0.6)
        itm1 <- phi**yL * (1-phi)**(nL-yL)
        itm2.all <- p2.sps**yC*(1-p2.sps)**(nC-yC)
        itm2 <- mean(itm2.all)
        return(itm1*itm2)
    }else if (model=="C"){
        p1.sps <- runif(10000, min=0, max=phi)
        itm1 <- phi**yC * (1-phi)**(nC-yC)
        itm2.all <- p1.sps**yL*(1-p1.sps)**(nL-yL)
        itm2 <- mean(itm2.all)
        return(itm1*itm2)
    }
}


BMS.prob2R <- function(phi, cys, cns, model="C"){
    yC <- cys[2]
    yR <- cys[3]
    
    nC <- cns[2]
    nR <- cns[3]
    
    if (model=="C"){
        p3.sps <- runif(10000, min=phi, max=0.6)
        itm1 <- phi**yC * (1-phi)**(nC-yC)
        itm2.all <- p3.sps**yR*(1-p3.sps)**(nR-yR)
        itm2 <- mean(itm2.all)
        return(itm1*itm2)
    }else if (model=="R"){
        p2.sps <- runif(10000, min=0, max=phi)
        itm1 <- phi**yR * (1-phi)**(nR-yR)
        itm2.all <- p1.sps**yL*(1-p1.sps)**(nL-yL) * p2.sps**yC*(1-p2.sps)**(nC-yC)
        itm2 <- mean(itm2.all)
        return(itm1*itm2)
    }
    
}

BMS.model.prob <- function(phi, cys, cns, cover.doses){
    if (cover.doses[2]==1){
        ps <- c(1, 0, 0)
    }else{
        if (is.na(cys[1]) & (cover.doses[3]==1)){
            ps <- c(0, 1, 0)
        }
        
        if (is.na(cys[1]) & (!(cover.doses[3]==1))){
            pC <- BMS.prob2R(phi, cys, cns, "C")    
            pR <- BMS.prob2R(phi, cys, cns, "R")    
            ups <- c(0, pC, pR)
            ps <- ups/sum(ups)
        }
        
        if (is.na(cys[3]) | (cover.doses[3]==1)){
            pL <- BMS.prob2L(phi, cys, cns, "L")    
            pC <- BMS.prob2L(phi, cys, cns, "C")    
            ups <- c(pL, pC, 0)
            ps <- ups/sum(ups)
        }
        
        if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
            pL <- BMS.prob3(phi, cys, cns, "L")    
            pC <- BMS.prob3(phi, cys, cns, "C")    
            pR <- BMS.prob3(phi, cys, cns, "R")    
            ups <- c(pL, pC, pR)
            ps <- ups/sum(ups)
        }
        
    }
    names(ps) <- c("D prob", "S prob", "E prob")
    ps
}



p.true <- c(0.1, 0.2, 0.3, 0.4, 0.5)
phi <- 0.2
ncohort <- 10
cohortsize <- 3
add.args <- list(alp.prior=phi, bet.prior=1-phi, p.prior=c(0.1, 0.2, 0.3, 0.4, 0.5), m=1)


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
    ps <- BMS.model.prob(phi, cys, cns, cover.doses)
    final.action <- make.move.fn(ps, m=add.args$m)
    idx.chg <- final.action - 2
    
    print(cidx)
    cidx <- idx.chg + cidx
    print(c(cidx, earlystop))
    print(ps)
}


if (earlystop==0){
    MTD <- select.mtd(phi, tns, tys)$MTD
}else{
    MTD <- 99
}
list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi)


