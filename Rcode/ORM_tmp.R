rm(list=ls())
#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
setwd("/Users/jinhuaqing/Documents/Projects_Code/phaseI/Rcode")
source("utilities.R")
source("butterfly_utils.R")
library(magrittr)



Odds.samples <- function(y1, n1, y2, n2, alp.prior, bet.prior){
    alp1 <- alp.prior + y1
    alp2 <- alp.prior + y2
    bet1 <- alp.prior + n1 - y1
    bet2 <- alp.prior + n2 - y2
    sps1 <- c()
    sps2 <- c()
    while (length(sps1)<10000){
        sp1 <- rbeta(1, alp1, bet1)
        sp2 <- rbeta(1, alp2, bet2)
        if (sp1 <= sp2){
            sps1 <- c(sp1, sps1)
            sps2 <- c(sp2, sps2)
        }
    }
    
    list(sps1=sps1, sps2=sps2)
}

prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
    alp1 <- alp.prior + y1
    alp2 <- alp.prior + y2
    bet1 <- alp.prior + n1 - y1
    bet2 <- alp.prior + n2 - y2
    fn.min <- function(x){
        dbeta(x, alp1, bet1)*(1-pbeta(x, alp2, bet2))
    }
    fn.max <- function(x){
        pbeta(x, alp1, bet1)*dbeta(x, alp2, bet2)
    }
    const.min <- integrate(fn.min, lower=0, upper=1)$value
    const.max <- integrate(fn.max, lower=0, upper=1)$value
    p1 <- integrate(fn.min, lower=0, upper=phi)$value/const.min
    p2 <- integrate(fn.max, lower=0, upper=phi)$value/const.max
    
    list(p1=p1, p2=p2)
}


OR.values <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior, type){
    ps <- prob.int(phi, y1, n1, y2, n2, alp.prior, bet.prior)
    if (type=="L"){
        pC <- 1 - ps$p2
        pL <- 1 - ps$p1
        oddsC <- pC/(1-pC)
        oddsL <- pL/(1-pL)
        OR <- oddsC*oddsL
        
    }else if (type=="R"){
        pC <- 1 - ps$p1
        pR <- 1 - ps$p2
        oddsC <- pC/(1-pC)
        oddsR <- pR/(1-pR)
        OR <- (1/oddsC)/oddsR
    }
    return(OR)
}

All.OR.table <- function(phi, n1, n2, type, alp.prior, bet.prior){
   ret.mat <- matrix(rep(0, (n1+1)*(n2+1)), nrow=n1+1)
   for (y1cur in 0:n1){
       for (y2cur in 0:n2){
           ret.mat[y1cur+1, y2cur+1] <- OR.values(phi, y1cur, n1, y2cur, n2, alp.prior, bet.prior, type)
       }
   }
   ret.mat
}

# compute the marginal prob when lower < phiL/phiC/phiR < upper
# i.e., Pr(Y=y|lower<phi<upper)
margin.phi <- function(y, n, lower, upper){
    C <- 1/(upper-lower)
    fn <- function(phi) {
       dbinom(y, n, phi)*C
    }
    integrate(fn, lower=lower, upper=upper)$value
}

# Obtain the table of marginal distribution of (y1, y2) 
# after intergrate out (phi1, phi2)
# under H0 and H1
# H0: phi1=phi, phi < phi2 < 2phi
# H1: phi2=phi, 0   < phi1 < phi
margin.ys.table <- function(n1, n2, phi, hyperthesis){
    if (hyperthesis=="H0"){
        p.y1s <- dbinom(0:n1, n1, phi)
        p.y2s <- sapply(0:n2, margin.phi, n=n2, lower=phi, upper=2*phi)
    }else if (hyperthesis=="H1"){
        p.y1s <- sapply(0:n1, margin.phi, n=n1, lower=0, upper=phi)
        p.y2s <- dbinom(0:n2, n2, phi)
    }
    p.y1s.mat <- matrix(rep(p.y1s, n2+1), nrow=n1+1)
    p.y2s.mat <- matrix(rep(p.y2s, n1+1), nrow=n1+1, byrow=TRUE)
    margin.ys <- p.y1s.mat * p.y2s.mat
    margin.ys
}

# Obtain the optimal gamma for the hypothesis test
optim.gamma.fn <- function(n1, n2, phi, type, alp.prior, bet.prior){
    OR.table <- All.OR.table(phi, n1, n2, type, alp.prior, bet.prior) 
    ys.table.H0 <- margin.ys.table(n1, n2, phi, "H0")
    ys.table.H1 <- margin.ys.table(n1, n2, phi, "H1")
    
    argidx <- order(OR.table)
    sort.OR.table <- OR.table[argidx]
    sort.ys.table.H0 <- ys.table.H0[argidx]
    sort.ys.table.H1 <- ys.table.H1[argidx]
    n.tol <- length(sort.OR.table)
    
    if (type=="L"){
        errs <- rep(0, n.tol-1)
        for (i in 1:(n.tol-1)){
            err1 <- sum(sort.ys.table.H0[1:i])
            err2 <- sum(sort.ys.table.H1[(i+1):n.tol])
            err <- err1 + err2
            errs[i] <- err
        }
        min.err <- min(errs)
        if (min.err > 1){
            gam <- 0
            min.err <- 1
        }else {
            minidx <- which.min(errs)
            gam <- sort.OR.table[minidx]
        }
    }else if (type=='R'){
        errs <- rep(0, n.tol-1)
        for (i in 1:(n.tol-1)){
            err1 <- sum(sort.ys.table.H1[1:i])
            err2 <- sum(sort.ys.table.H0[(i+1):n.tol])
            err <- err1 + err2
            errs[i] <- err
        }
        min.err <- min(errs)
        if (min.err > 1){
            gam <- 0
            min.err <- 1
        }else {
            minidx <- which.min(errs)
            gam <- sort.OR.table[minidx]
        }
        
    }
    list(gamma=gam, min.err=min.err)
}

phi <- 0.2
y1 <- 1
n1 <- 5
y2 <- 3
n2 <- 9
alp.prior <- phi
bet.prior <- 1 - phi

optim.gamma.fn(n1, n2, phi, "L", alp.prior, bet.prior)
OR.values(phi, y1, n1, y2, n2, alp.prior, bet.prior, type="L")
optim.gamma.fn(n1, n2, phi, "R", alp.prior, bet.prior)
OR.values(phi, y1, n1, y2, n2, alp.prior, bet.prior, type="R")

