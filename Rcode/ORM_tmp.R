rm(list=ls())
setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
source("utilities.R")
source("butterfly_utils.R")


phi <- 0.3
y1 <- 1
n1 <- 10
y2 <- 2
n2 <- 20
alp.prior <- phi
bet.prior <- 1 - phi

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

margin.phi <- function(y, n, lower, upper){
    C <- 1/(upper-lower)
    fn <- function(phi) {
       dbinom(y, n, phi)*C
    }
    integrate(fn, lower=lower, upper=upper)$value
}

vs <- c()
for (y in 0:5){
    v <- margin.phi(y, 5, 0, phi)
    vs <- c(v, vs)
}

OR.values(phi, y1, n1, y2, n2, alp.prior, bet.prior, type="R")
All.OR.table(phi, n1, n2, "R", alp.prior, bet.prior)
All.OR.table(phi, n1, n2, "L", alp.prior, bet.prior)
