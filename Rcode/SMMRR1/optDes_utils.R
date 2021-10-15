source("./ORM_utils.R")
## Optimal benchmark 

gen.comp.info.fn <- function(N, ps){
    ## args:
    ##  N: num of patients
    ##  ps: Reponse probs (can be DLT or Eff rates), vector of size K
    
    K <- length(ps)
    us <- runif(N)
    us.mat <- matrix(rep(us, K), nrow=N)
    ps.mat <- matrix(rep(ps, N), nrow=N, byrow=T)
    com.info.sps <- 1* (us.mat <= ps.mat)
    ps.hat <- colMeans(us.mat <= ps.mat)
    
    res <- list(ps.hat=ps.hat, sps=com.info.sps)
    res
}


# phase I
phaseI.opt.simu.fn <- function(nsim, N, phi, p.true){
    ## args:
    ##  nsim: num of simulation times
    ##  N: num of subjects
    ##  phi: target DLT rate
    
    K <- length(p.true)
    MTD.sels <- rep(0, K+1)
    for (i in 1:nsim){
        res <- gen.comp.info.fn(N, p.true)
        ps.hat <- res$ps.hat
        y <- colSums(res$sps)[1]
        add.args <- list(alp.prior=phi, bet.prior=1-phi, y=y, n=N)
        if (overdose.fn(phi, add.args)){
        #if (ps.hat[1] > phi + non.sel){
            MTD <- K+1
        }else{
            MTD <- which.min(abs(ps.hat-phi))
        }
        MTD.sels[MTD] <- MTD.sels[MTD]+1
    }
    MTD.sels.rate <- 100*(MTD.sels/nsim)
    names(MTD.sels.rate) <- c(1:K, "non.sel")
    MTD.sels.rate
}


phase12.opt.simu.fn <- function(nsim, N, phi, phiE, p.true, pE.true){
    ## args:
    ##  nsim: num of simulation times
    ##  N: num of subjects
    ##  phi: target DLT rate
    ##  phiE: minimal eff rate
    
    K <- length(p.true)
    OBD.sels <- rep(0, K+1)
    for (i in 1:nsim){
        res <- gen.comp.info.fn(N, p.true)
        ps.hat <- res$ps.hat
        y <- colSums(res$sps)[1]
        add.args <- list(alp.prior=phi, bet.prior=1-phi, y=y, n=N)
        psE.hat <- gen.comp.info.fn(N, pE.true)$ps.hat
        if (overdose.fn(phi, add.args)){
        #if (ps.hat[1] > phi + non.sel){
            OBD <- K+1
        }else{
            MTD <- which.min(abs(ps.hat-phi))
            GBDs <- 1:MTD
            if (sum(psE.hat[GBDs]<phiE)==MTD){
                OBD <- K+1
            }else{
                OBD <- which.max(psE.hat[GBDs])
            }
        }
        OBD.sels[OBD] <- OBD.sels[OBD]+1
    }
    OBD.sels.rate <- 100*(OBD.sels/nsim)
    names(OBD.sels.rate) <- c(1:K, "non.sel")
    OBD.sels.rate
}



