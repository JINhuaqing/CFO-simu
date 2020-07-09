library(TruncatedDistributions)
library(BOIN)
source("utilities.R")
#set.seed(10)
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

# Compute adaptive BOIN interval based on target DLT rate
aBOIN.int <- function(phi, n, N){
    # n: current number of cohorts
    # N: Total number of cohorts
    c <- 0.2
    cont <- 0.4 + (1/sqrt(n)-1)*(c-0.4)/(1/sqrt(N)- 1)
    
    phiL <- (1- cont)* phi
    phiU <- (1+ cont)* phi
    
    deltaL <- log((1-phiL)/(1-phi))/log(phi*(1-phiL)/phiL/(1-phi))
    deltaU <- log((1-phi)/(1-phiU))/log(phiU*(1-phi)/phi/(1-phiU))
    res <- c(deltaL, deltaU)
    names(res) <-  c("BOIN lb", "BOIN ub")
    res
}



# truncated Beta sampler
tbeta.sampler.low <- function(pcs, alp, bet, pU=0.9){
    n <- length(pcs)
    rtbeta(n, alp, bet, pcs, pU)
}

tbeta.sampler.up <- function(pcs, alp, bet, pL=0){
    n <- length(pcs)
    rtbeta(n, alp, bet, pL, pcs)
}


# probabilities of De-escalation, Stay or Escalation
move.dose.probs.fn <- function(ys, ns, alp.prior, bet.prior, BOINs, over.doses){
   # ys, ns: trial results for current 3 dose levels.
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
           p1.sps <- tbeta.sampler.up(p2.sps, alps[1], bets[1])
           p1 <- mean(p1.sps>=BOINs[1] & p1.sps<=BOINs[2])
       }else{
           p1 <- 0
       }
       if (!(is.na(ys[3]) | over.doses[3]==1) ){
           p3.sps <- tbeta.sampler.low(p2.sps, alps[3], bets[3])
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

# The posterior density of alpha for CRM model under power funtion model up to a constant
post.density.crm.powerfn <- function(alpha, p.prior, tys, tns) {
    sigma2 <- 2
    lik <- 1
    for(i in 1:length(tns))
    {
        pi <- p.prior[i]^(exp(alpha));
        lik <- lik*pi^tys[i]*(1-pi)^(tns[i]-tys[i]);
    }
    return(lik*exp(-0.5*alpha*alpha/sigma2));
}

crm.power.limfn <- function(cv, pj.prior){
    num <- log(cv)
    den <- log(pj.prior)
    return(log(num/den))
}

# probabilities of De-escalation, Stay or Escalation with CRM + interval
#move.dose.probs.crm.fn(tys=c(1, 1, 1), tns=c(3, 2, 6), p.prior=c(0.1, 0.2, 0.3), BOIN.int(0.5), 2, c(0, 0, 0))
move.dose.probs.crm.fn <- function(tys, tns, p.prior, BOINs, cidx, over.doses){
    # calculate posterior mean of toxicity probability at each dose leavel
    if (over.doses[2]==1){
        ps <- c(1, 0, 0)
    }else{
        #marginal <- integrate(post.density.crm.powerfn,lower=-Inf,upper=Inf,p.prior,tys,tns)$value
        U <- BOINs[2]
        L <- BOINs[1]
        p2.prior <- p.prior[cidx]
        p2.num <- integrate(post.density.crm.powerfn,
                            lower=crm.power.limfn(U, p2.prior),
                            upper=crm.power.limfn(L, p2.prior),
                            p.prior,tys,tns)$value
        if (cidx!=1){
            p1.prior <- p.prior[cidx-1]
            p1.num <- integrate(post.density.crm.powerfn,
                            lower=crm.power.limfn(U, p1.prior),
                            upper=crm.power.limfn(L, p1.prior),
                            p.prior,tys,tns)$value
        }else{
            p1.num <- 0
        }
        if (!(cidx==length(tys) | over.doses[3]==1) ){
            p3.prior <- p.prior[cidx+1]
            p3.num <- integrate(post.density.crm.powerfn,
                            lower=crm.power.limfn(U, p3.prior),
                            upper=crm.power.limfn(L, p3.prior),
                            p.prior,tys,tns)$value
        }else{
            p3.num <- 0
        }
        
        ups <- c(p1.num, p2.num, p3.num)
        ps <- ups/sum(ups)
    }
    names(ps) <- c("D prob", "S prob", "E prob")
    ps
}


# Make a decison among De-escalation, Stay and Escalation
make.move.fn <- function(ps, m=10){
    # ps: Output from move.dose.probs.fn 
    # m: number of samples to draw for majority vote
    # output: 
    #    action: D--1, S--2, E--3
    cps <- cumsum(ps)
    rvs <- runif(m)
    locs <- sapply(rvs, function(rv)rv<= cps)
    actions <- apply(locs, 2, which.max)
    if (m==1){
        final.action <- actions
    }else{
        res <- sapply(1:3, function(i)actions==i)
        res <- colSums(res)
        final.action <- which.max(res)
    }
    final.action
}

mul.make.move.fn <- function(ps1, ps2, ms=c(5, 5)){
    # ps: Output from move.dose.probs.fn 
    # m: number of samples to draw for majority vote
    # output: 
    #    action: D--1, S--2, E--3
    m1 <- ms[1]
    m2 <- ms[2]
    cps1 <- cumsum(ps1)
    cps2 <- cumsum(ps2)
    rvs1 <- runif(m1)
    rvs2 <- runif(m2)
    locs1 <- sapply(rvs1, function(rv)rv<= cps1)
    locs2 <- sapply(rvs2, function(rv)rv<= cps2)
    locs <- cbind(locs1, locs2)
    actions <- apply(locs, 2, which.max)
    res <- sapply(1:3, function(i)actions==i)
    res <- colSums(res)
    final.action <- which.max(res)
    final.action
}




# posterior probability of pj >= phi given data
post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
}

overdose.fn <- function(phi, type="BB", add.args=list()){
        args <- c(list(phi=phi), add.args)
    if (type=="BB"){
        y <- add.args$y
        n <- add.args$n
        alp.prior <- add.args$alp.prior
        bet.prior <- add.args$bet.prior
        pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
        if ((pp >= 0.95) & (add.args$n>=3)){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }else if (type=="CRM"){
        tys <- add.args$tys
        tns <- add.args$tns
        p.prior <- add.args$p.prior
        cidx <- add.args$cidx
        pc.prior <- p.prior[cidx]
        cv <- crm.power.limfn(phi, pc.prior)
        marginal <- integrate(post.density.crm.powerfn,lower=-Inf,upper=Inf,p.prior,tys,tns)$value
        num <- integrate(post.density.crm.powerfn,lower=-Inf,upper=cv,p.prior,tys,tns)$value
        pp <- num/marginal
        if ((pp >= 0.95) & (tns[cidx]>=3)){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }else if (type=="BB+CRM"){
        y <- add.args$y
        n <- add.args$n
        alp.prior <- add.args$alp.prior
        bet.prior <- add.args$bet.prior
        pp1 <- post.prob.fn(phi, y, n, alp.prior, bet.prior)

        tys <- add.args$tys
        tns <- add.args$tns
        p.prior <- add.args$p.prior
        cidx <- add.args$cidx
        pc.prior <- p.prior[cidx]
        cv <- crm.power.limfn(phi, pc.prior)
        marginal <- integrate(post.density.crm.powerfn,lower=-Inf,upper=Inf,p.prior,tys,tns)$value
        num <- integrate(post.density.crm.powerfn,lower=-Inf,upper=cv,p.prior,tys,tns)$value
        pp2 <- num/marginal
        if (((pp1 >= 0.95) | (pp2>=0.95)) & (add.args$n>=3)){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }else if (type=="BB+"){
        y <- add.args$y
        n <- add.args$n
        alp.prior <- add.args$alp.prior
        bet.prior <- add.args$bet.prior
        pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
        if ((pp >= 0.95) & (add.args$n>=3)){
            return(TRUE)
        }else{
            return(FALSE)
        }
    }
}

#phi <- 0.2
#p.true <- c(0.2, 0.3, 0.4)
#ncohort <- 12
#cohortsize <- 1
#alp.prior <- 0.1
#bet.prior <- 0.1

# Function to make a decision based on the odds close to 1 most
make.decision.odds.fn <- function(phi, ys, ns, alp.prior, bet.prior, over.doses, diag=FALSE){
   if (over.doses[2]==1){
       if (diag){
           rev <- list(final.action=1,
                       p1.sps=NULL,
                       p2.sps=NULL,
                       p3.sps=NULL,
                       oddss=c(NULL, NULL, NULL))
           return(rev)
       }else{
           return(1)
       }
   }else{
       alps <- ys + alp.prior
       bets <- ns - ys + bet.prior
       p2.under <- pbeta(phi, alps[2], bets[2])
       p2.odds <- (1-p2.under)/p2.under

       p2.sps <- rbeta(10000, alps[2], bets[2])
       p1.sps <- NULL
       p3.sps <- NULL
       if (!is.na(ys[1])){
           p1.sps <- tbeta.sampler.up(p2.sps, alps[1], bets[1])
           p1.over <- mean(p1.sps>=phi)
           p1.odds <- p1.over / (1 - p1.over)
       }else{
           p1.odds <- Inf
       }
       if (!(is.na(ys[3]) | over.doses[3]==1) ){
           p3.sps <- tbeta.sampler.low(p2.sps, alps[3], bets[3])
           p3.over <- mean(p3.sps>=phi)
           p3.odds <- p3.over / (1 - p3.over)
       }else{
           p3.odds <- Inf
       }
       
       oddss <- c(p1.odds, p2.odds, p3.odds)
       final.action <- which.min(abs(log(oddss)))
       if (diag){
           rev <- list(final.action=final.action, 
                       p1.sps=p1.sps, 
                       p2.sps=p2.sps,
                       p3.sps=p3.sps,
                       oddss=oddss)
           return(rev)
       }else{
           return(final.action)
       }
   }
}


# Function to make a decision based on the 
# Pr(\theta_R>\phi|D)/Pr(\theta_C>\phi|D) and 
# Pr(\theta_L<\phi|D)/Pr(\theta_C<\phi|D).
make.decision.BF.fn <- function(phi, ys, ns, alp.prior, bet.prior, over.doses, diag=FALSE){
   if (over.doses[2]==1){
       if (diag){
           rev <- list(final.action=1,
                       p1.sps=NULL,
                       p2.sps=NULL,
                       p3.sps=NULL,
                       BFs=c(NULL, NULL))
           return(rev)
       }else{
           return(1)
       }
   }else{
       alps <- ys + alp.prior
       bets <- ns - ys + bet.prior
       p2.under <- pbeta(phi, alps[2], bets[2])
       p2.over <- 1 - p2.under

       p2.sps <- rbeta(10000, alps[2], bets[2])
       p1.sps <- NULL
       p3.sps <- NULL
       if (!is.na(ys[1])){
           p1.sps <- tbeta.sampler.up(p2.sps, alps[1], bets[1])
           p1.under <- mean(p1.sps<=phi)
           p.under.BF <- p1.under/p2.under
       }else{
           p.under.BF <- 0.1
       }
       if (!(is.na(ys[3]) | over.doses[3]==1) ){
           p3.sps <- tbeta.sampler.low(p2.sps, alps[3], bets[3])
           p3.over <- mean(p3.sps>=phi)
           p.over.BF <- p3.over/p2.over
       }else{
           p.over.BF <- 0.1
       }
       CV1 <- 0.5
       CV2 <- 0.5
       dec.over <- (log10(p.over.BF)>CV1) 
       dec.under <-  (log10(p.under.BF)>CV2)
       if (dec.over & (!dec.under)){
           final.action <- 3
       }else if (dec.under & (!dec.over)){
           final.action <- 1
       }else{
           final.action <- 2
       }
       if (diag){
           rev <- list(final.action=final.action, 
                       p1.sps=p1.sps, 
                       p2.sps=p2.sps,
                       p3.sps=p3.sps,
                       BFs=c(p.under.BF, p.over.BF))
           return(rev)
       }else{
           return(final.action)
       }
   }
}
 
# Function to make a decision based on the odds close to 1 most
# with ditribution p1 < p2 < p3 
# Very bad performance
make.decision.odds.fn2 <- function(phi, ys, ns, alp.prior, bet.prior, over.doses){
   if (over.doses[2]==1){
       return(1)
   }else{
       alps <- ys + alp.prior
       bets <- ns - ys + bet.prior
       p2.under <- pbeta(phi, alps[2], bets[2])
       p2.odds <- (1-p2.under)/p2.under
       p2.sps <- rbeta(10000, alps[2], bets[2])
       if (is.na(ys[1])){
           p2.sps <- c() 
           p3.sps <- c()
           while (length(p2.sps)<=10000){
               p2.sp <- rbeta(1, alps[2], bets[2])
               p3.sp <- rbeta(1, alps[3], bets[3])
               if (p2.sp < p3.sp){
                   p2.sps <- c(p2.sps, p2.sp)
                   p3.sps <- c(p3.sps, p3.sp)
               }

           }
           p2.over <- mean(p2.sps>=phi)
           p2.odds <- p2.over / (1 - p2.over)
           p3.over <- mean(p3.sps>=phi)
           p3.odds <- p3.over / (1 - p3.over)
           p1.odds <- Inf
       }else if (is.na(ys[3])){
           p1.sps <- c() 
           p2.sps <- c()
           while (length(p2.sps)<=10000){
               p1.sp <- rbeta(1, alps[1], bets[1])
               p2.sp <- rbeta(1, alps[2], bets[2])
               if (p1.sp < p2.sp){
                   p1.sps <- c(p1.sps, p1.sp)
                   p2.sps <- c(p2.sps, p2.sp)
               }

           }
           p1.over <- mean(p1.sps>=phi)
           p1.odds <- p1.over / (1 - p1.over)
           p2.over <- mean(p2.sps>=phi)
           p2.odds <- p2.over / (1 - p2.over)
           p3.odds <- Inf
       }else{
           p1.sps <- c() 
           p2.sps <- c()
           p3.sps <- c()
           while (length(p2.sps)<=10000){
               p1.sp <- rbeta(1, alps[1], bets[1])
               p2.sp <- rbeta(1, alps[2], bets[2])
               p3.sp <- rbeta(1, alps[3], bets[3])
               if (p1.sp < p2.sp & p2.sp < p3.sp){
                   p1.sps <- c(p1.sps, p1.sp)
                   p2.sps <- c(p2.sps, p2.sp)
                   p3.sps <- c(p3.sps, p3.sp)
               }

           }
           p1.over <- mean(p1.sps>=phi)
           p1.odds <- p1.over / (1 - p1.over)
           p2.over <- mean(p2.sps>=phi)
           p2.odds <- p2.over / (1 - p2.over)
           p3.over <- mean(p3.sps>=phi)
           p3.odds <- p3.over / (1 - p3.over)
       }

       if (over.doses[3]==1){
           p3.odds <- Inf
       }
       
       oddss <- c(p1.odds, p2.odds, p3.odds)
       final.action <- which.min(abs(log(oddss)))
       return(final.action)
   }
}

# Simulation function
butterfly.simu.fn <- function(phi, p.true, ncohort=12,
                              cohortsize=1, type="BB", add.args=list()){
    # phi: Target DIL rate
    # p.true: True DIL rates under the different dose levels
    # ncohort: The number of cohorts
    # cohortsize: The sample size in each cohort
    # alp.prior, bet.prior: prior parameters
    # add.args, list of argments. 
    #    BB: list(alp.prior=, bet.prior, m)
    #    CRM: list(p.prior, m)
    #    BB+CRM: list(alp.prior=, bet.prior=, p.prior, m)
    #           m: number of samples to draw for majority vote
    earlystop <- 0
    ndose <- length(p.true)
    cidx <- ceiling(ndose/2)
     
    tys <- rep(0, ndose) # number of responses for different doses.
    tns <- rep(0, ndose) # number of subject for different doses.
    tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
    
    
    
    
    for (i in 1:ncohort){
        BOINs <- BOIN.int(phi)
        #BOINs <- aBOIN.int(phi, i, ncohort)
        pc <- p.true[cidx] 
        
        # sample from current dose
        cres <- rbinom(cohortsize, 1, pc)
        
        # update results
        tys[cidx] <- tys[cidx] + sum(cres)
        tns[cidx] <- tns[cidx] + cohortsize
        
        
        
        cy <- tys[cidx]
        cn <- tns[cidx]
        
        add.args.od <- c(list(y=cy, n=cn, tys=tys, tns=tns, cidx=cidx), add.args)
        if (type!="CRM" | type!="BB+CRM"){
            typeod <- "BB+"
        }else{
            typeod <- type
        }

        
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
        
        if (type=="BB"){
            m <- add.args$m
            alp.prior <- add.args$alp.prior
            bet.prior <- add.args$bet.prior
            ps <- move.dose.probs.fn(cys, cns, alp.prior, bet.prior, BOINs, cover.doses) 
            idx.chg <- make.move.fn(ps, m=m) - 2
        }else if (type=="CRM"){
            m <- add.args$m
            p.prior <- add.args$p.prior
            ps <- move.dose.probs.crm.fn(tys, tns, p.prior, BOINs, cidx, cover.doses)
            idx.chg <- make.move.fn(ps, m=m) - 2
        }else if (type=="BB+CRM"){
            m <- add.args$m
            alp.prior <- add.args$alp.prior
            bet.prior <- add.args$bet.prior
            p.prior <- add.args$p.prior
            ps1 <- move.dose.probs.fn(cys, cns, alp.prior, bet.prior, BOINs, cover.doses) 
            ps2 <- move.dose.probs.crm.fn(tys, tns, p.prior, BOINs, cidx, cover.doses)
            idx.chg <- mul.make.move.fn(ps1, ps2, ms=c(m/2, m/2)) - 2
        }else if (type=="Odds"){
            alp.prior <- add.args$alp.prior
            bet.prior <- add.args$bet.prior
            idx.chg <- make.decision.odds.fn(phi, cys, cns, alp.prior, bet.prior, cover.doses) - 2
        }else if (type=="BF"){
            alp.prior <- add.args$alp.prior
            bet.prior <- add.args$bet.prior
            idx.chg <- make.decision.BF.fn(phi, cys, cns, alp.prior, bet.prior, cover.doses) - 2

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
    list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi)
}

nsimu.fn <- function(phi, p.true, ncohort=12, cohortsize=1, nsimu=1000, m=10){
    MTDs <- list()
    dose.nss <- list()
    DLT.nss <- list()
    add.args <- list(p.prior=c(0.1, 0.2, 0.3))
    #add.args <- list(alp.prior=0.1, bet.prior=0.1)
    for (k in 1:nsimu){
        print(k)
        res <- butterfly.simu.fn(phi, p.true, type="CRM",
                                 ncohort=ncohort, cohortsize=cohortsize, m=m, add.args=add.args)
        dose.nss[[k]] <- res$dose.ns
        
        ncls <- length(res$dose.ns)
        
        MTDvec <- rep(0, ncls)
        if (res$MTD !=99){
            MTDvec[res$MTD] <- 1
            MTDs[[k]] <- MTDvec
        }
        DLT.nss[[k]] <- res$DLT.ns
    }
    
    list(MTDs=MTDs, dose=dose.nss, DLTs=DLT.nss, p.true=p.true, target=phi)
}

# simulation function for random scenarios
nsimu.rand.fn <- function(phi, ncohort=30, cohortsize=1, nsimu=1000, m=10){
    MTDs <- rep(0, nsimu)
    tMTDs <- rep(0, nsimu)
    dose.nss <- list()
    DLT.nss <- list()
    add.args <- list(p.prior=c(0.1, 0.2, 0.3))
    #add.args <- list(alp.prior=0.1, bet.prior=0.1)
    for (k in 1:nsimu){
        print(k)
        rand.case <- gen.rand.doses(3, phi, 0.65, 0.45)
        res <- butterfly.simu.fn(phi, rand.case$p.true, type="CRM",
                                 ncohort=ncohort, cohortsize=cohortsize, m=m, add.args=add.args)
        dose.nss[[k]] <- res$dose.ns
        
        ncls <- length(res$dose.ns)
        
        tMTDs[k] <- rand.case$mtd.level
        if (res$MTD !=99){
            MTDs[k] <- res$MTD
        }
        DLT.nss[[k]] <- res$DLT.ns
    }
    
    list(MTDs=MTDs, tMTDs=tMTDs, dose=dose.nss, DLTs=DLT.nss, target=phi)
}

post.process <- function(res){
    MTDs <- res$MTDs
    MTDs <- do.call(rbind, MTDs)
    MTDs.percent <- colSums(MTDs)/length(res$dose)
    
    dose <- res$dose
    dose <- do.call(rbind, dose)
    dose.mean <- colMeans(dose)
    
    DLTs <- res$DLTs
    DLTs <- do.call(rbind, DLTs)
    DLTs.mean <- colMeans(DLTs)
    
    list(MTDs.percent=MTDs.percent, av.dose=dose.mean, av.DLT=DLTs.mean, t.dose=sum(dose.mean), t.DLT=sum(DLTs.mean), p.true=res$p.true, target=res$target)
}
    
post.process.raw <- function(ress){
    MTDs <- list()
    dose <- list()
    DLTs <- list()
    K <- length(ress)
    for (k in 1:K){
        res <- ress[[k]]
        dose[[k]] <- res$dose.ns
        
        ncls <- length(res$dose.ns)
        
        MTDvec <- rep(0, ncls)
        if (res$MTD !=99){
            MTDvec[res$MTD] <- 1
            MTDs[[k]] <- MTDvec
        }
        DLTs[[k]] <- res$DLT.ns
    }

    MTDs <- do.call(rbind, MTDs)
    MTDs.percent <- colSums(MTDs)/K
    
    dose <- do.call(rbind, dose)
    dose.mean <- colMeans(dose)
    
    DLTs <- do.call(rbind, DLTs)
    DLTs.mean <- colMeans(DLTs)
    
    list(MTDs.percent=MTDs.percent, av.dose=dose.mean, av.DLT=DLTs.mean, t.dose=sum(dose.mean), t.DLT=sum(DLTs.mean), p.true=res$p.true, target=res$target)
}


# R code Related to BMS method

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
