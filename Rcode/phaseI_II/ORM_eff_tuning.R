setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
phi <- 0.30
phiE <- 0.4
p.true <- c(0.02, 0.05, 0.07, 0.1, 0.15)
pE.true <- c(0.05, 0.08, 0.15, 0.30, 0.45)

add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=1, bet.prior.eff=1)


earlystop <- 0
ndose <- length(p.true)
cidx <- 1

tys <- rep(0, ndose) # number of DLT responses for different doses.
txs <- rep(0, ndose) # number of efficacy responses for different doses.
tns <- rep(0, ndose) # number of subject for different doses.
tover.doses <- rep(0, ndose) # Whether each dose is too toxic or not, 1 yes.
tunder.effs <- rep(0, ndose) # Whether the dose is not efficacious or not, 1 yes

cohortsize <- 3

i <- 1
#for (i in 1:ncohort){
    pc <- p.true[cidx] ;pc
    pEc <- pE.true[cidx] ;pEc
    
    # sample from current dose
    cres <- rbinom(cohortsize, 1, pc);cres
    cEres <- rbinom(cohortsize, 1, pEc);cEres
    
    # update results
    tys[cidx] <- tys[cidx] + sum(cres);tys
    txs[cidx] <- txs[cidx] + sum(cEres);txs
    tns[cidx] <- tns[cidx] + cohortsize;tns
    
    
    
    cy <- tys[cidx];cy
    cx <- txs[cidx];cx
    cn <- tns[cidx];cn
    
    add.args <- c(list(y=cy, n=cn, x=cx, tys=tys, txs=txs, tns=tns, cidx=cidx), add.args)
    
    if (overdose.fn(phi, add.args)){
        tover.doses[cidx:ndose] <- 1
    }
    
    tover.doses
    
    if (under.eff.fn(phiE, add.args)){
        tunder.effs[cidx] <- 1
    }
    
    tunder.effs
    
    
    if ((tover.doses[1] == 1) | (sum(tunder.effs)==ndose)){
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
    
    cys
    cns
    cover.doses
    # The up.idx of the admissible set
    up.idx <- make.decision.ORM.fn(phi, cys, cns, add.args$alp.prior, add.args$bet.prior, cover.doses) - 2 + cidx
    up.idx
    if (up.idx == 1){
        cidx <- 1
    }else{
        if (cidx == 1){
            low.idx <- 1
        }else{
            low.idx <- cidx - 1
        }
        ad.xs <- txs[low.idx:up.idx]
        ad.ns <- tns[low.idx:up.idx]
        if (length(ad.xs)==1){
            cidx <- low.idx 
        }else{
            probs <- ORM.Eff.move.probs(ad.xs, ad.ns, add.args$alp.prior.eff, add.args$bet.prior.eff)
            cidx <- which.max(probs) + low.idx - 1
            #cidx <- make.move.fn(probs, m=10) + low.idx - 1
            print(probs)
            print(ad.xs)
            print(ad.ns)
        }
        
    }    
    cidx
    print(tns)
    print("--------------------------------------------------------------------------------")
    
#}

MTD <- select.mtd(phi, tns, tys)$MTD;MTD

OBD.probs <- ORM.Eff.move.probs(txs[1:MTD], tns[1:MTD], add.args$alp.prior.eff, add.args$bet.prior.eff)
OBD.probs
OBD <- which.max(OBD.probs)
OBD


phiE <- 0.5
post.prob.fn(phiE, 3, 12, 1, 1)
post.prob.fn(phiE, 1, 6, 1, 1)
post.prob.fn(phiE, 0, 3, 1, 1)
ad.xs <- c(3, 1, 0)
ad.ns <- c(12, 6, 3)
ORM.Eff.move.probs(ad.xs, ad.ns, add.args$alp.prior.eff, add.args$bet.prior.eff)

tns <- rep(6, 5)
tys <- rep(6, 5)

select.mtd(phi, tns,tys)$MTD
