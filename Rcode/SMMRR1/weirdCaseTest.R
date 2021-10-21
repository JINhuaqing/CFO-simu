# test whether the case that >gamL and > gamR will happens or not
#rm(list=ls())

library(reshape2)
library(ggplot2)
setwd("C:/Users/JINHU/Documents/ProjectCode/CFO/Rcode")
source("./utilities.R")
source("./ORM_utils.R")

phi <- 0.3
alp.prior <- phi
bet.prior <- 1-phi
nLs <- seq(3, 30, 3)
nCs <- seq(3, 30, 3)
nRs <- seq(0, 30, 3)
comBo.all <- c()
for (nC in nCs){
  for (nL in nLs){
    for (nR in nRs){
        print(c(nL, nC, nR))
        ORsL <- All.OR.table(phi, nL, nC, "L", alp.prior, bet.prior)
        ORsR <- All.OR.table(phi, nC, nR, "R", alp.prior, bet.prior)
        gamL <- optim.gamma.fn(nL, nC, phi, type="L", alp.prior, bet.prior)
        gamR <- optim.gamma.fn(nC, nR, phi, type="R", alp.prior, bet.prior)
        cases <- get.all.cases(ORsL, ORsR, gamL, gamR)
        if (length(cases) != 0){
            Nums <- matrix(rep(c(nL, nC, nR), dim(cases)[1]), nrow=dim(cases)[1], byrow=1)
            comBo <- cbind(cases, Nums)
            names(comBo) <- c("xL", "xC", "xR", "mL", "mC", "mR")
            comBo <- comBo[, c(1, 4, 2, 5, 3, 6)]
            comBo.all <- rbind(comBo.all, comBo)
        }
    }
  }
}

get.all.cases <- function(ORsL, ORsR, gamL, gamR){
    matL <- ORsL > gamL$gamma
    matR <- t(ORsR > gamR$gamma)
    kpidxs <- which(colSums(matL) !=0 & colSums(matR) !=0)
    if (length(kpidxs)==0){
        resAll <- c()
    }else{
      resAll <- c()
      for (idx in kpidxs){
        Ridxs <- which(matR[, idx]) - 1
        Lidxs <- which(matL[, idx]) - 1
        res <- expand.grid(Lidxs, idx-1, Ridxs)
        names(res) <- c("L", "C", "R")
        resAll <- rbind(resAll, res)
      }
    }
    resAll
}

# Test OR and gammas
phi <- 0.3
alp.prior <- phi
bet.prior <- 1-phi
mL <- 3
mC <- 3
mR <- 0
xL <- 0
xC <- 3
xR <- 0
ORL <- OR.values(phi, xL, mL, xC, mC, alp.prior, bet.prior, type="L");ORL
optim.gamma.fn(mL, mC, phi, type="L", alp.prior, bet.prior)

ORR <- OR.values(phi, xC, mC, xR, mR, alp.prior, bet.prior, type="R");ORR
optim.gamma.fn(mC, mR, phi, type="R", alp.prior, bet.prior)

  ORL <- OR.values(phi, 1, 3, 3, 9, alp.prior, bet.prior, type="L");ORL
optim.gamma.fn(nL, nC, phi, type="L", alp.prior, bet.prior)
optim.gamma.fn(nC, nR, phi, type="R", alp.prior, bet.prior)


