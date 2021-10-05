rm(list=ls())
setwd("C:/Users/JINHU/Documents/ProjectCode/CFO/")
source("Rcode/utilities.R")

phaseI_II.scenario.plot.fn <- function(p.true, pE.true, main="Plateau", OBD=NULL){
    ndose <- length(p.true)
    lvls <- 1:ndose
    cols <- c("black", "red")
    lwds <- c(2, 2)
    ltys <- c(1, 2)
    pchs <- c(15, 16)
    plot(1:ndose, p.true, ylab = "Probability", xlab="Dose level", ylim=c(0, 0.8), main=main, xaxt="n", 
         type="b", 
         lwd=lwds[1], col=cols[1], lty=ltys[1], pch=pchs[1])
    if (is.null(OBD)){
        axis(1, at=lvls, lvls)
    }else{
        axis(1, at=OBD, paste0(OBD, "*"), col.axis="red", font.axis=2)
        axis(1, at=lvls[-OBD], lvls[-OBD])
    }
    lines(1:ndose, pE.true, type="b", lwd=lwds[2], col=cols[2], lty=ltys[2], pch=pchs[2])
    #    legend("topleft", c("DLT rate", "Efficacy rate"), col=cols, lty=ltys, lwd=lwds, pch=pchs)
}

phi <- 0.3
psi <- 0.2
psi.U <- 0.7
mu1 <- 0.55
mu2 <- 0.55
ndose <- 5


gen.rand.doses.plateau <- function(ndose, phi, psi, psi.U, mu1, mu2){
    k.OBD <- sample.int(ndose, 1)
    k.MTD <- sample(k.OBD:ndose, 1)
    
    ps <- gen.rand.doses(ndose, phi, mu1, mu2, MTD=k.MTD)
    
    q.OBD <- runif(1, psi, psi.U)
    if (k.OBD == 1){
        qs <- rep(q.OBD, ndose)
    }else if (k.OBD == ndose){
        qs.l <- sort(runif(ndose-1, 0, q.OBD))
        qs <- c(qs.l, q.OBD)
    }else {
        qs.l <- sort(runif(k.OBD-1, 0, q.OBD))
        qs.u <- rep(q.OBD, ndose-k.OBD)
        qs <- c(qs.l, q.OBD, qs.u)
    }
    
    res <- list(
        qs=qs,
        ps=ps$p.trues, 
        k.MTD=k.MTD, 
        k.OBD=k.OBD
    )
    res
}

gen.rand.doses.umbrella <- function(ndose, phi, psi, psi.U, mu1, mu2){
    k.OBD <- sample.int(ndose, 1)
    k.MTD <- sample(k.OBD:ndose, 1)
    
    ps <- gen.rand.doses(ndose, phi, mu1, mu2, MTD=k.MTD)
    
    q.OBD <- runif(1, psi, psi.U)
    if (k.OBD == 1){
        qs.u <- sort(runif(ndose-1, 0, q.OBD), decreasing=TRUE)
        qs <- c(q.OBD, qs.u)
    }else if (k.OBD == ndose){
        qs.l <- sort(runif(ndose-1, 0, q.OBD))
        qs <- c(qs.l, q.OBD)
    }else {
        qs.l <- sort(runif(k.OBD-1, 0, q.OBD))
        qs.u <- sort(runif(ndose-k.OBD, 0, q.OBD), decreasing=TRUE)
        qs <- c(qs.l, q.OBD, qs.u)
    }
    
    res <- list(
        qs=qs,
        ps=ps$p.trues, 
        k.MTD=k.MTD, 
        k.OBD=k.OBD
    )
    res
}


res <- gen.rand.doses.umbrella(ndose, phi, psi, psi.U, mu1, mu2)
phaseI_II.scenario.plot.fn(res$ps, res$qs, OBD=res$k.OBD)
