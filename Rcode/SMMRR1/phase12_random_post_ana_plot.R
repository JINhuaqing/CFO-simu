rm(list=ls())
setwd("C:/Users/JINHU/Documents/ProjectCode/CFO/Rcode")
source("utilities.R")


phi <- 0.30
psi1 = phiE = 0.30
psi2 <- 0.65
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsim <- 5000
ncore <- 10
ndose <- 5

# Target = 0.3
# dose 5, mu1=mu2=0.23, 0.05
# dose 5, mu1=mu2=0.38, 0.07 
# dose 5, mu1=mu2=0.53, 0.1
# dose 5, mu1=mu2=0.71, 0.15

mu1 <- 0.23
typ <- "um"
res.all <- list()
for (mu1 in c(0.23, 0.38, 0.53, 0.71)){
for (typ in c("um", "pl")){
mu2 <- mu1
res.name <- paste0(typ, 100*mu1)
# random scenario results
prefix <- paste0("../results/SMMR-R1/phase12_", typ)
save.f.Name.STEIN <- paste0(prefix, "_STEIN2_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.CFO <- paste0(prefix, "_CFO_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.CFO2 <- paste0(prefix, "_CFO2_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.WT <- paste0(prefix, "_WT_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
save.f.Name.MADA <- paste0(prefix, "_MADA_nsim", nsim, "_ndose", ndose, "_phi", 
         phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
load(save.f.Name.STEIN)
load(save.f.Name.WT)
load(save.f.Name.CFO)
load(save.f.Name.MADA)
load(save.f.Name.CFO2)
tb <- phase12.post.process.random.all(ress.WT, ress.STEIN, ress.MADA, 
                                ress.CFO, ress.CFO2, names=c("WT", "STEIN", "MADA", "CFO", "CFO2"))
tb <- tb[c(5, 3, 2, 1), ]
rownames(tb) <- c("CFO", "MADA", "STEIN", "WT")
        
res.all[[res.name]] <- tb
}
}

# OBD Sel
typ <- "pl"
png(filename=paste0("../plots/SMMRR1/", typ, "Sel.png"), unit="in", height=6, width=6, res=300)
flag <- 1
for (mu1 in c(0.23, 0.38, 0.53, 0.71)){
    res.name <- paste0(typ, 100*mu1)
    res.cur <- res.all[[res.name]]
    if (mu1 == 0.23){
        if (typ == "pl")
        plot(res.cur$OBD.Sel, type="b", ylim=c(0.28, 0.40), lty = 1, col=flag, pch=flag, lwd=2,
             xlab = "", xaxt="n", ylab="OBD Selection (%)")
        else 
        plot(res.cur$OBD.Sel, type="b", ylim=c(0.32, 0.48), lty = 1, col=flag, pch=flag, lwd=2,
             xlab = "", xaxt="n", ylab="OBD Selection (%)")
        axis(side=1, at=1:4, labels=rownames(res.cur))
    }else{
        lines(res.cur$OBD.Sel, type="b", lty=1, col=flag, pch=flag, lwd=2)
    }
    abline(h=mean(res.cur$OBD.Sel[1]), col=flag, lty=2)
    flag <- flag + 1
}

legend("bottomright", legend = c("0.05", "0.07", "0.10", "0.15"), lty=1, col=1:4, pch=1:4, lwd=2)
dev.off()

# OBD Allo
typ <- "um"
png(filename=paste0("../plots/SMMRR1/", typ, "Allo.png"), unit="in", height=6, width=6, res=300)
flag <- 1
for (mu1 in c(0.23, 0.38, 0.53, 0.71)){
    res.name <- paste0(typ, 100*mu1)
    res.cur <- res.all[[res.name]]
    if (mu1 == 0.23){
        if (typ == "pl")
        plot(res.cur$OBD.Allo, type="b", ylim=c(0.24, 0.40), lty = 1, col=flag, pch=flag, lwd=2,
             xlab = "", xaxt="n", ylab="OBD Allocation (%)")
        else 
        plot(res.cur$OBD.Allo, type="b", ylim=c(0.25, 0.42), lty = 1, col=flag, pch=flag, lwd=2,
             xlab = "", xaxt="n", ylab="OBD Allocation (%)")
        axis(side=1, at=1:4, labels=rownames(res.cur))
    }else{
        lines(res.cur$OBD.Allo, type="b", lty=1, col=flag, pch=flag, lwd=2)
    }
    abline(h=mean(res.cur$OBD.Allo[1]), col=flag, lty=2)
    flag <- flag + 1
}
legend("bottomright", legend = c("0.05", "0.07", "0.10", "0.15"), lty=1, col=1:4, pch=1:4, lwd=2)
dev.off()



res.both.all <- list()
for (mu1 in c(0.23, 0.38, 0.53, 0.71)){
mu2 <- mu1
res.name <- paste0("Delta", mu1*100)
# random scenario results
prefix <- paste0("../results/SMMR-R1/phase12_", "um")
save.f.Name.STEIN <- paste0(prefix, "_STEIN2_nsim", nsim, "_ndose", ndose, "_phi",  phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
load(save.f.Name.STEIN)
ress.STEIN1 <- ress.STEIN
prefix <- paste0("../results/SMMR-R1/phase12_", "pl")
save.f.Name.STEIN <- paste0(prefix, "_STEIN2_nsim", nsim, "_ndose", ndose, "_phi",  phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
load(save.f.Name.STEIN)
ress.STEIN.b <- c(ress.STEIN1, ress.STEIN)

prefix <- paste0("../results/SMMR-R1/phase12_", "um")
save.f.Name.CFO <- paste0(prefix, "_CFO_nsim", nsim, "_ndose", ndose, "_phi",  phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
load(save.f.Name.CFO)
ress.CFO1 <- ress.CFO
prefix <- paste0("../results/SMMR-R1/phase12_", "pl")
save.f.Name.CFO <- paste0(prefix, "_CFO_nsim", nsim, "_ndose", ndose, "_phi",  phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
load(save.f.Name.CFO)
ress.CFO.b <- c(ress.CFO1, ress.CFO)

prefix <- paste0("../results/SMMR-R1/phase12_", "um")
save.f.Name.CFO2 <- paste0(prefix, "_CFO2_nsim", nsim, "_ndose", ndose, "_phi",  phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
load(save.f.Name.CFO2)
ress.CFO21 <- ress.CFO2
prefix <- paste0("../results/SMMR-R1/phase12_", "pl")
save.f.Name.CFO2 <- paste0(prefix, "_CFO2_nsim", nsim, "_ndose", ndose, "_phi",  phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
load(save.f.Name.CFO2)
ress.CFO2.b <- c(ress.CFO21, ress.CFO2)

prefix <- paste0("../results/SMMR-R1/phase12_", "um")
save.f.Name.MADA <- paste0(prefix, "_MADA_nsim", nsim, "_ndose", ndose, "_phi",  phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
load(save.f.Name.MADA)
ress.MADA1 <- ress.MADA
prefix <- paste0("../results/SMMR-R1/phase12_", "pl")
save.f.Name.MADA <- paste0(prefix, "_MADA_nsim", nsim, "_ndose", ndose, "_phi",  phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
load(save.f.Name.MADA)
ress.MADA.b <- c(ress.MADA1, ress.MADA)

prefix <- paste0("../results/SMMR-R1/phase12_", "um")
save.f.Name.WT <- paste0(prefix, "_WT_nsim", nsim, "_ndose", ndose, "_phi",  phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
load(save.f.Name.WT)
ress.WT1 <- ress.WT
prefix <- paste0("../results/SMMR-R1/phase12_", "pl")
save.f.Name.WT <- paste0(prefix, "_WT_nsim", nsim, "_ndose", ndose, "_phi",  phi*100, "_phiE", phiE*100, "_mu", mu1*100, "_ncoh", ncohort, "_cohSize", cohortsize, ".RData")
load(save.f.Name.WT)
ress.WT.b <- c(ress.WT1, ress.WT)

tb <- phase12.post.process.random.all(ress.WT.b, ress.STEIN.b, ress.MADA.b,
                                ress.CFO.b, ress.CFO2.b, names=c("WT", "STEIN", "MADA", "CFO", "CFO2"))
        
tb <- tb[c(5, 3, 2, 1), ]
rownames(tb) <- c("CFO", "MADA", "STEIN", "WT")
res.both.all[[res.name]] <- tb
}


# OBD Sel
png(filename="../plots/SMMRR1/bothSel.png", unit="in", height=6, width=6, res=300)
flag <- 1
for (mu1 in c(0.23, 0.38, 0.53, 0.71)){
    res.name <- paste0("Delta", 100*mu1)
    res.cur <- res.both.all[[res.name]]
    if (mu1 == 0.23){
        plot(res.cur$OBD.Sel, type="b", ylim=c(0.28, 0.45), lty = 1, col=flag, pch=flag, lwd=2,
             xlab = "", xaxt="n", ylab="OBD Selection (%)")
        axis(side=1, at=1:4, labels=rownames(res.cur))
    }else{
        lines(res.cur$OBD.Sel, type="b", lty=1, col=flag, pch=flag, lwd=2)
    }
    abline(h=mean(res.cur$OBD.Sel[1]), col=flag, lty=2)
    flag <- flag + 1
}

legend("bottomright", legend = c("0.05", "0.07", "0.10", "0.15"), lty=1, col=1:4, pch=1:4, lwd=2)
dev.off()

# OBD Allo
png(filename="../plots/SMMRR1/bothAllo.png", unit="in", height=6, width=6, res=300)
flag <- 1
for (mu1 in c(0.23, 0.38, 0.53, 0.71)){
    res.name <- paste0("Delta", 100*mu1)
    res.cur <- res.both.all[[res.name]]
    if (mu1 == 0.23){
        plot(res.cur$OBD.Allo, type="b", ylim=c(0.25, 0.42), lty = 1, col=flag, pch=flag, lwd=2,
             xlab = "", xaxt="n", ylab="OBD Allocation(%)")
        axis(side=1, at=1:4, labels=rownames(res.cur))
    }else{
        lines(res.cur$OBD.Allo, type="b", lty=1, col=flag, pch=flag, lwd=2)
    }
    abline(h=mean(res.cur$OBD.Allo[1]), col=flag, lty=2)
    flag <- flag + 1
}
legend("bottomright", legend = c("0.05", "0.07", "0.10", "0.15"), lty=1, col=1:4, pch=1:4, lwd=2)
dev.off()



