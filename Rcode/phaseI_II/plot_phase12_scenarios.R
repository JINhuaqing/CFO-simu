setwd("C:/Users/JINHU/Documents/ProjectCode/CFO/Rcode")
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


scenarios <- list()
scenarios[[1]] <- list(p.true=c(0.05, 0.1, 0.3, 0.5, 0.6),
                      pE.true=c(0.2, 0.3, 0.5, 0.5, 0.5))
scenarios[[2]] <- list(p.true=c(0.15, 0.25, 0.3, 0.35, 0.4),
                      pE.true=c(0.2, 0.5, 0.5, 0.5, 0.5))

scenarios[[3]] <- list(p.true=c(0.10, 0.22, 0.25, 0.3, 0.4),
                       pE.true=c(0.3, 0.6, 0.55, 0.35, 0.2))
scenarios[[4]] <- list(p.true=c(0.05, 0.15, 0.25, 0.40, 0.45),
                        pE.true=c(0.08, 0.17, 0.45, 0.30, 0.25))

scenarios[[5]] <- list(p.true=c(0.05, 0.07, 0.1, 0.12, 0.16),
                     pE.true=c(0.35, 0.45, 0.5, 0.55, 0.75))

scenarios[[6]] <- list(p.true=c(0.40, 0.5, 0.55, 0.6, 0.7),
                       pE.true=c(0.15, 0.25, 0.5, 0.5, 0.5))
# for (i in 1:6){
#     sce <- scenarios[[i]]
#     outstr <- paste0("$(", sce$p.true, ",", " ", sce$pE.true, ")$")
#     outstr <- paste(outstr, collapse = " & ")
#     outstr <- paste0("&", outstr, " \\\\", "\n")
#     cat(outstr)
# }


pdf("../plots/SMMRR1/eff_tox_sc6new.pdf", width=6, height=8)
de.par <- par()
mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 7), nrow=4, byrow=TRUE)
layout(mat=mat, heights=c(0.5, 0.5, 0.5, 0.2))
par(mar=de.par$mar)

phaseI_II.scenario.plot.fn(scenarios[[1]]$p.true, scenarios[[1]]$pE.true, main="scenario 1", OBD=3)
phaseI_II.scenario.plot.fn(scenarios[[2]]$p.true, scenarios[[2]]$pE.true, main="scenario 2", OBD=2)
phaseI_II.scenario.plot.fn(scenarios[[3]]$p.true, scenarios[[3]]$pE.true, main="scenario 3", OBD=2)
phaseI_II.scenario.plot.fn(scenarios[[4]]$p.true, scenarios[[4]]$pE.true, main="scenario 4", OBD=3)
phaseI_II.scenario.plot.fn(scenarios[[5]]$p.true, scenarios[[5]]$pE.true, main="scenario 5", OBD=5)
phaseI_II.scenario.plot.fn(scenarios[[6]]$p.true, scenarios[[6]]$pE.true, main="scenario 6", OBD=NULL)

par(mar = c(0, 2, 0, 2))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
cols <- c("black", "red")
lwds <- c(2, 2)
ltys <- c(1, 2)
pchs <- c(15, 16)
legend(x="top", inset=0, legend=c("DLT rate", "Efficacy rate"), 
       col=cols, lty=ltys, lwd=lwds, pch=pchs, horiz=TRUE)
dev.off()
