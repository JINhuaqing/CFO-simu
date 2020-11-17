phaseI_II.scenario.plot.fn <- function(p.true, pE.true, main="Plateau"){
    ndose <- length(p.true)
    cols <- c("black", "red")
    lwds <- c(2, 2)
    ltys <- c(1, 2)
    pchs <- c(15, 16)
    plot(1:ndose, p.true, ylab = "Probability", xlab="Dose level", main=main, 
         type="b", 
         lwd=lwds[1], col=cols[1], lty=ltys[1], pch=pchs[1])
    lines(1:ndose, pE.true, type="b", lwd=lwds[2], col=cols[2], lty=ltys[2], pch=pchs[2])
#    legend("topleft", c("DLT rate", "Efficacy rate"), col=cols, lty=ltys, lwd=lwds, pch=pchs)
}


p.true <- c(0.1, 0.2, 0.3, 0.4, 0.5)
pE.true <- c(0.1, 0.2, 0.3, 0.2, 0.1)

dev.off()
de.par <- par()
mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 7), nrow=4, byrow=TRUE)
layout(mat=mat, heights=c(0.5, 0.5, 0.5, 0.2))
par(mar=de.par$mar)

phaseI_II.scenario.plot.fn(p.true, pE.true)
phaseI_II.scenario.plot.fn(p.true, pE.true)
phaseI_II.scenario.plot.fn(p.true, pE.true)
phaseI_II.scenario.plot.fn(p.true, pE.true)
phaseI_II.scenario.plot.fn(p.true, pE.true)
phaseI_II.scenario.plot.fn(p.true, pE.true)

par(mar = c(0,2,0,2))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
cols <- c("black", "red")
lwds <- c(2, 2)
ltys <- c(1, 2)
pchs <- c(15, 16)
legend(x="top", inset=0, legend=c("DLT rate", "Efficacy rate"), col=cols, lty=ltys, lwd=lwds, pch=pchs, horiz=TRUE)
