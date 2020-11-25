setwd("C:/Users/Dell/Documents/ProjectCode/phaseI")
phi <- 0.3

xs <- seq(0.01, 0.99, length.out=100)


cols <- rainbow(3)
ltys <- c(2, 1, 4)

pdf("./plots/odds_L.pdf", width=8, height=4)
par(mfrow=c(1,2))
par(mar = c(2, 4, 2, 2))
alps <- c(3.5, 5, 7)
bets <- c(7, 5, 3)
plot(xs, dbeta(xs, alps[1], bets[1]), type="l", lwd=2, lty=ltys[1], 
     xlab="", ylim=c(0, 4), xaxt="n", ylab="Posterior density")
lines(xs, dbeta(xs, alps[2], bets[2]), type="l", col=cols[1], lty=ltys[2], lwd=2)
lines(xs, dbeta(xs, alps[3], bets[3]), type="l", col=cols[2], lty=ltys[3], lwd=2)
axis(1, phi, expression(phi))
abline(v=phi, lty=3)
legend.labs <- c(expression(p[L]), expression(p[C]), expression(p[R])) 
legend("topright", legend.labs, lty=ltys, col=c("black", cols[1:2]))

alps <- c(2, 3.5, 7)
bets <- c(8, 7, 3)
plot(xs, dbeta(xs, alps[1], bets[1]), type="l", lwd=2, lty=ltys[1], 
     xlab="", ylim=c(0, 4), xaxt="n", ylab="Posterior density")
lines(xs, dbeta(xs, alps[2], bets[2]), type="l", col=cols[1], lty=ltys[2], lwd=2)
lines(xs, dbeta(xs, alps[3], bets[3]), type="l", col=cols[2], lty=ltys[3], lwd=2)
axis(1, phi, expression(phi))
abline(v=phi, lty=3)
legend.labs <- c(expression(p[L]), expression(p[C]), expression(p[R])) 
legend("topright", legend.labs, lty=ltys, col=c("black", cols[1:2]))

dev.off()

