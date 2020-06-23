x1s <- rbeta(100000, 1, 3)
x2s <- rbeta(100000, 2, 2)
x3s <- rbeta(100000, 3, 1)
xs <- rbind(x1s, x2s, x3s)
xs.sorted <- apply(xs, 2, sort)
x1s <- xs.sorted[1, ]
x2s <- xs.sorted[2, ]
x3s <- xs.sorted[3, ]

y1s <- c()
y2s <- c()
y3s <- c()
while (length(y1s)<=100000){
    y1 <- rbeta(1, 1, 3)
    y2 <- rbeta(1, 2, 2)
    y3 <- rbeta(1, 3, 1)
    if (y1 < y2 & y2 < y3){
        y1s <- c(y1s, y1)
        y2s <- c(y2s, y2)
        y3s <- c(y3s, y3)
    }
}

plot(density(x1s), main="X1", col=2, lty=2, lwd=2, ylim=c(0, 4.5), xlim=c(0, 1))
lines(density(y1s), col=3, lty=3, lwd=2)
lines(density(x2s), col=2, lty=2, lwd=2)
lines(density(y2s), col=3, lty=3, lwd=2)
lines(density(x3s), col=2, lty=2, lwd=2, ylim=c(0, 4))
lines(density(y3s), col=3, lty=3, lwd=2)
legend("topleft", legend=c("Method 1", "Method 2"), col=2:3, lty=2:3, lwd=rep(2, 2))
