library(TruncatedDistributions)
library(magrittr)

alp <- 0.1 + 1
bet <- 0.1 + 1
alp1 <- 0.1 + 0
bet1 <- 0.1 + 1
a <- 0.2 
b <- 0.9

sps.beta <- rbeta(10000, alp, bet)
plot(density(sps.beta), main="Density plot", ylim=c(0, 1))
sps.trbeta <- rtbeta(10000, alp, bet, sps.beta, b)



tbeta.sampler.low <- function(pcs, alp, bet, pU=0.9){
       n <- length(pcs)
       rtbeta(n, alp, bet, pcs, pU)
}

tbeta.sampler.up <- function(pcs, alp, bet, pL=0){
  n <- length(pcs)
  rtbeta(n, alp, bet, pL, pcs)
}

sps.trbeta2 <- c()
for (i in 1:length(sps.beta)){
  sp <- rtbeta(1, 0.1, 1.1, 0, sps.beta[i])
  sps.trbeta2 <- c(sps.trbeta2, sp)
}

# 10000 vs 1000000
sps.beta <- rbeta(10000, alp, bet)
sps.trbeta <- tbeta.sampler.up(sps.beta, 0.1, 1.1)
plot(density(sps.trbeta), main="Density plot", ylim=c(0, 1))
lines(density(sps.trbeta2), main="Density plot", col=2)
