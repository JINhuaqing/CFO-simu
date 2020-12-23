rm(list=ls())

library(reshape2)
library(ggplot2)
setwd("C:/Users/JINHU/Documents/ProjectCode/CFO/Rcode")
setwd("C:/Users/Dell/Documents/ProjectCode/CFO/Rcode")
source("./utilities.R")
source("./ORM_utils.R")


# 1. test the value of gamma and OR
y1 <- 1
n1 <- 6
y2 <- 2
n2 <- 6
y3 <- 3
n3 <- 6
phi <- 0.3
alp.prior <- phi
bet.prior <- 1-phi
ORL <- OR.values(phi, y1, n1, y2, n2, alp.prior, bet.prior, type="L");ORL
optim.gamma.fn(n1, n2, phi, type="L", alp.prior, bet.prior)
ORR <- OR.values(phi, y2, n2, y3, n3, alp.prior, bet.prior, type="R");ORR
optim.gamma.fn(n2, n3, phi, type="R", alp.prior, bet.prior)
ORL <- OR.values(phi, y1, n1, y2+1, n2+1, alp.prior, bet.prior, type="L");ORL
optim.gamma.fn(n1, n2+1, phi, type="L", alp.prior, bet.prior)


# 2 draw the plots for gammaL
n1s <- 1:30
n2s <- 1:30
gams.raw.L <- matrix(NA, nrow=30, ncol=30)
for (n1 in n1s){
  for (n2 in n2s){
    print(c(n1, n2))
    gams.raw.L[n1, n2] <- optim.gamma.fn(n1, n2, phi, type="L", alp.prior, bet.prior)$gam
  }
}



gams <- as.data.frame(gams.raw.L)
rownames(gams) <- n1s
colnames(gams) <- n2s
gams$nL <- rownames(gams)

data <- melt(gams, id.vars=c("nL"))
data$nL <- as.numeric(data$nL)
data$nL <- as.factor(data$nL)

names(data) <- c("nL", "nC", "gammaL")

p <- ggplot(data=data, mapping=aes(x=nL, y=nC)) + 
  geom_tile(aes(fill=gammaL)) + 
  scale_fill_gradient2(expression(gamma[L]), low="#5bc5f4", mid="white", high="red",  
                       midpoint=mean(range(data$gammaL))-0.1, limits=range(data$gammaL)) +
  theme(legend.position = "right") + 
  ylab(expression(m[C])) +  xlab(expression(m[L])) +
  scale_x_discrete(labels=seq(5, 30, 5), breaks=seq(5, 30, 5)) + 
  scale_y_discrete(labels=seq(5, 30, 5), breaks=seq(5, 30, 5)); p
ggsave("../plots/gammaL.jpg", height=5, width = 6, unit="in")


# 3 draw the plots for gammaR
n1s <- 1:30
n2s <- 1:30
gams.raw.R <- matrix(NA, nrow=30, ncol=30)
for (n1 in n1s){
  for (n2 in n2s){
    print(c(n1, n2))
    gams.raw.R[n1, n2] <- optim.gamma.fn(n1, n2, phi, type="R", alp.prior, bet.prior)$gam
  }
}


gams <- as.data.frame(gams.raw.R)
rownames(gams) <- n1s
colnames(gams) <- n2s
gams$nC <- rownames(gams)

data <- melt(gams, id.vars=c("nC"))
data$nC <- as.numeric(data$nC)
data$nC <- as.factor(data$nC)

names(data) <- c("nC", "nR", "gammaR")

p <- ggplot(data=data, mapping=aes(x=nR, y=nC)) + 
  geom_tile(aes(fill=gammaR)) + 
  scale_fill_gradient2(expression(gamma[R]), low="#5bc5f4", mid="white", high="red",  
                       midpoint=mean(range(data$gammaR)), limits=range(data$gammaR)) +
  theme(legend.position = "right") + 
  ylab(expression(m[C])) +  xlab(expression(m[R])) +
  scale_x_discrete(labels=seq(5, 30, 5), breaks=seq(5, 30, 5)) + 
  scale_y_discrete(labels=seq(5, 30, 5), breaks=seq(5, 30, 5))
p
ggsave("../plots/gammaR.jpg", height=5, width = 6, unit="in")

