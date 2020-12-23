rm(list=ls())

setwd("C:/Users/JINHU/Documents/ProjectCode/CFO/Rcode")
source("./utilities.R")
source("./ORM_utils.R")


cover.doses <- c(0, 0, 0)
phi <- 0.3
alp.prior <- phi
bet.prior <- 1-phi

# Short-term memory coherence counter example
cys <- c(0, 0, 0)
cns <- c(4, 4, 3)
make.decision.ORM.fn(phi, cys, cns, alp.prior, bet.prior, cover.doses, diag=FALSE)
# move farword

cys2 <- c(0, 1, 0)
cns2 <- c(4, 4, 3)
make.decision.ORM.fn(phi, cys2, cns2, alp.prior, bet.prior, cover.doses, diag=FALSE)



# Long-term memory coherence counter example
cys <- c(0, 4, 1)
cns <- c(6, 10, 10)
make.decision.ORM.fn(phi, cys, cns, alp.prior, bet.prior, cover.doses, diag=FALSE)
