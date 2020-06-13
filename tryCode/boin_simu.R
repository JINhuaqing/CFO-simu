library(BOIN)

target <- 0.2

p.true1 <- c(0.1, 0.2, 0.3)
p.true2 <- c(0.05, 0.22, 0.38)

p.true3 <- c(0.2, 0.3, 0.4)
p.true4 <- c(0.18, 0.3, 0.45)

p.true5 <- c(0.07, 0.13, 0.21)
p.true6 <- c(0.04, 0.1, 0.2)

ncohort <- 12
cohortsize <- 1

simu <- get.oc(target=target, 
       p.true=p.true6,
       ncohort=ncohort,
       cohortsize=cohortsize, startdose=2,
       ntrial=1000)

summary(simu)

bdry <- get.boundary(target=0.4, ncohort=4, cohortsize=3)
summary(bdry)
