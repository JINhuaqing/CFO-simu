library(BOIN)

target <- 0.2

p.true <- c(0.01, 0.05, 0.1)
#p.true <- c(0.05, 0.1, 0.2)
#p.true <- c(0.1, 0.2, 0.3)
ncohort <- 4
cohortsize <- 3

simu <- get.oc(target=target, 
       p.true=p.true,
       ncohort=ncohort,
       cohortsize=cohortsize,
       ntrial=1000)

summary(simu)
