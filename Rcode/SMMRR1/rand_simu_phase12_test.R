setwd("/home/r6user2/MyResearch/CFO/Rcode/")
source("utilities.R")
library(parallel)

phi <- 0.30
psi1 = phiE = 0.30
psi2 <- 0.65
cohortsize = 3 # cohort size
ncohort = 20 # number of cohorts
nsim <- 5000
ncore <- 40
ndose <- 5
seeds <- 1:nsim
add.args <- list(alp.prior=phi, bet.prior=1-phi, alp.prior.eff=0.5, bet.prior.eff=1-0.5)

# Target = 0.3
# dose 5, mu1=mu2=0.23, 0.05
# dose 5, mu1=mu2=0.38, 0.07 
# dose 5, mu1=mu2=0.53, 0.1
# dose 5, mu1=mu2=0.71, 0.15
mu1 <- 0.23
typ <- "um"

util.fn <- function(sc, phi=0.3){
    u<-sc$qs-0.33*sc$ps-1.09*(sc$ps>phi) 	
    uMax <- which.max(u)
    return(uMax==sc$k.OBD)
}

for (mu1 in c(0.23, 0.38, 0.53, 0.71)){
for (typ in c("pl", "um")){
mu2 <- mu1



f.Name.um <- paste0("../results/SMMR-R1/phase12_rc_umbrellaU_nsim", nsim, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")
f.Name.pl <- paste0("../results/SMMR-R1/phase12_rc_plateauU_nsim", nsim, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")
load(f.Name.um)
load(f.Name.pl)



if (typ=="um"){
    scs <- scs.um
}else {
    scs <- scs.pl
}

rates <- sapply(scs, util.fn)
print(mean(rates))


}
}
