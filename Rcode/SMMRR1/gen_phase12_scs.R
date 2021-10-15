# To generate the phase I/II scenarios
setwd("/home/r6user2/MyResearch/CFO/Rcode/")
source("utilities.R")
set.seed(1)


# Target = 0.3
# dose 5, mu1=mu2=0.23, 0.05
# dose 5, mu1=mu2=0.38, 0.07 
# dose 5, mu1=mu2=0.53, 0.1
# dose 5, mu1=mu2=0.71, 0.15

phi <- 0.30 # target DLT rate
phiE <- 0.30 # the minimal acceptable efficate rate
phiE.U <- 0.80 # The upper bound of the efficacy rate
ndose <- 5
mu1 = mu2 = 0.71  # Delta is 0.1 for ndose = 5
nsim <- 5000

# with U, we make the sc consistent with STEIN
f.Name.um <- paste0("../results/SMMR-R1/phase12_rc_umbrellaU_nsim", nsim, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")
f.Name.pl <- paste0("../results/SMMR-R1/phase12_rc_plateauU_nsim", nsim, "_ndose", ndose, "_phi", phi*100, "_phiE", phiE*100, "_mu", mu1*100, ".RData")

scs.um <- lapply(1:nsim, function(i)gen.rand.doses.umbrella.U(ndose, phi, phiE, phiE.U, mu1, mu2))
save(scs.um, file=f.Name.um)
scs.pl <- lapply(1:nsim, function(i)gen.rand.doses.plateau.U(ndose, phi, phiE, phiE.U, mu1, mu2))
save(scs.pl, file=f.Name.pl)
