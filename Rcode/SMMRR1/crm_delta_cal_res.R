setwd("/home/r6user2/MyResearch/CFO/")

phi <- 0.33
psi <- 2
ndose <- 5
nu <- 3
ncohort <- 10
deltas <- seq(0.01, min(phi-0.005, 0.15), 0.01)
f.Name <- paste0("./results/SMMR-R1/CRM_cal_phi", phi*100, "_ndose_", ndose, "_ncohort_", ncohort, ".RData")
load(f.Name)
print(res)
print(deltas[which.max(res)])


fasdfas
for (ndose in c(3, 5, 7, 9)){
    for (ncohort in c(7, 10, 16, 20)){
        for (phi in c(0.25, 0.30, 0.33)){
            nu <- ceiling(ndose/2)
            deltas <- seq(0.01, min(phi-0.005, 0.15), 0.01)
            f.Name <- paste0("./results/SMMR-R1/CRM_cal_phi", phi*100, "_ndose_", ndose, "_ncohort_", ncohort, ".RData")
            load(f.Name)
            print(deltas[which.max(res)])
        }
    }
}
