source("../Rcode/utilities.R")
# plateau: 2, 4
# umbrella: 7, 12, 11
# increasing: 1, 4
# on OBD: 2, 3, 5
idx <- 5
fName <- "plateau"
fName <- "umbrella"
fName <- "increa"
fName <- "noOBD"

fName.all <- paste0("res_", fName, "_1000_", idx, ".RData")
#fName.efftox <- paste0("efftox_res_", fName, "_1000_", idx, ".RData")
#load(fName.efftox)
load(fName.all)
#sum.all[["efftox"]] <- sum.res.efftox
#ls()
phase.I.II.pretty.tb(sum.all)
