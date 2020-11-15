source("../Rcode/utilities.R")
# plateau: 2, 4
# umbrella: 7
idx <- 9
fName <- "plateau"
fName <- "umbrella"

fName.all <- paste0("res_", fName, "_5000_", idx, ".RData")
fName.efftox <- paste0("efftox_res_", fName, "_", idx, ".RData")
load(fName.efftox)
load(fName.all)
sum.all[["efftox"]] <- sum.res.efftox
ls()
phase.I.II.pretty.tb(sum.all)
