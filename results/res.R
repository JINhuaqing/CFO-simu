#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/results")
source("../Rcode/utilities.R")
# plateau: 2, 4 / 2, 4
# umbrella: 7, 12, 11 / 7, 12
# increasing: 1, 4 / 1
# on OBD: 2 / 2
idx <- 2
fName <- "plateau"
fName <- "umbrella"
fName <- "increa"
fName <- "noOBD"

fName.all <- paste0("res_", fName, "_5000_", idx, ".RData")
fName.efftox <- paste0("efftox_res_", fName, "_5000_", idx, ".RData")
#fName.all <- paste0("res_", fName, "_5000_", idx, ".RData")
#fName.efftox <- paste0("efftox_res_", fName, "_", idx, ".RData")

load(fName.efftox)
load(fName.all)
sum.all[["efftox"]] <- sum.res.efftox
res <- phase.I.II.pretty.tb(sum.all);res
fsadf

nams <- c("STEIN", "CF", "MADA", "EffTox")
#nams <- c("STEIN", "CF", "CF.alter", "MADA", "EffTox")
res.trans <- data.frame(lapply(res, as.character), stringsAsFactors = FALSE)
row.names(res.trans) <- nams
res.trans <- res.trans[, -8]


trans.f <- function(eachRow){
    level5 <- sapply(eachRow[1:5], function(i)paste(strsplit(i, "(", fixed=TRUE)[[1]], collapse = " ("))
    itm1 <- paste(level5, "& ")
    itm1 <- paste(itm1, collapse = "")
    itm2 <- paste0(eachRow[6], "/", eachRow[7])
    res.str <- paste(itm1, itm2, "&", eachRow[8])
    res.str <- paste(row.names(eachRow), "&", res.str, "\\\\")
    res.str
}

for (i in c(2, 4, 3, 1)){
#for (i in c(2, 5, 4, 1)){
    cat(trans.f(res.trans[i, ]), "\n")
}



