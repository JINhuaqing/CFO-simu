source("utilities.R")
load("results4.RData")

sum.fn <- function(ress){
    nams <- names(ress)
    ndose <- length(ress[[1]]$Allocation)
    res.df <- data.frame(levels=1:ndose)
    for (nam in nams){
        res.df[paste0(nam, ".Sel")] <- ress[[nam]]$Selection
    }
    for (nam in nams){
        res.df[paste0(nam, ".Allo")] <- ress[[nam]]$Allocation
    }
    res.df
}


sum.fn(results)
