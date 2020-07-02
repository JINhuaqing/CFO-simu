rm(list=ls())
setwd("C:/Users/Dell/Documents/ProjectCode/phaseI")
source("Rcode/utilities.R")


load("results/Allmethods10000TSS_30_Level_5.RData")
res <- post.process.random(results)
res <- res[c(1, 2, 6, 7), ]
res
M.names <- c("Butterfly-BF", "Butterfly-Odds", "CRM", "BOIN")
res.ggplot.fn(res, filename="results/Level5.jpg", main="Cohortsize=3, ncohort=10", M.names=M.names, angle=45)

