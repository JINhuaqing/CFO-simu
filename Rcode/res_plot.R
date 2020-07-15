rm(list=ls())
setwd("C:/Users/Dell/Documents/ProjectCode/phaseI")
source("Rcode/utilities.R")


load("results/Simu10000Level_3_10_ms.RData")
res <- post.process.random(results)
res
M.names <- c("BF", "BMS-m1", "BMS-m10", "BMS-m50", "BMS-mInf", "CRM", "BOIN")
res.ggplot.fn(res, filename="results/Level3_10_BMS.jpg", main="Cohortsize=3, ncohort=10, 3 Levels, 0.10", M.names=M.names, angle=45)


