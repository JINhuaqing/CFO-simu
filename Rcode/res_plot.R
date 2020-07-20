rm(list=ls())
#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI")
setwd("/Users/jinhuaqing/Documents/Projects_Code/phaseI/")
source("Rcode/utilities.R")


load("results/Simu10000Level_7_10_ms_corhortsize1.RData")
res <- post.process.random(results)
res
M.names <- c("BF", "BMS-m1", "BMS-m10", "BMS-m50", "BMS-mInf", "CRM", "BOIN")
res.ggplot.fn(res, filename="results/Level7_10_BMS_corhortsize1.jpg", main="Cohortsize=1, ncohort=30, 7 Levels, 0.10", M.names=M.names, angle=45)

