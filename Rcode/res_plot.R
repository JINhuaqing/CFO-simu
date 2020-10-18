rm(list=ls())
setwd("C:/Users/Dell/Documents/ProjectCode/phaseI")
#setwd("/Users/jinhuaqing/Documents/Projects_Code/phaseI/")
source("Rcode/utilities.R")


load("results/Simu10000Level_5_10_CS3.RData")
res <- post.process.random(results)
res
res <- res[-1, ]
M.names <- c("ODDs", "BMS-N1", "BMS-N10", "BMS-N50", "BMS-NInf", "CRM", "BOIN")
res.ggplot.fn(res, filename="results/Level5_10_ms_corhortsize1.jpg", main="Cohortsize=1, ncohort=30, 5 Levels, 0.10", 
              M.names=M.names, angle=45)

