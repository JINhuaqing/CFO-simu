rm(list=ls())
setwd("C:/Users/Dell/Documents/ProjectCode/phaseI")
source("Rcode/utilities.R")


load("results/Simu10000Level_7_15.RData")
res <- post.process.random(results)
res
M.names <- c("Butterfly-BF", "Butterfly-Odds", "CRM", "BOIN")
res.ggplot.fn(res, filename="results/Level7_15.jpg", main="Cohortsize=3, ncohort=10, 7 Levels, 0.15", M.names=M.names, angle=45)


