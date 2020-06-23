rm(list=ls())
setwd("C:/Users/Dell/Documents/ProjectCode/phaseI")
source("tryCode/utilities.R")


load("results/Odds_rand10000.RData")
res <- post.process.random(results)
M.names <- c("Butterfly-Odds", "CRM", "BOIN")
#M.names <- c("Butterfly-BB", "Butterfly-CRM", "Butterfly-BB-CRM", "CRM", "BOIN")
res.ggplot.fn(res, filename="Odds_1_30.jpg", main="Odds, cohortsize=1, ncohort=30", M.names=M.names, angle=45)

