setwd("../../Rcode/")
source("utilities.R")
idx <- 5
nsimu <- 5000
file.name <- paste0("../results/JRSSC-R/", "Hard_MTDSimu_", nsimu, "fix_",  idx, ".RData")
load(file.name)
crm.ress <- lapply(1:nsimu, function(i)results[[i]]$crm)
boin.ress <- lapply(1:nsimu, function(i)results[[i]]$boin)
orm.ress <- lapply(1:nsimu, function(i)results[[i]]$orm)
sum.all <- list(
                BOIN = phase1.post.fn(boin.ress),
                ORM = phase1.post.fn(orm.ress),
                CRM = phase1.post.fn(crm.ress)
                )
print(phase.I.pretty.tb(sum.all)[, -c(6, 7)])
print(phase.I.pretty.tb(sum.all))
