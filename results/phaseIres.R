setwd("../../results/")
fil <- dir(pattern="5000fix_4")
fil
load(fil)

setwd("../Rcode")
source("utilities.R")


nsimu <- 5000
crm.ress <- lapply(1:nsimu, function(i)results[[i]]$crm)
boin.ress <- lapply(1:nsimu, function(i)results[[i]]$boin)
orm.ress <- lapply(1:nsimu, function(i)results[[i]]$orm)
sum.all <- list(
                BOIN = phase1.post.fn(boin.ress),
                ORM = phase1.post.fn(orm.ress),
                CRM = phase1.post.fn(crm.ress)
                )
phase.I.pretty.tb(sum.all)
