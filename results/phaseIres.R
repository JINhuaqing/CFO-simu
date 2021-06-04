setwd("../results/")
fil <- dir(pattern="5000random_0.15")
fil
load(fil)

setwd("../Rcode")
source("utilities.R")
# post.process.random(results) # for random cases


nsimu <- 5000
crm.ress <- lapply(1:nsimu, function(i)results[[i]]$crm)
boin.ress <- lapply(1:nsimu, function(i)results[[i]]$boin)
orm.ress <- lapply(1:nsimu, function(i)results[[i]]$orm)
sum.all <- list(
                BOIN = phase1.post.fn(boin.ress),
                ORM = phase1.post.fn(orm.ress),
                CRM = phase1.post.fn(crm.ress)
                )
tb <- phase.I.pretty.tb(sum.all); tb
tb <- tb[, -7];tb

tb[c(2, 1, 3), ]

