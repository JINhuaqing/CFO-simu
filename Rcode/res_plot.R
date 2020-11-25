rm(list=ls())
library(ggplot2)
setwd("C:/Users/Dell/Documents/ProjectCode/phaseI")
#setwd("/Users/jinhuaqing/Documents/Projects_Code/phaseI/")
source("Rcode/utilities.R")

fils <- dir("results", pattern="MTD.+random.+15");fils

fil <- paste0("results/", fils[1])
load(fil)

grp.names <- c("MTD selection", "MTD allocation", "Overdose selection", "Overdose allocation" 
               , "Risk of high toxicity",  "Average DLT rate")
m.names <- c("BOIN", "CFOD", "CRM")

g.var <- rep(grp.names, times=3)
m.var <- rep(m.names, each=6)

tb <- post.process.random(results);tb
#tb <- tb[, 1:4];tb
v.var <- c(as.vector(unlist(tb["boin", ])),as.vector(unlist(tb["orm", ])), as.vector(unlist(tb["crm", ]))) * 100

data <- data.frame(g=factor(g.var, levels=grp.names), m=factor(m.var, levels=m.names), v=v.var)
ggplot(data = data, mapping = aes(x = g, y = v, fill = m)) + geom_bar(stat = 'identity', position = 'dodge') +
    theme(legend.position = "bottom", plot.title = element_text(hjust=0.5)) +  xlab("") + ylab("Percentage (%)") + 
    guides(fill=guide_legend(title='Methods')) + ggtitle("Average probability difference around the target = 0.15")
    
ggsave("plots/MTD_random_15.jpg", width=10, height = 4.5, units="in")
