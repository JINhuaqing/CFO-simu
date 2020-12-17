library(dplyr)
setwd("C:/Users/Dell/Documents/ProjectCode/CFO")
source("./Rcode/utilities.R")
anova.ress <- dir("./results/Anovas", pattern="*.RData", full.names = TRUE)

sep.paths <- strsplit(anova.ress, "_")
raw.labs <- sapply(sep.paths, function(i)i[3])
labs <- lapply(strsplit(raw.labs, ""), function(i)as.numeric(i))
mat.labs <- do.call(rbind, labs)

anova.data.list <- list()
for (i in 1:length(anova.ress)){
    anova.fil <- anova.ress[i]
    s.labs <- mat.labs[i, ]
    m.labs <- matrix(rep(s.labs, 3), nrow=3, byrow=TRUE)
    load(anova.fil)
    res <- post.process.random(results)
    com.res <- cbind(m.labs, c("CFO", "CRM", "BOIN"), res)
    
    anova.data.list[[i]] <- com.res
}


anova.data <- do.call(rbind, anova.data.list)
names(anova.data) <- c(c("nlevels", "CohortSize", "SampleSizes", "DiffProbs", "Targets", "Methods"), names(anova.data)[7:12])
anova.data.df <- as.data.frame(anova.data)
rownames(anova.data.df) <- NULL
anova.data.df


fit <- aov(MTD.Sel~nlevels+CohortSize+SampleSizes+DiffProbs+Targets+Methods, data=anova.data.df)
summary(fit)

dat <- group_by(anova.data.df, Methods)
dat <- group_by(anova.data.df,  SampleSizes, Methods)
summarise(dat, MTD.Sel=mean(MTD.Sel), MTD.Allo=mean(MTD.Allo))
summarise(dat, Over.Sel=mean(Over.Sel), Over.Allo=mean(Over.Allo), 
          Risk.of.HT=mean(Risk.of.HT), PerDLT=mean(PerDLT))


# Targets is less important
# nlevel has some effect
# DiffProbs
sub.data <- filter(anova.data.df, CohortSize==1)
sub.data %>% group_by(SampleSizes, Methods) %>%
        summarise(MTD.Sel=mean(MTD.Sel), MTD.Allo=mean(MTD.Allo))
