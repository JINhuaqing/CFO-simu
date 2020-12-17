library(dplyr)
setwd("C:/Users/Dell/Documents/ProjectCode/CFO")
setwd("C:/Users/JINHU/Documents/ProjectCode/phaseI")
source("./Rcode/utilities.R")
source("./Rcode/phaseI/anova_settings.R")
anova.ress <- dir("./results/Anovas", pattern="*.RData", full.names = TRUE)

sep.paths <- strsplit(anova.ress, "_")
raw.labs <- sapply(sep.paths, function(i)i[3])
labs <- lapply(strsplit(raw.labs, ""), function(i)as.numeric(i))
mat.labs <- do.call(rbind, labs)

anova.data.list <- list()
for (i in 1:length(anova.ress)){
    print(i)
    anova.fil <- anova.ress[i]
    s.labs <- mat.labs[i, ]
    s.true.labs <- rep(0, length(s.labs))
    s.true.labs[1] <- nlevels[s.labs[1]]
    s.true.labs[2] <- csizes[s.labs[2]]
    s.true.labs[3] <- sample.sizes[s.labs[3]]
    s.true.labs[4] <- diff.probs[s.labs[4]]
    s.true.labs[5] <- targets[s.labs[5]]
    
    m.labs <- matrix(rep(s.true.labs, 3), nrow=3, byrow=TRUE)
    load(anova.fil)
    res <- post.process.random(results)
    com.res <- cbind(m.labs, c("CFO", "CRM", "BOIN"), res)
    
    anova.data.list[[i]] <- com.res
}


anova.data <- do.call(rbind, anova.data.list)
names(anova.data) <- c(c("nLevels", "CohortSize", "SampleSizes", 
                         "DiffProbs", "Targets", "Methods"), 
                       names(anova.data)[7:12])
anova.data.df <- as.data.frame(anova.data)
rownames(anova.data.df) <- NULL
anova.data.df

anova.data.df$nLevels <- as.factor(anova.data.df$nLevels)
anova.data.df$CohortSize <- as.factor(anova.data.df$CohortSize)
anova.data.df$SampleSizes <- as.factor(anova.data.df$SampleSizes)
anova.data.df$DiffProbs <- as.factor(anova.data.df$DiffProbs)
anova.data.df$Targets <- as.factor(anova.data.df$Targets)


fit <- aov(MTD.Sel~nLevels+CohortSize+SampleSizes+DiffProbs+Targets+Methods, data=anova.data.df)
summary(fit)

# Sample sizes 30, 60
dat <- group_by(anova.data.df, Methods)
dat <- group_by(anova.data.df, SampleSizes, Methods)
summarise(dat, MTD.Sel=mean(MTD.Sel), MTD.Allo=mean(MTD.Allo))
summarise(dat, Over.Sel=mean(Over.Sel), Over.Allo=mean(Over.Allo), 
          Risk.of.HT=mean(Risk.of.HT), PerDLT=mean(PerDLT))


# Targets is less important
# nlevel has some effect
# DiffProbs
#sample sizes 30, 42
sub.data <- filter(anova.data.df, 
                   SampleSizes==48|SampleSizes==30|SampleSizes==42, CohortSize==1|CohortSize==3)
sub.data %>% group_by(Methods) %>%
        summarise(MTD.Sel=mean(MTD.Sel), MTD.Allo=mean(MTD.Allo))
