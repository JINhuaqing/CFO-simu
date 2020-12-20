#library(bestNormalize)
library(dplyr)
setwd("C:/Users/Dell/Documents/ProjectCode/CFO")
setwd("C:/Users/JINHU/Documents/ProjectCode/CFO")
source("./Rcode/utilities.R")
source("./Rcode/phaseI/anova_settings.R")
anova.ress <- dir("./results/Anovas", pattern="*.RData", full.names = TRUE)

sep.paths <- strsplit(anova.ress, "_")
raw.labs <- sapply(sep.paths, function(i)i[3])
tmp.fn <- function(i){
    if (length(i)==5){
        return(as.numeric(i))
    }else{
       num.i <- as.numeric(i) 
       rv <- c(num.i[1:2], 10*num.i[3]+num.i[4], num.i[5:6])
       return(rv)
    }
}
labs <- lapply(strsplit(raw.labs, ""), tmp.fn)
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

anova.data.df <- filter(anova.data.df, CohortSize==1|CohortSize==3,
                        SampleSizes==21|SampleSizes==30|SampleSizes==48|SampleSizes==60,
                        DiffProbs!=0.15)
fit <- aov(MTD.Sel~(nLevels+CohortSize+SampleSizes+DiffProbs+Targets+Methods)*(nLevels+CohortSize+SampleSizes+DiffProbs+Targets+Methods), data=anova.data.df)
#fit <- aov(MTD.Sel~nLevels+CohortSize+SampleSizes+DiffProbs+Targets+Methods, data=anova.data.df)
summary(fit)
# main factors, DiffProbs, nLevels, SampleSize

# fit diagnostic
errs <- fit$residuals
norm.errs <- errs/sd(errs)
ks.test(norm.errs, pnorm)
#plot(norm.errs)
#norm.errs %>% density %>% plot
#rnorm(10000) %>% density %>% lines(col="red")

fit.method <- aov(MTD.Allo~Methods, data=anova.data.df)
summary(fit.method)



dat <- group_by(anova.data.df, Methods)
dat <- group_by(anova.data.df, Methods, SampleSizes)
dat <- group_by(anova.data.df, Methods, DiffProbs)
dat <- group_by(anova.data.df, Methods, nLevels)
summarise(dat, MTD.Sel=mean(MTD.Sel), MTD.Allo=mean(MTD.Allo))
summarise(dat, Over.Sel=mean(Over.Sel), Over.Allo=mean(Over.Allo), 
          Risk.of.HT=mean(Risk.of.HT), PerDLT=mean(PerDLT))


# Targets is less important
# nlevel has some effect
# DiffProbs
#sample sizes 21, 30, 48, 60
sub.data <- filter(anova.data.df, 
                   SampleSizes==42)
sub.data %>% group_by(Methods) %>%
        summarise(MTD.Sel=mean(MTD.Sel), MTD.Allo=mean(MTD.Allo))

MTD.Sel <- anova.data.df$MTD.Sel
SST <- sum((MTD.Sel - mean(MTD.Sel))**2)
sum.fit <- summary(fit)
SSE.main <- sum(sum.fit[[1]]$`Sum Sq`[c(1, 3, 4)])
SSE.main/SST
