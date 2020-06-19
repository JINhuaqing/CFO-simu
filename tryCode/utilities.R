library(magrittr)
# return latex scr code for output
latex.out.fn <- function(res, prefix, tidx, suffix=NULL){
  MTDs <- res$MTDs.percent  
  MTDs <- format(round(MTDs*100, 1), nsmall=1)
  if (!missing(tidx)){
      MTDs[tidx] <- paste0("\\bf{", MTDs[tidx], "}")
  }
  av.dose <- res$av.dose
  tss <- res$t.dose
  tdlt <- res$t.DLT
  row1 <- paste("&", MTDs, collapse=" ") 
  row1 <- paste(prefix, "& MTD \\%", row1, "\\\\"); 
  
  row2 <- paste("&", format(round(av.dose, 1), nsmall=1), collapse=" ") 
  row2 <- paste("& ASN", row2, "&", format(round(tdlt, 1), nsmall=1), "&",
                format(round(tss, 1), nsmall=1), "\\\\")
  out <- paste(row1, "\n", row2, suffix)
  cat(out)
}


logistic <- function(x){
    num <- 1
    den <- 1+exp(-x)
    return(num/den)
}

# Gnerate dose toxicity relations randomly
gen.rand.doses <- function(ndose, target, eps=0.05, diff=0.15){
    raw.p.trues <- rep(0, ndose)
    mtd.level <- sample.int(ndose, 1)
    raw.p.trues[mtd.level] <- rnorm(1, logit(target), 1)
    
    if (mtd.level != 1){
        eps.minus <- rnorm(1, 0, 1)
        if (raw.p.trues[mtd.level] > logit(target)){
            dlt <- raw.p.trues[mtd.level] - logit(target)
        }else{
            dlt <- 0 
        }
        raw.p.trues[mtd.level-1] <- raw.p.trues[mtd.level] - (dlt +  eps.minus**2)
    }
    
    if (mtd.level != ndose){
        eps.plus <- rnorm(1, 0, 1)
        if (raw.p.trues[mtd.level] < logit(target)){
            dlt <- logit(target) - raw.p.trues[mtd.level]
        }else{
            dlt <- 0 
        }
        raw.p.trues[mtd.level+1] <- raw.p.trues[mtd.level] + (dlt +  eps.plus**2)
    }
    
    if ((mtd.level-2)>0){
        for (i in (mtd.level-2):1){
           eps.minus <- rnorm(1, 0, 1)
           raw.p.trues[i]  <- raw.p.trues[i+1] - eps.minus**2
        }
    }
    if ((mtd.level+2)<=ndose){
        for (j in (mtd.level+2):ndose){
           eps.plus <- rnorm(1, 0, 1)
           raw.p.trues[j] <- raw.p.trues[j-1] + eps.plus**2
        }
    }
    p.trues <- logistic(raw.p.trues)
    list(p.trues=p.trues, mtd.level=mtd.level)
}

# Gnerate dose toxicity relations randomly
gen.rand.doses <- function(ndose, target, eps=0.05, diff=0.125){
    p.trues <- rep(0, ndose)
    mtd.level <- sample.int(ndose, 1)
    p.trues[mtd.level] <- runif(1, target-eps, target+eps)
    
    if ((mtd.level-1)>=1){
        for (i in (mtd.level-1):1){
           p.trues[i]  <- p.trues[i+1] - runif(1, diff-0.025, diff+0.025)
        }
    }
    if ((mtd.level+1)<=ndose){
        for (j in (mtd.level+1):ndose){
           p.trues[j] <- p.trues[j-1] + runif(1, diff-0.025, diff+0.025)
        }
    }
    list(p.trues=p.trues, mtd.level=mtd.level)
}


gen.rand.doses(5, 0.3)
