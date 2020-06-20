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

# Gnerate dose toxicity relations randomly
gen.rand.doses <- function(ndose, target, mu1=0.55, mu2=0.55, inv.fn=qnorm, fn=pnorm){
    sigma0 <- 0.05
    sigma1 <- 0.35
    sigma2 <- 0.35
    raw.p.trues <- rep(0, ndose)
    mtd.level <- sample.int(ndose, 1)
    raw.p.trues[mtd.level] <- rnorm(1, inv.fn(target), sigma0)
    
    if (mtd.level != 1){
        eps.minus <- rnorm(1, mu1, sigma1)
        if (raw.p.trues[mtd.level] > inv.fn(target)){
            dlt <- raw.p.trues[mtd.level] - inv.fn(2*target-fn(raw.p.trues[mtd.level]))
        }else{
            dlt <- 0 
        }
        raw.p.trues[mtd.level-1] <- raw.p.trues[mtd.level] - (dlt +  eps.minus**2)
    }
    
    if (mtd.level != ndose){
        eps.plus <- rnorm(1, mu2, sigma2)
        if (raw.p.trues[mtd.level] < inv.fn(target)){
            dlt <- inv.fn(2*target-fn(raw.p.trues[mtd.level])) - raw.p.trues[mtd.level]
        }else{
            dlt <- 0 
        }
        raw.p.trues[mtd.level+1] <- raw.p.trues[mtd.level] + (dlt +  eps.plus**2)
    }
    
    if ((mtd.level-2)>0){
        for (i in (mtd.level-2):1){
           eps.minus <- rnorm(1, mu1, sigma1)
           raw.p.trues[i]  <- raw.p.trues[i+1] - eps.minus**2
        }
    }
    if ((mtd.level+2)<=ndose){
        for (j in (mtd.level+2):ndose){
           eps.plus <- rnorm(1, mu2, sigma2)
           raw.p.trues[j] <- raw.p.trues[j-1] + eps.plus**2
        }
    }
    p.trues <- fn(raw.p.trues)
    list(p.trues=p.trues, mtd.level=mtd.level)
}



# compute the prob diff around the target level
prob.diff.fn <- function(res, target=0.3){
    p.trues <- res$p.trues
    mtd <- res$mtd.level
    ndose <- length(p.trues)
    
    diffs <- c()
    if (mtd!=1){
        diffs <- c(abs(p.trues[mtd-1]-target), diffs)
    }

    if (mtd!=ndose){
        diffs <- c(abs(p.trues[mtd+1]-target), diffs)
    }
    min(diffs)
}

# Compute the prob diff around the target level, separately
prob.diff.fn.sep <- function(res, target=0.3){
    p.trues <- res$p.trues
    mtd <- res$mtd.level
    ndose <- length(p.trues)
    
    diffsL <- NULL
    diffsU <- NULL
    if (mtd!=1){
        diffsL <- abs(p.trues[mtd-1]-target)
    }

    if (mtd!=ndose){
        diffsU <- abs(p.trues[mtd+1]-target)
    }
    list(L=diffsL, U=diffsU) 
}



# ress <- lapply(1:10000, function(i)gen.rand.doses(3, 0.3, mu1=0.55, mu2=0.40))
# sapply(ress, prob.diff.fn, target=0.3) %>% mean
# tmp <- lapply(ress, prob.diff.fn.sep, target=0.3)
# sapply(tmp, function(i)i$L) %>% unlist %>% mean
# sapply(tmp, function(i)i$U) %>% unlist %>% mean
# gen.rand.doses(5, 0.3, 0.65, 0.45)









