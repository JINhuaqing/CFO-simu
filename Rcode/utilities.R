library(magrittr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

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



#ress <- lapply(1:10000, function(i)gen.rand.doses(7, 0.3, mu1=0.70, mu2=0.50))
#sapply(ress, prob.diff.fn, target=0.3) %>% mean
#tmp <- lapply(ress, prob.diff.fn.sep, target=0.3)
#sapply(tmp, function(i)i$L) %>% unlist %>% mean
#sapply(tmp, function(i)i$U) %>% unlist %>% mean




post.process.onemethod <- function(res, paras){
    tmtd <- paras$mtd
    target <- paras$target
    ndose <- length(res$p.true)
    rv <- rep(0, 7)
    rv[1] <- sum(res$MTD==tmtd)
    rv[2] <- res$dose.ns[tmtd]
    if (res$MTD == 99){
        rv[3] <- 0
    }else{
        rv[3] <- sum(res$MTD>tmtd)
    }
    if (tmtd==ndose){
        rv[4] <- 0
    }else{
        rv[4] <- sum(res$dose.ns[(tmtd+1):ndose])
    }
    rv[5] <- sum((sum(res$DLT.ns)/sum(res$dose.ns))>target)
    rv[6] <- sum(res$DLT.ns)
    rv[7] <- sum(res$dose.ns)
    names(rv) <- c(
        "MTD Sel", "MTD Allo", "Over Sel", 
        "Over Allo", "Risk of HT", "No DLT",
        "No Subject")
    rv
    
}
post.process.single <- function(result){
    nMethod <- length(result)-1
    paras <- result$paras
    tmtd <- paras$mtd
    methods <- names(result)[1:nMethod]
    res.v <- list()
    for (name in methods){
        res.v[[name]] <- post.process.onemethod(result[[name]], paras)
    }
    do.call(rbind, res.v)
}

# Function to handle the results under random scenario
post.process.random <- function(results){
    nsimu <- length(results)
    res.all <- 0
    for (result in results){
        res.all <- post.process.single(result) + res.all
    }
    res.all.df <- data.frame(res.all)
    final.res <- transmute(res.all.df, 
                           MTD.Sel=MTD.Sel/nsimu,
                           MTD.Allo=MTD.Allo/No.Subject,
                           Over.Sel=Over.Sel/nsimu,
                           Over.Allo=Over.Allo/No.Subject,
                           Risk.of.HT=Risk.of.HT/nsimu,
                           PerDLT=No.DLT/No.Subject
    )
    rownames(final.res) <- rownames(res.all.df)
    final.res
}


# function to plot the results under random scenarios
res.plot.fn <- function(res, filename, M.names, is.save=TRUE, angle=45){
  if (missing(filename)){
    is.save <- FALSE
  }
  myBarplot <- function(data, labels_vec, rot_angle, main, ylim=NULL) {
    nM <- dim(res)[1]
    plt <- barplot(data, col=2:(nM+1), xaxt="n", main=main, ylim=ylim)
    text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), 
         xpd = TRUE, cex=0.8, font=1) 
  }
  if (missing(M.names)){
      M.names <- c("Butterfly-BB", "Butterfly-CRM", "Butterfly-BB-CRM", "CRM", "BOIN")
  }
  if (is.save){
    png(filename, width=360*2, height=240*2)
  }
  par(mfrow=c(2, 3))
  myBarplot(res$MTD.Sel, main="MTD selection %", labels_vec=M.names, rot_angle=angle)
  myBarplot(res$MTD.Allo, main="MTD Allocation %", labels_vec=M.names, rot_angle=angle)
  myBarplot(res$Over.Sel, main="Overdose selection %", labels_vec=M.names, rot_angle=angle)
  myBarplot(res$Over.Allo, main="Overdose Allocation %", labels_vec=M.names, rot_angle=angle)
  myBarplot(res$Risk.of.HT, main="Risk of high toxicity %", labels_vec=M.names, rot_angle=angle)
  myBarplot(res$PerDLT, main="Avergae DLT rate %", labels_vec=M.names, rot_angle=angle)
  par(mfrow=c(1, 1))
  if (is.save){
    dev.off()
  }
}

res.ggplot.fn <- function(res, filename, main, M.names, is.save=TRUE, angle=NULL){
    if (missing(filename)){
        is.save <- FALSE
    }
    myggBarplot <- function(data, angle, main, theylim) {
         g <- ggplot(mapping=aes(x=M.names, y=data, fill=M.names))
         g <- g + geom_bar(stat="identity") + guides(fill=FALSE) + ggtitle(main)
         g <- g + theme(axis.title = element_blank())
         if (!is.null(angle)){
             g <- g + theme(axis.text.x = element_text(angle=angle, vjust=1, hjust=1))
         }
         if (!missing(theylim)){
             g <- g + ylim(theylim)
         }
             
         g
    }
    if (missing(M.names)){
        M.names <- c("Butterfly-BB", "Butterfly-CRM", "Butterfly-BB-CRM", "CRM", "BOIN")
    }
    g1 <- myggBarplot(res$MTD.Sel, main="MTD selection %", angle=angle)
    g2 <- myggBarplot(res$MTD.Allo, main="MTD Allocation %", angle=angle)
    g3 <- myggBarplot(res$Over.Sel, main="Overdose selection %", angle=angle)
    g4 <- myggBarplot(res$Over.Allo, main="Overdose Allocation %", angle=angle)
    g5 <- myggBarplot(res$Risk.of.HT, main="Risk of high toxicity %", angle=angle)
    g6 <- myggBarplot(res$PerDLT, main="Avergae DLT rate %", angle=angle)
    gg <- grid.arrange(g1, g2, g3, g4, g5, g6, nrow=3, top=textGrob(main, gp=gpar(fontsize=20, fontface=2)))
    if (is.save){
        ggsave(filename, gg, width=8, height=12, unit="in")
    }
}


# Plot the randomly generated dose levels
gen.rand.doses.plot <- function(cases, phi){
    ncase <- length(cases)
    case <- cases[[1]]
    ndose <- length(case$p.trues)
    plot(1:ndose, case$p.trues, type="b", ylim=c(0, 1), ylab="Prob of toxicity", xlab="Dose levels", xaxt="n")
    points(case$mtd.level, case$p.trues[case$mtd.level], pch=19, col=1)
    axis(1, at=1:ndose, labels=1:ndose)
    for (i in 2:ncase){
        ccase <- cases[[i]]
        lines(1:ndose, ccase$p.trues, type="b", ylim=c(0, 1), lty=i)
        points(ccase$mtd.level, ccase$p.trues[ccase$mtd.level], pch=19, col=1)
        i <- i + 1
    }
    abline(h=phi, col=2, lwd=2)
    legend("topleft", legend=c("Target prob", "MTD level"), col=c(2, 1), pch=c(NA, 19), 
           lty=c(1, NA), lwd=c(2, 1))
    
}

