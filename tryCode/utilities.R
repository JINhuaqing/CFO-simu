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


#latex.out.fn(res, prefix="m=1", 1, "[1pt]")

