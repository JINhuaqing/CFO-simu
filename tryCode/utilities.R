library(magrittr)
# return latex scr code for output
latex.out.fn <- function(res, prefix, tidx){
  MTDs <- res$MTDs.percent  
  MTDs <- round(MTDs*100, 1)
  if (!missing(tidx)){
      MTDs[tidx] <- paste0("\bf{", MTDs[tidx], "}")
  }
  av.dose <- res$av.dose
  tss <- res$t.dose
  tdlt <- res$t.DLT
  row1 <- paste("&", MTDs, collapse=" ") 
  row1 <- paste(prefix, "& MTD %", row1); 
  
  row2 <- paste("&", round(av.dose, 1), collapse=" ") 
  row2 <- paste("& ASN", row2, "&", round(tdlt, 1), "&", round(tss, 1))
  print(row1)
  print(row2)
}


latex.out.fn(res, prefix="m=1", 1)
