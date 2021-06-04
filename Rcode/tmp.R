#source("../utilities.R")
#load("../res_plateau_4.RData")
#phase.I.II.pretty.tb(sum.all)


# Table 4 MTD comparison
DLTs <- c(9.6, 9.4, 9.7, 8.6, 8.7, 8.4, 7.4, 7.5, 6.9, 5.0, 4.9, 4.2)
numSubs <- list()
numSubs[[1]] <- list()
numSubs[[2]] <- list()
numSubs[[3]] <- list()
numSubs[[4]] <- list()

numSubs[[1]]$CFO <- c(19.5, 6.8, 1.0, 0.1, 0.0)
numSubs[[1]]$BOIN <- c(18.8, 6.8, 1.1, 0.1, 0.0)
numSubs[[1]]$CRM <- c(20.0, 7.4, 0.8, 0.0, 0.0)

numSubs[[2]]$CFO <- c(11.1, 14.0, 4.1, 0.4, 0.0)
numSubs[[2]]$BOIN <- c(11.2, 13.3, 4.4, 0.5, 0.0)
numSubs[[2]]$CRM <- c(11.1, 15.1, 3.3, 0.2, 0.0)

numSubs[[3]]$CFO <- c(6.4, 9.8, 9.0, 3.8, 1.0)
numSubs[[3]]$BOIN <- c(6.3, 9.5, 8.7, 4.1, 1.1)
numSubs[[3]]$CRM <- c(6.3, 12.5, 8.3, 2.6, 0.2)

numSubs[[4]]$CFO <- c(3.1, 3.6, 4.3, 6.1, 12.8)
numSubs[[4]]$BOIN <- c(3.1, 3.7, 4.2, 6.7, 12.3)
numSubs[[4]]$CRM <- c(3.1, 4.0, 6.0, 8.0, 8.9)

numSubs.list <- lapply(1:4, function(i)do.call(rbind, numSubs[[i]]))
numSubs.mat <- do.call(rbind, numSubs.list)
Table4.data <- cbind(numSubs.mat, DLTs)
tol.numSubs <- rowSums(Table4.data[, 1:5])
rep.tol.mat <- matrix(rep(tol.numSubs, 6), ncol=6)

round(100*Table4.data/rep.tol.mat, 1)[10:12, ]


# Table 5 OBD comparison
Table5.list <- list()

Table5.list[[1]] <- c(10.8, 15.6, 29.8, 3.5, 0.2, 12.9, 23.6)
Table5.list[[2]] <- c(6.5 , 15.5, 28.2, 4.7, 1.5, 13.5, 22.9)
Table5.list[[3]] <- c(9.7 , 21.0, 23.4, 4.9, 1.0, 12.7, 22.9)
Table5.list[[4]] <- c(9.2 , 18.3, 28.6, 3.3, 0.3, 12.8, 23.4)
Table5.list[[5]] <- c(15.5, 31.1, 9.5 , 3.0 , 0.6 , 14.2, 25.0) 
Table5.list[[6]] <- c(7.4 , 16.8, 16.9, 6.7 , 10.8, 17.0, 27.0) 
Table5.list[[7]] <- c(19.1, 24.3, 10.9, 3.7 , 1.4 , 14.2, 24.2) 
Table5.list[[8]] <- c(12.5, 31.6, 10.7, 3.1 , 0.3 , 14.1, 25.2) 
Table5.list[[9]] <- c(11.0, 36.8, 10.2, 1.6 , 0.5 , 12.2, 31.7) 
Table5.list[[10]] <- c(12.6, 22.4, 18.1, 3.5 , 2.8 , 12.9, 28.8) 
Table5.list[[11]] <- c(15.5, 24.2, 13.4, 4.6 , 2.2 , 12.5, 28.7) 
Table5.list[[12]] <- c(11.4, 37.5, 9.3 , 1.4 , 0.3 , 12.2, 31.4) 
Table5.list[[13]] <- c(10.2, 12.9, 31.4, 4.6 , 0.9 , 12.5, 18.8) 
Table5.list[[14]] <- c(4.4 , 7.2 , 26.5, 7.6 , 5.4 , 13.3, 17.1) 
Table5.list[[15]] <- c(9.1 , 17.9, 24.5, 6.7 , 1.9 , 12.8, 17.2) 
Table5.list[[16]] <- c(6.9 , 12.7, 31.7, 6.0 , 1.1 , 13.1, 18.8) 
Table5.list[[17]] <- c(8.8 , 11.5, 12.1, 12.8, 14.7, 6.4 , 32.5) 
Table5.list[[18]] <- c(8.2 , 8.8 , 14.7, 11.4, 16.9, 6.5 , 33.0) 
Table5.list[[19]] <- c(4.8 , 6.9 , 9.2 , 12.3, 26.8, 7.4 , 36.3) 
Table5.list[[20]] <- c(8.0 , 14.4, 15.5, 12.1, 10.0, 6.0 , 30.9) 
Table5.list[[21]] <- c(35.2, 12.7, 2.3 , 0.2 , 0.0 , 15.0, 9.7 ) 
Table5.list[[22]] <- c(7.7 , 12.3, 7.9 , 3.0 , 3.8 , 14.9, 11.7) 
Table5.list[[23]] <- c(39.9, 13.9, 2.8 , 0.5 , 0.1 , 17.0, 11.1) 
Table5.list[[24]] <- c(26.4, 14.2, 3.3 , 0.3 , 0.0 , 14.0, 9.3 ) 

Table5.mat <- do.call(rbind, Table5.list)
colnames(Table5.mat) <- c(1:5, "DLT", "Eff")
rownames(Table5.mat) <- rep(c("CFO", "efftox", "MADA", "STEIN"), 6)
tol.subjs5 <- rowSums(Table5.mat[, 1:5])
tol.mat5 <- matrix(rep(tol.subjs5, 7), ncol=7)
i <- 1
round(100*Table5.mat/tol.mat5, 1)[21:24, ]

# Table 2
MTD.sel.list <- list()

MTD.sel.list[[1]] <- c(35.9, 43.1,  51.8,  62.2) 
MTD.sel.list[[2]] <- c(36.3, 42.8,  50.9,  60.6) 
MTD.sel.list[[3]] <- c(34.0, 41.7,  51.5,  65.3)
MTD.sel.list[[4]] <- c(43.4, 47.6,  50.6,  51.5) 
MTD.sel.list[[5]] <- c(43.0, 47.0,  49.8,  50.7) 
MTD.sel.list[[6]] <- c(43.6, 47.1,  50.2,  51.5)
MTD.sel.list[[7]] <- c(52.5, 48.8,  46.9,  44.9) 
MTD.sel.list[[8]] <- c(51.3, 48.3,  46.6,  44.4) 
MTD.sel.list[[9]] <- c(52.5, 48.8,  46.8,  44.4) 

MTD.sel <- do.call(rbind, MTD.sel.list)


MTD.allo.list <- list()
MTD.allo.list[[1]] <- c(31.7, 36.7, 43.0, 51.1)
MTD.allo.list[[2]] <- c(31.4, 35.4, 40.6, 47.2)
MTD.allo.list[[3]] <- c(30.4, 35.6, 42.2, 51.8)
MTD.allo.list[[4]] <- c(34.7, 38.7, 43.5, 45.4)
MTD.allo.list[[5]] <- c(33.1, 36.8, 41.4, 43.2)
MTD.allo.list[[6]] <- c(34.2, 38.1, 42.8, 44.9)
MTD.allo.list[[7]] <- c(50.6, 41.6, 36.9, 32.4)
MTD.allo.list[[8]] <- c(48.8, 39.6, 34.8, 31.3)
MTD.allo.list[[9]] <- c(50.1, 41.1, 36.3, 32.4)

MTD.allo <- do.call(rbind, MTD.allo.list)

rowMeans(MTD.sel)
rowMeans(MTD.allo)
