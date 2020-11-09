rm(list=ls(all=TRUE))
gc()

library(MTS)
library(edfReader)
library(forecast)
library("mvtnorm")
library("xtable")
library("lattice")
library("factoextra")
library("Rcpp")
library("RcppArmadillo")
library("GGMselect")
library("network")
library('Matrix')
library('glmnet')
library('glasso')
######## FUNCTIONS ##########################
source("functions_TBFL.R")
sourceCpp("functions_TBFL.cpp")



####################################################
detrending <- function(x){
    td <- ma(x, order = 256)
    return(x - td)
}

channels <- c('Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'T7', 'C3', 'Cz', 'C4', 'T8',
              'P7', 'P3', 'Pz', 'P4', 'P8', 'O1', 'Oz', 'O2')
prefix <- 'database/EEG_Cat_Study4_Resting_S'

i <- 8
name <- paste(prefix, i, '.bdf', sep = '')
print(name)
bdf_header <- readEdfHeader(name)
signal <- readEdfSignals(bdf_header)

ts_dataset <- c()
channels <- c()
# for(ch in channels){
for(ch in 1:(length(signal)-1) ){
    temp <- signal[[ch]]$signal
    ts_dataset <- cbind(ts_dataset, temp)
    channels <- c(channels, signal[[ch]]$name)
    print(paste('Loading... channel:', ch, sep = ' '))
}
colnames(ts_dataset) <- channels

ts_dataset <- ts_dataset[(65000:95000),]
# choose the observation every 10 points
ts_dataset <- ts_dataset[seq(1,nrow(ts_dataset), by = 10),]
ts_dataset <- apply(ts_dataset, 2, FUN = detrending)
ts_dataset <- scale(ts_dataset)
data1 <- na.omit(ts_dataset) # remove NA values (after detrending, we have NA for the begining and ending time series)

MTSplot(data1)
print('======================================')


method <- c("GGM")
b_n <- 50
ptm.temp <- proc.time()
temp <- tbfl(method, data_y = data1, block.size =  b_n)
time.temp <- proc.time() - ptm.temp

cp.final <- temp$cp.final
MTSplot(data1)
abline(v = cp.final, lwd = 2)


selectFast(data1[1:1500,], family="C01")

s <-  var(data1[1:1000, ])
a1 <-glassopath(s)
plot.omega.matrix(a1$wi[, , 1], p = 1)
# plot.omega.matrix(a1$wi[, , 10], p = 1)

s <-  var(data1[1500:2745, ])
a2 <-glassopath(s)
plot.omega.matrix(a2$wi[, , 1], p = 1)
# plot.omega.matrix(a1$wi[, , 10], p = 1)
