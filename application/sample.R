#### This sample is for subject 8 EEG data

rm(list=ls())
library("MTS")
library("sparsevar")
library("factoextra")
library("edfReader")
library("forecast")
library("Rcpp")
library("RcppArmadillo")
source("functions_BFL.R")
sourceCpp("functions_BFL.cpp")

#### loading data
bdf_header <- readEdfHeader("database/EEG_Cat_Study4_Resting_S8.bdf")
signal <- readEdfSignals(bdf_header)

### 20 channels case
channels <- c('RHEOG', 'LHEOG', 'POz', 'PO8', 'PO7', 'P8', 'P7',
              'T8', 'C4', 'Cz', 'C3', 'T7', 'F8', 'F4', 'Fz', 'F3', 'F7', 'Fp1', 'Fp2', 'Fpz')
ts_dataset <- c()
for(ch in channels){
    temp <- signal[[ch]]$signal
    ts_dataset <- cbind(ts_dataset, temp)
    print(paste('Loading... channel:', ch, sep = ' '))
}
colnames(ts_dataset) <- channels

###################
ts_dataset <- scale(ts_dataset)
data <- na.omit(ts_dataset)   # remove NA values (after detrending, we have NA for the begining and ending time series)
# data <- data1[30000:109999,]   # remove the first redundant 8 resting seconds
# MTSplot(data)

################################# Robust detrending #######################################
detrended_ts <- c()
for(j in 1:dim(data)[2]){
    ### initial step, set weight equal to 1
    feature_ <- data[,j]
    weight <- rep(1, length(feature_))
    new_weight <- 0
    
    ### check if the weight vector will be changed or not
    iter <- 1
    while(sqrt(sum((new_weight - weight)^2)) > 1e-3 && iter < 10){
        new_weight <- weight
        smooth_fit <- smooth.spline(feature_, w = new_weight)
        smoothed_vals <- smooth_fit$y
        d <- (feature_ - smoothed_vals)
        
        ### update weight
        temp_weight <- abs(d) / sd(abs(d))
        weight <- as.numeric(temp_weight > 0.1)
        iter <- iter + 1
        print(sqrt(sum((new_weight - weight)^2)))
    }
    detrended_ts <- cbind(detrended_ts, feature_ - smoothed_vals)
    print(paste("Feature", j, "is done!", sep = " "))
}

### down-sampling step
period <- dim(detrended_ts)[1]/(256*8)
start <- period*256*2
end <- period*256*7
period*256/16
true_cp <- c(1:4)* (period*256/16)
true_cp
detrended_ts <- detrended_ts[start:end,]
# detrended_ts <- detrended_ts[30000:109999,]
freq <- 16
sub_data <- c()
for(t in 1:dim(detrended_ts)[1]){
    if(t %% freq == 0){
        sub_data <- rbind(sub_data, detrended_ts[t,])
    }
}
sub_data <- scale(sub_data)

### remove outliers
for(t in 1:dim(sub_data)[1]){
    sub_data[t, sub_data[t,] > 5] <- 5
    sub_data[t, sub_data[t,] < -5] <- -5
}
MTSplot(sub_data)
abline(v = true_cp, col ='red', lwd = 2)

########################### use detrended data to run VAR model ############################
T <- dim(sub_data)[1]
k <- dim(sub_data)[2]
lambda.1.cv <- c(20, 10, 5, 1, 0.5, 0.1, 0.01, 0.001, 0.0001)
lambda.2.cv <- c(0.1*sqrt(log(k)/T))
b_n <- floor(0.8*sqrt(T))
b_n_bound = 200 #block size for boundary
blocks <- c(seq(1,b_n_bound*2+1 , b_n_bound),
            seq(b_n_bound*2+b_n+1, T+1-2*b_n_bound, b_n),
            seq(T+1-b_n_bound ,T+1,  b_n_bound));
if(blocks[length(blocks)] < T+1){
    blocks <- c(blocks[-length(blocks)], T+1)
}
print(blocks)
fit <- tbfl(method = "VAR", data_y = sub_data, max.iteration = 100, tol = 5e-3, 
            blocks = blocks, 
            lambda.1.cv = lambda.1.cv, 
            lambda.2.cv = lambda.2.cv
)
ret <- fit$cp.final

############################################################################################
############################ estimate segments for the data ################################
est_mats <- vector('list', length(true_cp)+1)
T <- dim(sub_data)[1]
k <- dim(sub_data)[2]
breaks <- c(1, true_cp, T)
for(i in 1:(length(breaks)-1)){
    s <- breaks[i] + 20
    e <- breaks[i+1] - 20
    
    ### segments
    temp_data <- sub_data[s:e, ]
    fit <- fitVAR(temp_data, p=1, alpha = 1, intercept = FALSE)
    mat <- fit$A[[1]]
    
    ### select the edges with threshold 0.1
    for(r in 1:k){
        for(c in 1:k){
            if(abs(mat[r,c]) < 0.1){
                mat[r,c] <- 0
            }
        }
    }
    est_mats[[i]] <- mat
}



