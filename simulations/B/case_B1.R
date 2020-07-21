setwd("d:/Work/Dropbox (UFL)/NeuroImage/simulations/new-simus/B")
#################################################
rm(list=ls())
gc()

######## Loading Packages #######################
library("mvtnorm")
library("lattice")
library("Rcpp")
library("RcppArmadillo")
library("MTS")
library("sparsevar")
library("ggplot2")

######## Call Functions #########################
source("Functions_BSS.R")
#cpp code for block fused lasso and block lasso
sourceCpp("Functions_BSS.cpp")

# evaluate function
evaluate <- function(true_mat, est_mat, method){
    method.full <- c('SEN', 'FNR', 'ACC')
    nr <- nrow(true_mat)
    nc <- ncol(true_mat)
    
    ### calculate confusion matrix
    tp = fn = tn = fp <- 0
    for(i in 1:nr){
        for(j in 1:nc){
            if(true_mat[i,j] != 0){
                if(est_mat[i,j] != 0){
                    tp <- tp + 1
                }else if(est_mat[i,j] == 0){
                    fn <- fn + 1
                }
            }else if(true_mat[i,j] == 0){
                if(est_mat[i,j] != 0){
                    fp <- fp + 1
                }else if(est_mat[i,j] == 0){
                    tn <- tn + 1
                }
            }
        }
    }
    if(method == 'SEN'){
        return(tp / (tp + fn))
    }else if(method == 'FNR'){
        return(fp / (fp + tn))
    }else if(method == 'ACC'){
        return((tp + tn) / (tp + tn + fp + fn))
    }
}

#################################################
### universal parameter settings
T <- 20000    ### the number of observations
p <- 20       ### dimension of time series
brk <- NULL
for(i in 1:19){
    brk <- c(brk, floor(T*i/20))
}
brk <- c(brk, T+1)                            ### true interval end points
m0 <- length(brk) - 1   ### true number of break points
q.t <- 1      ### true time lag

### generate time series
m <- m0 + 1
phi.full <- matrix(0, p, p*q.t*m)
for(mm in 1:m){
    phi.full[,((mm-1)*q.t*p+1):((mm)*q.t*p)] <- 0
    if(mm %% 2 == 0 ){
        temp <- matrix(0, p, p)
        for(j in 1:(p-1)){
            temp[j,j+1] <- 0.8
        }
        phi.full[,((mm-1)*q.t*p+1):((mm)*q.t*p)]  <- temp
    }else{
        temp <- matrix(0, p, p)
        for(j in 1:(p-1)){
            temp[j,j+1] <- -0.8
        }
        phi.full[,((mm-1)*q.t*p+1):((mm)*q.t*p)]  <- temp
    }
}  

#display the true AR coefficients matrice
print(plot.matrix((phi.full), m))

#################################################
max_iter <- 50
cases_ret <- list()
estimated_times <- NULL
estimated_cps <- vector('list', max_iter)
estimated_mats <- vector('list', max_iter)
for(epoch in 1:max_iter){
    ### generate time series
    set.seed(123456*epoch)
    e.sigma <- as.matrix(0.1*diag(p));
    try <- var.sim.break(T, arlags = q.t, malags = NULL, phi = phi.full, sigma = e.sigma, brk = brk)
    data <- try$series
    data <- as.matrix(data)
    
    ### run the bss method, get the estimated change points and a_n
    begin_time <- Sys.time()
    ret <- bss(data)
    end_time <- Sys.time()
    print(paste("running time is:", end_time - begin_time))
    
    ### get estimated transition matrices
    final.estimated.intervals <- c(1, ret$final.selected.points, T)
    an <- ret$an
    estimated_blockmats <- NULL
    for(j in 2:length(final.estimated.intervals)){
        s <- final.estimated.intervals[j-1] + floor(an/2)
        e <- final.estimated.intervals[j] - floor(an/2)
        
        interval_data <- data[s:e, ]
        var.fit <- fitVAR(interval_data, penalty = 'ENET', alpha = 1, intercept = FALSE)
        temp_mat <- var.fit$A[[1]]
        estimated_blockmats <- cbind(estimated_blockmats, temp_mat)
    }
    
    ### record results into lists 
    estimated_times <- c(estimated_times, (end_time - begin_time))
    estimated_cps[[epoch]] <- ret$final.selected.points
    estimated_mats[[epoch]] <- estimated_blockmats
    
    print(epoch)
}
### record results for cases lists
cases_ret$times <- estimated_times
cases_ret$CPs <- estimated_cps
cases_ret$Matrices <- estimated_mats
save(cases_ret, file = "results_B1.RData")


### evaluate performance
# 1. calculate average running time
for(i in 1:max_iter){
    if(estimated_times[i] < 50){
        estimated_times[i] <- estimated_times[i] * 60
    }
}
mean(estimated_times)

# evaluate the mean, sd, and selection rate of estimated CPs
true_brk <- brk[1:(length(brk)-1)]
cp_loc <- vector('list', length(true_brk))
for(i in 1:max_iter){
    iter_cps <- estimated_cps[[i]]
    if(length(iter_cps) > 0){
        for(j in 1:length(iter_cps)){
            current_cp <- iter_cps[j]
            if(current_cp <= 1200 && current_cp >= 800){
                cp_loc[[1]] <- c(cp_loc[[1]], current_cp/T)
            }else if(current_cp <= 2200 && current_cp >= 1800){
                cp_loc[[2]] <- c(cp_loc[[2]], current_cp/T)
            }else if(current_cp <= 3200 && current_cp >= 2800){
                cp_loc[[3]] <- c(cp_loc[[3]], current_cp/T)
            }else if(current_cp <= 4200 && current_cp >= 3800){
                cp_loc[[4]] <- c(cp_loc[[4]], current_cp/T)
            }else if(current_cp <= 5200 && current_cp >= 4800){
                cp_loc[[5]] <- c(cp_loc[[5]], current_cp/T)
            }else if(current_cp <= 6200 && current_cp >= 5800){
                cp_loc[[6]] <- c(cp_loc[[6]], current_cp/T)
            }else if(current_cp <= 7200 && current_cp >= 6800){
                cp_loc[[7]] <- c(cp_loc[[7]], current_cp/T)
            }else if(current_cp <= 8200 && current_cp >= 7800){
                cp_loc[[8]] <- c(cp_loc[[8]], current_cp/T)
            }else if(current_cp <= 9200 && current_cp >= 8800){
                cp_loc[[9]] <- c(cp_loc[[9]], current_cp/T)
            }else if(current_cp <= 10200 && current_cp >= 9800){
                cp_loc[[10]] <- c(cp_loc[[10]], current_cp/T)
            }else if(current_cp <= 11200 && current_cp >= 10800){
                cp_loc[[11]] <- c(cp_loc[[11]], current_cp/T)
            }else if(current_cp <= 12200 && current_cp >= 11800){
                cp_loc[[12]] <- c(cp_loc[[12]], current_cp/T)
            }else if(current_cp <= 13200 && current_cp >= 12800){
                cp_loc[[13]] <- c(cp_loc[[13]], current_cp/T)
            }else if(current_cp <= 14200 && current_cp >= 13800){
                cp_loc[[14]] <- c(cp_loc[[14]], current_cp/T)
            }else if(current_cp <= 15200 && current_cp >= 14800){
                cp_loc[[15]] <- c(cp_loc[[15]], current_cp/T)
            }else if(current_cp <= 16200 && current_cp >= 15800){
                cp_loc[[16]] <- c(cp_loc[[16]], current_cp/T)
            }else if(current_cp <= 17200 && current_cp >= 16800){
                cp_loc[[17]] <- c(cp_loc[[17]], current_cp/T)
            }else if(current_cp <= 18200 && current_cp >= 17800){
                cp_loc[[18]] <- c(cp_loc[[18]], current_cp/T)
            }else if(current_cp <= 19200 && current_cp >= 18800){
                cp_loc[[19]] <- c(cp_loc[[19]], current_cp/T)
            }
        }
    }else{
        next
    }
}
lapply(cp_loc, FUN = 'mean')
lapply(cp_loc, FUN = 'sd')
for(i in 1:length(cp_loc)){
    print(length(cp_loc[[i]])/max_iter)
}

# 2. evaluate the transition matrices
mats <- cases_ret$Matrices
sen = spc = acc <- NULL
re <- NULL
for(i in 1:max_iter){
    est_mats <- mats[[i]]
    if(ncol(est_mats) / nrow(est_mats) == 20){
        sen <- c(sen, evaluate(phi.full, est_mats, method = 'SEN'))
        spc <- c(spc, 1-evaluate(phi.full, est_mats, method = 'FNR'))
        acc <- c(acc, evaluate(phi.full, est_mats, method = 'ACC'))
        re <- c(re, norm(phi.full - est_mats, "F")/norm(phi.full, "F"))
    }
}
mean(sen)
mean(spc)
mean(re)