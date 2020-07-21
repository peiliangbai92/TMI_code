
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


#################################################
### universal parameter settings
T <- 6000     ### the number of observations
p <- 20       ### dimension of time series
brk <- c(floor(T/6), floor(T/3), floor(T/2), floor(2*T/3), floor(5*T/6), T+1)   ### true interval end points
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
max_iter <- 1
block_sizes <- c(40, 60, 80, 100, 150, 200)
# block_sizes <- c(100, 150, 200, 500, 600, 700, 1000, 1200, 1400, 1600)
cases_ret <- vector('list', length(block_sizes))
for(i in 1:length(block_sizes)){
    block_size <- block_sizes[i]
    running_times <- rep(0, max_iter)
    estimated_cps <- vector('list', max_iter)
    estimated_mats <- vector('list', max_iter)
    for(epoch in 1:max_iter){
        ### generate time series
        set.seed(123456*epoch)
        e.sigma <- as.matrix(1*diag(p));
        try <- var.sim.break(T, arlags=q.t, malags=NULL, phi = phi.full, sigma = e.sigma, brk = brk)
        data <- try$series
        data <- as.matrix(data)
        
        ### run the bss method, get the estimated change points and a_n
        begin_time <- Sys.time()
        ret <- bss(data, block.size = block_size)
        # ret <- bss(data)
        end_time <- Sys.time()
        
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
        running_times[epoch] <- end_time - begin_time
        estimated_cps[[epoch]] <- ret$final.selected.points
        estimated_mats[[epoch]] <- estimated_blockmats
        
        print(epoch)
        print(end_time - begin_time)
    }
    ### record results for cases lists
    cases_ret[[i]]$runningTime <- running_times
    cases_ret[[i]]$CPs <- estimated_cps
    cases_ret[[i]]$Matrices <- estimated_mats
    
    print(paste("Finished block size case:", block_size))
}
# save(cases_ret, file = "results.RData")

### evaluate performance
# 1. calculate average running time
run_times <- NULL
for(i in 2:length(block_sizes)){
    rt <- mean(cases_ret[[i]]$runningTime)
    run_times <- c(run_times, rt)
    print(paste("average running times is:", rt))
}

# run_times <- c(177.07, 85.28, 49.11, 43.04, 17.18, 44.85, 59.86, 89.88, 107.64, 137.94)

# plot averaged running time in line plot
plot(block_sizes, (as.numeric(run_times)), type = "o", col = 1, pch = 21, bg = 1,
     lty = 1,lwd = 1.2 ,cex = 1.2,
     ylim = c(min((as.numeric(run_times))), max(as.numeric(run_times))),
     ylab ='avg. comp. time(sec)', xlab = 'block sizes' ,cex.lab = 1.5, cex.axis = 1.5)

# evaluate the mean, sd, and selection rate of estimated CPs
true_brk <- brk[1:5]
cps <- cases_ret[[6]]$CPs
cp1 = cp2 = cp3 = cp4 = cp5 <- NULL
for(i in 1:max_iter){
    iter_cps <- cps[[i]]
    if(length(iter_cps) > 0){
        for(j in 1:length(iter_cps)){
            current_cp <- iter_cps[j]
            if(current_cp <= 1200 && current_cp >= 800){
                cp1 <- c(cp1, current_cp)
            }else if(current_cp <= 2200 && current_cp >= 1800){
                cp2 <- c(cp2, current_cp)
            }else if(current_cp <= 3200 && current_cp >= 2800){
                cp3 <- c(cp3, current_cp)
            }else if(current_cp <= 4200 && current_cp >= 3800){
                cp4 <- c(cp4, current_cp)
            }else if(current_cp <= 5200 && current_cp >= 4800){
                cp5 <- c(cp5, current_cp)
            }
        }
    }else{
        next
    }
}
mean(cp1/T); sd(cp1/T); length(cp1)/max_iter
mean(cp2/T); sd(cp2/T); length(cp2)/max_iter 
mean(cp3/T); sd(cp3/T); length(cp3)/max_iter
mean(cp4/T); sd(cp4/T); length(cp4)/max_iter
mean(cp5/T); sd(cp5/T); length(cp5)/max_iter

# 2. evaluate the transition matrices
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

mats <- cases_ret[[6]]$Matrices
sen = spc = acc <- NULL
re <- NULL
for(i in 1:max_iter){
    est_mats <- mats[[i]]
    if(ncol(est_mats) / nrow(est_mats) == 6){
        sen <- c(sen, evaluate(phi.full, est_mats, method = 'SEN'))
        spc <- c(spc, 1-evaluate(phi.full, est_mats, method = 'FNR'))
        acc <- c(acc, evaluate(phi.full, est_mats, method = 'ACC'))
        re <- c(re, norm(phi.full - est_mats, "F")/norm(phi.full, "F"))
    }
}
mean(sen)
mean(spc)
mean(re)