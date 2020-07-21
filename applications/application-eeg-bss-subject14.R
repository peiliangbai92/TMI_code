library(MTS)
library(edfReader)
library(sparsevar)
library(forecast)
library(doParallel)
library(lattice)


detrending <- function(x){
    td <- ma(x, order = 256)
    return(x - td)
}

bdf_header <- readEdfHeader("database/EEG_Cat_Study4_Resting_S14.bdf")
signal <- readEdfSignals(bdf_header)

idx <- 1:72
ts_dataset <- c()
for(ch in 1:72){
    temp <- signal[[ch]]$signal
    ts_dataset <- cbind(ts_dataset, temp)
    print(paste('Loading... channel:', ch, sep = ' '))
}
# colnames(ts_dataset) <- channels
colnames(ts_dataset) <- names(signal)[1:72]

ts_dataset <- ts_dataset[-(1:2560),]
ts_dataset <- apply(ts_dataset, 2, FUN = detrending)
ts_dataset <- scale(ts_dataset)
data1 <- na.omit(ts_dataset)   # remove NA values (after detrending, we have NA for the begining and ending time series)

data <- data1[30000:109999,]   # remove the first redundant 8 resting seconds
MTSplot(data)

######### Use BSS algorithm to estimate change-point ##########
library(Rcpp)
library(RcppArmadillo)
source("functions_BSS_C.R")
sourceCpp("var_break_fit_block.cpp")
sourceCpp("var_lasso_brk.cpp")


### Parameters settings
T <- dim(data)[1]
k <- dim(data)[2]

tol <- 2*10^(-2)                ### tolerance
step.size <- 2*10^(-4)          ### step size
p = p.t <- 1                    ### the selected AR order
max.iteration <- 200            ### max iteration for blocked detection code.

### Tuning parameters settings ###
lambda.1.cv <- c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
lambda.2.cv <-  c(0.01*sqrt(log(k)/T))
# lambda.2.cv <- 0

set.seed(123456)
b_n <- floor(1.75*sqrt(T))
blocks <- seq(0, T, b_n)

ptm <- proc.time()
temp <- bss("LASSO", data, lambda.1.cv, lambda.2.cv, p, max.iteration = max.iteration,
            tol = tol, step.size = (1/2)*10^(-4), blocks, cpp = TRUE)
print(proc.time() - ptm)


MTSplot(data)
pts.final <- temp$pts.3[[1]]
abline(v = pts.final, col = 'red', lwd = 2)

####################################################################################
###### estimate segment Granger Causal networks
T <- dim(data)[1]
k <- dim(data)[2]
pts.final <- temp$pts.3[[1]]
an <- temp$an
m <- length(pts.final)

segments <- vector('list', m+1)
segments[[1]] <- data[3000:3199,]
segments[[2]] <- data[26000:26199,]
segments[[3]] <- data[42000:42199,]
segments[[4]] <- data[53000:53199,]
segments[[5]] <- data[70000:70199,]

library(NetworkDistance)
library(ggplot2)
library(reshape2)

# lambda.1 <- c(1e-4, 1e-4, 1e-4, 1e-4, 1e-4)
lambda.1 <- c(3e-4, .75e-4, 3e-4, .75e-4, 3e-4)
networks <- vector('list', m+1)
for(i in 1:(m+1)){
    seg_data <- segments[[i]]
    seg.fit <- var_lasso_brk(seg_data, lambda = lambda.1[i], p = 1, max_iteration = 500, tol = 1e-4)
    networks[[i]] <- seg.fit$phi.hat
    # seg.fit <- fitVAR(seg_data, alpha = 1)
    # networks[[i]] <- seg.fit$A[[1]]
    print(i)
}

adj_mats <- vector('list', m+1)
for(i in 1:(m+1)){
    temp <- networks[[i]]
    for(r in 1:nrow(temp)){
        for(c in 1:ncol(temp)){
            if(abs(temp[r,c]) > 0){
                temp[r,c] <- 1
            }else{
                temp[r,c] <- 0
            }
        }
    }
    adj_mats[[i]] <- temp
}

output_networks <- nd.hamming(adj_mats)
dist_ <- as.matrix(output_networks$D)

### plot Hamming distance matrix as heatmap
rownames(dist_) = colnames(dist_) <- c('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5')
melt_dist_ <- melt(dist_)
ggheatmap <- ggplot(data = melt_dist_, aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = 'white') + xlab("") + ylab("") + 
    geom_text(aes(label = round(value, 3))) + scale_fill_gradient(low = 'white', high = 'red') + guides(fill = FALSE)
print(ggheatmap)

for(i in 1:(m+1)){
    print(networks[[i]][1:20,1:20])
}

save(networks, file = "estimated_network_s14.RData")

####### calculate the basic statistics for each networks #######
for(i in 1:length(networks)){
    colnames(networks[[i]]) = rownames(networks[[i]]) <- rownames(eegcoords)
    temp <- networks[[i]]
    edges <- 0
    for(r in 1:nrow(temp)){
        for(c in 1:ncol(temp)){
            if(abs(temp[r,c]) > 0.2){
                edges <- edges + 1
                temp[r,c] <- 1
            }else{
                temp[r,c] <- 0
            }
        }
    }
    print(edges)
    
    nconn_vertex <- apply(temp, 2, FUN = 'sum')
    print(sort(nconn_vertex, decreasing = TRUE)[1:5])
}
