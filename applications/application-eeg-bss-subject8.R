library(MTS)
library(edfReader)
library(sparsevar)
library(forecast)
library(doParallel)
library(lattice)
library(changepoint)


detrending <- function(x){
    td <- ma(x, order = 256)
    return(x - td)
}

bdf_header <- readEdfHeader("database/EEG_Cat_Study4_Resting_S8.bdf")
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
# MTSplot(data)

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
# lambda.1.cv <- c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)   ### good one
lambda.1.cv <- c(1e-1, 1e-2, 1e-3, 1e-4, 5e-4)
lambda.2.cv <-  c(0.01*sqrt(log(k)/T))
# lambda.2.cv <- 0

set.seed(123456)
# b_n <- floor(1.75*sqrt(T))    ### good one
b_n <- floor(5*sqrt(T))
blocks <- seq(0, T, b_n)

ptm <- Sys.time()
temp <- bss("LASSO", data, lambda.1.cv, lambda.2.cv, p, max.iteration = max.iteration,
             tol = tol, step.size = (1/2)*10^(-4), blocks, cpp = TRUE)
print(paste("total time:", Sys.time() - ptm))


MTSplot(data)
pts.final <- temp$pts.3[[1]]
abline(v = pts.final, col = 'red', lwd = 2)
save(temp, file = 'estimated_cps_s8_var1_allchannels.RData')

####################################################################################
###### estimate segment Granger Causal networks
T <- dim(data)[1]
k <- dim(data)[2]
pts.final <- temp$pts.3[[1]]
an <- temp$an
m <- length(pts.final)

####### use the segments to estimate Granger causal networks #######
segments <- vector('list', m+1)
segments[[1]] <- data[3000:3199,]
segments[[2]] <- data[26000:26199,]
segments[[3]] <- data[42000:42199,]
segments[[4]] <- data[53000:53199,]
segments[[5]] <- data[70000:70199,]

library(NetworkDistance)
library(ggplot2)
library(reshape2)

lambda.1 <- c(3e-4, .75e-4, 3e-4, .75e-4, 3e-4)
networks <- vector('list', m+1)
for(i in 1:(m+1)){
    seg_data <- segments[[i]]
    seg.fit <- var_lasso_brk(seg_data, lambda = lambda.1[i], p = 1, max_iteration = 500, tol = 1e-4)
    networks[[i]] <- seg.fit$phi.hat
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

save(networks, file = "estimated_network_s8.RData")

###### plot estimated networks
est.phi <- NULL
for(i in 1:(m+1)){
    est.phi <- cbind(est.phi, networks[[i]])
}
print(plot.matrix(abs(est.phi), m+1))

### truncate the estimated networks, only keeps the magnitudes larger than 0.05
### construct the adjacency matrices and plot networks
library(eegkit)
library(igraph)
data("eegcoord")
channels <- c(toupper(names(signal)[1:64]), "NZ", "F10", "F9", "TP9", "TP10")

### plot EEG cap
eegcap(channels, type = '2d', plotaxes = FALSE, ears = TRUE, cex.label = .75, 
       cex.point = 3.25, col.label = 'black', col.point = 'green')

### create EEG channels coordinates
eegcoords <- eegcoord[toupper(names(signal)[1:64]), c("xproj", "yproj")]

### add additional eeg channels and their coordinates
eegcoords <- rbind(eegcoords, eegcoord[c('TP9', 'TP10', 'NZ'), c("xproj", "yproj")])
VEOG <- matrix(c(-5.5263, 5.5164,11.11425, 11.11425), 2, 2)
rownames(VEOG) <- c('LVEOG', 'RVEOG')
colnames(VEOG) <- c('xproj', 'yproj')
eegcoords <- rbind(eegcoords, VEOG)
eegcoords <- rbind(eegcoords, eegcoord[c('F9', 'F10'), c('xproj', 'yproj')])
NFPZ <- matrix(c(0.012599879, 14.15510), nrow = 1, ncol = 2)
rownames(NFPZ) <- c('NFPZ')
colnames(NFPZ) <- c('xproj', 'yproj')
eegcoords <- rbind(eegcoords, NFPZ)

n <- dim(eegcoords)[1]
k <- 72
nodes.coord <- as.data.frame(eegcoords)
for(i in 1:(m+1)){
    adjmat <- networks[[i]]
    rownames(adjmat) <- rownames(eegcoords)
    colnames(adjmat) <- rownames(eegcoords)
    for(r in 1:k){
        for(c in 1:k){
            if(abs(adjmat[r,c]) < 0.1){
                adjmat[r,c] <- 0
            }else{
                adjmat[r,c] <- 1
            }
        }
    }
    gg1 <- graph_from_adjacency_matrix(adjmat, mode = 'directed', diag = FALSE)
    plot.igraph(gg1, layout = as.matrix(nodes.coord[, c("xproj","yproj")]), edge.arrow.size = 0.25,
         vertex.color = 'lightblue', vertex.size = 14, vertex.label.dist = 0, vertex.label.cex = 0.75,
         edge.color = 'darkgrey')
}

####### plot the differences between two consecutive segments #######
adj_mats <- vector('list', m+1)
for(i in 1:(m+1)){
    temp <- networks[[i]]
    for(r in 1:k){
        for(c in 1:k){
            if(abs(temp[r,c]) < 0.2){
                temp[r,c] <- 0
            }else{
                temp[r,c] <- 1
            }
        }
    }
    adj_mats[[i]] <- temp
    colnames(adj_mats[[i]]) = rownames(adj_mats[[i]]) <- rownames(eegcoords)
}

remaining_edges <- matrix(0, k, k)
increment_edges = decrease_edges <- matrix(0, k, k)
colnames(remaining_edges) = colnames(increment_edges) = colnames(decrease_edges) <- rownames(eegcoords)
rownames(remaining_edges) = rownames(increment_edges) = rownames(decrease_edges) <- rownames(eegcoords)

for(i in 1:m){
    for(r in 1:k){
        for(c in 1:k){
            if(adj_mats[[i+1]][r,c] == 1 && adj_mats[[i]][r,c] == 1){
                remaining_edges[r,c] <- 1
            }else if(adj_mats[[i+1]][r,c] == 1 && adj_mats[[i]][r,c] == 0){
                increment_edges[r,c] <- 1
            }else if(adj_mats[[i+1]][r,c] == 0 && adj_mats[[i]][r,c] == 1){
                decrease_edges[r,c] <- 1
            }
        }
    }
    g_increase <- graph_from_adjacency_matrix(increment_edges, mode = 'directed', diag = FALSE) %>% set_edge_attr("color", value = "red")
    g_decrease <- graph_from_adjacency_matrix(decrease_edges, mode = 'directed', diag = FALSE) %>% set_edge_attr("color", value = "blue")
    
    g_union <- g_increase %u% g_decrease
    
    color_union <- edge_attr(g_union, 'color_1')
    E2 <- which(edge_attr(g_union, 'color_2') == 'blue')
    color_union[E2] <- edge_attr(g_union, 'color_2')[E2]
    g_union2 <- set_edge_attr(g_union, "color", value = color_union)
    plot.igraph(g_union2, layout = as.matrix(nodes.coord[, c("xproj","yproj")]), edge.arrow.size = 0.25, 
         vertex.color = 'lightblue', vertex.size = 14, vertex.label.dist = 0, vertex.label.cex = 0.75)
}


####### plot Hamming distance for all estimated segments #######
library(NetworkDistance)
library(ggplot2)
library(reshape2)

### calculate Hamming distance
output_networks <- nd.hamming(networks)
dist_ <- as.matrix(output_networks$D)

### plot Hamming distance matrix as heatmap
rownames(dist_) = colnames(dist_) <- c('Segment 1', 'Segment 2', 'Segment 3', 'Segment 4', 'Segment 5')
melt_dist_ <- melt(dist_)
ggheatmap <- ggplot(data = melt_dist_, aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = 'white') + xlab("") + ylab("") + 
    geom_text(aes(label = round(value, 3))) + scale_fill_gradient(low = 'white', high = 'red') + guides(fill = FALSE)
print(ggheatmap)


####### calculate the basic statistics for each networks #######
for(i in 1:length(networks)){
    colnames(networks[[i]]) = rownames(networks[[i]]) <- rownames(eegcoords)
    temp <- networks[[i]]
    edges <- 0
    for(r in 1:nrow(temp)){
        for(c in 1:ncol(temp)){
            if(abs(temp[r,c]) > 0.1){
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
