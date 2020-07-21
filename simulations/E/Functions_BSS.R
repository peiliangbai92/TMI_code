#' Generating non-stationary ARMA data.
#' 
#' @param nobs number of time points
#' @param arlags the true AR order
#' @param malags the true MA order
#' @param cnst the constant
#' @param phi parameter matrix of the AR model
#' @param theta parameter matrix of the MA model
#' @param sigma covariance matrix of the white noise
#' @return Matrice of time series data and white noise data 

var.sim.break <- function (nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL,theta = NULL, skip = 200, sigma, brk = nobs+1) {
  if (!is.matrix(sigma)) 
    sigma = as.matrix(sigma)
  k <- nrow(sigma); m <- length(brk); nT <- nobs + skip
  
  #generate multivariate normal distributed data as the white noise data
  at <- rmvnorm(nT, rep(0, k), sigma)
  
  #generate the ARMA time series data
  nar <- length(arlags); p <- 0
  if (nar > 0) {
    arlags <- sort(arlags)
    p <- arlags[nar]
  }
  
  nma <- length(malags); q <- 0
  if (nma > 0) {
    malags <- sort(malags)
    q <- malags[nma]
  }
  
  ist = max(p, q) + 1
  zt = matrix(0, nT, k)
  if (length(cnst) == 0) 
    cnst = rep(0, k)
  if (m == 1){
    for (it in ist:nT) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
  }
  
  #if there are some break points
  if (m > 1){
    for (it in ist:(skip+brk[1]-1)) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
    for ( mm in 1:(m-1)){
      for (it in (skip+brk[mm]):(skip+brk[mm+1]-1) ) {
        tmp = matrix(at[it, ], 1, k)
        if (nma > 0) {
          for (j in 1:nma) {
            jdx = (j - 1) * k
            thej = theta[, (jdx + 1):(jdx + k)]
            atm = matrix(at[it - malags[j], ], 1, k)
            tmp = tmp - atm %*% t(thej)
          }
        }
        if (nar > 0) {
          for (i in 1:nar) {
            idx = (i - 1) * k
            phj = phi[, ((mm)*p*k+idx + 1):((mm)*p*k+idx + k)]
            ztm = matrix(zt[it - arlags[i], ], 1, k)
            tmp = tmp + ztm %*% t(phj)
          }
        }
        zt[it, ] = cnst + tmp
      }
    }
  }
  
  zt = zt[(1 + skip):nT, ]
  at = at[(1 + skip):nT, ]
  VARMAsim <- list(series = zt, noises = at)
}


#' block segmentation scheme (BSS).
#' 
#' @description Perform the block segmentation scheme (BSS) algorithm to detect the structural breaks 
#' in large scale high-dimensional non-stationary VAR models.
#' 
#' @param data input data matrix, with each column representing the time series component 
#' @param lambda.1.cv tuning parmaeter lambda_1 for fused lasso
#' @param lambda.2.cv tuning parmaeter lambda_2 for fused lasso
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso 
#' @param block.size the block size
#' @param blocks the blocks
#' @return A list oject, which contains the followings
#' \describe{
#'   \item{pts.1}{a set of selected break point after the first block fused lasso step}
#'   \item{pts.2}{a set of selected break point after the second local screening step}
#'   \item{pts.3}{a set of selected break point after the thrid exhaustive search step}
#'   \item{an}{the selected neighborhood size a_n after the grid search}
#' }

bss <- function( data, lambda.1.cv= NULL, lambda.2.cv = NULL, q=1, max.iteration = 100, tol = 10^(-2), block.size = NULL, blocks=NULL){
  
  T <- length(data[,1]); p <- length(data[1,]); 
  second.brk.points <- c(); pts.final <- c();
  
  ############# block size and blocks ###########
  if(is.null(block.size) && is.null(blocks) ){
    block.size = floor(sqrt(T));
    blocks <- seq(0,T,block.size);
  }else if( !is.null(block.size) && is.null(blocks)){
    blocks <- seq(0,T,block.size);
  }else if(!is.null(block.size) && !is.null(blocks)){
    #check if the block.size and blocks match
    n.new <- length(blocks) - 1;
    blocks.size.check <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
    if( sum(blocks.size.check[1: (length(blocks.size.check)-1 )] != block.size ) >0 ){
      stop("Error: The block.size and blocks can't match!")
    }
  }
  
  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  
  
  #sample the cv index for cross-validation
  bbb <- floor(n.new/5);
  aaa <- sample(1:5, 1);
  cv.index <- seq(aaa,n.new,floor(n.new/bbb));
  
  ############# Tuning parameter ################
  if(is.null(lambda.1.cv)){
    lambda.1.max <- lambda_warm_up(data, q, blocks, cv.index)$lambda_1_max
    if(blocks[2] <= 2*p ){
      epsilon <-  10^(-3)
    }
    if(blocks[2] >= 2*p ){
      epsilon <-  10^(-4)
    }
    nlam <- 10 
    lambda.1.min <-  lambda.1.max*epsilon
    delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
    lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
  }
  
  if(is.null(lambda.2.cv)){
    lambda.2.cv <-  c(10*sqrt(log(p)/T),1*sqrt(log(p)/T),0.10*sqrt(log(p)/T))
  }

  ######################################################################
  ######## First Step: Initial Break Points Selection ##################
  ######################################################################
  
  #run the first block fused lasso step 
  temp.first <- first.step.blocks( data, lambda.1.cv, lambda.2.cv, q, max.iteration = max.iteration, tol = tol, cv.index, blocks=blocks)
  first.brk.points <- temp.first$brk.points;

  #construct the grid values of neighborhood size a_n
  n <- T - q;
  an.lb <- max(floor(mean(blocks.size)),floor( (log(n)*log(p))^1 ));
  an.ub <-  min(10*an.lb,0.95*(min(first.brk.points)-1-q),0.95*(n - max(first.brk.points)-1))
  an.grid <- seq(an.lb,  an.ub, length.out = 5);
  an.idx.final <- length(an.grid)
  an.grid <- floor(an.grid);
  final.pts.res <- vector("list",length(an.grid));
  flag <- c("FALSE");
  an.idx <- 0;
  
  #for each a_n, run the second and thrid step
  while(an.idx < length(an.grid)  ){
    an.idx <- an.idx + 1;
    an <- an.grid[an.idx]
    
    #remove the boundary points
    remove.ind <- c();
    if(length(first.brk.points) != 0){
      for(i in 1:length(first.brk.points)){
        if ( first.brk.points[i] < (an-1-q)   ){remove.ind <- c(remove.ind,i);}
        if ( (T-first.brk.points[i]) < (an-1-q)   ){remove.ind <- c(remove.ind,i);}
      }
    }
    if( length(remove.ind) > 0  ){first.brk.points <- first.brk.points[-remove.ind];}
    
    #if there are selected break points after the first step
    if( length(first.brk.points) != 0){
      
      #####################################################
      ######## Second Step: Local Screening      ##########
      #####################################################
      eta <- (1/1)*(log(2*an)*log(p))/(2*an); # the tuning parameter for second and third steps.
      #run the second local screening step
      temp <- second.step.local(data, eta = eta, q, max.iteration = 1000, tol = tol, first.brk.points, an)
      
      #record the selected break points after local screening step
      second.brk.points <- temp$pts;
      
      ######################################################
      ######## Thrid Step: Exhaustive Search      ##########
      ######################################################
      pts.final <- second.brk.points;
      #keep running until none of the selected break points close to any other selected break points 
      while( min(abs(diff(pts.final)),3*an) <  2*an  ){
        if( length(pts.final) != 0){
          #cluster the selected break points by size 2a_n
          pts.list <- block.finder(pts.final,2*an)
          # run the third exhaustive search step for each cluster
          pts.final <- third.step.exhaustive(data, q, max.iteration = 1000, tol = tol, pts.list, an , eta )
        }
      }
      
      #record the final selected break points for each given a_n
      final.pts.res[[an.idx]] <- pts.final
      
      #terminate the grid search of an if the number of final selected break points is stable
      if(an.idx > 2){
        if( length(final.pts.res[[an.idx]]) == length(final.pts.res[[an.idx-1]]) && length(final.pts.res[[an.idx-1]]) == length(final.pts.res[[an.idx-2]]) ){
          flag <- c("TRUE");
          an.idx.final <- an.idx;
          an.sel <- an.grid[an.idx];
          break;
        }
      }
      
    }
  }
  
  #if the stable criterion hasn't been met
  #find the length that happen the most
  #if there are multiple lengths with same occurrence, find the longest one
  if(flag == FALSE){
    loc.final <- rep(0,length(an.grid));
    for(i in 1:length(an.grid)){
      loc.final[i] <- length(final.pts.res[[i]]);
    }
    loc.table <- table(loc.final)
    counts.final <- sort(loc.table,decreasing=TRUE)[1]
    len.final <- max(as.integer(names(loc.table)[loc.table == counts.final]))
    an.idx.final <- max(c(1:length(loc.final))[loc.final == len.final])
    an.sel <- an.grid[an.idx.final];
  }
  
  return(list(first.selected.points = first.brk.points, second.selected.points = second.brk.points, final.selected.points = final.pts.res[[an.idx.final]], final.selected.points.grid = final.pts.res, an = an.sel)) 
}


#' block fused lasso step (first step).
#' 
#' @description Perform the block fused lasso to detect candidate break points.
#' 
#' @param data.temp input data matrix, with each column representing the time series component 
#' @param lambda.1.cv tuning parmaeter lambda_1 for fused lasso
#' @param lambda.2.cv tuning parmaeter lambda_2 for fused lasso
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso 
#' @param cv.index the index of time points for cross-validation
#' @param blocks the blocks
#' @return A list oject, which contains the followings
#' \describe{
#'   \item{brk.points}{a set of selected break point after the first block fused lasso step}
#'   \item{cv}{the cross validation values for tuning parmeter selection}
#'   \item{cv1.final}{the selected lambda_1}
#'   \item{cv2.final}{the selected lambda_2}
#' }

first.step.blocks <- function(data.temp, lambda.1.cv, lambda.2.cv, q, max.iteration = max.iteration, tol = tol,cv.index, blocks){
  
  cv.l <- length(cv.index); data.org <- data.temp; T.org <- length(data.temp[,1]); p <- length(data.temp[1,]); n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  
  #create the tuning parmaeter combination of lambda1 and lambda2
  lambda.full <- expand.grid(lambda.1.cv,lambda.2.cv)
  kk <- length(lambda.full[,1]);
  
  cv <- rep(NA,kk); 
  phi.final <- vector("list",kk);
  T <- length(data.temp[,1]); p <- length(data.temp[1,]);
  brk.points.final <- vector("list",kk);
  flag.full <- rep(0,kk);
  
  #cross-validation for each values of lambda1 and lambda2
  nlam1 <- length(lambda.1.cv)
  nlam2 <- length(lambda.2.cv)
  kk <- nlam1*nlam2
  i = 1
  while(i <= kk) {
    i.lam1 <- i%% nlam1
    if(i.lam1 == 0 ){
      i.lam1 = nlam1
    }
    i.lam2 <- floor((i-1)/nlam1)+1
    if ( i == 1){
      test <- var_break_fit_block_cpp(data.temp, lambda.full[i,1],lambda.full[i,2], q, max.iteration, tol = tol, initial_phi = 0.0+matrix(0.0,p,p*q*n.new), blocks, cv.index)
      flag.full[i] <- test$flag;
    }else if(is.na(cv[i-1]) ){
      test <- var_break_fit_block_cpp(data.temp, lambda.full[i,1],lambda.full[i,2], q, max.iteration, tol = tol, initial_phi = 0.0+matrix(0.0,p,p*q*n.new), blocks, cv.index)
      flag.full[i] <- test$flag;
    }else{
      initial.phi <- phi.final[[(i-1)]]
      if(max(abs(phi.final[[(i-1)]])) > 10^3  ){initial.phi <- 0*phi.final[[(i-1)]];}
      test <- var_break_fit_block_cpp(data.temp, lambda.full[i,1],lambda.full[i,2], q, max.iteration, tol = tol, initial_phi = initial.phi, blocks, cv.index)
      flag.full[i] <- test$flag;
    }
    phi.hat.full <- test$phi.hat;
    phi.final[[i]] <- phi.hat.full;
    
    ll <- c(0);
    brk.points.list <- vector("list",length(ll));
    
    for(j in 1:length(ll)){
      phi.hat <- phi.hat.full;
      n <- T - q;
      m.hat <- 0; brk.points <- rep(0,n.new);
      
      for (iii in 1:(n.new-1))
      {
        if ( sum((phi.hat[,((iii-1)*p*q+1):(iii*p*q)] )^2 ) > tol   ){
          m.hat <- m.hat + 1; brk.points[m.hat] <- blocks[iii+1];
        }
      }
      
      
      loc <- rep(0,m.hat);
      brk.points <- brk.points[1:m.hat];
      
      #remove the bounary points and clean up 
      brk.points <- brk.points[which(brk.points > 3*mean(blocks.size))]; brk.points <- brk.points[which(brk.points < (n-3*mean(blocks.size)))];
      m.hat <- length(brk.points);
      del <- 0;
      if(m.hat >= 2){
        while(del < m.hat){
          if(length(brk.points) <= 1){break;}
          del <- del + 1;
          deleted <- 0;
          for (i.3 in 2:length(brk.points)) {
            if(deleted == 0 &&  abs(brk.points[i.3] - brk.points[i.3-1]) <= (1)*max(q,min(blocks.size))  ){
              brk.points <- brk.points[-i.3]; deleted <- 1;
            }
          }
        }
      }
      
      brk.points.list[[j]] <- brk.points;
    }
    
    brk.points.final[[i]] <- brk.points;
    m.hat <- length(brk.points);
    
    #forecast the time series based on the estimated matrix Phi
    #and compute the forecast error
    phi.full.all <- vector("list",n.new);
    forecast <- matrix(0,p,T);
    phi.full.all[[1]] <- phi.hat[,(1):(p*q)];
    for(i.1 in 2:n.new){
      phi.full.all[[i.1]] <- phi.full.all[[i.1-1]] + phi.hat[,((i.1-1)*p*q+1):(i.1*p*q)];
      forecast[,(blocks[i.1]+1):(blocks[i.1+1])] <- pred.block(t(data.org),phi.full.all[[i.1-1]],q,blocks[i.1],p,blocks[i.1+1]-blocks[i.1]);
    }
    forecast.new <- matrix(0,p,cv.l);
    for(j in (1):cv.l){
      forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j])]],q,blocks[cv.index[j]+1]-1,p,1)
    }
    temp.index <- rep(0,cv.l);
    for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1];}
    cv[i] <- (1/(p*cv.l))*sum( (forecast.new - t(data.org[temp.index,])  )^2 );
    
    #break condition  
    if( !(i %in% seq(1, kk, nlam1))  && cv[i] > cv[i-1]){
      i.lam2 <- i.lam2 + 1
      i <- (i.lam2-1)*nlam1 + 1
    }
    else{
      i <- i+1
    }
  }
  
  #select the tuning parmaete that has the small cross-validation value
  lll <- min(which(cv == min(cv, na.rm = TRUE)));
  phi.hat.full <- phi.final[[lll]];

  return(list(brk.points = brk.points.final[[lll]], cv = cv, cv1.final = lambda.full[lll,1], cv2.final = lambda.full[lll,2]))
}


#' local sreening step (second step).
#' 
#' @description Perform the local sreening to "thin out" redundant break points. 
#' 
#' @param data input data matrix, with each column representing the time series component 
#' @param eta tuning parmaeter eta for lasso
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso 
#' @param pts the selected break points after the first step
#' @param an the neighborhood size a_n
#' @return A list oject, which contains the followings
#' \describe{
#'   \item{pts}{a set of selected break point after the second local screening step}
#'   \item{omega}{the selected Omega value}
#' }

second.step.local <- function(data, eta, q, max.iteration = 1000, tol = 10^(-4), pts, an){
  m <- length(pts); if( m == 0){break;}
  
  #compute the local loss functions for each selected break points 
  try <- break.var.local.new(data, eta, q, max.iteration = 1000, tol = tol, pts, an);
  # record the local loss function that include or exclude some break point 
  L.n.1 = try$L.n.1; L.n.2 = try$L.n.2; 
  
  #OMEGA is selected by data-driven method
  #first, compute the V value as the difference of loss functions that include and exclude some break point 
  V = rep(0, m)
  for(i in 1:m ){
    V[i] = L.n.2[i] - (L.n.1[2*i-1] +  L.n.1[2*i])
  }
  
  #add two bounary points as reference points (by assumption, their V values shoule be extremly small)
  T <- length(data[,1]);
  pts.redundant <- c(an+q, T-an)
  try.redundant <- break.var.local.new(data, eta, q, max.iteration = 1000, tol = tol,  pts.redundant, an);
  L.n.1.redundant = try.redundant$L.n.1; L.n.2.redundant = try.redundant$L.n.2;
  V.redundant <- rep(0, 2)
  for(i in 1:2 ){
    V.redundant[i] = L.n.2.redundant[i] - (L.n.1.redundant[2*i-1] +  L.n.1.redundant[2*i])
  }

  #use the maximum value of V.redundant as the reference V value
  V <- c(V,rep(max(V.redundant),floor(2*length(V))))
  if( length(unique(V)) <= 2 ){
    omega <- max(V);
  }
  if( length(unique(V)) > 2 ){
    #use kmeans to cluster the V 
    clus.2 <- kmeans(V, centers = 2); fit.2 <- clus.2$betweenss/clus.2$totss; 
    if(fit.2 < 0.20){
      omega <- max(V);
    }
    if( fit.2 >= 0.20 ){
      #if the reference point is in the subset with larger center, this means no change points: set omeage = max(V)
      #otherwise, set omega = min(V) -1 
      loc <- clus.2$cluster;
      if( clus.2$centers[1] > clus.2$centers[2]  ){
        omega <- min(V[which(loc==1)]) -1 ;
        if(loc[length(loc)] == 1){
          omega <- max(V);
        }
      }
      if( clus.2$centers[1] < clus.2$centers[2]  ){
        omega <- min(V[which(loc==2)]) -1 ;
        if(loc[length(loc)] == 2){
          omega <- max(V);
        }
      }
    }
  }
  
  #select the break points by localized information criterion (LIC)
  L.n.1.temp <- L.n.1; L.n.2.temp <- L.n.2; L.n.plot <- rep(0,m+1); L.n.plot[1] <- sum(L.n.1) + m*omega; 
  mm <- 0; ic <- 0; add.temp <- 0; pts.full <- vector("list",m+1); pts.full[[1]] <- pts; ind.pts <- rep(0,m);
  while(mm < m){
    mm <- mm + 1;
    L.n.temp <- rep(0,length(pts));
    for(i in 1:length(pts)){
      L.n.temp[i] <- sum(L.n.1.temp) - L.n.1.temp[(2*i-1)] - L.n.1.temp[(2*i)] + L.n.2.temp[i] + 1*add.temp;
    }
    ll <- min(which.min(L.n.temp)); ind.pts[mm] <- ll;
    pts <- pts[-ll]; 
    L.n.1.temp <- L.n.1.temp[-c(2*ll-1,2*ll)]; add.temp <- add.temp + 1*L.n.2.temp[ll]; 
    L.n.2.temp <- L.n.2.temp[-ll]; 
    L.n.plot[mm+1] <- L.n.temp[ll] + (m - mm)*omega;
    pts.full[[mm+1]] <- pts;
  }
  
  ind <- 0;
  ind <- min(which.min(L.n.plot))
  
  return(list(pts = pts.full[[ind]], omega = omega  ))
}


#' Compute local loss function.
#' 
#' @param data input data matrix, with each column representing the time series component 
#' @param eta tuning parmaeter eta for lasso
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso 
#' @param pts the selected break points after the first step
#' @param an the neighborhood size a_n
#' @return A list oject, which contains the followings
#' \describe{
#'   \item{L.n.1}{A vector of loss functions that include some break point}
#'   \item{L.n.2}{A vector of loss functions that exclude some break point}
#' }

break.var.local.new <- function(data, eta, q, max.iteration = 1000, tol = 10^(-4),  pts, an){
  p <- length(data[1,]); T <- length(data[,1]); m <- length(pts);
  
  #construct the local interval for computing the loss function
  bounds.1 <- vector("list",2*m); bounds.2 <- vector("list",m);
  for(i in 1:m){
    bounds.1[[(2*i-1)]] <- c(pts[i] - an, pts[i] - 1 );
    bounds.1[[(2*i)]] <- c(pts[i], pts[i] + an );
    bounds.2[[(i)]] <- c(pts[i] - an, pts[i] + an );
  }
  
  #compute the local loss function that include the given break point
  L.n.1 <- c()
  for(mm in 1:(2*m)){
    data.temp <- data[(bounds.1[[mm]][1]):(bounds.1[[mm]][2]),];
    try <- var_lasso_brk(data = data.temp, eta, q, 1000, tol = tol)
    L.n.1 <- c(L.n.1 , try$pred.error)
  }
  
  #compute the local loss function that include the given break point
  L.n.2 <- c()
  for(mm in 1:m){
    data.temp <- data[(bounds.2[[mm]][1]):(bounds.2[[mm]][2]),];
    try <- var_lasso_brk(data = data.temp, eta, q, 1000, tol = tol)
    L.n.2 <- c(L.n.2 , try$pred.error)
  }
  
  return(list(L.n.1 = L.n.1, L.n.2 = L.n.2))
}


#' exhaustive search step (third step).
#' 
#' @description Perform the exhaustive search to select the break point for each cluster. 
#' 
#' @param data input data matrix, with each column representing the time series component 
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso 
#' @param pts.list the selected break points clustered by a_n after the second step
#' @param an the neighborhood size a_n
#' @param eta tuning parmaeter eta for lasso
#' @return A list oject, which contains the followings
#' \describe{
#'   \item{pts}{a set of final selected break point after the third exhaustive search step}
#' }

third.step.exhaustive <- function(data, q, max.iteration = 1000, tol = tol,  pts.list, an, eta ){
  N <- length(data[,1]); p <- length(data[1,]);
  n <- length(pts.list);  #number of cluster
  final.pts <- rep(0,n);
  pts.list.full <- pts.list
  pts.list.full <- c(1, pts.list.full , N)
  
  #construct the interval for performing the lasso and computing the loss function
  for(i in 1:n){
    pts.temp <- pts.list.full[[i+1]];
    m <- length(pts.temp); 
    if( m <= 1  ) {
      final.pts[i] <- pts.temp;
    }
    if( m > 1  ){
      bounds.1 <- vector("list",2*m);
      lb <- max(pts.list.full[[i]]) + an;
      ub <- min(pts.list.full[[i+2]]) - an;
      for(ii in 1:m){
        bounds.1[[(2*ii-1)]] <- c(pts.temp[ii] - an , pts.temp[ii] - 1 );
        bounds.1[[(2*ii)]] <- c(pts.temp[ii], pts.temp[ii] + an );
      }

      L.n <- c()
      for(mm in 1:(2*m)){
        data.temp <- data[(bounds.1[[mm]][1]):(bounds.1[[mm]][2]),];
        try <- var_lasso_brk(data = data.temp, eta, q, 1000, tol = tol)
        L.n <- c(L.n , try$pred.error)
      }
      
      sse.full <- rep(0,m)
      for(ii in 1:m ){
        sse.full[ii] = abs(L.n[2*ii-1]  +  L.n[2*ii])
      }
      
      #select the point that has the smallest SSE among the cluster
      final.pts[i] <- pts.list.full[[i+1]][min(which(sse.full == min(sse.full)))];
    }
  }
  
  return(pts = final.pts) 
}


#' cluster the points by neighborhood size a_n
block.finder <- function(pts,an){
  nn <- length(pts);
  if( nn == 1){b <- pts;}
  if( nn > 1){
    b <- vector("list",nn);
    i.ind <- 1;
    jj <- 0;
    while (i.ind < nn) {
      ct <- 1;
      jj <- jj + 1;
      for (j in (i.ind+1):nn) {
        if( abs(pts[i.ind] - pts[j]  ) <= an   ){ct <- ct + 1;}
      }
      b[[jj]] <- pts[(i.ind):(i.ind+ct-1)];
      i.ind <- i.ind + ct;
    }
    l <- length(b[[jj]]);
    if(b[[jj]][l] != pts[nn]  ){
      jj <- jj + 1;
      b[[(jj)]] <- c(pts[nn])   
    }
    b <- b[(1):(jj)];
  }
  
  return(b = b)
}


#' Prediction function 1
pred.block <- function(Y,phi,q,T,p,h){
  concat.Y <- matrix(0,p,q+h); concat.Y[,1:q] <- Y[,(T-q+1):T];
  for ( j in 1:h){
    temp <- matrix(0,p,1);
    for (i in 1:q){temp <- temp +  phi[,((i-1)*p+1):(i*p)]%*%concat.Y[,q+j-i];}
    concat.Y[,q+j] <- temp; 
  }
  return(as.matrix(concat.Y[,(q+1):(q+h)]))
}

#' Prediction function 2
pred <- function(Y,phi,q,T,p,h){
  concat.Y <- matrix(0,p,q+h); concat.Y[,1:q] <- Y[,(T-q+1):T];
  for ( j in 1:h){
    temp <- matrix(0,p,1);
    for (i in 1:q){temp <- temp +  phi[,((i-1)*p+1):(i*p)]%*%concat.Y[,q+j-i];}
    concat.Y[,q+j] <- temp; 
  }
  return(as.matrix(concat.Y[,q+h]))
}

#' Plot the AR coefficient matrix
plot.matrix <- function (phi,p,name = NULL) {
  B <- phi
  if (nrow(B) == 1) {
    B <- matrix(B[, 1:ncol(B)], nrow = 1)
  }
  else {
    B <- B[, 1:ncol(B)]
  }
  k <- nrow(B)
  s1 <- 0
  m <- 0
  s <- 0
  s <- s + s1
  text <- c()
  for (i in 1:p) {
    text1 <- as.expression(bquote(bold(Phi)^(.(i))))
    text <- append(text, text1)
  }
  if (m > 0) {
    for (i in (p + 1):(p + s + 1)) {
      text1 <- as.expression(bquote(bold(beta)^(.(i - p - 
                                                    s1))))
      text <- append(text, text1)
    }
  }
  f <- function(m) t(m)[, nrow(m):1]
  rgb.palette <- colorRampPalette(c("blue", "white", "red"), space = "Lab")
  at <- seq(k/2 + 0.5, p * (k) + 0.5, by = k)
  if (m > 0) {
    at2 <- seq(p * k + m/2 + 0.5, p * k + s * m + 0.5, by = m)
  }
  else {
    at2 = c()
  }
  at <- c(at, at2)
  se2 = seq(1.75, by = k, length = k)
  L2 <- levelplot(as.matrix(f(B)), at =seq( -1, 1, length=101),col.regions = rgb.palette(100), 
                  colorkey = NULL, xlab = NULL, ylab = NULL, main = list(label = name, 
                                                                         cex = 1), panel = function(...) {
                                                                           panel.levelplot(...)
                                                                           panel.abline(a = NULL, b = 1, h = seq(1.5, m * s + 
                                                                                                                   p * k + 0.5, by = 1), v = seq(1.5, by = 1, length = p * 
                                                                                                                                                   k + m * s), lwd = 0.5)
                                                                           bl1 <- seq(k + 0.5, p * k + 0.5, by = k)
                                                                           b23 <- seq(p * k + 0.5, p * k + 0.5 + s * m, by = m)
                                                                           b1 <- c(bl1, b23)
                                                                           panel.abline(a = NULL, b = 1, v = p * k + 0.5, lwd = 3)
                                                                           panel.abline(a = NULL, b = 1, v = b1, lwd = 2)
                                                                         }, scales = list(x = list(alternating = 1, labels = text, 
                                                                                                   cex = 2, at = at, tck = c(0, 0)), y = list(alternating = 0, 
                                                                                                                                              tck = c(0, 0))))
  return(L2)
}







