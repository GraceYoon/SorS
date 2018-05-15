# All functions for SorS analysis and all other comparisons.
# Including functions for generating data, performance measures and CV for choosing tuning parameters.

library(sandwich)
library(BhGLM)
library(glmnet)
library(caret)
library(bindata)

GenerateData <- function(n, p, heterolevel = 1, corrtype = "corr", bset = c(2, 2, -2, 2, -2), s = 4){
  # heterolevel 1 = moderate heteroscedasticity (var = mu)
  # heterolevel 2 = severe heteroscedasticity (var = mu^2)
  # corrtype = "corr" and "uncr"
  
  trueb <- matrix(0, ncol = 1, nrow = (p + 1))
  trueb[1:(s+1)] <- bset
  
  if(corrtype == "corr"){
    rho <- 0.5
    times <- 1:p
    H <- abs(outer(times, times, "-"))
    m <- rho^H
    X <- rmvbin(n, margprob=rep(0.5, p), bincorr = m)
  } else if(corrtype == "uncr"){
    X <- matrix(rbinom(n*p,size = 1,prob = 0.5), ncol = (p), nrow = (n), byrow = T)
  }
  
  mu <- exp(cbind(1,X)%*%trueb)
  
  if( heterolevel == 0 ){
    y <- rgamma(n, shape = mu^2, rate = mu)
  } else if( heterolevel == 1 ){
    y <- rgamma(n, shape = mu, rate = 1)
  } else if( heterolevel == 2 ){
    y <- rgamma(n, shape = 1, rate = 1/mu)
  }
  
  return (list( y = y, X = X, truecoef = trueb, heterolevel = heterolevel, corrtype = corrtype))
}

brmse <- function(bhat, trueb){
  out <- sqrt(mean((bhat-trueb)^2))
  return(out)
}

rmse <- function(yhat, realy){
  out <- sqrt(mean((yhat-realy)^2))
  return(out)
}

size <- function(coeffest, tol = 1e-6){
  out <- sum(abs(coeffest)>tol)
  return(out)
}


full <- function(y,X,family="poisson"){
  ptm <- proc.time()
  full <- glm(y ~ X, family="poisson")
  bhat0 <- matrix(full$coefficients, ncol=1)
  Comptime <- (proc.time()-ptm)[3]
  return(list( fit=full, coef=bhat0, comptime=Comptime ))
}

SorS <- function(y, X, a, A=nrow(X)){
  n <- nrow(X); p <- ncol(X)
  full <- glm(y ~ X, family="poisson")
  bhat0 <- matrix(full$coefficients, ncol=1)
  cov.m1 <- vcovHC(full, type="HC0")
  std.err <- sqrt(diag(cov.m1))[-1]
  z.val <- bhat0[-1]/std.err
  ix.sel_SorS <- sort(abs(z.val), decreasing=T, index.return=T)$ix
  pU <- c() # log-posterior probability
  sel <- rep(a,(p+1)); sel[1] <- A
  for ( kk in 0:p ){
    if(kk>0){sel[ix.sel_SorS[1:kk]+1] <- A}
    U <- diag(diag(cov.m1)*sel)
    pU[kk+1] <- -0.5*n*t(bhat0)%*%solve(n*U+n*cov.m1)%*%bhat0+p*log(n)/2-0.5*log(det(2*pi*(n*U+n*cov.m1)))
  }
  gam <- which.max(pU)-1
  ixx <- sort(ix.sel_SorS[1:gam])
  fit <- glm(y ~ X[,ixx], family="poisson")
  bhat <- fit$coefficients
  return(list( fit=fit, ix=ixx, coef=bhat, a=a, A=A, post.prob=pU, rank.ix=ix.sel_SorS))
}

cv.SorS <- function(y,X, nfolds=10,
                    a.cand=NULL){
  xx <- createFolds(y,k=nfolds, list=TRUE, returnTrain = TRUE)
  ptm <- proc.time()
  
  if( is.null(a.cand)==TRUE ){ a.cand <- c(0.001, 0.005, seq(0.01,0.09,by=0.01),seq(0.1,2,by=0.05))}
  # a.cand <- c(0.01, seq(0.05, 2, by=0.05)) } #length=41
  # a.cand <- c(0.001, 0.005, seq(0.01,0.09,by=0.01),seq(0.1,2,by=0.05)) } # length=50
  err <- c()
  for ( j in 1:length(a.cand)){
    
    err.folds <- c()
    for (fold in 1:nfolds){
      x.train <- X[xx[[fold]],]
      x.test <- X[-xx[[fold]],]
      y.train <- y[xx[[fold]]]
      y.test <- y[-xx[[fold]]]
      
      result <- SorS(y.train, x.train, a=a.cand[j])
      ixx <- result$ix
      yhat <- exp(cbind(1,x.test[,ixx])%*%result$coef)
      err.folds[fold] <- sqrt(mean((y.test-yhat)^2))
    }
    err[j] <- mean(err.folds)
    cat("Done with a=", a.cand[j], "\n")
  }
  
  #select.a <- a.cand[which.min(err)]
  select.a <- a.cand[which(err == min(err))]
  
  result.entire <- SorS(y,X,a=select.a[1])
  
  Comptime_SorS <- (proc.time()-ptm)[3]
  
  ixx <- result.entire$ix
  yhat <- exp(cbind(1,X[,ixx])%*%result.entire$coef)
  names(err) <- a.cand
  return(list( fit=result.entire, ix=ixx, coef=result.entire$coef, comptime=Comptime_SorS, cv.err=err, a.cand=a.cand))
}


cv.bmlasso <- function(y,X, nfolds=10,
                       s0.cand=NULL){
  xx <- createFolds(y,k=nfolds, list=TRUE, returnTrain = TRUE)
  ptm <- proc.time()
  if( is.null(s0.cand)==TRUE ){ s0.cand <- seq(0.005, 0.2, by = 0.005) }
  err <- c()
  for ( j in 1:length(s0.cand)){
    
    err.folds <- c()
    for (fold in 1:nfolds){
      x.train <- X[xx[[fold]],]
      x.test <- X[-xx[[fold]],]
      y.train <- y[xx[[fold]]]
      y.test <- y[-xx[[fold]]]
      
      result <- bmlasso(x=x.train, y=y.train, family="poisson", prior="mde", ss=c(s0.cand[j], 1))
      bhat.bmlasso <- as.vector(result$coefficients)
      ixx <- which(bhat.bmlasso[-1]!=0)
      
      yhat <- exp(cbind(1,x.test)%*%bhat.bmlasso)
      err.folds[fold] <- sqrt(mean((y.test-yhat)^2))
    }
    err[j] <- mean(err.folds)
    cat("Done with ", s0.cand[j], "\n")
  }
  
  #select.s0 <- s0.cand[which.min(err)]
  select.s0 <- s0.cand[which(err == min(err))]
  
  result.entire <- bmlasso(x=X, y=y, family="poisson", prior="mde", ss=c(select.s0[1], 1))
  
  Comptime_bmlasso <- (proc.time()-ptm)[3]
  bhat <- as.vector(result.entire$coefficients)
  ixx <- which(bhat[-1]!=0)
  yhat <- exp(cbind(1,X)%*%bhat)
  names(err) <- s0.cand
  return(list( fit=result.entire, ix=ixx, coef=bhat, comptime=Comptime_bmlasso, cv.err=err, s0.cand=s0.cand)) 
}



cv.lasso <- function(y,X,alpha=1,nfolds=10){
  xx <- createFolds(y,k=nfolds, list=TRUE, returnTrain = TRUE)
  ptm <- proc.time()
  lamb <- cv.glmnet(X, y, family="poisson", alpha=1, nfolds=nfolds, type.measure="mse")$lambda
  
  err <- c()
  
  err.folds <- matrix(NA, nrow=nfolds, ncol=length(lamb))
  for (fold in 1:nfolds){
    x.train <- X[xx[[fold]],]
    x.test <- X[-xx[[fold]],]
    y.train <- y[xx[[fold]]]
    y.test <- y[-xx[[fold]]]
    
    result <- glmnet(x=x.train, y=y.train, family="poisson", alpha=1, lambda=lamb)
    bhat <- rbind(result$a0,result$beta)
    pred <- function(coef){
      coef <- matrix(coef,ncol=1)
      exp(cbind(1,x.test)%*%coef)
    }
    rmse.lamb <- function(yhat){
      sqrt(mean((y.test-yhat)^2))
    }
    yhat <- apply(bhat,MARGIN = 2,pred)
    err.folds[fold,] <- apply(yhat,MARGIN=2,rmse.lamb)
    cat("Done with fold: ", fold , "\n")
  }
  err <- colMeans(err.folds)
  
  select.lamb <- lamb[which(err == min(err))]
  
  result.entire <- glmnet(x=X, y=y, family="poisson", alpha=1)
  bhat.lasso <- coef(result.entire)[, which(result.entire$lambda == select.lamb)]
  ixx <- which(bhat.lasso[-1]!=0)
  
  if(length(ixx) > 0){
    fit_lasso <- glm(y ~ X[,ixx], family="poisson")
  } else if(length(ixx) == 0){
    fit_lasso <- glm(y ~ 1, family="poisson")
    ixx <- 0
  }
  
  Comptime_lasso <- (proc.time()-ptm)[3]
  yhat <- fit_lasso$fitted.values
  
  return(list( ix=ixx, coef=fit_lasso$coefficients))
}

pred <- function(coef){
  coef <- matrix(coef,ncol=1)
  exp(cbind(1,x.test)%*%coef)
}
rmse.lamb <- function(yhat){
  sqrt(mean((y.test-yhat)^2))
}



extcnt <- function(model, trueb){
  crtcnt <- 0
  for(i in 1:length(model)){
    if( length(model[[i]])==sum(trueb[-1]!=0) & sum(model[[i]] %in% which(trueb[-1]!=0))==sum(trueb[-1]!=0) ){
      crtcnt <- crtcnt+1
    }
  }
  return(crtcnt/length(model))
}
cov.prob <- function(model, trueb){
  coverage <- 0
  for(i in 1:length(model)){
    if( sum( model[[i]] %in% which(trueb[-1]!=0) ) == sum(trueb[-1]!=0) ){
      coverage <- coverage+1
    }
  }
  return(coverage/length(model))
}
corzeros <- function(model, trueb){
  correctzeros <- c()
  for(i in 1:length(model)){
    correctzeros[i] <- sum( setdiff(1:p,model[[i]]) %in% which(trueb[-1]==0) )/sum(trueb[-1]==0)
  }
  return(mean(correctzeros))
}

incorzeros <- function(model, trueb){
  incorrectzeros <- c()
  for(i in 1:length(model)){
    incorrectzeros[i] <- sum( setdiff(1:p,model[[i]]) %in% which(trueb[-1]!=0) )/sum(trueb[-1]!=0)
  }
  return(mean(incorrectzeros))
}


avg.accu <- function(model, trueb){
  avg.accu <- c()
  correct <- as.numeric(trueb[-1]!=0)
  p <- length(correct)
  for ( i in 1:length(model)){
    binary <- rep(0,p)
    binary[model[[i]]] <- 1
    avg.accu[i] <- mean(correct==binary)
  }
  return(mean(avg.accu))
}
