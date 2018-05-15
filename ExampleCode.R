# Example code for Case 1 with correlated predictors case.

rm(list=ls())

source("helpers.R")
n <- 1000; p <- 50
set.seed(1)

dat <- GenerateData(n, p, heterolevel = 1, corrtype = "corr")
testdat <- GenerateData(n, p, heterolevel = 1, corrtype = "corr")

Full <- LASSO <- SSLASSO <- mySorS <- list()

  fit_SorS <- cv.SorS(y=dat$y, X=dat$X)
  bhat_SorS <- fit_SorS$coef
  ixx_SorS <- fit_SorS$ix
  yhat_SorS <- exp(cbind(1,testdat$X[,ixx_SorS])%*%bhat_SorS)
  mySorS$rmse <- sqrt(mean((testdat$y-yhat_SorS)^2))
  mySorS$coef <- rep(0, p+1); mySorS$coef[c(1, 1+ixx_SorS)] <- bhat_SorS
  
  fit_lasso <- cv.lasso(y=dat$y, X=dat$X)
  bhat_lasso <- fit_lasso$coef
  ixx_lasso <- fit_lasso$ix
  yhat_lasso <- exp(cbind(1,testdat$X[,ixx_lasso])%*%bhat_lasso)
  LASSO$rmse <- sqrt(mean((testdat$y-yhat_lasso)^2))
  LASSO$coef <- rep(0, p+1); LASSO$coef[c(1, 1+ixx_lasso)] <- bhat_lasso

  fit_bmlasso <- cv.bmlasso(dat$y, dat$X)
  bhat_bmlasso <- fit_bmlasso$coef
  ixx_bmlasso <- fit_bmlasso$ix
  yhat_bmlasso <- exp(cbind(1,testdat$X)%*%bhat_bmlasso)
  SSLASSO$rmse <- sqrt(mean((testdat$y-yhat_bmlasso)^2))
  SSLASSO$coef <- bhat_bmlasso
  
  fit_full <- glm(dat$y~dat$X, family="poisson")
  bhat_full <- matrix(fit_full$coefficients, ncol=1)
  yhat_full <- exp(cbind(1,testdat$X)%*%bhat_full)
  Full$rmse <- sqrt(mean((testdat$y-yhat_full)^2))
  Full$coef <- bhat_full


  output <- list( Full = Full,
                  LASSO = LASSO,
                  SSLASSO = SSLASSO,
                  SorS = mySorS
                  )


