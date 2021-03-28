library(ranger)
library(Matrix)
# library(Rcpp)
# library(RcppArmadillo)

### data split(2) random forest
rf.2split <- function(X,y,num.trees=200,mtry=NULL,max.depth=0,min.node.size=5,MSE.thol=1e6,forest.save=F) {
  n <- nrow(X); p <- ncol(X)
  if (is.null(mtry)) mtry <- round(p/3)

  Data <- data.frame(cbind(y, X))
  names(Data) <- c("y", paste("X", 1:p, sep = ""))
  ### grid search
  params.grid <- expand.grid(
    num.trees = num.trees,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )

  # split the data into two parts A1 and A2
  # use A2 to build the random forest and use
  # A1 to do inference
  ind.sep <- cut(1:n, breaks = 2, labels = FALSE)
  A1.ind <- which(ind.sep==1, arr.ind = TRUE)
  A2.ind <- which(ind.sep==2, arr.ind = TRUE)
  Data.A1 <- Data[A1.ind, ]
  Data.A2 <- Data[A2.ind, ]

  forest.A2 <- NULL;
  MSE.oob.A2 <- MSE.thol
  params.A2 <- NULL
  ### use oob error to do hyper-parameter tuning
  for (i in 1:nrow(params.grid)) {
    temp.A2 <- ranger(y~., data = Data.A2,
                      num.trees=params.grid$num.trees[i],
                      mtry=params.grid$mtry[i],
                      max.depth = params.grid$max.depth[i],
                      min.node.size = params.grid$min.node.size[i]
    )
    if (temp.A2$prediction.error <= MSE.oob.A2) {
      forest.A2 <- temp.A2
      params.A2 <- params.grid[i,]
      MSE.oob.A2 <- temp.A2$prediction.error
    }
  }
  # nodes information of A1 on the random forest built on A2
  # nodes is a matrix of n.A1 by num.trees
  # the information of random forest built on A2 has been used here
  nodes.A1 <- predict(forest.A2, data = Data.A1, type = "terminalNodes")$predictions

  returnList <- list(forest.A2 = forest.A2,
                     params.A2 = params.A2,
                     A1.ind = A1.ind,
                     A2.ind = A2.ind,
                     nodes.A1 = nodes.A1,
                     MSE.oob.A2 = MSE.oob.A2)
  if (!forest.save) returnList <- returnList[-1]
  returnList
}


### compute the weight matrix for data split(2) random forest
### nodes.A1 is the n.A1 by n.A1 nodes matrix by rf.2split
weight.2split <- function(nodes.A1) {
  n.A1 <- nrow(nodes.A1); num.trees <- ncol(nodes.A1)
  w.list <- rep(list(NA), num.trees)
  for (j in 1:num.trees) {
    w.mat <- matrix(0, n.A1, n.A1)
    for (i in 1:n.A1) {
      ind <- nodes.A1[,j]==nodes.A1[i, j]
      ind[i] <- FALSE # to remove self-prediction
      w <- 1/sum(ind)
      w.vec <- ifelse(ind,yes=w,no=0)
      w.mat[i, ] <- w.vec/num.trees
    }
    w.list[[j]] <- Matrix(w.mat, sparse = TRUE)
  }
  Reduce("+", w.list)
}


### k split cross-fitting random forest
rf.kcross <- function(X,y,k=2,num.trees=200,mtry=NULL,max.depth=0,min.node.size=5,MSE.thol=1e6,forest.save=F) {
  n <- nrow(X);p<-ncol(X)
  if (is.null(mtry)) mtry <- round(p/3)

  Data <- data.frame(cbind(y, X))
  names(Data) <- c("y", paste("X", 1:p, sep = ""))
  ### grid search
  params.grid <- expand.grid(
    num.trees = num.trees,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )

  ind.sep <- cut(1:n, breaks = k, labels = FALSE)
  k.ind <- rep(list(NA), k)
  for (j in 1:k) {
    k.ind[[j]] <- which(ind.sep==j)
  }

  forest <- params <- rep(list(NA), k)
  MSE.oob <- rep(MSE.thol, k)
  for (j in 1:k) {
    for (i in 1:nrow(params.grid)) {
      temp.k <- ranger(y~., data = Data[-k.ind[[j]],],
                       num.trees=params.grid$num.trees[i],
                       mtry=params.grid$mtry[i],
                       max.depth = params.grid$max.depth[i],
                       min.node.size = params.grid$min.node.size[i]
      )
      if (temp.k$prediction.error < MSE.oob[j]) {
        forest[[j]] <- temp.k
        params[[j]] <- params.grid[i,]
        MSE.oob[j] <- temp.k$prediction.error
      }
    }
  }

  ### nodes contains k matrices
  nodes <- list(rep(NA), k)
  for (j in 1:k) {
    nodes[[j]] <- predict(forest[[j]],data=Data[k.ind[[j]],],type="terminalNodes")$predictions
  }

  returnList <- list(forest = forest,
                     params = params,
                     k.ind = k.ind,
                     nodes = nodes,
                     MSE.oob = MSE.oob)
  if (!forest.save) returnList <- returnList[-1]
  returnList
}


### compute the weight matrix for k split cross-fitting random forest
### nodes is the list of nodes matrices by rf.kcross
weight.kcross <- function(nodes) {
  w.list <- lapply(nodes, weight.2split)
  w.out <- bdiag(w.list)
}



### random forest using full data
### directly use the weight.2split() for weight matrix computation
rf.full <- function(X,y,num.trees=200,mtry=NULL,max.depth=0,min.node.size=5,MSE.thol=1e6,forest.save=F) {
  n <- nrow(X);p<-ncol(X)
  if (is.null(mtry)) mtry <- round(p/3)

  Data <- data.frame(cbind(y, X))
  names(Data) <- c("y", paste("X", 1:p, sep = ""))
  ### grid search
  params.grid <- expand.grid(
    num.trees = num.trees,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )

  forest <- params <- NA
  MSE.oob <- MSE.thol
  for (i in 1:nrow(params.grid)) {
    temp<- ranger(y~., data = Data,
                  num.trees=params.grid$num.trees[i],
                  mtry=params.grid$mtry[i],
                  max.depth = params.grid$max.depth[i],
                  min.node.size = params.grid$min.node.size[i]
    )
    if (temp$prediction.error < MSE.oob) {
      forest <- temp
      params <- params.grid[i,]
      MSE.oob <- temp$prediction.error
    }
  }

  predicted.values <- predict(forest,data = Data, type="response")$predictions
  nodes <- predict(forest,data = Data, type="terminalNodes")$predictions

  returnList = list(forest = forest,
                    params = params,
                    predicted.values = predicted.values,
                    nodes = nodes,
                    MSE.oob = MSE.oob)
  if (!forest.save) returnList <- returnList[-1]
  returnList
}


### helper function for check the invertability of t(W.rep)%*%W.rep
### if not invertible, add ridge type of penalty in the stats.rf() function
try.inverse <- function(m) class(try(solve(m),silent=T))[1]=="matrix"


### calculate the standard error and conduct IV strength test
statRF.2split <- function(Y.rep, D.rep, Y, D, VW.rep, VW, betaHat, weight, n, SigmaSqY, SigmaSqD, L=500, boot=TRUE, c0=0.01, C1=1, lam=0.05) {
  n.A1 <- length(D.rep)
  p <- ncol(VW.rep)
  ### check the inverse of t(VW.rep)%*%VW.rep
  try.mat <- t(VW.rep)%*%VW.rep
  if (try.inverse(try.mat)) {
    invert <- FALSE
    Portho.VW.rep <- diag(1,n.A1,n.A1) - VW.rep %*% solve(t(VW.rep)%*%VW.rep) %*% t(VW.rep)
    TRF.V <- t(as.matrix(weight))%*%Portho.VW.rep%*%as.matrix(weight)
  } else {
    invert <- TRUE
    ### ridge penalty
    Portho.VW.rep <- diag(1,n.A1,n.A1) - VW.rep %*% solve(t(VW.rep)%*%VW.rep+lam*diag(1,p,p)) %*% t(VW.rep)
    TRF.V <- t(as.matrix(weight))%*%Portho.VW.rep%*%as.matrix(weight)
  }

  ### the standard error of the estimator and iv strength test
  tr.TRF <- sum(diag(TRF.V))
  tr.TRF2 <- sum(diag(TRF.V%*%TRF.V))
  LHS <- t(D)%*%TRF.V%*%D
  Sd <- sqrt(SigmaSqY*(t(D)%*%TRF.V%*%TRF.V%*%D)/(LHS^2))


  ### standard errors of  bias-corrected estimator
  ### the normal method, assuming gaussian
  try.mat <- t(VW)%*%VW
  if (try.inverse(try.mat)) {
    Portho.VW <- diag(1,n.A1,n.A1) - VW %*% solve(t(VW)%*%VW) %*% t(VW)
  } else {
    ### ridge penalty
    Portho.VW <- diag(1,n.A1,n.A1) - VW %*% solve(t(VW)%*%VW+lam*diag(1,p,p)) %*% t(VW)
  }
  SigmaYD <- t(D-D.rep)%*%Portho.VW%*%(Y-D*betaHat)/(n.A1-p)
  beta.cor <- betaHat - SigmaYD*tr.TRF/LHS
  Sd.cor <- sqrt((SigmaSqY*(t(D)%*%TRF.V%*%TRF.V%*%D)+(SigmaSqD*SigmaSqY+SigmaYD^2)*(tr.TRF^2)/(n.A1-p))/(LHS^2))

  ### using bootstrap
  # SigmaHat <- matrix(c(SigmaSqY,SigmaYD,SigmaYD,SigmaSqD),2,2)
  # T.boot <- rep(NA,L)
  # Sd.cor.boot <- NA
  # for (l in 1:L) {
  #   Error.boot <- mvrnorm(n.A1, rep(0,2), SigmaHat)
  #   T.boot[l] <- t(Error.boot[,1])%*%TRF.V%*%D - (t(Error.boot[,2])%*%Portho.VW%*%Error.boot[,1]*tr.TRF)/(n.A1-p)
  # }
  # Sd.cor.boot <- sqrt(var(T.boot)/(LHS^2))


  ### the IV strength test
  tau.n <- log(log(n))
  # cn <- sqrt(tau.n/tr.TRF)/tau.n*(2+sqrt(tr.TRF2/tr.TRF)/tau.n)+1/tau.n^2
  cn <- 0
  Cn.V <- sqrt(tr.TRF2) + 2*sqrt(tr.TRF)*max(tau.n, sqrt(LHS/((1-cn)*SigmaSqD*tr.TRF)))
  RHS <- (1+c0)*tr.TRF*SigmaSqD + sqrt(tau.n)*SigmaSqD*Cn.V


  ### the Signal Strength Test for splitting estimator
  SignalTest.LHS <- sqrt(SigmaSqY*t(D)%*%TRF.V%*%TRF.V%*%D)
  SignalTest.RHS <- C1*(SigmaYD*tr.TRF+sqrt(tr.TRF2)*sqrt(tau.n))

  ### the Signal Strength Test for bias-corrected estimator
  SignalTest.cor.LHS <- SignalTest.LHS
  SignalTest.cor.RHS <- C1*sqrt(tr.TRF2)*sqrt(tau.n) + tau.n


  ### calculate the decomposition
  term1 <- (t(Error[,2][forest.2$A1.ind])%*%TRF.V%*%D)/LHS
  term2 <- -(SigmaYD-0.5)*tr.TRF/LHS
  term3 <- (t(Error[,2][forest.2$A1.ind])%*%TRF.V%*%Error[,1][forest.2$A1.ind]-0.5*tr.TRF)/LHS


  returnList <- list(Sd = Sd,
                     Sd.cor = 1.1*Sd.cor,
                     SigmaYD = SigmaYD,
                     beta.cor = beta.cor,
                     # Sd.cor.boot = 1.1*Sd.cor.boot,
                     tr.TRF = tr.TRF,
                     tr.TRF2 = tr.TRF2,
                     LHS = LHS,
                     RHS = RHS,
                     SignalTest.LHS = SignalTest.LHS,
                     SignalTest.RHS = SignalTest.RHS,
                     SignalTest.cor.LHS = SignalTest.cor.LHS,
                     SignalTest.cor.RHS = SignalTest.cor.RHS,
                     invert = invert,
                     term1 = term1,
                     term2 = term2,
                     term3 = term3)
  returnList
}


### the same function for kcross
statRF.kcross <- function(D.rep, D, VW.rep, weight, SigmaSqY, SigmaSqD, c0=0.01, lam=0.05) {
  n <- length(D.rep)
  p <- ncol(VW.rep)
  ### check the inverse of t(VW.rep)%*%VW.rep
  try.mat <- t(VW.rep)%*%VW.rep
  if (try.inverse(try.mat)) {
    invert <- FALSE
    Portho.VW <- diag(1,n,n) - VW.rep %*% solve(t(VW.rep)%*%VW.rep) %*% t(VW.rep)
    TRF.V <- t(as.matrix(weight))%*%Portho.VW%*%as.matrix(weight)
  } else {
    invert <- TRUE
    ### ridge penalty
    Portho.VW <- diag(1,n,n) - VW.rep %*% solve(t(VW.rep)%*%VW.rep+lam*diag(1,p,p)) %*% t(VW.rep)
    TRF.V <- t(as.matrix(weight))%*%Portho.VW%*%as.matrix(weight)
  }

  ### the standard error of the estimator and iv strength test
  LHS <- t(D)%*%TRF.V%*%D
  Sd <- sqrt(SigmaSqY*(t(D)%*%TRF.V%*%TRF.V%*%D)/(LHS)^2)

  ### the IV strength test
  # tr.TRF <- sum(diag(TRF.V))
  # taun <- log(n)
  # cn <- sqrt(2*log(n)/tr.TRF)/taun*(2+1/taun)+1/(taun)^2
  # Cn.V <- sqrt(tr.TRF)*sqrt(log(n))*SigmaSqD*(1+max(2*taun,2*sqrt(LHS/((1-cn)*tr.TRF*SigmaSqD))))
  # RHS <- (1+c0)*tr.TRF*SigmaSqD + Cn.V

  returnList <- list(Sd = Sd,
                     # LHS = LHS,
                     # RHS = RHS,
                     invert = invert)
  returnList
}


### sd function for TSLS random forest
### D.rep and VW contain only samples from A1
statRF.TSLS <- function(D.rep, VW, SigmaSqY, lam=0.05) {
  n.A1 <- length(D.rep)
  p <- ncol(VW)
  ### check the inverse of t(VW.rep)%*%VW.rep
  try.mat <- t(VW)%*%VW
  if (try.inverse(try.mat)) {
    invert <- FALSE
    Portho.VW <- diag(1,n.A1,n.A1) - VW %*% solve(t(VW)%*%VW) %*% t(VW)
  } else {
    invert <- TRUE
    ### ridge penalty
    Portho.VW <- diag(1,n.A1,n.A1) - VW %*% solve(t(VW)%*%VW+lam*diag(1,p,p)) %*% t(VW)
  }
  Sd <- sqrt(SigmaSqY/(t(D.rep)%*%Portho.VW%*%D.rep))

  returnList <- list(Sd = Sd,
                     invert = invert)
  returnList
}


### sd function for random forest with full data
statRF.full <- function(D.rep, D, VW.rep, weight, SigmaSqY, lam=0.05) {
  n <- length(D.rep)
  p <- ncol(VW.rep)
  ### check the inverse of t(VW.rep)%*%VW.rep
  try.mat <- t(VW.rep)%*%VW.rep
  if (try.inverse(try.mat)) {
    invert <- FALSE
    Portho.VW <- diag(1,n,n) - VW.rep %*% solve(t(VW.rep)%*%VW.rep) %*% t(VW.rep)
    TRF.V <- t(as.matrix(weight))%*%Portho.VW%*%as.matrix(weight)
  } else {
    invert <- TRUE
    ### ridge penalty
    Portho.VW <- diag(1,n,n) - VW.rep %*% solve(t(VW.rep)%*%VW.rep+lam*diag(1,p,p)) %*% t(VW.rep)
    TRF.V <- t(as.matrix(weight))%*%Portho.VW%*%as.matrix(weight)
  }

  ### the standard error of the estimator and iv strength test
  LHS <- t(D)%*%TRF.V%*%D
  Sd <- sqrt(SigmaSqY*(t(D)%*%TRF.V%*%TRF.V%*%D)/(LHS)^2)

  returnList <- list(Sd = Sd,
                     invert = invert)
  returnList
}
