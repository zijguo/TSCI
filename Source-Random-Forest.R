### Source-Random-Forest.R
### Functions: Two Stage Curvature Identification method using Random Forest
###            Estimate the treatment effect using the instrumental variable
###            approach with violation space selection


### required packages
library(ranger)
library(Matrix)
### for weight matrix computation in C++, not in use yet
# library(Rcpp)
# library(RcppArmadillo)


### TSRF.fit
### Function: Model fitting function for TWo Stage Random Forest
###           using data splitting approach
###           Use out-of-bag error for hyper-parameter tuning
### Input: X: continuous or binary, n by p_x covariates
###        Y: continuous or binary, n by 1 outcome vector
###        num.trees: integer, the number of trees in random forest
###        mtry: integer, the number of covariates to split at each node
###        max.depth: integer, the maximal depth of each tree, 0 refers to unlimited depth
###        min.node.size: integer, the minimal size(# samples in it) of each leaf node
###        MSE.thol: numeric, a large value of MSE, used for the start of hyper-parameter selection
###        forest.save: logic, to save the random forest object or not, default by FALSE to save memory
### Output: forest.A2: random forest built on subsample A2, available if forest.save=TRUE
###         params.A2: best hyper-parameters selected by out-of-bag error
###         A1.ind: index of data in subsample A1
###         A2.ind: index of data in subsample A2
###         nodes.A1: a n_A1 by num.trees matrix containing the leaf nodes' indices of subsample A1
###                   in each tree of forest.A2
###         MSE.oob: the minimal out-of-bag error using the best hyper-parameters
TSRF.fit <- function(X,Y,num.trees=200,mtry=NULL,max.depth=0,min.node.size=5,MSE.thol=1e6,forest.save=F) {
  n <- nrow(X); p <- ncol(X)
  if (is.null(mtry)) mtry <- round(p/3)

  Data <- data.frame(cbind(Y, X))
  names(Data) <- c("Y", paste("X", 1:p, sep = ""))
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
    temp.A2 <- ranger(Y~., data = Data.A2,
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



### TSRF.weight
### Function: Compute the weight matrix of data splitting random forest and random forest with full data
### Input: nodes: a n by num.trees matrix containing the nodes information, the (ith, jth) entry refers
###               to the leaf node index of the jth tree where the ith sample falls into
### Output: w.mat: the n_A1 by n_A1 symmetric sparse weight matrix(class dgCMatrix), the ith row represents
###                the weights of each sample on the prediction of the ith outcome
TSRF.weight <- function(nodes) {
  n.A1 <- nrow(nodes); num.trees <- ncol(nodes)
  w.list <- rep(list(NA), num.trees)
  for (j in 1:num.trees) {
    w.mat <- matrix(0, n.A1, n.A1)
    for (i in 1:n.A1) {
      ind <- nodes[,j]==nodes[i, j]
      ind[i] <- FALSE # to remove self-prediction
      w <- 1/sum(ind)
      w.vec <- ifelse(ind,yes=w,no=0)
      w.mat[i, ] <- w.vec/num.trees
    }
    w.list[[j]] <- Matrix(w.mat, sparse = TRUE)
  }
  w.mat <- Reduce("+", w.list)
  return(w.mat)
}


### TSRF.crossfit
### Function: Model fitting function for TWo Stage Random Forest
###           using cross fitting approach
###           Use out-of-bag error for hyper-parameter tuning
### Input: X: continuous or binary, n by p_x covariates
###        Y: continuous or binary, n by 1 outcome vector
###        k: integer, number of subsamples for cross fitting
###        num.trees: integer, the number of trees in random forest
###        mtry: integer, the number of covariates to split at each node
###        max.depth: integer, the maximal depth of each tree, 0 refers to unlimited depth
###        min.node.size: integer, the minimal size(# samples in it) of each leaf node
###        MSE.thol: numeric, a large value of MSE, used for the start of hyper-parameter selection
###        forest.save: logic, to save the random forest object or not, default by FALSE to save memory
### Output: forest: a list of k random forest objects, available if forest.save=TRUE
###         params: a list of k sets of best hyper-parameters selected by out-of-bag error
###         k.ind: a list of vectors indicating the index of the each sample in k subsamples
###         nodes: a list of nodes information of the k subsamples, similar to nodes.A1 in TSRF.fit
###         MSE.oob: the minimal out-of-bag error using the best hyper-parameters
TSRF.crossfit <- function(X,Y,k=2,num.trees=200,mtry=NULL,max.depth=0,min.node.size=5,MSE.thol=1e6,forest.save=F) {
  n <- nrow(X);p<-ncol(X)
  if (is.null(mtry)) mtry <- round(p/3)

  Data <- data.frame(cbind(Y, X))
  names(Data) <- c("Y", paste("X", 1:p, sep = ""))
  ### grid search
  params.grid <- expand.grid(
    num.trees = num.trees,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )

  # split the data into k subsamples
  # use cross-ftting procedures
  ind.sep <- cut(1:n, breaks = k, labels = FALSE)
  k.ind <- rep(list(NA), k)
  for (j in 1:k) {
    k.ind[[j]] <- which(ind.sep==j)
  }

  forest <- params <- rep(list(NA), k)
  ### use oob error to do hyper-parameter tuning
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


### TSRF.weight.crossfit
### Function: Compute the weight matrix of cross-fitting random forest, this function
###           inherits the TSRF.weight() function and augment it to cross-fitting
### Input: nodes: a list of nodes information, typically from the output of TSRF.crossfit()
### Output: w.mat: a n by n symmetric sparse weight matrix(class dgCMatrix), the ith row
###         represents the weights of each sample on the prediction of the ith outcome
TSRF.weight.crossfit <- function(nodes) {
  w.list <- lapply(nodes, weight.2split)
  w.mat <- bdiag(w.list)
  return(w.mat)
}



### TSRF.full
### Function: Model fitting function for TWo Stage Random Forest using
###           full data. Use out-of-bag error for hyper-parameter tuning
### Input: X: continuous or binary, n by p_x covariates
###        Y: continuous or binary, n by 1 outcome vector
###        k: integer, number of subsamples for cross fitting
###        num.trees: integer, the number of trees in random forest
###        mtry: integer, the number of covariates to split at each node
###        max.depth: integer, the maximal depth of each tree, 0 refers to unlimited depth
###        min.node.size: integer, the minimal size(# samples in it) of each leaf node
###        MSE.thol: numeric, a large value of MSE, used for the start of hyper-parameter selection
###        forest.save: logic, to save the random forest object or not, default by FALSE to save memory
### Output: forest: random forest object using full data, available if forest.save=TRUE
###         params: a list of best hyper-parameters selected by the out-of-bag error
###         predicted.values: the predicted values of outcome Y using full data
###         nodes: a n by num.trees nodes information matrix, similar to nodes.A1 in TSRF.fit
###         MSE.oob: the minimal out-of-bag error using the best hyper-parameters
TSRF.full <- function(X,Y,num.trees=200,mtry=NULL,max.depth=0,min.node.size=5,MSE.thol=1e6,forest.save=F) {
  n <- nrow(X);p<-ncol(X)
  if (is.null(mtry)) mtry <- round(p/3)

  Data <- data.frame(cbind(Y, X))
  names(Data) <- c("Y", paste("X", 1:p, sep = ""))
  ### grid search
  params.grid <- expand.grid(
    num.trees = num.trees,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )

  # use oob error to do hyper-parameter tuning
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

  # conduct prediction
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


### try.inverse
### Function: a helper function to identify the singularity of a matrix
### Input: m: a square matrix
### Output: a logical value whether the input matrix m is invertible
try.inverse <- function(m) class(try(solve(m),silent=T))[1]=="matrix"


### TSRF.stat
### Function: Compute the basic statistics for Two Stage Random Forest using the
###           data splitting approach, including the data splitting estimator
###           and its standard deviation, the bias corrected estimator and its
###           standard deviation, the IV strength test and the signal strength test
### Input: Y.A1: continuous, a n.A1 by 1 vector denoting the outcome in subsample A1
###        D: continuous or binary, a n.A1 by 1 vector denoting the treatment in subsample A1
###        VW.A1: continuous or binary, the instruments-covariates matrix expanded by
###               corresponding basis in subsample A1, first column is intercept
###        betaHat: the estimated treatment effect by data splitting estimator
###        weight: the n.A1 by n.A1 weight matrix, usually a sparse matrix with class dgCMatrix
###        n: full sample size
###        SigmaSqY: estimated noise level of outcome model
###        SigmaSqD: estimated noise level of treatment model
###        c0: positive, the inflation of the threshold in IV strength test, default by 0.01
###        C1: positive, the inflation of the threshold in signal strength test, default by 1.25
###        tau.n: the tuning parameter in IV strength test and signal strength test, default by log(log(n))
###        lam: the tuning parameter of ridge regression if t(VW.rep)%*%VW.rep is singular, default by 0.05
### Output: Sd: the estimated standard deviation of data splitting estimator
###         betaHat.cor: the bias corrected estimator
###         Sd.cor: the estimated standard deviation of bias corrected estimator
###         SigmaYD: the estimated covariance between \sigma ans \delta
###         T.V: the T^{RF}(V) defined by V, saved for violation space selection
###         trace.T: the trace of T^{RF}(V)
###         trace.T2: the trace of T^{RF}(V)%*%T^{RF}(V)
###         iv.str: the left hand side of IV strength test
###         iv.thol: the right hand side of IV strength test
###         DT.sq: t(D.A1)%*%T.V%*%T.V%*%D.A1, saved for violation space selection
###         signal.str: the left hand side of signal strength test for data splitting estimator
###         signal.str: the right hand side of signal strength test for data splitting estimator
###         signal.str: the left hand side of signal strength test for bias corrected estimator
###         signal.str: the right hand side of signal strength test for bias corrected estimator
###         Singularity: logical, whether t(VW.rep)%*%VW.rep is singular
TSRF.stat <- function(Y.A1, D.A1, VW.A1, betaHat, weight, n, SigmaSqY, SigmaSqD, c0=0.01, C1=2, tau.n=NULL, lam=0.05) {

  ### Constants
  if (is.null(tau.n)) {tau.n <- log(log(n))}
  n.A1 <- length(A1.ind); r.VW <- ncol(VW.A1) # the rank of (V, W)
  ### Compute the representations
  Y.rep <- as.matrix(weight)%*%Y.A1; D.rep <- as.matrix(weight)%*%D.A1
  VW.rep <- as.matrix(weight)%*%VW.A1


  ### check the inverse of t(VW.rep)%*%VW.rep
  try.mat <- t(VW.rep)%*%VW.rep
  if (try.inverse(try.mat)) {
    sigularity <- FALSE
    Portho.VW.rep <- diag(1,n.A1,n.A1) - VW.rep %*% solve(try.mat) %*% t(VW.rep)
    T.V <- t(as.matrix(weight))%*%Portho.VW.rep%*%as.matrix(weight)
  } else {
    sigularity <- TRUE
    ### ridge penalty
    Portho.VW.rep <- diag(1,n.A1,n.A1) - VW.rep %*% solve(try.mat+lam*diag(1,r.VW,r.VW)) %*% t(VW.rep)
    T.V <- t(as.matrix(weight))%*%Portho.VW.rep%*%as.matrix(weight)
  }

  ### the standard error of the estimator and iv strength test
  trace.T <- sum(diag(T.V))
  trace.T2 <- sum(diag(T.V%*%T.V))
  iv.str <- t(D.A1)%*%T.V%*%D.A1
  # this is the numerator of the variance of betaHat
  DT.Sq <- t(D.A1)%*%T.V%*%T.V%*%D.A1
  Sd <- sqrt(SigmaSqY*DT.Sq/(iv.str^2))


  ### standard errors of  bias-corrected estimator
  ### the normal method, assuming gaussian
  try.mat <- t(VW.A1)%*%VW.A1
  if (try.inverse(try.mat)) {
    Portho.VW <- diag(1,n.A1,n.A1) - VW.A1 %*% solve(t(VW.A1)%*%VW.A1) %*% t(VW.A1)
  } else {
    ### ridge penalty
    Portho.VW <- diag(1,n.A1,n.A1) - VW.A1 %*% solve(t(VW.A1)%*%VW.A1+lam*diag(1,r.VW,r.VW)) %*% t(VW.A1)
  }
  SigmaYD <- t(D.A1-D.rep)%*%Portho.VW%*%(Y.A1-D.A1*betaHat)/(n.A1-r.VW)
  betaHat.cor <- betaHat - SigmaYD*trace.T/iv.str
  Sd.cor <- sqrt((SigmaSqY*DT.Sq+(SigmaSqD*SigmaSqY+SigmaYD^2)*(trace.T^2)/(n.A1-r.VW))/(iv.str^2))


  ### the IV strength test
  # cn <- sqrt(tau.n/trace.T)/tau.n*(2+sqrt(trace.T2/trace.T)/tau.n)+1/tau.n^2
  cn <- 0
  Cn.V <- sqrt(trace.T2) + 2*sqrt(trace.T)*max(tau.n, sqrt(iv.str/((1-cn)*SigmaSqD*trace.T)))
  iv.thol <- (1+c0)*trace.T*SigmaSqD + sqrt(tau.n)*SigmaSqD*Cn.V


  ### the Signal Strength Test for splitting estimator
  signal.str <- sqrt(SigmaSqY*DT.Sq)
  signal.thol <- C1*(SigmaYD*trace.T+sqrt(trace.T2)*sqrt(tau.n))

  ### the Signal Strength Test for bias-corrected estimator
  signal.str.cor <- signal.str
  signal.thol.cor <- C1*sqrt(trace.T2)*sqrt(tau.n) + tau.n


  returnList <- list(Sd = Sd,
                     betaHat.cor = betaHat.cor,
                     Sd.cor = 1.1*Sd.cor,
                     SigmaYD = SigmaYD,
                     T.V = T.V,
                     trace.T = trace.T,
                     trace.T2 = trace.T2,
                     iv.str = iv.str,
                     iv.thol = iv.thol,
                     DT.Sq = DT.Sq,
                     signal.str = signal.str,
                     signal.thol = signal.thol,
                     signal.str.cor = signal.str.cor,
                     signal.thol.cor  = signal.thol.cor,
                     sigularity = sigularity)
  returnList
}


### get.sigma
### Function: a helper function for the estimate of noise level
### Input: betaHat: the estimated treatment effect by data splitting estimator
###        Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        VW: the instruments-covariance matrix expanded by corresponding basis,
###            first column is intercept
### Output: the estimate of noise level, typically SigmaSqY
get.sigma <- function(betaHat, Y, D, VW) {
  n <- length(Y)
  Portho.VW <- diag(1,n,n) - VW%*%solve(t(VW)%*%VW)%*%t(VW)
  sum((Portho.VW%*%(Y-betaHat*D))^2)/(n-ncol(VW))
}


### TSRF.Selection
### Function: Conduct the violation space selection of Two Stage Random Forest
###           using data splitting approach
### Input: Y: continuous, n by 1 outcome vector
###        D: Y: continuous or binary, n by 1 treatment vector
###        Cov.aug: the augmented instruments-covariates matrix used for selection,
###                 first column is intercept
###        A1.ind: a set of indices indicating the samples in A1
###        weight: the n.A1 by n.A1 weight matrix, usually a sparse matrix with class dgCMatrix
###        Q: number of violation space to consider
###        alpha: significance level, default by 0.05
###        tuning: tuning parameter for C.alpha, default by log(log(n))
### Output: Coef.vec: a vector of length 2*Q, the fixed violation space estimators
###         Sd.vec: the estimated standard deviation of Coef.vec
###         Coef.robust: a vector of length 4, the data splitting and bias corrected
###                      estimator using qhat.c(selection) or qhat.r(robust method)
###         Sd.robust: the estimated standard deviation of Coef.robust
###         SigmaSqY: the noise level of outcome model
###         SigmaSqD: the noise level of treatment model
###         SigmaYD: the estimated covariance between \epsilon ans \delta
###         trace.T: the trace of T^{RF}(V)
###         trace.T2: the trace of T^{RF}(V)%*%T^{RF}(V)
###         iv.str: a vector of length Q, the left hand side of IV strength test
###         iv.thol: a vector of length Q,, the right hand side of IV strength test
###         DT.sq: t(D.A1)%*%T.V%*%T.V%*%D.A1, saved for violation space selection
###         signal.str: a vector of length Q, the left hand side of signal strength test for data splitting estimator
###         signal.str: a vector of length Q, the right hand side of signal strength test for data splitting estimator
###         signal.str: a vector of length Q, the left hand side of signal strength test for bias corrected estimator
###         signal.str: a vector of length Q, the right hand side of signal strength test for bias corrected estimator
###         H: a vector of length Q-1, the variance of betaHat(q)-betaHat(q+1)
###         C.alpha: a 0-1 vector of length Q-1 indicating whether there is difference between betaHat(q) and betaHat(q+1)
###         Q.max; Qmax chosen by IV strength test
###         SigmaSqY.Qmax: the noise level of outcome model using violation space V_{Qmax}
###         qhat.c: violation space indicator selected
###         qhat.r: violation space indicator selected using robust method
###         validity: the validity test of TSLS
TSRF.Selection <- function(Y, D, Cov.aug, A1.ind, weight, Q, alpha=0.05, tuning=NULL) {

  ### constants
  n <- length(Y); n.A1 <- length(A1.ind)
  if (is.null(tuning)) {tuning <- log(log(n))}
  Y.A1 <- Y[A1.ind]; D.A1 <- D[A1.ind]; Cov.aug.A1 <- Cov.aug[A1.ind,]
  ### compute the representations
  Y.rep <- as.matrix(weight)%*%Y.A1; D.rep <- as.matrix(weight)%*%D.A1
  Cov.rep <- as.matrix(weight)%*%Cov.aug.A1

  ### the noise level of treatment model
  SigmaSqD <- mean((D.rep-D.A1)^2)

  ### save estimates for selection part
  names <- c(paste("RF-q",0:(Q-1),sep=""),paste("RF-Cor-q",0:(Q-1),sep=""))
  Coef.vec <- sd.vec <- rep(NA,2*Q)
  names(Coef.vec) <- names(sd.vec) <- names

  # IV strength test and signal strength test
  iv.str <- iv.thol <- signal.str <- signal.thol <- signal.str.cor <- signal.thol.cor <- rep(NA,Q)
  names(iv.str) <- names(iv.thol) <- names(signal.str) <- names(signal.thol) <-
    names(signal.str.cor) <- names(signal.thol.cor) <- paste("q",0:(Q-1),sep="")

  # the noise level of outcome model and covariance of \epsilon and \delta
  SigmaSqY <- SigmaYD <- rep(NA,Q)
  names(SigmaSqY) <- names(SigmaYD) <- names[1:Q]

  trace.T2 <- trace.T <- DT.Sq <- rep(NA,Q)
  names(trace.T) <- names(trace.T2) <- names(DT.Sq) <- paste("q",0:(Q-1),sep="")


  ### fixed violation space, compute necessary inputs of selection part
  # save T.V for the computation of H
  T.V <- rep(list(NA),Q)
  for (q in 0:(Q-1)) {
    if (q==Q-1) {
      reg.rf <- lm(Y.rep~D.rep+Cov.rep) # remove intercept, Cov.rep includes intercept
      betaHat <- coef(reg.rf)[2]
      Coef.vec[q+1] <- betaHat
      SigmaSqY[q+1] <- get.sigma(betaHat, Y, D, cbind(1,Cov.aug))
      stat.inputs <- TSRF.stat(Y.A1, D.A1, cbind(1,Cov.aug.A1), betaHat, weight, n, SigmaSqY[q+1], SigmaSqD)
    } else {
      reg.rf <- lm(Y.rep~D.rep+Cov.rep[,-(1:(Q-1-q))]) # remove intercept, Cov.rep includes intercept
      betaHat <- coef(reg.rf)[2]
      Coef.vec[q+1] <- betaHat
      SigmaSqY[q+1] <- get.sigma(betaHat, Y, D, cbind(1,Cov.aug[,-(1:(Q-1-q))]))
      stat.inputs <- TSRF.stat(Y.A1, D.A1, cbind(1,Cov.aug.A1[,-(1:(Q-1-q))]), betaHat, weight, n, SigmaSqY[q+1], SigmaSqD)
    }
    ### the statistics
    Coef.vec[q+Q+1] <- stat.inputs$betaHat.cor
    SigmaYD[q+1] <- stat.inputs$SigmaYD
    sd.vec[q+1] <- stat.inputs$Sd
    sd.vec[q+1+Q] <- stat.inputs$Sd.cor
    iv.str[q+1] <- stat.inputs$iv.str; iv.thol[q+1] <- stat.inputs$iv.thol;
    signal.str[q+1] <- stat.inputs$signal.str; signal.thol[q+1] <- stat.inputs$signal.thol;
    signal.str.cor[q+1] <- stat.inputs$signal.str.cor; signal.thol.cor[q+1] <- stat.inputs$signal.thol.cor;
    trace.T[q+1] <- stat.inputs$trace.T
    trace.T2[q+1] <- stat.inputs$trace.T2
    DT.Sq[q+1] <- stat.inputs$DT.Sq
    T.V[[q+1]] <- stat.inputs$T.V
  }


  ### violation space selection
  ### all of the q here are from 0 to 4, so use q+1 to index the columns
  # H is the variance of betaHat(q)-betaHat(q+1)
  # C.alpha is a 0-1 vector indicating the difference of betaHat(q) and betaHat(q+1)
  H <- C.alpha <- rep(NA,Q-1)
  ### robust estimator
  Coef.robust <- sd.robust <- rep(NA,4)
  names(Coef.robust) <- names(sd.robust) <- c("RF-c","RF-Cor-c","RF-r","RF-Cor-r")

  ivtest.vec <- iv.str>=iv.thol
  if (sum(ivtest.vec)==0) {
    Q.max <- 0
  } else {
    Q.max <- sum(ivtest.vec)-1
  }
  # plus one if Q.max==0
  if (Q.max==0) {Q.max <- Q.max+1}
  ### noise level defined by V_{Q.max}
  SigmaSqY.Qmax <- SigmaSqY[Q.max+1]
  ### compute C.alpha
  # difference between betaHat.cor, use bias corrected estimator for selection
  beta.diff <- rep(NA,Q-1)
  # the threshold of estimator difference
  thol <- rep(NA,Q-1)
  ### selection
  for (q in 0:(Q-2)) {
    H[q+1] <- DT.Sq[q+1]/(iv.str[q+1]^2)+DT.Sq[q+2]/(iv.str[q+2]^2)-2*t(D.A1)%*%T.V[[q+2]]%*%T.V[[q+1]]%*%D.A1/(iv.str[q+1]*iv.str[q+2])
    beta.diff[q+1] <- abs(Coef.vec[q+Q+2]-Coef.vec[q+Q+1])
    thol[q+1] <- qnorm(1-alpha/tuning)*SigmaSqY.Qmax*sqrt(H[q+1])
  }
  C.alpha <- ifelse(beta.diff<=thol,0,1)
  if (sum(C.alpha)==Q-1) {
    qhat.c <- Q.max
  } else {
    qhat.c <- min(which(C.alpha==0))-1
  }
  ### validity of TSLS
  if (qhat.c>=1) {
    validity <- 1
  } else {
    validity <- 0
  }
  qhat.r <- min(qhat.c+1, Q.max)
  Coef.robust[1] <- Coef.vec[qhat.c+1]
  Coef.robust[2] <- Coef.vec[qhat.c+Q+1]
  Coef.robust[3] <- Coef.vec[qhat.r+1]
  Coef.robust[4] <- Coef.vec[qhat.r+Q+1]
  sd.robust[1] <- sd.vec[qhat.c+1]
  sd.robust[2] <- sd.vec[qhat.c+Q+1]
  sd.robust[3] <- sd.vec[qhat.r+1]
  sd.robust[4] <- sd.vec[qhat.r+Q+1]

  returnList = list(Coef.vec = Coef.vec,
                    sd.vec = sd.vec,
                    Coef.robust = Coef.robust,
                    sd.robust = sd.robust,
                    SigmaSqY = SigmaSqY,
                    SigmaSqD = SigmaSqD,
                    SigmaYD = SigmaYD,
                    trace.T = trace.T, trace.T2 = trace.T2,
                    iv.str = iv.str, iv.thol = iv.thol, DT.Sq = DT.Sq,
                    signal.str = signal.str, signal.thol = signal.thol,
                    signal.str.cor = signal.str.cor,
                    signal.thol.cor = signal.thol.cor,
                    H = H, C.alpha = C.alpha,
                    Q.max = Q.max,
                    SigmaSqY.Qmax = SigmaSqY.Qmax,
                    qhat.c =qhat.c, qhat.r = qhat.r,
                    validity = validity)

}



### TSRF.stat.crossfit
### Function: Compute the point estimator and  standard deviation of
###           Two Stage Random Forest using cross-fitting approach
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        VW: the instruments-covariates matrix, first column is intercept
###        weight: a n by n symmetric sparse weight matrix(class dgCMatrix), the ith row
###                represents the weights of each sample on the prediction of the ith outcome
###        lam: the tuning parameter of ridge regression if t(VW.rep)%*%VW.rep is singular,
###             default by 0.05
### Output: betaHat: the estimated treatment effect
###         Sd: the standard deviation of betaHat
###         Singularity: logical, whether t(VW.rep)%*%VW.rep is singular
TSRF.stat.crossfit <- function(Y, D, VW, weight, lam=0.05) {
  n <- length(Y);  r.VW <- ncol(VW)
  VW.rep <- as.matrix(weight)%*%VW; Y.rep <- as.matrix(weight)%*%Y
  D.rep <- as.matrix(weight)%*%Y

  # point estimator
  reg.rf <- lm(Y.rep~D.rep+VW.rep[,-1])
  betaHat <- coef(reg.rf)[2]

  ### noise level of outcome model
  SigmaSqY <- get.sigma(betaHat, Y, D, VW)

  ### check the inverse of t(VW.rep)%*%VW.rep
  try.mat <- t(VW.rep)%*%VW.rep
  if (try.inverse(try.mat)) {
    singularity <- FALSE
    Portho.VW <- diag(1,n,n) - VW.rep %*% solve(try.mat) %*% t(VW.rep)
    T.V <- t(as.matrix(weight))%*%Portho.VW%*%as.matrix(weight)
  } else {
    singularity <- TRUE
    ### ridge penalty
    Portho.VW <- diag(1,n,n) - VW.rep %*% solve(try.mat+lam*diag(1,r.VW,r.VW)) %*% t(VW.rep)
    T.V <- t(as.matrix(weight))%*%Portho.VW%*%as.matrix(weight)
  }

  ### the standard error of the estimator and iv strength test
  iv.str <- t(D)%*%T.V%*%D
  Sd <- sqrt(SigmaSqY*(t(D)%*%T.V%*%T.V%*%D)/(iv.str)^2)

  returnList <- list(betaHat = betaHat,
                     Sd = Sd,
                     singularity = singularity)
  returnList
}


### naiveRF.stat
### Function: Compute the plug-in TSLS using random forest with
###           data splitting as the first stage
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        VW: the instruments-covariates matrix, first column is intercept
###        A1.ind: a set of indices indicating the samples in A1
###        weight: a n by n symmetric sparse weight matrix(class dgCMatrix), the ith row
###                represents the weights of each sample on the prediction of the ith outcome
### Output: betaHat: the estimated treatment effect
###         Sd: the standard deviation of betaHat
###         Singularity: logical, whether t(VW.rep)%*%VW.rep is singular
naiveRF.stat <- function(Y, D, VW, A1.ind, weight, lam=0.05) {
  n.A1 <- length(A1.ind); r.VW <- ncol(VW)
  Y.A1 <- Y[A1.ind]; VW.A1 <- VW[A1.ind,]; D.A1 <- D[A1.ind]
  D.rep <- as.matrix(weight)%*%D.A1

  ### point estimator
  reg.rf <- lm(Y.A1~D.rep+VW.A1[,-1])
  betaHat <- coef(reg.rf)[2]

  ### noise level of outcome model
  SigmaSqY <- get.sigma(betaHat, Y, D, VW)

  ### check the inverse of t(VW.rep)%*%VW.rep
  try.mat <- t(VW.A1)%*%VW.A1
  if (try.inverse(try.mat)) {
    singularity <- FALSE
    Portho.VW <- diag(1,n.A1,n.A1) - VW.A1 %*% solve(try.mat) %*% t(VW.A1)
  } else {
    singularity <- TRUE
    ### ridge penalty
    Portho.VW <- diag(1,n.A1,n.A1) - VW.A1 %*% solve(try.mat+lam*diag(1,r.VW,r.VW)) %*% t(VW.A1)
  }
  Sd <- sqrt(SigmaSqY/(t(D.rep)%*%Portho.VW%*%D.rep))

  returnList <- list(betaHat = betaHat,
                     Sd = Sd,
                     singularity = singularity)
  returnList
}


### TSRF.stat.full
### Function: Compute the Two Stage Random Forest using full data
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        VW: the instruments-covariates matrix, first column is intercept
###        weight: a n by n symmetric sparse weight matrix(class dgCMatrix), the ith row
###                represents the weights of each sample on the prediction of the ith outcome
### Output: betaHat: the estimated treatment effect
###         Sd: the standard deviation of betaHat
###         Singularity: logical, whether t(VW.rep)%*%VW.rep is singular
TSRF.stat.full <- function(Y, D, VW, weight, lam=0.05) {
  n <- length(Y); r.VW <- ncol(VW)
  D.rep <- as.matrix(weight)%*%D; Y.rep <- as.matrix(weight)%*%Y
  VW.rep <- as.matrix(weight)%*%VW

  ### point estimator
  reg.rf <- lm(Y.rep~D.rep+VW.rep[,-1])
  betaHat <- coef(reg.rf)[2]

  ### noise level of outcome model
  SigmaSqY <- get.sigma(betaHat, Y, D, VW)

  ### check the inverse of t(VW.rep)%*%VW.rep
  try.mat <- t(VW.rep)%*%VW.rep
  if (try.inverse(try.mat)) {
    singularity <- FALSE
    Portho.VW <- diag(1,n,n) - VW.rep %*% solve(try.mat) %*% t(VW.rep)
    T.V <- t(as.matrix(weight))%*%Portho.VW%*%as.matrix(weight)
  } else {
    singularity <- TRUE
    ### ridge penalty
    Portho.VW <- diag(1,n,n) - VW.rep %*% solve(try.mat+lam*diag(1,r.VW,r.VW)) %*% t(VW.rep)
    T.V <- t(as.matrix(weight))%*%Portho.VW%*%as.matrix(weight)
  }

  ### the standard error of the estimator
  iv.str <- t(D)%*%T.V%*%D
  Sd <- sqrt(SigmaSqY*(t(D)%*%T.V%*%T.V%*%D)/(iv.str)^2)

  returnList <- list(betaHat = betaHat,
                     Sd = Sd,
                     singularity = singularity)
  returnList
}

