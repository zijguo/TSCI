### Source-Random-Forest.R
### Functions: Two Stage Curvature Identification method using Random Forest
###            Estimate the treatment effect using the instrumental variable
###            approach with violation space selection


### required packages
library(ranger)
library(Matrix)


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
###        split.prop: the proportion of data to be used in A1 for inference
###        forest.save: logic, to save the random forest object or not, default by FALSE to save memory
### Output: forest.A2: random forest built on subsample A2, available if forest.save=TRUE
###         params.A2: best hyper-parameters selected by out-of-bag error
###         A1.ind: index of data in subsample A1
###         A2.ind: index of data in subsample A2
###         nodes.A1: a n_A1 by num.trees matrix containing the leaf nodes' indices of subsample A1
###                   in each tree of forest.A2
###         MSE.oob: the minimal out-of-bag error using the best hyper-parameters
TSRF.fit <- function(X,Y,num.trees=200,mtry=NULL,max.depth=0,min.node.size=5,MSE.thol=1e6,split.prop=2/3,forest.save=F) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  n <- NROW(X); p <- NCOL(X)
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
  n.A1 <- round(split.prop*n)
  A1.ind <- 1:n.A1
  A2.ind <- setdiff(1:n,A1.ind)
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
### Output: out.w: the n_A1 by n_A1 symmetric sparse weight matrix(class dgCMatrix), the ith row represents
###                the weights of each sample on the prediction of the ith outcome
### speed up version
TSRF.weight <- function(nodes) {
  n.A1 <- nrow(nodes); num.trees <- ncol(nodes)
  out.w <- matrix(0,n.A1,n.A1)
  for (j in 1:num.trees) {
    w.mat <- matrix(0, n.A1, n.A1)
    unique.nodes <- unique(nodes[,j])
    for (i in 1:length(unique.nodes)) {
      ind <- nodes[,j]==unique.nodes[i]
      inverse.w <- sum(ind)
      w <- 1/(inverse.w-1)  # to remove self-prediction
      w.vec <- ifelse(ind,yes=w,no=0)
      w.mat[ind,] <- matrix(rep(w.vec,inverse.w),inverse.w,byrow=T)/num.trees
    }
    diag(w.mat) <- 0
    out.w <- out.w + w.mat
  }
  # out.w <- out.w/rowSums(out.w)
  out.w <- Matrix(out.w, sparse = T)
  return(out.w)
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
###               corresponding basis in subsample A1
###        betaHat: the estimated treatment effect by data splitting estimator
###        weight: the n.A1 by n.A1 weight matrix, usually a sparse matrix with class dgCMatrix
###        n: full sample size
###        SigmaSqY: estimated noise level of outcome model
###        SigmaSqD: estimated noise level of treatment model
### Output: sd: the estimated standard deviation of data splitting estimator
###         betaHat.cor: the bias corrected estimator
###         sd.cor: the estimated standard deviation of bias corrected estimator
###         SigmaYD: the estimated covariance between \sigma ans \delta
###         D.resid: residuals of second stage regression
###         iv.str: the left hand side of IV strength test
###         iv.thol: the right hand side of IV strength test
###         DT.sq: t(D.A1)%*%T.V%*%T.V%*%D.A1, saved for violation space selection
TSRF.stat <- function(Y.A1, D.A1, VW.A1, betaHat, weight, n, SigmaSqY, SigmaSqD) {
  
  n.A1 <- length(Y.A1); r.VW <- NCOL(VW.A1) # the rank of (V, W)
  ### Compute the representations
  Y.rep <- as.matrix(weight%*%Y.A1); D.rep <- as.matrix(weight%*%D.A1)
  VW.rep <- as.matrix(weight%*%VW.A1)
  
  ### the standard error of the estimator and iv strength test
  trace.T <- 0
  for (j in 1:n.A1) {
    trace.T <- trace.T + sum(resid(lm(weight[,j]~VW.rep))^2)
  }
  D.resid <- resid(lm(D.rep~VW.rep))
  iv.str <- sum(D.resid^2)/(SigmaSqD)
  # this is the numerator of the variance of betaHat
  DT.Sq <- as.numeric(t(D.resid)%*%weight%*%weight%*%D.resid)
  sd <- sqrt(SigmaSqY*DT.Sq/(iv.str^2))
  
  
  ### standard errors of bias-corrected estimator
  res <- resid(lm(Y.A1-D.A1*betaHat~VW.A1))
  SigmaYD <- sum((D.A1-D.rep)*res)/(n.A1-r.VW-1)
  betaHat.cor <- betaHat - SigmaYD*trace.T/iv.str
  sd.cor <- sqrt((SigmaSqY*DT.Sq+(SigmaSqD*SigmaSqY+SigmaYD^2)*(trace.T^2)/(n.A1-r.VW-1))/(iv.str^2))
  
  
  boot.vec <- rep(NA,300)
  for (i in 1:300) {
    delta <- rnorm(n.A1, 0, sqrt(SigmaSqD))
    delta.rep <- weight%*%delta
    delta.resid <- resid(lm(as.matrix(delta.rep)~VW.rep))
    D.rep2 <- weight%*%D.rep
    boot.vec[i] <- sum(delta.resid^2) + 2*sum(D.rep2*delta.resid)
  }
  iv.thol <- max(quantile(boot.vec,0.975),20)/(SigmaSqD) # use max or 0.99
  # iv.thol <- quantile(boot.vec,0.975)
  
  returnList <- list(sd = sd,
                     betaHat.cor = betaHat.cor,
                     sd.cor = 1.05*sd.cor,
                     SigmaYD = SigmaYD,
                     D.resid = D.resid,
                     iv.str = iv.str,
                     iv.thol = iv.thol,
                     DT.Sq = DT.Sq)
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
  Y.D <- Y-betaHat*D
  Y.resid <- resid(lm(Y.D~VW))
  sum(Y.resid^2)/(n-NCOL(VW)-1)
}


### TSRF.Selection
### Function: Conduct the violation space selection of Two Stage Random Forest
###           using data splitting approach
### Input: Y: continuous, n by 1 outcome vector
###        D: Y: continuous or binary, n by 1 treatment vector
###        Cov.aug: the augmented instruments-covariates matrix used for selection
###        A1.ind: a set of indices indicating the samples in A1
###        weight: the n.A1 by n.A1 weight matrix, usually a sparse matrix with class dgCMatrix
###        Q: number of violation space to consider
### Output: Coef.vec: a vector of length 2*Q, the fixed violation space estimators
###         sd.vec: the estimated standard deviation of Coef.vec
###         Coef.robust: a vector of length 4, the data splitting and bias corrected
###                      estimator using q.comp(selection) or q.robust(robust method)
###         sd.robust: the estimated standard deviation of Coef.robust
###         SigmaSqY: the noise level of outcome model
###         SigmaSqD: the noise level of treatment model
###         SigmaYD: the estimated covariance between \epsilon ans \delta
###         iv.str: a vector of length Q, the left hand side of IV strength test
###         iv.thol: a vector of length Q,, the right hand side of IV strength test
###         C.alpha: a 0-1 vector of length Q-1 indicating whether there is difference between betaHat(q) and betaHat(q+1)
###         Q.max; Qmax chosen by IV strength test
###         SigmaSqY.Qmax: the noise level of outcome model using violation space V_{Qmax}
###         q.comp: violation space indicator selected
###         q.robust: violation space indicator selected using robust method
###         validity: the validity test of TSLS
TSRF.Selection <- function(Y, D, Cov.aug, A1.ind, weight, Q, intercept=TRUE) {
  Y <- as.matrix(Y); D <- as.matrix(D); Cov.aug <- as.matrix(Cov.aug)
  ### constants
  n <- length(Y); n.A1 <- length(A1.ind)
  Y.A1 <- Y[A1.ind]; D.A1 <- D[A1.ind]; Cov.aug.A1 <- Cov.aug[A1.ind,]
  ### compute the representations
  Y.rep <- as.matrix(weight%*%Y.A1); D.rep <- as.matrix(weight%*%D.A1)
  Cov.rep <- as.matrix(weight%*%Cov.aug.A1)
  
  ### the noise level of treatment model
  SigmaSqD <- mean((D.rep-D.A1)^2)
  
  ### save estimates for selection part
  names <- c(paste("RF-q",0:(Q-1),sep=""),paste("RF-Cor-q",0:(Q-1),sep=""))
  Coef.vec <- sd.vec <- rep(NA,2*Q)
  names(Coef.vec) <- names(sd.vec) <- names
  
  # IV strength test and signal strength test
  iv.str <- iv.thol <- rep(NA,Q)
  names(iv.str) <- names(iv.thol) <- paste("q",0:(Q-1),sep="")
  
  # the noise level of outcome model and covariance of \epsilon and \delta
  SigmaSqY <- SigmaYD <- rep(NA,Q)
  names(SigmaSqY) <- names(SigmaYD) <- names[1:Q]
  
  ### the numerator of variance
  DT.Sq <- rep(NA,Q)
  names(DT.Sq) <- paste("q",0:(Q-1),sep="")
  
  ### fixed violation space, compute necessary inputs of selection part
  # save D.resid for the computation of H and z.alpha
  D.resid <- rep(list(NA),Q)
  for (q in 0:(Q-1)) {
    if (q==Q-1) {
      if (intercept) {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep)
        betaHat <- coef(reg.rf)[2]
      } else {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep-1)
        betaHat <- coef(reg.rf)[1]
      }
      Coef.vec[q+1] <- betaHat
      SigmaSqY[q+1] <- get.sigma(betaHat, Y, D, Cov.aug)
      stat.inputs <- TSRF.stat(Y.A1, D.A1, Cov.aug.A1, betaHat, weight, n, SigmaSqY[q+1], SigmaSqD)
    } else {
      if (intercept) {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep[,-(1:(Q-1-q))])
        betaHat <- coef(reg.rf)[2]
      } else {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep[,-(1:(Q-1-q))]-1)
        betaHat <- coef(reg.rf)[1]
      }
      Coef.vec[q+1] <- betaHat
      SigmaSqY[q+1] <- get.sigma(betaHat, Y, D, Cov.aug[,-(1:(Q-1-q))])
      stat.inputs <- TSRF.stat(Y.A1, D.A1, Cov.aug.A1[,-(1:(Q-1-q))], betaHat, weight, n, SigmaSqY[q+1], SigmaSqD)
    }
    ### the statistics
    Coef.vec[q+Q+1] <- stat.inputs$betaHat.cor
    SigmaYD[q+1] <- stat.inputs$SigmaYD
    sd.vec[q+1] <- stat.inputs$sd
    sd.vec[q+1+Q] <- stat.inputs$sd.cor
    iv.str[q+1] <- stat.inputs$iv.str; iv.thol[q+1] <- stat.inputs$iv.thol;
    DT.Sq[q+1] <- stat.inputs$DT.Sq
    D.resid[[q+1]] <- stat.inputs$D.resid
  }
  
  
  ### violation space selection
  ### all of the q here are from 0 to 4, so use q+1 to index the columns
  ### robust estimator
  Coef.robust <- sd.robust <- rep(NA,4)
  names(Coef.robust) <- names(sd.robust) <- c("RF-comp","RF-Cor-comp","RF-robust","RF-Cor-robust")
  
  ivtest.vec <- iv.str>=iv.thol
  run.OLS <- weak.iv <- FALSE
  if (sum(ivtest.vec)==0) {
    warning("Weak IV: Even if the IV is assumed to be valid, run OLS") # stop, output results of OLS
    run.OLS <- TRUE
    Q.max <- 1
  } else {
    Q.max <- sum(ivtest.vec)-1
    if (Q.max==0) {
      ### if Q.max==0, redefine Qmax by log(log(n)), ignore at this stage
      warning("Weak IV: If the IV is linearly invalid. Set Qmax as 1.")
      Q.max <- 1
      weak.iv = TRUE
    }
  }
  
  ### noise level defined by V_{Q.max}
  SigmaSqY.Qmax <- SigmaSqY[Q.max+1]
  
  ### selection
  ### define comparison matrix
  H <- beta.diff <- matrix(0,Q.max,Q.max)
  ### compute H matrix
  for (q1 in 0:(Q.max-1)) {
    for (q2 in (q1+1):(Q.max)) {
      H[q1+1,q2] <- as.numeric(DT.Sq[q1+1]/(iv.str[q1+1]^2)+DT.Sq[q2+1]/(iv.str[q2+1]^2)-
                                 2*t(D.resid[[q1+1]])%*%weight%*%weight%*%D.resid[[q2+1]]/(iv.str[q1+1]*iv.str[q2+1]))
      
    }
  }
  
  ### compute beta difference matrix
  for (q in 0:(Q.max-1)) {
    beta.diff[q+1,(q+1):(Q.max)] <- abs(Coef.vec[q+Q+1]-Coef.vec[(q+Q+2):(Q.max+Q+1)]) # use original or bias-corrected?
  }
  
  
  ## bootstrap for the quantile of the difference
  max.val <- rep(NA,300)
  for (i in 1:300) {
    diff.mat <- matrix(0,Q.max,Q.max)
    eps <- rnorm(n.A1, 0, sqrt(SigmaSqY.Qmax))
    eps.rep <- weight%*%eps
    for (q1 in 0:(Q.max-1)) {
      for (q2 in (q1+1):(Q.max)) {
        diff.mat[q1+1, q2] <- sum(D.resid[[q2+1]]*eps.rep)/(iv.str[q2+1])-sum(D.resid[[q1+1]]*eps.rep)/(iv.str[q1+1])
      }
    }
    diff.mat <- abs(diff.mat)/sqrt(SigmaSqY.Qmax*H)
    max.val[i] <- max(diff.mat,na.rm = TRUE)
  }
  # z.alpha <- quantile(max.val,0.99) # use max or 0.99
  z.alpha <- quantile(max.val,0.975)
  
  diff.thol <- z.alpha*sqrt(SigmaSqY.Qmax)*sqrt(H)
  C.alpha <- ifelse(beta.diff<=diff.thol,0,1)
  ### a vector indicating the selection of each layer
  sel.vec <- apply(C.alpha,1,sum)
  if (all(sel.vec != 0)) {
    q.comp = Q.max
  } else {
    q.comp = min(which(sel.vec==0))-1
  }
  
  
  ### validity of TSLS
  if (q.comp>=1) {
    validity <- 1
  } else {
    validity <- 0
  }
  q.robust <- min(q.comp+1, Q.max)
  Coef.robust[1] <- Coef.vec[q.comp+1]
  Coef.robust[2] <- Coef.vec[q.comp+Q+1]
  Coef.robust[3] <- Coef.vec[q.robust+1]
  Coef.robust[4] <- Coef.vec[q.robust+Q+1]
  sd.robust[1] <- sd.vec[q.comp+1]
  sd.robust[2] <- sd.vec[q.comp+Q+1]
  sd.robust[3] <- sd.vec[q.robust+1]
  sd.robust[4] <- sd.vec[q.robust+Q+1]
  
  returnList = list(Coef.vec = Coef.vec,
                    sd.vec = sd.vec,
                    Coef.robust = Coef.robust,
                    sd.robust = sd.robust,
                    SigmaSqY = SigmaSqY,
                    SigmaSqD = SigmaSqD,
                    SigmaYD = SigmaYD,
                    iv.str = iv.str, iv.thol = iv.thol,
                    Q.max = Q.max,
                    SigmaSqY.Qmax = SigmaSqY.Qmax,
                    q.comp =q.comp, q.robust = q.robust,
                    validity = validity,
                    run.OLS = run.OLS,
                    weak.iv = weak.iv)
  
}


### naiveRF.stat
### Function: Compute the plug-in TSLS using random forest with
###           data splitting as the first stage
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        VW: the instruments-covariates matrix
###        A1.ind: a set of indices indicating the samples in A1
###        weight: a n by n symmetric sparse weight matrix(class dgCMatrix), the ith row
###                represents the weights of each sample on the prediction of the ith outcome
### Output: betaHat: the estimated treatment effect
###         sd: the standard deviation of betaHat
###         singularity: logical, whether t(VW.rep)%*%VW.rep is singular
naiveRF.stat <- function(Y, D, VW, A1.ind, weight, lam=0.05) {
  n.A1 <- length(A1.ind); r.VW <- ncol(VW)
  Y.A1 <- Y[A1.ind]; VW.A1 <- VW[A1.ind,]; D.A1 <- D[A1.ind]
  D.rep <- as.matrix(weight%*%D.A1)
  
  ### point estimator
  reg.rf <- lm(Y.A1~D.rep+VW.A1)
  betaHat <- coef(reg.rf)[2]
  
  ### noise level of outcome model
  SigmaSqY <- get.sigma(betaHat, Y, D, VW)
  
  D.resid <- resid(lm(D.rep~VW.A1))
  sd <- sqrt(SigmaSqY/sum(D.resid^2))
  
  returnList <- list(betaHat = betaHat,
                     sd = sd,
                     SigmaSqY = SigmaSqY)
  returnList
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
    temp<- ranger(Y~., data = Data,
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


### TSRF.stat.full
### Function: Compute the Two Stage Random Forest using full data
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        VW: the instruments-covariates matrix
###        weight: a n by n symmetric sparse weight matrix(class dgCMatrix), the ith row
###                represents the weights of each sample on the prediction of the ith outcome
### Output: betaHat: the estimated treatment effect
###         sd: the standard deviation of betaHat
###         singularity: logical, whether t(VW.rep)%*%VW.rep is singular
TSRF.stat.full <- function(Y, D, VW, weight) {
  n <- length(Y); r.VW <- ncol(VW)
  D.rep <- as.matrix(weight%*%D); Y.rep <- as.matrix(weight%*%Y)
  VW.rep <- as.matrix(weight%*%VW)
  
  ### point estimator
  reg.rf <- lm(Y.rep~D.rep+VW.rep)
  betaHat <- coef(reg.rf)[2]
  
  ### noise level of outcome model
  SigmaSqY <- get.sigma(betaHat, Y, D, VW)
  
  ### the standard error of the estimator
  D.resid <- resid(lm(D.rep~VW.rep))
  iv.str <- sum(D.resid^2)
  sd <- sqrt(SigmaSqY*(t(D.resid)%*%as.matrix(weight)%*%as.matrix(weight)%*%D.resid)/(iv.str^2))
  
  returnList <- list(betaHat = betaHat,
                     sd = sd,
                     SigmaSqY = SigmaSqY)
  returnList
}


# ### TSRF.crossfit
# ### Function: Model fitting function for TWo Stage Random Forest
# ###           using cross fitting approach
# ###           Use out-of-bag error for hyper-parameter tuning
# ### Input: X: continuous or binary, n by p_x covariates
# ###        Y: continuous or binary, n by 1 outcome vector
# ###        k: integer, number of subsamples for cross fitting
# ###        num.trees: integer, the number of trees in random forest
# ###        mtry: integer, the number of covariates to split at each node
# ###        max.depth: integer, the maximal depth of each tree, 0 refers to unlimited depth
# ###        min.node.size: integer, the minimal size(# samples in it) of each leaf node
# ###        MSE.thol: numeric, a large value of MSE, used for the start of hyper-parameter selection
# ###        forest.save: logic, to save the random forest object or not, default by FALSE to save memory
# ### Output: forest: a list of k random forest objects, available if forest.save=TRUE
# ###         params: a list of k sets of best hyper-parameters selected by out-of-bag error
# ###         k.ind: a list of vectors indicating the index of the each sample in k subsamples
# ###         nodes: a list of nodes information of the k subsamples, similar to nodes.A1 in TSRF.fit
# ###         MSE.oob: the minimal out-of-bag error using the best hyper-parameters
# TSRF.crossfit <- function(X,Y,k=2,num.trees=200,mtry=NULL,max.depth=0,min.node.size=5,MSE.thol=1e6,forest.save=F) {
#   n <- nrow(X);p<-ncol(X)
#   if (is.null(mtry)) mtry <- round(p/3)
#   
#   Data <- data.frame(cbind(Y, X))
#   names(Data) <- c("Y", paste("X", 1:p, sep = ""))
#   ### grid search
#   params.grid <- expand.grid(
#     num.trees = num.trees,
#     mtry = mtry,
#     max.depth = max.depth,
#     min.node.size = min.node.size
#   )
#   
#   # split the data into k subsamples
#   # use cross-ftting procedures
#   ind.sep <- cut(1:n, breaks = k, labels = FALSE)
#   k.ind <- rep(list(NA), k)
#   for (j in 1:k) {
#     k.ind[[j]] <- which(ind.sep==j)
#   }
#   
#   forest <- params <- rep(list(NA), k)
#   ### use oob error to do hyper-parameter tuning
#   MSE.oob <- rep(MSE.thol, k)
#   for (j in 1:k) {
#     for (i in 1:nrow(params.grid)) {
#       temp.k <- ranger(y~., data = Data[-k.ind[[j]],],
#                        num.trees=params.grid$num.trees[i],
#                        mtry=params.grid$mtry[i],
#                        max.depth = params.grid$max.depth[i],
#                        min.node.size = params.grid$min.node.size[i]
#       )
#       if (temp.k$prediction.error < MSE.oob[j]) {
#         forest[[j]] <- temp.k
#         params[[j]] <- params.grid[i,]
#         MSE.oob[j] <- temp.k$prediction.error
#       }
#     }
#   }
#   
#   ### nodes contains k matrices
#   nodes <- list(rep(NA), k)
#   for (j in 1:k) {
#     nodes[[j]] <- predict(forest[[j]],data=Data[k.ind[[j]],],type="terminalNodes")$predictions
#   }
#   
#   returnList <- list(forest = forest,
#                      params = params,
#                      k.ind = k.ind,
#                      nodes = nodes,
#                      MSE.oob = MSE.oob)
#   if (!forest.save) returnList <- returnList[-1]
#   returnList
# }
# 
# 
# ### TSRF.weight.crossfit
# ### Function: Compute the weight matrix of cross-fitting random forest, this function
# ###           inherits the TSRF.weight() function and augment it to cross-fitting
# ### Input: nodes: a list of nodes information, typically from the output of TSRF.crossfit()
# ### Output: out.w: a n by n symmetric sparse weight matrix(class dgCMatrix), the ith row
# ###         represents the weights of each sample on the prediction of the ith outcome
# TSRF.weight.crossfit <- function(nodes) {
#   w.list <- lapply(nodes, weight.2split)
#   out.w <- bdiag(w.list)
#   return(out.w)
# }
# 
# 
# ### TSRF.stat.crossfit
# ### Function: Compute the point estimator and  standard deviation of
# ###           Two Stage Random Forest using cross-fitting approach
# ### Input: Y: continuous, n by 1 outcome vector
# ###        D: continuous or binary, n by 1 treatment vector
# ###        VW: the instruments-covariates matrix, first column is intercept
# ###        weight: a n by n symmetric sparse weight matrix(class dgCMatrix), the ith row
# ###                represents the weights of each sample on the prediction of the ith outcome
# ###        lam: the tuning parameter of ridge regression if t(VW.rep)%*%VW.rep is singular,
# ###             default by 0.05
# ### Output: betaHat: the estimated treatment effect
# ###         sd: the standard deviation of betaHat
# ###         singularity: logical, whether t(VW.rep)%*%VW.rep is singular
# TSRF.stat.crossfit <- function(Y, D, VW, weight, lam=0.05) {
#   n <- length(Y);  r.VW <- ncol(VW)
#   VW.rep <- as.matrix(weight)%*%VW; Y.rep <- as.matrix(weight)%*%Y
#   D.rep <- as.matrix(weight)%*%Y
#   
#   # point estimator
#   reg.rf <- lm(Y.rep~D.rep+VW.rep[,-1])
#   betaHat <- coef(reg.rf)[2]
#   
#   ### noise level of outcome model
#   SigmaSqY <- get.sigma(betaHat, Y, D, VW)
#   
#   ### check the inverse of t(VW.rep)%*%VW.rep
#   try.mat <- t(VW.rep)%*%VW.rep
#   if (try.inverse(try.mat)) {
#     singularity <- FALSE
#     Portho.VW <- diag(1,n,n) - VW.rep %*% solve(try.mat) %*% t(VW.rep)
#     T.V <- t(as.matrix(weight))%*%Portho.VW%*%as.matrix(weight)
#   } else {
#     singularity <- TRUE
#     ### ridge penalty
#     Portho.VW <- diag(1,n,n) - VW.rep %*% solve(try.mat+lam*diag(1,r.VW,r.VW)) %*% t(VW.rep)
#     T.V <- t(as.matrix(weight))%*%Portho.VW%*%as.matrix(weight)
#   }
#   
#   ### the standard error of the estimator and iv strength test
#   iv.str <- t(D)%*%T.V%*%D
#   sd <- sqrt(SigmaSqY*(t(D)%*%T.V%*%T.V%*%D)/(iv.str)^2)
#   
#   returnList <- list(betaHat = betaHat,
#                      sd = sd,
#                      singularity = singularity)
#   returnList
# }

