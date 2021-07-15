### Source-Random-Forest.R
### Functions: Two Stage Curvature Identification using Random Forest.
###            Estimate the treatment effect using the instrumental variable
###            approach with violation space selection.


### required packages
library(ranger)
library(Matrix)


### TSCI.RF
### Function: Two Stage Curvature Identification using Random Forest with
###           violation space selection
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        Z: continuous or binary, n by 1 IV vector
###        X: continuous or binary, n by p_x baseline covariates matrix
###        intercept: Whether to include intercept or not in the outcome model, default by TRUE
###        vio.space: n by (Q-1) matrix, each column refers to a violation form of Z,
###                   default by NULL assumes linear and quadratic violation (Z^2,Z),
###                   choose between quadratic, linear and no violation(Q=3)
###        layer: logic, do layer selection of violation space or not, default by TRUE
###        split.prop: A value between 0 and 1, the proportion of samples we use in A1
###                    default by 2/3
###        num.trees: Number of trees, default by 200
###        mtry: Number of covariates to possibly split at in each node, default by p/3 to 2p/3
###        max.depth: Maximal depth of each tree, default by sequence (0,2,4,6) with 0 refers to unlimited depth
###        min.node.size: Minimal size of each leaf node, default by sequence (5,10,15)
###        str.thol: the minimal value of the threshold of IV strength test, default by 20
### Output: The same output of TSCI.RF.Selection
TSCI.RF <- function(Y,D,Z,X,intercept=TRUE,vio.space=NULL,layer=TRUE,split.prop=2/3,
                    num.trees=NULL,mtry=NULL,max.depth=NULL,min.node.size=NULL,str.thol=20) {
  Y <- as.matrix(Y); D <- as.matrix(D); Z <- as.matrix(Z); X <- as.matrix(X);
  # constants
  n <- NROW(X); p <- NCOL(X)
  # default value for hyper-parameters
  if (is.null(num.trees)) num.trees <- 200
  if (is.null(mtry)) mtry <- seq(round(p/3),round(2*p/3),by=1)
  if (is.null(max.depth)) max.depth <- c(0,2,4,6)
  if(is.null(min.node.size)) min.node.size <- c(5,10,20)
  # define the augmentation of covariates,
  # which is the combination of violation space and baseline covariates
  if (is.null(vio.space)) {
    Q = 3
    vio.space <- poly(Z,degree = Q-1,raw = TRUE)
  } else {
    Q = NCOL(vio.space) + 1
  }
  Cov.aug <- cbind(vio.space,X)

  # Treatment model fitting
  forest <- TSCI.RF.fit(D,Z,X,num.trees=num.trees,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size,split.prop=split.prop)
  A1.ind <- forest$A1.ind
  # Compute weight matrix
  weight <- TSCI.RF.weight(forest$nodes.A1)

  # Selection
  outputs <- TSCI.RF.Selection(Y,D,Cov.aug,A1.ind,weight=weight,Q=Q,intercept=intercept,layer=layer,str.thol=str.thol)
  return(outputs)
}


### TSCI.RF.fit
### Function: First Stage Model Fitting for TSCI using Random Forest with data splitting.
###           Use out-of-bag error for hyper-parameter tuning.
### Input: D: continuous or binary, n by 1 treatment vector
###        Z: continuous or binary, n by 1 IV vector
###        X: continuous or binary, n by p_x baseline covariates matrix
###        num.trees: Number of trees, default by 200
###        mtry: Number of covariates to possibly split at in each node, default by p/3 to 2p/3
###        max.depth: Maximal depth of each tree, default by sequence (0,2,4,6) with 0 refers to unlimited depth
###        min.node.size: Minimal size of each leaf node, default by sequence (5,10,15)
###        split.prop: Proportion of data in subsample A1 for inference
###        MSE.thol: A very large value of MSE, default by 1e12, used for the start of hyper-parameter selection
###        forest.save: To save the Random Forest output or not, default by FALSE to save memory
### Output: forest.A2: Random Forest built on subsample A2, available if forest.save=TRUE
###         params.A2: Best hyper-parameters for forest.A2 selected by out-of-bag error
###         A1.ind: Indices of subsample A1
###         nodes.A1: A n_A1 by num.trees matrix, rows refer to different samples, columns refer to
###                   different trees, the entrees are leaf node indices of each sample in each tree
###         MSE.oob: Minimal out-of-bag error using the best hyper-parameters
TSCI.RF.fit <- function(D,Z,X,num.trees,mtry,max.depth,min.node.size,split.prop,MSE.thol=1e12,forest.save=T) {
  W <- as.matrix(cbind(Z,X)); D <- as.matrix(D)
  n <- NROW(W); p <- NCOL(W)
  Data <- data.frame(cbind(D,W))
  names(Data) <- c("D", paste("W", 1:p, sep = ""))
  # grid search
  params.grid <- expand.grid(
    num.trees = num.trees,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )
  # split the data into two parts A1 and A2
  # use A2 to build the random forest and use A1 to do inference
  n.A1 <- round(split.prop*n)
  A1.ind <- 1:n.A1
  Data.A1 <- Data[A1.ind,]
  Data.A2 <- Data[-A1.ind,]
  forest.A2 <- NULL;
  MSE.oob.A2 <- MSE.thol
  params.A2 <- NULL
  ### use oob error to do hyper-parameter tuning
  for (i in 1:nrow(params.grid)) {
    temp.A2 <- ranger(D~., data = Data.A2,
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
  # leaf nodes information of A1 on the Random Forest built on A2
  nodes.A1 <- predict(forest.A2, data = Data.A1, type = "terminalNodes")$predictions
  returnList <- list(forest.A2 = forest.A2,
                     params.A2 = params.A2,
                     A1.ind = A1.ind,
                     nodes.A1 = nodes.A1,
                     MSE.oob.A2 = MSE.oob.A2)
  if (!forest.save) returnList <- returnList[-1]
  returnList
}


### TSCI.RF.weight
### Function: Compute the weight matrix of Random Forest with data splitting and Random Forest with full data
### Input: nodes: A n_A1 by num.trees matrix, rows refer to different samples, columns refer to
###               different trees, the entrees are leaf node indices of each sample in each tree
### Output: out.weight: A n_A1 by n_A1 symmetric sparse weight matrix(class dgCMatrix), the ith row represents
###                the weights of outcome of each sample on the prediction of the ith outcome
TSCI.RF.weight <- function(nodes) {
  n.A1 <- NROW(nodes); num.trees <- NCOL(nodes)
  out.weight <- matrix(0,n.A1,n.A1)
  for (j in 1:num.trees) {
    weight.mat <- matrix(0,n.A1,n.A1) # weight matrix for single tree
    unique.nodes <- unique(nodes[,j])
    for (i in 1:length(unique.nodes)) {
      ind <- nodes[,j]==unique.nodes[i] # indices of samples in the node
      num.samples <- sum(ind) # number of samples in the node
      w <- 1/(num.samples-1)  # weight, to remove self-prediction
      weight.vec <- ifelse(ind,yes=w,no=0)
      weight.mat[ind,] <- matrix(rep(weight.vec,num.samples),num.samples,byrow=T)/num.trees
    }
    diag(weight.mat) <- 0 # remove self prediction
    out.weight <- out.weight + weight.mat
  }
  out.weight <- Matrix(out.weight, sparse = T) # sparse matrix to save memory
  return(out.weight)
}


### TSCI.RF.stat
### Function: Compute the necessary statistics for TSCI with Random Forest
### Input: Y.rep: continuous, a n.A1 by 1 vector denoting the outcome representation in subsample A1
###        D.rep: continuous or binary, a n.A1 by 1 vector denoting the treatment representation in subsample A1
###        Cov.rep: continuous or binary, the augmented instruments-covariates matrix representation in subsample A1
###        betaHat: Estimated treatment effect by data splitting estimator
###        weight: n.A1 by n.A1 weight matrix
###        n: full sample size
###        SigmaSqY: Estimated noise level of outcome model
###        SigmaSqD: Estimated noise level of treatment model
###        SigmaYD: Estimated covariance between error terms in outcome and treatment models
###        str.thol: the minimal value of the threshold of IV strength test, default by 20
### Output: betaHat: Estimated treatment effect by data splitting estimator
###         sd: the estimated standard deviation of data splitting estimator
###         betaHat.cor: the bias corrected estimator
###         sd.cor: the estimated standard deviation of bias corrected estimator
###         D.resid: Residuals of second stage regression
###         iv.str: IV strength
###         iv.thol: IV Strength Test threshold
###         explained.iv: t(D.A1)%*%T.V%*%T.V%*%D.A1, saved for violation space selection
TSCI.RF.stat <- function(Y.rep, D.rep, Cov.rep, betaHat, weight, n, SigmaSqY, SigmaSqD, SigmaYD, str.thol) {
  n.A1 <- length(Y.rep); r.aug <- NCOL(Cov.rep)
  # compute the trace of T(V)
  trace.T <- 0 
  for (j in 1:n.A1) {
    trace.T <- trace.T + sum(resid(lm(weight[,j]~Cov.rep))^2)
  }
  D.resid <- resid(lm(D.rep~Cov.rep))
  D.RSS <- sum(D.resid^2)
  iv.str <- D.RSS/(SigmaSqD)
  # this is the numerator of the variance of betaHat
  explained.iv <- as.numeric(t(D.resid)%*%weight%*%weight%*%D.resid)
  sd <- sqrt(SigmaSqY*explained.iv/(D.RSS^2))

  ### standard errors of bias-corrected estimator
  betaHat.cor <- betaHat - SigmaYD*trace.T/D.RSS
  sd.cor <- sqrt((SigmaSqY*explained.iv+(SigmaSqD*SigmaSqY+SigmaYD^2)*(trace.T^2)/(n.A1-r.aug-1))/(D.RSS^2))

  # bootstrap for the threshold of IV strength test
  boot.vec <- rep(NA,300)
  for (i in 1:300) {
    delta <- rnorm(n.A1,0,sqrt(SigmaSqD))
    delta.rep <- weight%*%delta
    delta.resid <- resid(lm(as.matrix(delta.rep)~Cov.rep))
    D.rep2 <- weight%*%D.rep
    boot.vec[i] <- sum(delta.resid^2) + 2*sum(D.rep2*delta.resid)
  }
  iv.thol <- max(quantile(boot.vec,0.975),str.thol)/(SigmaSqD)
  scale <- 1
  returnList <- list(betaHat = betaHat,
                     sd = sd,
                     betaHat.cor = betaHat.cor,
                     sd.cor = scale*sd.cor,
                     D.resid = D.resid,
                     iv.str = iv.str,
                     iv.thol = iv.thol,
                     explained.iv = explained.iv)
  returnList
}


### TSCI.RF.Selection
### Function: Violation space selection of Two Stage Curvature Identification using Random Forest
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        Cov.aug: Augmented combiantion of violation space and baseline covariates
###        A1.ind: Indices of samples in A1
###        weight: n.A1 by n.A1 weight matrix
###        Q: Number of violation space, including no violation
###        intercept: Include intercept in the outcome model or not
###        layer: logic, layer selection or not
###        str.thol: the minimal value of the threshold of IV strength test
### Output: Coef.vec: A vector of length 2*Q, the point and bias corrected estimators with fixed violation spaces
###         sd.vec: Estimated standard deviation of Coef.vec
###         Coef.robust: A vector of length 4, the point and bias corrected estimators in the violation space using
###                      comparison method and robust method
###         sd.robust: Estimated standard deviation of Coef.robust
###         SigmaSqY: Estimated noise level of outcome model
###         SigmaSqD: Estimated noise level of treatment model
###         SigmaYD: Estimated covariance between error terms in outcome and treatment models
###         iv.str: A vector of length Q, the left hand side of IV strength test
###         iv.thol: A vector of length Q,, the right hand side of IV strength test
###         Q.max; Maximal violation space chosen by IV strength test
###         SigmaSqY.Qmax: Noise level of outcome model in violation space Q.max
###         q.comp: Violation space using comparison method
###         q.robust: Violation space using robust method
###         invalidity: Invalidity of TSLS, means there is at least linear violation
###         weak.iv: logic, whether the IV is weak when assuming linearly invalid, IV strength test only passes with no violation
###         run.OLS: logic, whether the IV is weak when assuming no violation, IV strength test does not pass even with no violation
###                         We need to run OLS if TRUE since TSLS result is not reliable
TSCI.RF.Selection <- function(Y, D, Cov.aug, A1.ind, weight, Q, intercept, layer, str.thol) {
  Y <- as.matrix(Y); D <- as.matrix(D); Cov.aug <- as.matrix(Cov.aug)
  # constants
  n <- length(Y); n.A1 <- length(A1.ind); r.aug <- NCOL(Cov.aug)
  Y.A1 <- Y[A1.ind]; D.A1 <- D[A1.ind]; Cov.aug.A1 <- Cov.aug[A1.ind,]
  # compute the representations
  Y.rep <- as.matrix(weight%*%Y.A1); D.rep <- as.matrix(weight%*%D.A1)
  Cov.rep <- as.matrix(weight%*%Cov.aug.A1)
  # the noise level of treatment model
  SigmaSqD <- mean((D.rep-D.A1)^2)
  # save estimates for selection part
  names <- c(paste("RF-q",0:(Q-1),sep=""),paste("RF-Cor-q",0:(Q-1),sep=""))
  Coef.vec <- sd.vec <- rep(NA,2*Q)
  names(Coef.vec) <- names(sd.vec) <- names
  # IV strength test and signal strength test
  iv.str <- iv.thol <- rep(NA,Q)
  names(iv.str) <- names(iv.thol) <- paste("q",0:(Q-1),sep="")
  # the noise level of outcome model and covariance of epsilon and delta
  SigmaSqY <- SigmaYD <- rep(NA,Q)
  names(SigmaSqY) <- names(SigmaYD) <- names[1:Q]
  # the numerator of variance
  explained.iv <- rep(NA,Q)
  names(explained.iv) <- paste("q",0:(Q-1),sep="")

  ### fixed violation space, compute necessary inputs of selection part
  # save D.resid for the computation of H and z.alpha
  D.resid <- rep(list(NA),Q)
  for (index in 1:Q) {
    q <- index-1
    if (q==Q-1) {
      if (intercept) {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep)
        betaHat <- coef(reg.rf)[2]
      } else {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep-1)
        betaHat <- coef(reg.rf)[1]
      }
      Coef.vec[index] <- betaHat
      SigmaSqY[index] <- get.sigma(betaHat, Y, D, Cov.aug)
      SigmaYD[index] <- sum((D.A1-D.rep)*resid(lm(Y.A1-D.A1*betaHat~Cov.aug.A1)))/(n.A1-r.aug-1)
      stat.outputs <- TSCI.RF.stat(Y.rep,D.rep,Cov.rep,betaHat,weight,n,SigmaSqY[index],SigmaSqD,SigmaYD[index],str.thol=str.thol)
    } else {
      if (intercept) {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep[,-(1:(Q-1-q))])
        betaHat <- coef(reg.rf)[2]
      } else {
        reg.rf <- lm(Y.rep~D.rep+Cov.rep[,-(1:(Q-1-q))]-1)
        betaHat <- coef(reg.rf)[1]
      }
      Coef.vec[index] <- betaHat
      SigmaSqY[index] <- get.sigma(betaHat, Y, D, Cov.aug[,-(1:(Q-1-q))])
      SigmaYD[index] <- sum((D.A1-D.rep)*resid(lm(Y.A1-D.A1*betaHat~Cov.aug.A1[,-(1:(Q-1-q))])))/(n.A1-r.aug-1)
      stat.outputs <- TSCI.RF.stat(Y.rep,D.rep,Cov.rep[,-(1:(Q-1-q))],betaHat,weight,n,
                                   SigmaSqY[index],SigmaSqD,SigmaYD[index],str.thol=str.thol)
    }
    # the statistics
    Coef.vec[index+Q] <- stat.outputs$betaHat.cor
    sd.vec[index] <- stat.outputs$sd
    sd.vec[index+Q] <- stat.outputs$sd.cor
    iv.str[index] <- stat.outputs$iv.str; iv.thol[index] <- stat.outputs$iv.thol;
    explained.iv[index] <- stat.outputs$explained.iv
    D.resid[[index]] <- stat.outputs$D.resid
  }
  # Residual sum of squares of D.rep~Cov.rep
  D.RSS <- iv.str*SigmaSqD


  # violation space selection
  # all of the q below are from 0 to (Q-1), so use q+1 to index the columns
  # Comparison and robust estimators
  Coef.robust <- sd.robust <- rep(NA,4)
  names(Coef.robust) <- names(sd.robust) <- c("RF-comp","RF-Cor-comp","RF-robust","RF-Cor-robust")
  ivtest.vec <- (iv.str>=iv.thol)
  run.OLS <- weak.iv <- FALSE
  if (sum(ivtest.vec)==0) {
    warning("Weak IV, if the IV is assumed to be valid, run OLS") # stop, output results of OLS
    run.OLS <- TRUE
    Q.max <- 1
  } else {
    Q.max <- sum(ivtest.vec)-1
    if (Q.max==0) {
      ### if Q.max==0, redefine Qmax by log(log(n)), ignore at this stage
      warning("Weak IV, if the IV is invalid. We still test the IV invalidity.")
      Q.max <- 1
      weak.iv = TRUE
    }
  }
  # noise level estimated in violation space Q.max
  SigmaSqY.Qmax <- SigmaSqY[Q.max+1]

  ### Selection
  # define comparison matrix
  H <- beta.diff <- matrix(0,Q.max,Q.max)
  # compute H matrix
  for (q1 in 0:(Q.max-1)) {
    for (q2 in (q1+1):(Q.max)) {
      H[q1+1,q2] <- as.numeric(explained.iv[q1+1]/(D.RSS[q1+1]^2)+explained.iv[q2+1]/(D.RSS[q2+1]^2)-
                                 2*t(D.resid[[q1+1]])%*%weight%*%weight%*%D.resid[[q2+1]]/(D.RSS[q1+1]*D.RSS[q2+1]))

    }
  }
  # compute beta difference matrix
  for (q in 0:(Q.max-1)) {
    beta.diff[q+1,(q+1):(Q.max)] <- abs(Coef.vec[q+Q+1]-Coef.vec[(q+Q+2):(Q.max+Q+1)]) # use original or bias-corrected?
  }
  # bootstrap for the quantile of the differences
  max.val <- rep(NA,300)
  for (i in 1:300) {
    diff.mat <- matrix(0,Q.max,Q.max)
    eps <- rnorm(n.A1, 0, sqrt(SigmaSqY.Qmax))
    eps.rep <- weight%*%eps
    for (q1 in 0:(Q.max-1)) {
      for (q2 in (q1+1):(Q.max)) {
        diff.mat[q1+1, q2] <- sum(D.resid[[q2+1]]*eps.rep)/(D.RSS[q2+1])-sum(D.resid[[q1+1]]*eps.rep)/(D.RSS[q1+1])
      }
    }
    diff.mat <- abs(diff.mat)/sqrt(SigmaSqY.Qmax*H)
    max.val[i] <- max(diff.mat,na.rm = TRUE)
  }
  z.alpha <- quantile(max.val,0.975)
  diff.thol <- z.alpha*sqrt(SigmaSqY.Qmax)*sqrt(H)
  # comparison matrix
  C.alpha <- ifelse(beta.diff<=diff.thol,0,1)

  # layer selection or not
  if (layer==TRUE) {
    # a vector indicating the selection of each layer
    sel.vec <- apply(C.alpha,1,sum)
    if (all(sel.vec != 0)) {
      q.comp = Q.max
    } else {
      q.comp = min(which(sel.vec==0))-1
    }
  } else {
    ### What if Q = 3(q2) and Q.max = 1?
    sel.val <- C.alpha[1,Q.max]
    if (sel.val==1) {
      q.comp = Q.max
    } else {
      q.comp = 0
    }
  }

  ### invalidity of TSLS
  if (q.comp>=1) {
    invalidity <- 1
  } else {
    invalidity <- 0
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
                    invalidity = invalidity,
                    run.OLS = run.OLS,
                    weak.iv = weak.iv)
  returnList
}


### get.sigma
### Function: A helper function for the estimate of noise level
### Input: betaHat: the estimated treatment effect by data splitting estimator
###        Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        Cov.aug: the instruments-covariance matrix expanded by corresponding basis,
###            first column is intercept
### Output: the estimate of noise level, typically SigmaSqY
get.sigma <- function(betaHat, Y, D, Cov.aug) {
  n <- length(Y)
  Y.D <- Y-betaHat*D
  Y.resid <- resid(lm(Y.D~Cov.aug))
  sum(Y.resid^2)/(n-NCOL(Cov.aug)-1)
}

