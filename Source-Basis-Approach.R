library(MASS)
library(fda)
library(AER)

### get.design
### Function: obtain the design matrix for the first column
###           of covaraites, prepared for semi-parametric model
### Input: x: n by p covariated matrix
###        knots:
### Output: a list of
###         m: the design matrix
###         basis: the basis object of package fda
getDesign <- function(x, knots) {
  p <- NCOL(x)
  
  if(is.null(colnames(x)))
    colnames(x) <- paste("x", 1:p, sep = "")
  
  x.current <- x[,1]
  
  ux.current <- unique(x.current)
  
  knots.use <- quantile(ux.current, seq(0, 1, length = knots))
  
  
  
  stopifnot(all.equal(range(knots.use), range(x.current)))
  
  basis <- create.bspline.basis(rangeval = range(knots.use), 
                                breaks = knots.use, norder = 4)
  
  m <- eval.basis(x.current, basis)
  list(m = m, basis = basis)
}


### create_knots
### Function: A helper function to create knots
### Input: P: the number of covariates
###        N: sample size
###        folds:
###        begin:
###        end:
###        min_knots:
###        max_knots:
###        num:
###        exe_max:
###        p_max:
### Output: knots:
create_knots <- function(P, N, folds, begin, end, l, min_knots = 2, max_knots = 50, num = 10, exe_max = FALSE, p_max = 0.5) {
  l_P <- length(P)
  l_N <- length(N)
  knots <- NULL
  for(i in 1:l_P)
  {
    knots_i <- NULL
    for(j in 1:l_N)
    {
      n <- ((P[i] * N[j] * (1 - 1 / folds)))^0.8
      if (max_knots >= n) max_knots <- round(n * p_max0)  
      knot_p <- seq(begin, end, length.out = l)
      knot <- n * knot_p
      knot <- sapply(knot, floor)
      avai_k <- knot_p[which(knot >= min_knots & knot <= max_knots)]
      len_knot <- length(knot)
      
      if (length(avai_k) == 0)
      {
        a <- min_knots / n
        b <- max_knots / n
        avai_k <- seq(a, b, length.out = l)
      }else{
        
        if (knot[1] > min_knots)
        {
          portion_min <- min_knots / n 
          portion_plus <- seq(portion_min, knot_p[1], length.out = num + 2)
          avai_k <- c(portion_plus, avai_k[-1])
        }
        
        if (exe_max == TRUE)
        {
          if (knot[len_knot] < max_knots)
          {
            portion_max <- max_knots / n 
            portion_plus <- seq(knot_p[len_knot], portion_max, length.out = num + 2)
            avai_k <- c(avai_k[-len_knot], portion_plus)
          }
        }
        
      }
      
      knots_i[[j]] <- avai_k
      
    }
    knots[[i]] <- knots_i
  }  
  return (knots)
}


### ESTIMATE
### Function: Semi-parametric model fitting function, nonliner
###           in the first covariates
### Input: X: the n by p covariates matrix
###        Y: the n by 1 outcome
###        knots:
### Output: an object with
###         m: the design matrix
###         basis: the the basis object of package fda
###         coefs: the estimated coefficients of the model
ESTIMATE <- function(X, Y, knots) {
  n_row <- length(X[,1])
  D_X <- X
  D_Y <- Y
  obj <- getDesign(D_X, knots)
  
  D_X <- cbind(obj$m, D_X[,-1])
  lmod <- lm(D_Y~ D_X - 1)
  coefs <- coef(lmod)
  obj$coefs <- coefs
  return(obj)
}


### pred
### Function: Predict the representation of the outcome Y 
###           defined in ESTIMATE()
### Input: object: an object obtained by ESTIMATE()
###        newdata: new covariates to use
### Output: pred.pre: the predicted representation of outcome
pred <- function(object, newdata) {
  x <- newdata
  n <- NROW(x)
  p <- NCOL(x)
  
  lowdiff <- highdiff <-
    matrix(FALSE, nrow = NROW(newdata), ncol = 1)
  
  x.cut  <- matrix(0, nrow = n, ncol = 1)
  
  x.current     <- x[,1]
  x.cut[,1]     <- x.current
  
  bas <- object$basis
  
  lower.end <- bas$rangeval[1]
  upper.end <- bas$rangeval[2]
  
  
  ind.lower <- x.current < lower.end
  ind.upper <- x.current > upper.end
  
  lowdiff[ind.lower,1]  <- (x.current - lower.end)[ind.lower]
  highdiff[ind.upper,1] <- (x.current - upper.end)[ind.upper]
  
  x.cut[ind.lower,1] <- lower.end
  x.cut[ind.upper,1] <- upper.end
  
  ## Get the slopes at the boundaries
  m <- eval.basis(x.cut[,1], bas)
  deriv.info  <- eval.basis(c(lower.end, upper.end), bas, Lfdobj = 1)
  
  df <- NCOL(m)
  
  lower.slopes <- deriv.info[1,]
  upper.slopes <- deriv.info[2,]
  
  
  beta <- object$coefs
  
  pred.pre <- cbind(m, x[, -1]) %*% beta
  
  ## Put the design matrix of the first derivates (lower.slopes,
  ## upper.slopes) into one long vector each (same length as index) and
  ## multiply with beta vector and take the sum. I.e. perform the matrix
  ## operation in a bit a special form.
  ## The result are the derivatives at the left- and the right-hand side
  ## boundaries of the training range (of the fitted object with the
  ## current coefficients)
  slopes.left  <- rowsum(lower.slopes * beta[2:(df+1)], group = rep(1, df))
  slopes.right <- rowsum(upper.slopes * beta[2:(df+1)], group = rep(1, df))
  
  ## Now we have to multiply the derivatives with the difference
  ## in the x-values (contained in lowdiff and highdiff)
  ## lowdiff and highdiff are matrices with dimensions n x p, i.e. the
  ## dimension of the newdata object.
  ## Each column of lowdiff and highdiff is multiplied with the slope
  ## value. The result will be what we have to add beyond the boundaries.
  ## add.left and add.right will also have dimension n x p.
  
  ## 'as.array' is here to force a warning message if recycling would
  ## take place (see help file of sweep)
  add.left  <- sweep(lowdiff, MARGIN = 2, STATS = as.array(slopes.left), FUN = "*")
  add.right <- sweep(highdiff, MARGIN = 2, STATS = as.array(slopes.right), FUN = "*")
  
  ## Calculate the final prediction:
  ## Take the prediction of the 'cut-down' matrix and add the linear
  ## extrapolation part (add.left + add.right). We have to take the sum
  ## in each row of the linear extrapolation part (add.left + add.right)
  pred.pre <- pred.pre + rowSums(add.left + add.right)
  return(pred.pre)
}


### Cross_Validation
### Function: Use cross validation to get the best knots for 
###           the semi-parametric model
### Input: DATA: first column is the n by 1 outcome, the rest is 
###        the covariate matrix
###        folds: number of folds for cross validation
###        knots:
### Output: a matrix of MSE in different folds
Cross_Validation <- function(DATA, folds, knots) {
  
  l_knots <- length(knots)
  obs <- length(DATA[,1])
  outcome <- rep(0, folds*l_knots)
  dim(outcome) <- c(folds, l_knots)
  Resample <- sample(obs)
  sub_obs <- floor(obs/folds)
  for (i in 1:(folds - 1))
  {
    sub_test_index <- seq((i - 1) * sub_obs + 1, i * sub_obs)
    sub_test_index <- Resample[sub_test_index] 
    
    sub_training <- DATA[-sub_test_index,]
    sub_testing <- DATA[sub_test_index,]
    knots_inner <- knots * (length(sub_training[,1])^0.8)
    for (j in 1:l_knots)
    {
      knot <-as.integer(round(knots_inner[j]))
      MODEL <- ESTIMATE(sub_training[,-1], sub_training[,1], knot)
      pred.resp <- pred(MODEL, sub_testing[,-1])
      outcome[i, j] <- sum((pred.resp - sub_testing[,1])^2) / length(sub_testing[,1])
      
    }
  }
  sub_train_index <- Resample[1:((folds - 1)*sub_obs)]
  sub_training <- DATA[sub_train_index,]
  sub_testing <- DATA[-sub_train_index,]
  knots_inner <- knots * length(sub_train_index^0.8)
  for (j in 1:l_knots)
  {
    knot <-as.integer(round(knots_inner[j]))
    MODEL <- ESTIMATE(sub_training[,-1], sub_training[,1], knot)
    pred.resp <- pred(MODEL, sub_testing[,-1])
    outcome[folds, j] <- sum((pred.resp - sub_testing[,1])^2) / length(sub_testing[,1])
  }
  return(outcome)
}


### TSCI.basis.fit
### Function: Model fitting part of TSCI using basis approach
### Input: W: covariates with the first column being an instrument
###        D: outcome, the treatment
### folds: number of folds for cross validation
### knots: 
### Output: D.rep: prediction of D
###         sd.D: estimate of noise level in treatment model
###         M:
TSCI.basis.fit <- function(W, D, folds=5, knots=NULL) {
  n <- length(D)
  Data <- cbind(D,W)
  if (is.null(knots)) {
    knots <- unlist(create_knots(1, n, folds, 0.01, 0.1, 20, min_knots = 2, max_knots = 100, num = 10, exe_max = FALSE)) 
  }
  ### cross validation to choose the best knot values
  knot.values<-apply(Cross_Validation(Data, folds, knots),2,mean)
  opt.index<-which.min(knot.values)
  knot<-round(knots[opt.index] * n^0.8)+1
  
  ### model fitting
  MODEL <- ESTIMATE(W, D, knot)
  D.rep<- pred(MODEL, W)
  resid<-D-D.rep
  ### the standard error of D-model
  sd.D=sqrt(sum(resid^2)/length(resid))
  # what is the name of M?
  M <- knot + 2
  
  returnList <- list(D.rep = D.rep,
                     sd.D = sd.D,
                     knot = knot,
                     M = M)
  returnList
}



### TSCI.basis.selection
### Function: Violation space selection for TSCI
###           using basis approach
### IV strength vector: str.vec
### Threshold vector for IV strength: thre.vec
### Proposed point estimator for a given Q: prop.vec.all
### Standard error of the proposed point estimator: se.vec
### Estimation of regression error: noise.vec
### inver.design
TSCI.basis.selection <- function(Y, D, W, D.rep, knot, M, Q=5) {
  Z <- W[,1]; X <- W[,-1]
  sd.D <- sqrt(mean((D-D.rep)^2))
  Q<-min(Q, M-1)

  str.vec<-rep(NA,Q)
  thre.vec<-rep(NA,Q)
  prop.vec.all<-rep(NA,Q)
  se.vec<-rep(NA,Q)
  inver.design<-rep(NA,Q)
  noise.vec<-rep(NA,Q)
  Coef.vec <- sd.vec <- rep(NA,2)
  names(Coef.vec) <- names(sd.vec) <- c("Basis-comp","Basis-robust")
  ####### do the computation from violation space 1 (poly 0) to Q (poly Q-1)
  for(index in 1:Q){
    q=index-1
    if(q==0){
      V=NULL
    }else{
      V=poly(Z,q) # raw=TRUE?
    }
    if(is.null(V)){
      Cov.total<-W[,-1]
    }else{
      V.rep=V
      q=dim(V)[2]
      for(l in 1:q){
        MODEL.V<-ESTIMATE(W,V[,l],knot)
        V.rep[,l] <- pred(MODEL.V, W)
      }
      Cov.total<-cbind(V.rep, W[,-1])
    }
    D.resid<-resid(lm(D.rep~Cov.total))
    str.vec[index]<-sum(D.resid^2)/(sd.D^2)
    taun<-1/log(M-q)
    thre.vec[index]<-(M-q)+sqrt((M-q)*log(M-q))*(1+2*max(taun,sqrt(str.vec[index]/(M-q))))
    thre.vec[index]<-(sd.D^2)*thre.vec[index]
    MODEL.Y <- ESTIMATE(W, Y, knot)
    Y.rep<- pred(MODEL.Y, W)
    #D.resid<-resid(lm(D.rep~Cov.total,-1))
    Y.resid<-resid(lm(Y.rep~Cov.total))
    ### estimate the point estimator
    prop.vec.all[index]<-sum(Y.resid*D.resid)/sum(D.resid^2)
    ### estimate the standard error
    D.res<-D-D.rep
    Y.res<-Y-pred(MODEL.Y, W)
    sigsq.hat<-mean((Y.res-prop.vec.all[index]*D.res)^2)
    noise.vec[index]<-sigsq.hat
    inver.design[index]<-1/sum(D.resid^2)
    scale<-1
    se.vec[index]<-scale*sqrt(sigsq.hat/sum(D.resid^2))
  }
  ### selection
  ### test iv strength of order q polynomial in TSCI
  str.iv.set<-str.vec >= thre.vec
  # str.iv.set<-(str.vec>thre.vec)*(str.vec>10)
  if (sum(str.iv.set)==0) {
    warning("Weak IV: Even if the IV is assumed to be valid, run OLS") # stop, output results of OLS
    Q.max <- 0
  } else {
    Q.max <- sum(str.iv.set)-1
    if (Q.max==0) {
      ### if Q.max==0, redefine Qmax by log(log(n)), ignore at this stage
      warning("Weak IV: IV Strenght Test Failed. Set Qmax as 1.")
      Q.max <- 1
    }
  }
  prop.vec<-as.vector(prop.vec.all[1:(Q.max+1)])
  inver.design<-as.vector(inver.design[1:(Q.max+1)])
  ### conduct the selection
  test.threshold <- test.diff <- matrix(0,Q.max,Q.max)
  for (q in 0:(Q.max-1)) {
    test.diff[q+1,(q+1):(Q.max)] <- abs(prop.vec[q+1]-prop.vec[(q+2):(Q.max+1)])
    H.vec <- inver.design[(q+2):(Q.max+1)] - inver.design[q+1]
    test.threshold[q+1,(q+1):(Q.max)] <- qnorm(1-0.05/log(M))*sqrt(noise.vec[Q.max+1])*sqrt(H.vec)
  }
  C.alpha <- ifelse(test.diff<=test.threshold,0,1)
  # layer selection
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
  
  Coef.vec[1] <- prop.vec[q.comp+1]
  Coef.vec[2] <- prop.vec[q.robust+1]
  sd.vec[1] <- se.vec[q.comp+1]
  sd.vec[2] <- se.vec[q.robust+1]
  
  
  returnList <- list(prop.vec.all = prop.vec.all,
                     se.vec = se.vec,
                     inver.design = inver.design,
                     noise.vec = noise.vec,
                     SigmaSqY.Qmax = noise.vec[Q.max+1],
                     str.vec = str.vec,
                     thre.vec = thre.vec,
                     Q.max = Q.max,
                     q.comp = q.comp,
                     q.robust = q.robust,
                     Coef.vec = Coef.vec,
                     sd.vec = sd.vec,
                     validity = validity)
  returnList
}

