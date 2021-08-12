#' @title TSCI with Basis Approach
#' @description Implement basis approach using data splitting and do inference with violation space selection
#'
#' @param Y continuous, n by 1 outcome vector.
#' @param D continuous or binary, n by 1 treatment vector.
#' @param Z continuous or binary, n by 1 Instrumental Variable, only one instrument is implemented for violation space selection now.
#' @param X continuous or binary, n by p_x baseline covariates matrix.
#' @param vio.space n by Q-1 matrix, each column refers to a violation form of Z, default by NULL assumes linear and quadratic violation.
#' @param layer logic, do layer selection of violation space or not, default by TRUE.
#' @param intercept logic, whether to include intercept or not in the outcome model, default by TRUE.
#' @param str.thol minimal value of the threshold of IV strength test, default by 20.
#'
#' @return
#'     \item{\code{Coef.vec}}{a vector of length 2*Q, the point and bias corrected estimators with fixed violation spaces.}
#'     \item{\code{sd.vec}}{estimated standard deviation of Coef.vec.}
#'     \item{\code{Coef.robust}}{a vector of length 4, the point and bias corrected estimators in the violation space using comparison method and robust method.}
#'     \item{\code{sd.robust}}{estimated standard deviation of Coef.robust.}
#'     \item{\code{SigmaSqY}}{estimated noise level of outcome model.}
#'     \item{\code{SigmaSqD}}{estimated noise level of treatment model.}
#'     \item{\code{SigmaSqY.Qmax}}{noise level of outcome model in violation space Q.max.}
#'     \item{\code{iv.str}}{a vector of length Q, the IV strength.}
#'     \item{\code{iv.thol}}{a vector of length Q,, IV strength threshold.}
#'     \item{\code{Q.max}}{maximal violation space chosen by IV strength test.}
#'     \item{\code{q.comp}}{estimated violation space using comparison method.}
#'     \item{\code{q.robust}}{estimated violation space using robust method.}
#'     \item{\code{invalidity}}{invalidity of TSLS, TRUE means there is at least linear violation.}
#'     \item{\code{run.OLS}}{ogic, whether the IV is weak when assuming linearly invalid, IV strength test only passes with no violation.}
#'     \item{\code{weak.iv}}{ogic, whether the IV is weak when assuming no violation, IV strength test does not pass even with no violation We need to run OLS if TRUE since TSLS result is not reliable.}
#' @export
#'
#' @examples
#' \donttest{
#' # function to generate covariance matrix
#' A1gen<-function(rho,p){
#'   A1=matrix(0,p,p)
#'   for(i in 1:p){
#'     for(j in 1:p){
#'       A1[i,j]<-rho^(abs(i-j))
#'     }
#'   }
#'   A1
#' }
#'
#'
#' ### dimension change the dimension 10,20
#' p = 10
#' ### please change this n = 1000, 3000, 5000
#' n = 1000
#' ### setting, change across 1(for previous 1 and 4), 2(for previous 3 and 6)
#' f.index = 2
#' ### change the interaction 0, 0.5, 1, 1.5
#' inter.val = 0
#' ### a denotes the IV strength, set as 1
#' a = 1
#' ### violation index, set as 1 or 2
#' vio.index = 1
#' ### tau denotes the violation strength
#' ### set tau as 1
#' tau = 1
#'
#'
#' ##############################
#' f_1 <- function(x){x+a*(x^2+0.5*x^4) -25/12}
#' # f_2 <- function(x){exp(2*a+x+0.5*x^3)+a*(x^2+0.5*x^4)}
#' f_2 <- function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
#' rho1=0.5
#' Cov<-(A1gen(rho1,p+1))
#' mu<-rep(0,p+1)
#' beta=1
#' alpha<-as.matrix(rep(-0.3,p))
#' gamma<-as.matrix(rep(0.2,p))
#' inter<-as.matrix(c(rep(inter.val,5),rep(0,p-5)))
#'
#'
#' #### generate the data
#' mu.error<-rep(0,2)
#' Cov.error<-matrix(c(1,0.5,0.5,1),2,2)
#' Error<-mvrnorm(n, mu.error, Cov.error)
#' W.original<-mvrnorm(n, mu, Cov)
#' W<-pnorm(W.original)
#' Z<-W[,1]
#' X<-W[,-1]
#' ###### generate the data for the treatment variable D
#' if(f.index==1){
#'   D=f_1(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
#' }
#' if(f.index==2){
#'   D=f_2(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
#' }
#' ####### generate the outcome variable
#' # if(vio.index==0){
#' #   Y=D*beta+ X%*% gamma+Error[,2]
#' # }
#' if(vio.index==1){
#'   Y=D*beta+ tau*Z+ X%*% gamma+Error[,2]
#' }
#' if(vio.index==2){
#'   Y=D*beta+ tau*(Z^2+Z-1)+ X%*% gamma+Error[,2] # difficult if no Z
#' }
#'
#'
#' ### basis
#' output.basis <- TSCI.basis(Y,D,Z,X)
#' # point estimate
#' output.basis$Coef.robust["Basis-robust"]
#' # standard error
#' output.basis$sd.robust["Basis-robust"]
#' }
#'
#'
TSCI.basis <- function(Y,D,Z,X,vio.space=NULL,layer=TRUE,intercept=TRUE,str.thol=20) {
  Y <- as.matrix(Y); D <- as.matrix(D); Z <- as.matrix(Z); X <- as.matrix(X);
  W <- as.matrix(cbind(Z,X))
  # constants
  n <- NROW(X); p <- NCOL(X)
  # define the augmentation of covariates,
  # which is the combination of violation space and baseline covariates
  if (is.null(vio.space)) {
    Q = 3
    vio.space <- poly(Z,degree = Q-1,raw = TRUE)
  } else {
    Q = NCOL(vio.space) + 1
  }
  Cov.aug <- cbind(vio.space,X)

  # model fitting
  basis.fit <- TSCI.basis.fit(D,W)
  D.rep <- basis.fit$D.rep
  M <- basis.fit$M
  knot <- basis.fit$knot

  outputs <- TSCI.basis.Selection(Y, D, W, D.rep, knot, M, Q, vio.space, intercept, str.thol)
  return(outputs)
}


#' Title
#'
#' @param D continuous or binary, n by 1 treatment vector.
#' @param W continuous or binary, instrument-covariates matrix.
#' @param folds number of folds for cross validation,default by 5.
#'
#' @return
#'     \item{\code{D.rep}}{representation of D.}
#'     \item{\code{SigmaSqD}}{estimated noise level in treatment model.}
#'     \item{\code{knot}}{number of knots.}
#'     \item{\code{M}}{number of basis functions, equals to #knot + 2.}
#' @export
#'
#' @examples
#' \donttest{
#'
#' }
TSCI.basis.fit <- function(D, W, folds=5) {
  n <- length(D)
  Data <- cbind(D,W)
  knots <- unlist(create_knots(1, n, folds, 0.01, 0.1, 20, min_knots = 2, max_knots = 100, num = 10, exe_max = FALSE))
  # cross validation to choose the best knot values
  knot.values<-apply(Cross_Validation(Data, folds, knots),2,mean)
  opt.index<-which.min(knot.values)
  knot<-round(knots[opt.index] * n^0.8)+1

  ### model fitting
  MODEL <- ESTIMATE(W, D, knot)
  D.rep<- pred(MODEL, W)
  resid<-D-D.rep
  # the standard error of D-model
  SigmaSqD=sum(resid^2)/length(resid)
  M <- knot + 2

  returnList <- list(D.rep = D.rep,
                     SigmaSqD = SigmaSqD,
                     knot = knot,
                     M = M)
  returnList
}


###### helper functions

### TSCI.basis.selection
### Function: Violation space selection for TSCI
###           using basis approach
### Input: Y: continuous, n by 1 outcome vector
###        D: continuous or binary, n by 1 treatment vector
###        W: Instrument-covariates matrix
###        D.rep: Representation of D
###        knot: Number of knots
###        M: Number of basis functions, equals to #knot+2
###        vio.space: n by (Q-1) matrix, each column refers to a violation form of Z,
###                   default by NULL assumes linear and quadratic violation (Z^2,Z),
###                   choose between quadratic, linear and no violation(Q=3)
###        intercept: Include intercept in the outcome model or not
###        str.thol: the minimal value of the threshold of IV strength test
### Output: Coef.vec: A vector of length 2*Q, the point and bias corrected estimators with fixed violation spaces
###         sd.vec: Estimated standard deviation of Coef.vec
###         Coef.robust: A vector of length 4, the point and bias corrected estimators in the violation space using
###                      comparison method and robust method
###         sd.robust: Estimated standard deviation of Coef.robust
###         SigmaSqY: Estimated noise level of outcome model
###         SigmaSqD: Estimated noise level of treatment model
###         SigmaSqY.Qmax: Noise level of outcome model in violation space Q.max
###         iv.str: A vector of length Q, the left hand side of IV strength test
###         iv.thol: A vector of length Q,, the right hand side of IV strength test
###         Q.max; Maximal violation space chosen by IV strength test
###         q.comp: Violation space using comparison method
###         q.robust: Violation space using robust method
###         invalidity: Invalidity of TSLS, means there is at least linear violation
###         weak.iv: logic, whether the IV is weak when assuming linearly invalid, IV strength test only passes with no violation
###         run.OLS: logic, whether the IV is weak when assuming no violation, IV strength test does not pass even with no violation
###                         We need to run OLS if TRUE since TSLS result is not reliable
TSCI.basis.Selection <- function(Y, D, W, D.rep, knot, M, Q, vio.space, intercept, str.thol) {
  Y <- as.matrix(Y); D <- as.matrix(D);
  W <- as.matrix(W); D.rep <- as.matrix(D.rep)
  Z <- W[,1]
  n <- length(Y)
  SigmaSqD <- mean((D-D.rep)^2)
  Q<-min(Q, M-1)

  iv.str<-rep(NA,Q)
  iv.thol<-rep(NA,Q)
  Coef.vec<-rep(NA,Q)
  sd.vec<-rep(NA,Q)
  inver.design<-rep(NA,Q)
  SigmaSqY<-rep(NA,Q)
  Coef.robust <- sd.robust <- rep(NA,2)
  names(Coef.robust) <- names(sd.robust) <- c("Basis-comp","Basis-robust")
  D.resid.list <- rep(list(NA),Q)
  ####### do the computation from violation space 1 (poly 0) to Q (poly Q-1)
  for(index in 1:Q){
    q=index-1
    if (is.null(vio.space)) {
      if(q==0){
        V=NULL
      }else{
        V=poly(Z,q,raw = TRUE) # raw=TRUE?
      }
      if(is.null(V)){
        Cov.total <- W[,-1]
      }else{
        V.rep=V
        q=dim(V)[2]
        for(l in 1:q){
          MODEL.V<-ESTIMATE(W,V[,l],knot)
          V.rep[,l] <- pred(MODEL.V, W)
        }
        Cov.total<-cbind(V.rep, W[,-1])
      }
    } else {
      if (q==0) {
        Cov.total <- W[,-1]
      } else if (q==Q-1) {
        V <- vio.space
        V.rep <- V
        for (l in 1:q) {
          MODEL.V<-ESTIMATE(W,V[,l],knot)
          V.rep[,l] <- pred(MODEL.V, W)
        }
        Cov.total<-cbind(V.rep, W[,-1])
      } else {
        V <- as.matrix(vio.space[,-(1:(Q-1-q))])
        V.rep <- V
        for (l in 1:q) {
          MODEL.V<-ESTIMATE(W,V[,l],knot)
          V.rep[,l] <- pred(MODEL.V, W)
        }
        Cov.total<-cbind(V.rep, W[,-1])
      }
    }


    D.resid<-resid(lm(D.rep~Cov.total))
    D.resid.list[[index]] <- D.resid
    iv.str[index]<-sum(D.resid^2)/(SigmaSqD)
    taun<-1/log(M-q)
    boot.vec <- rep(NA,300)
    for (i in 1:300) {
      delta <- rnorm(n,0,sqrt(SigmaSqD))
      MODEL.delta <- ESTIMATE(W,delta,knot)
      delta.rep <- pred(MODEL.delta,W)
      delta.resid <- resid(lm(delta.rep~Cov.total))
      boot.vec[i] <- sum(delta.resid^2) + 2*sum(D.resid*delta.resid) # abs of second term?
    }
    iv.thol[index] <- max(quantile(boot.vec,0.975),20)/(SigmaSqD)
    # iv.thol[index] <- quantile(boot.vec,0.975)
    MODEL.Y <- ESTIMATE(W, Y, knot)
    Y.rep<- pred(MODEL.Y, W)
    #D.resid<-resid(lm(D.rep~Cov.total,-1))
    if (intercept) {
      Y.resid<-resid(lm(Y.rep~Cov.total))
    } else {
      Y.resid<-resid(lm(Y.rep~Cov.total-1))
    }
    ### estimate the point estimator
    Coef.vec[index]<-sum(Y.resid*D.resid)/sum(D.resid^2)
    ### estimate the standard error
    D.res<-D-D.rep
    Y.res<-Y-pred(MODEL.Y, W)

    SigmaSqY[index]<-mean((Y.res-Coef.vec[index]*D.res)^2)
    inver.design[index]<-1/sum(D.resid^2)
    scale<-1
    sd.vec[index]<-scale*sqrt(SigmaSqY[index]/sum(D.resid^2))
  } # index
  D.RSS <- iv.str*SigmaSqD # residual sum of squares of D-model

  ### selection
  ### test iv strength of order q polynomial in TSCI
  ivtest.vec<-iv.str >= iv.thol
  run.OLS <- weak.iv <- FALSE
  if (sum(ivtest.vec)==0) {
    warning("Weak IV: Even if the IV is assumed to be valid, run OLS") # stop, output results of OLS
    Q.max <- 1
    run.OLS <- TRUE
  } else {
    Q.max <- sum(ivtest.vec)-1
    if (Q.max==0) {
      ### if Q.max==0, redefine Qmax by log(log(n)), ignore at this stage
      warning("Weak IV, if the IV is invalid. We still test the IV invalidity.")
      weak.iv <- TRUE
      Q.max <- 1
    }
  }
  Coef.vec.Qmax<-as.vector(Coef.vec[1:(Q.max+1)])
  inver.design<-as.vector(inver.design[1:(Q.max+1)])
  ### conduct the selection
  diff.thol <- beta.diff <- matrix(0,Q.max,Q.max)
  for (q in 0:(Q.max-1)) {
    beta.diff[q+1,(q+1):(Q.max)] <- abs(Coef.vec.Qmax[q+1]-Coef.vec.Qmax[(q+2):(Q.max+1)])
    H.vec <- inver.design[(q+2):(Q.max+1)] - inver.design[q+1]
    diff.thol[q+1,(q+1):(Q.max)] <- sqrt(SigmaSqY[Q.max+1])*sqrt(H.vec) # need to multiply by z.alpha
  }


  ### bootstrap to get z.alpha
  max.val <- rep(NA,300)
  for (i in 1:300) {
    diff.mat <- matrix(0,Q.max,Q.max)
    eps <- rnorm(n, 0, sqrt(SigmaSqY[Q.max+1]))
    MODEL.eps <- ESTIMATE(W,eps,knot)
    eps.rep <- pred(MODEL.eps,W)
    for (q1 in 0:(Q.max-1)) {
      for (q2 in (q1+1):(Q.max)) {
        diff.mat[q1+1, q2] <- sum(D.resid.list[[q2+1]]*eps.rep)/D.RSS[q2+1]-sum(D.resid.list[[q1+1]]*eps.rep)/D.RSS[q1+1]
      }
    }
    diff.mat <- abs(diff.mat)/diff.thol
    max.val[i] <- max(diff.mat,na.rm = TRUE)
  }
  # z.alpha <- quantile(max.val,0.99)
  z.alpha <- quantile(max.val,0.975)
  diff.thol <- z.alpha*diff.thol


  C.alpha <- ifelse(beta.diff<=diff.thol,0,1)
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

  Coef.robust[1] <- Coef.vec.Qmax[q.comp+1]
  Coef.robust[2] <- Coef.vec.Qmax[q.robust+1]
  sd.robust[1] <- sd.vec[q.comp+1]
  sd.robust[2] <- sd.vec[q.robust+1]


  returnList <- list(Coef.vec = Coef.vec,
                     sd.vec = sd.vec,
                     Coef.robust = Coef.robust,
                     sd.robust = sd.robust,
                     SigmaSqY = SigmaSqY,
                     SigmaSqD = SigmaSqD,
                     SigmaSqY.Qmax = SigmaSqY[Q.max+1],
                     iv.str = iv.str, iv.thol = iv.thol,
                     Q.max = Q.max, q.comp = q.comp, q.robust = q.robust,
                     validity = validity,
                     run.OLS = run.OLS,
                     weak.iv = weak.iv)
  returnList
}


### get.design
### Function: obtain the design matrix for the first column
###           of covaraites, prepared for semi-parametric model
### Input: x: n by p covariate matrix
###        knots: Number of knots
### Output: m: Design matrix
###         basis: Basis object of package fda
getDesign <- function(x, knots) {
  p <- NCOL(x)
  if(is.null(colnames(x))) colnames(x) <- paste("x", 1:p, sep = "")
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
### Input: P: Number of covariates
###        N: Sample size
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
      if (max_knots >= n) max_knots <- round(n * p_max)
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
### Input: X: n by p covariates matrix
###        Y: n by 1 outcome
###        knots: Knots created by create_knots function
### Output: m: Design matrix
###         basis:Basis object of package fda
###         coefs: Estimated coefficients of the model
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
### Input: object: An object obtained by ESTIMATE()
###        newdata: New covariates to use
### Output: pred.pre: Predicted representation of outcome
pred <- function(object, newdata) {
  x <- newdata
  n <- NROW(x)
  p <- NCOL(x)
  lowdiff <- highdiff <- matrix(FALSE, nrow = NROW(newdata), ncol = 1)
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
### Input: DATA: First column is the n by 1 outcome, the rest is
###        the covariate matrix
###        folds: Number of folds for cross validation
###        knots: Knots created by create_knots function
### Output: A matrix of MSE in different folds
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
