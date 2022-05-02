#' @title Two Stage Curvature Identification with Basis Approach
#' @description This functions implements Two Stage Curvature Identification with the Basis Approach. It tests the IV strength and chooses the violation space among a series of candidate spaces, and also constructs the confidence interval for the treatment effect with the selected violation space.
#'
#' @param Y outcome with dimension n by 1
#' @param D treatment with dimension n by 1
#' @param Z instrument variable with dimension n by 1
#' @param X baseline covariates with dimension n by p
#' @param vio.space a matrix with each column corresponding to a violation form of Z; If NULL, then default by the n by 3 matrix (Z^3, Z^2, Z). Violation space selection will be performed according to provided violation space, for example, null violation space vs Z vs (Z^2, Z) vs (Z^3, Z^2, Z) in the default case
#' @param intercept logic, including the intercept or not in the outcome model, default by TRUE
#' @param str.thol the minimal value of the threshold of IV strength test, default by 10
#' @param alpha the significance level, default by 0.05
#'
#' @return
#'     \item{\code{Coef.all}}{a series of point estimators of treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{sd.all}}{standard errors of Coef.all}
#'     \item{\code{CI.all}}{confidence intervals for the treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{Coef.robust}}{the point estimators corresponding to the violation space selected by the robust comparison}
#'     \item{\code{sd.robust}}{the standard errors of Coef.robust}
#'     \item{\code{CI.robust}}{confidence intervals for the treatment effect with the violation space selected by the robust comparison}
#'     \item{\code{iv.str}}{IV strength corresponding to different violation spaces}
#'     \item{\code{iv.thol}}{the threshold of IV strength test corresponding to different violation spaces}
#'     \item{\code{Qmax}}{the index of largest violation space selected by IV strength test. If -1, the IV strength test fails for null violation space and run OLS. If 0, the IV Strength test fails for the first violation space and run TSBA only for null violation space. In other cases, violation space selection is performed}
#'     \item{\code{q.hat}}{the index of estimated violation space corresponding to Qmax}
#'     \item{\code{invalidity}}{invalidity of TSLS. If TRUE, the IV is invalid; Otherwise, the IV is valid}
#' @export
#'
#' @examples
#' \dontrun{
#'
#'
#' # dimension
#' p = 10
#' #
#' n = 100
#' # interaction value
#' inter.val = 1
#' # the IV strength
#' a = 1
#' # violation strength
#' tau = 1
#' f = function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
#' rho1=0.5
#' # function to generate covariance matrix
#' A1gen=function(rho,p){
#'   A1=matrix(0,p,p)
#'   for(i in 1:p){
#'     for(j in 1:p){
#'       A1[i,j]=rho^(abs(i-j))
#'     }
#'   }
#'   A1
#' }
#' Cov=(A1gen(rho1,p+1))
#' mu=rep(0,p+1)
#' # true effect
#' beta=1
#' alpha=as.matrix(rep(-0.3,p))
#' gamma=as.matrix(rep(0.2,p))
#' inter=as.matrix(c(rep(inter.val,5),rep(0,p-5)))
#'
#'
#' # generate the data
#' mu.error=rep(0,2)
#' Cov.error=matrix(c(1,0.5,0.5,1),2,2)
#' Error=mvrnorm(n, mu.error, Cov.error)
#' W.original=mvrnorm(n, mu, Cov)
#' W=pnorm(W.original)
#' Z=W[,1]
#' X=W[,-1]
#' # generate the treatment variable D
#' D=f(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
#' # generate the outcome variable Y
#' Y=D*beta+ tau*Z+ X%*% gamma+Error[,2]
#'
#'
#' # Two Stage Basis Approach
#' output.BA = TSBA(Y,D,Z,X)
#' # point estimates
#' output.BA$Coef.robust
#' # standard errors
#' output.BA$sd.robust
#' # confidence intervals
#' output.BA$CI.robust
#' }
#'
#' @importFrom fda create.bspline.basis eval.basis
#' @importFrom stats coef lm poly predict qnorm quantile resid rnorm
#'
TSBA = function(Y,D,Z,X,vio.space=NULL,intercept=TRUE,str.thol=10,alpha=0.05) {
  Y = as.matrix(Y); D = as.matrix(D); Z = as.matrix(Z); X = as.matrix(X);
  # constants
  n = nrow(X); p = ncol(X)
  # define the augmentation of covariates, which is the combination of violation space and baseline covariates
  if (is.null(vio.space)) {
    Q = 4
    vio.space = matrix(NA,nrow(Z),0)
    for (q in 1:(Q-1)) {
      vio.space = cbind(Z^q,vio.space)
    }
  } else {
    Q = ncol(vio.space) + 1
  }
  Cov.aug = cbind(vio.space,X)

  # model fitting
  basis.fit = TSBA.fit(D,Z,X)
  D.rep = basis.fit$D.rep
  M = basis.fit$M
  knot = basis.fit$knot

  outputs = TSBA.Selection(Y,D,Z,X,vio.space,Q,D.rep,knot,M,intercept,str.thol,alpha)
  return(outputs)
}



#' @title Regression with B-spline basis
#' @description Implement regression for treatment model with B-spline basis
#'
#' @param D treatment with dimension n by 1
#' @param Z instrument variable with dimension n by 1
#' @param X baseline covariates with dimension n by p
#' @param folds number of folds for cross validation
#'
#' @return:
#'     \item{\code{D.rep}}{transformed treatment of dimension n_A1 by 1 corresponding to a violation space}
#'     \item{\code{knot}}{number of knots for B-spline regression}
#'     \item{\code{M}}{number of basis functions, equals number of basis + 2}
#' @noRd
#'
TSBA.fit = function(D,Z,X,folds=5) {
  n = length(D)
  W = cbind(Z,X)
  Data = cbind(D,W)
  knots = unlist(create_knots(1, n, folds, 0.01, 0.1, 20, min_knots = 2, max_knots = 100, num = 10, exe_max = FALSE))
  # cross validation to choose the best knot values
  knot.values=apply(Cross_Validation(Data, folds, knots),2,mean)
  opt.index=which.min(knot.values)
  knot=round(knots[opt.index] * n^0.8)+1

  ### model fitting
  MODEL = ESTIMATE(W, D, knot)
  D.rep= pred(MODEL, W)
  resid=D-D.rep
  M = knot + 2

  returnList = list(D.rep = D.rep,
                     knot = knot,
                     M = M)
  returnList
}



#' @title Violation Space Selection of Basis Approach
#' @description Select violation space for Basis Approach and construct confidence intervals
#'
#' @param Y outcome with dimension n by 1
#' @param D treatment with dimension n by 1
#' @param Z instrument variable with dimension n by 1
#' @param X baseline covariates with dimension n by p
#' @param vio.space a matrix with each column corresponding to a violation form of Z
#' @param Q the number of violation spaces, including the null violation space, but the final possible violation space is defined as min(Q,M-1)
#' @param D.rep transformed treatment of dimension n_A1 by 1 corresponding to a violation space
#' @param knot the number of knots for B-spline regression
#' @param M the number of basis functions, equals number of basis + 2
#' @param intercept logic, including the intercept or not in the outcome model
#' @param str.thol the minimal value of the threshold of IV strength test
#' @param alpha the significance level
#'
#' @return
#'     \item{\code{Coef.all}}{a series of point estimators of treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{sd.all}}{standard errors of Coef.all}
#'     \item{\code{CI.all}}{confidence intervals for the treatment effect corresponding to different violation spaces and the OLS}
#'     \item{\code{Coef.robust}}{the point estimators corresponding to the violation space selected by the robust comparison}
#'     \item{\code{sd.robust}}{the standard errors of Coef.robust}
#'     \item{\code{CI.robust}}{confidence intervals for the treatment effect with the violation space selected by the robust comparison}
#'     \item{\code{iv.str}}{IV strength corresponding to different violation spaces}
#'     \item{\code{iv.thol}}{the threshold of IV strength test corresponding to different violation spaces}
#'     \item{\code{Qmax}}{the index of largest violation space selected by IV strength test. If -1, the IV strength test fails for null violation space and run OLS. If 0, the IV Strength test fails for the first violation space and run TSBA only for null violation space. In other cases, violation space selection is performed}
#'     \item{\code{q.hat}}{the index of estimated violation space corresponding to Qmax}
#'     \item{\code{invalidity}}{invalidity of TSLS. If TRUE, the IV is invalid; Otherwise, the IV is valid}
#' @noRd
#'
TSBA.Selection = function(Y,D,Z,X,vio.space,Q,D.rep,knot,M,intercept,str.thol,alpha) {
  Y = as.matrix(Y); D = as.matrix(D);
  W = as.matrix(cbind(Z,X)); D.rep = as.matrix(D.rep)
  n = length(Y)
  SigmaSqD = mean((D-D.rep)^2)
  Q=min(Q, M-1)
  iv.str=rep(NA,Q)
  iv.thol=rep(NA,Q)
  Coef.all=rep(NA,Q)
  sd.all=rep(NA,Q)
  Coef.names = paste("TSBA-q",0:(Q-1),sep="")
  names(Coef.all) = names(sd.all) = Coef.names

  inver.design=rep(NA,Q)
  SigmaSqY=rep(NA,Q)
  Coef.robust = sd.robust = rep(NA,2)
  names(Coef.robust) = names(sd.robust) = c("TSBA-comp","TSBA-robust")
  D.resid.list = rep(list(NA),Q)
  ### do the computation from violation space 1 (poly 0) to Q (poly Q-1)
  for(index in 1:Q) {
    q=index-1
    if (q==0) {
      Cov.total = W[,-1]
    } else if (q==Q-1) {
      V = vio.space
      V.rep = V
      for (l in 1:q) {
        MODEL.V=ESTIMATE(W,V[,l],knot)
        V.rep[,l] = pred(MODEL.V, W)
      }
      Cov.total=cbind(V.rep, W[,-1])
    } else {
      V = as.matrix(vio.space[,-(1:(Q-1-q))])
      V.rep = V
      for (l in 1:q) {
        MODEL.V=ESTIMATE(W,V[,l],knot)
        V.rep[,l] = pred(MODEL.V, W)
      }
      Cov.total=cbind(V.rep, W[,-1])
    }


    D.resid=resid(lm(D.rep~Cov.total))
    D.resid.list[[index]] = D.resid
    iv.str[index]=sum(D.resid^2)/(SigmaSqD)
    boot.vec = rep(NA,300)
    for (i in 1:300) {
      delta = rnorm(n,0,sqrt(SigmaSqD))
      MODEL.delta = ESTIMATE(W,delta,knot)
      delta.rep = pred(MODEL.delta,W)
      delta.resid = resid(lm(delta.rep~Cov.total))
      boot.vec[i] = sum(delta.resid^2) + 2*sum(D.resid*delta.resid)
    }
    iv.thol[index] = quantile(boot.vec,0.975)/SigmaSqD + max(2*(M-q),str.thol)
    # iv.thol[index] = quantile(boot.vec,0.975)
    MODEL.Y = ESTIMATE(W, Y, knot)
    Y.rep= pred(MODEL.Y, W)
    if (intercept) {
      Y.resid=resid(lm(Y.rep~Cov.total))
    } else {
      Y.resid=resid(lm(Y.rep~Cov.total-1))
    }
    ### point estimator
    Coef.all[index]=sum(Y.resid*D.resid)/sum(D.resid^2)
    ### the standard error
    D.res=D-D.rep
    Y.res=Y-pred(MODEL.Y, W)

    SigmaSqY[index]=mean((Y.res-Coef.all[index]*D.res)^2)
    inver.design[index]=1/sum(D.resid^2)
    scale=1
    sd.all[index]=scale*sqrt(SigmaSqY[index]/sum(D.resid^2))
  } # index
  D.RSS = iv.str*SigmaSqD # residual sum of squares of the treatment model


  ### iv strength test
  ivtest.vec=iv.str >= iv.thol
  if (sum(ivtest.vec)==0) {
    warning("Weak IV: Even if the IV is assumed to be valid, run OLS") # stop, output results of OLS
    Qmax = -1
  } else {
    Qmax = sum(ivtest.vec)-1
    if (Qmax==0) {
      warning("Weak IV, if the IV is invalid. We still test the IV invalidity.")
    }
  }


  if (Qmax>=1) { # selection
    Coef.all.Qmax=as.vector(Coef.all[1:(Qmax+1)])
    inver.design=as.vector(inver.design[1:(Qmax+1)])
    ### conduct the selection
    diff.thol = beta.diff = matrix(0,Qmax,Qmax)
    for (q in 0:(Qmax-1)) {
      beta.diff[q+1,(q+1):(Qmax)] = abs(Coef.all.Qmax[q+1]-Coef.all.Qmax[(q+2):(Qmax+1)])
      H.vec = inver.design[(q+2):(Qmax+1)] - inver.design[q+1]
      diff.thol[q+1,(q+1):(Qmax)] = sqrt(SigmaSqY[Qmax+1])*sqrt(H.vec) # need to multiply by z.alpha
    }


    ### bootstrap to get z.alpha
    max.val = rep(NA,300)
    for (i in 1:300) {
      diff.mat = matrix(0,Qmax,Qmax)
      eps = rnorm(n, 0, sqrt(SigmaSqY[Qmax+1]))
      MODEL.eps = ESTIMATE(W,eps,knot)
      eps.rep = pred(MODEL.eps,W)
      for (q1 in 0:(Qmax-1)) {
        for (q2 in (q1+1):(Qmax)) {
          diff.mat[q1+1, q2] = sum(D.resid.list[[q2+1]]*eps.rep)/D.RSS[q2+1]-sum(D.resid.list[[q1+1]]*eps.rep)/D.RSS[q1+1]
        }
      }
      diff.mat = abs(diff.mat)/diff.thol
      max.val[i] = max(diff.mat,na.rm = TRUE)
    }
    z.alpha = quantile(max.val,0.975)
    diff.thol = z.alpha*diff.thol


    C.alpha = ifelse(beta.diff<=diff.thol,0,1)
    # layer selection
    sel.vec = apply(C.alpha,1,sum)
    if (all(sel.vec != 0)) {
      q.comp = Qmax
    } else {
      q.comp = min(which(sel.vec==0))-1
    }
  } # selection


  ### invalidity of TSLS
  if (Qmax>=1) {
    if (q.comp>=1) {
      invalidity = 1
    } else {
      invalidity = 0
    }
  } else {
    invalidity = 0
  }


  # OLS estimator
  OLS = summary(lm(Y~D+X))$coefficients
  Coef.OLS = OLS[2,1]
  sd.OLS = OLS[2,2]
  CI.OLS = c(Coef.OLS + qnorm(alpha/2)*sd.OLS,Coef.OLS + qnorm(1-alpha/2)*sd.OLS)
  names(CI.OLS) = c("lower","upper")
  # add OLS to Coef.all
  Coef.all = c(Coef.all,Coef.OLS)
  names(Coef.all) = c(Coef.names,"OLS")
  sd.all = c(sd.all,sd.OLS)
  names(sd.all) = c(Coef.names,"OLS")


  # confidence intervals
  CI.all = rbind(Coef.all + qnorm(alpha/2)*sd.all,Coef.all + qnorm(1-alpha/2)*sd.all)
  rownames(CI.all) = c("lower","upper")


  if (Qmax>=1) {
    q.robust = min(q.comp+1, Qmax)
    q.hat = c(q.comp,q.robust)
    names(q.hat) = c("comp","robust")
    Coef.robust[1] = Coef.all[q.comp+1]
    Coef.robust[2] = Coef.all[q.robust+1]
    sd.robust[1] = sd.all[q.comp+1]
    sd.robust[2] = sd.all[q.robust+1]
  }
  if (Qmax==0) {
    if (Qmax==0) {
      q.hat = c(0)
      names(q.hat) = "q0"
      Coef.robust = Coef.all[1]
      sd.robust = sd.all[1]
    }
    if(Qmax==-1) {
      q.hat = c(-1)
      names(q.hat) = "OLS"
      Coef.robust = Coef.all[Q+1]
      sd.robust = sd.all[Q+1]
    }
  }

  # confidence intervals
  CI.robust = rbind(Coef.robust + qnorm(alpha/2)*sd.robust,Coef.robust + qnorm(1-alpha/2)*sd.robust)
  rownames(CI.robust) = c("lower","upper")


  returnList = list(Coef.all = Coef.all,
                     sd.all = sd.all,
                     CI.all = CI.all,
                     Coef.robust = Coef.robust,
                     sd.robust = sd.robust,
                     CI.robust = CI.robust,
                     iv.str = iv.str, iv.thol = iv.thol,
                     Qmax = Qmax, q.hat = q.hat,
                     invalidity = invalidity)
  returnList
}


##### helper functions for basis approach #####

### get.design
### Function: obtain the design matrix for the first column of covaraites, prepared for semi-parametric model
### Input: x: n by p covariate matrix
###        knots: Number of knots
### Output: m: Design matrix
###         basis: Basis object of package fda
getDesign = function(x, knots) {
  p = ncol(x)
  if(is.null(colnames(x))) colnames(x) = paste("x", 1:p, sep = "")
  x.current = x[,1]
  ux.current = unique(x.current)
  knots.use = quantile(ux.current, seq(0, 1, length = knots))
  stopifnot(all.equal(range(knots.use), range(x.current)))
  basis = create.bspline.basis(rangeval = range(knots.use),
                                breaks = knots.use, norder = 4)
  m = eval.basis(x.current, basis)
  list(m = m, basis = basis)
}


### create_knots
### Function: A helper function to create knots
create_knots = function(P, N, folds, begin, end, l, min_knots = 2, max_knots = 50, num = 10, exe_max = FALSE, p_max = 0.5) {
  l_P = length(P)
  l_N = length(N)
  knots = NULL
  for(i in 1:l_P)
  {
    knots_i = NULL
    for(j in 1:l_N)
    {
      n = ((P[i] * N[j] * (1 - 1 / folds)))^0.8
      if (max_knots >= n) max_knots = round(n * p_max)
      knot_p = seq(begin, end, length.out = l)
      knot = n * knot_p
      knot = sapply(knot, floor)
      avai_k = knot_p[which(knot >= min_knots & knot <= max_knots)]
      len_knot = length(knot)

      if (length(avai_k) == 0)
      {
        a = min_knots / n
        b = max_knots / n
        avai_k = seq(a, b, length.out = l)
      }else{

        if (knot[1] > min_knots)
        {
          portion_min = min_knots / n
          portion_plus = seq(portion_min, knot_p[1], length.out = num + 2)
          avai_k = c(portion_plus, avai_k[-1])
        }

        if (exe_max == TRUE)
        {
          if (knot[len_knot] < max_knots)
          {
            portion_max = max_knots / n
            portion_plus = seq(knot_p[len_knot], portion_max, length.out = num + 2)
            avai_k = c(avai_k[-len_knot], portion_plus)
          }
        }

      }

      knots_i[[j]] = avai_k

    }
    knots[[i]] = knots_i
  }
  return (knots)
}


### ESTIMATE
### Function: Semi-parametric model fitting function, nonliner in the first covariates
### Input: X: n by p covariates matrix
###        Y: n by 1 outcome
###        knots: Knots created by create_knots function
### Output: m: Design matrix
###         basis:Basis object of package fda
###         coefs: Estimated coefficients of the model
ESTIMATE = function(X, Y, knots) {
  n_row = length(X[,1])
  D_X = X
  D_Y = Y
  obj = getDesign(D_X, knots)

  D_X = cbind(obj$m, D_X[,-1])
  lmod = lm(D_Y~ D_X - 1)
  coefs = coef(lmod)
  obj$coefs = coefs
  return(obj)
}


### pred
### Function: Predict the representation of the outcome Y defined in ESTIMATE()
### Input: object: An object obtained by ESTIMATE()
###        newdata: New covariates to use
### Output: pred.pre: Predicted representation of outcome
pred = function(object, newdata) {
  x = newdata
  n = nrow(x)
  p = ncol(x)
  lowdiff = highdiff = matrix(FALSE, nrow = nrow(newdata), ncol = 1)
  x.cut  = matrix(0, nrow = n, ncol = 1)
  x.current     = x[,1]
  x.cut[,1]     = x.current
  bas = object$basis
  lower.end = bas$rangeval[1]
  upper.end = bas$rangeval[2]
  ind.lower = x.current < lower.end
  ind.upper = x.current > upper.end

  lowdiff[ind.lower,1]  = (x.current - lower.end)[ind.lower]
  highdiff[ind.upper,1] = (x.current - upper.end)[ind.upper]

  x.cut[ind.lower,1] = lower.end
  x.cut[ind.upper,1] = upper.end

  ## Get the slopes at the boundaries
  m = eval.basis(x.cut[,1], bas)
  deriv.info  = eval.basis(c(lower.end, upper.end), bas, Lfdobj = 1)
  df = ncol(m)
  lower.slopes = deriv.info[1,]
  upper.slopes = deriv.info[2,]

  beta = object$coefs

  pred.pre = cbind(m, x[, -1]) %*% beta

  ## Put the design matrix of the first derivates (lower.slopes,
  ## upper.slopes) into one long vector each (same length as index) and
  ## multiply with beta vector and take the sum. I.e. perform the matrix
  ## operation in a bit a special form.
  ## The result are the derivatives at the left- and the right-hand side
  ## boundaries of the training range (of the fitted object with the
  ## current coefficients)
  slopes.left  = rowsum(lower.slopes * beta[2:(df+1)], group = rep(1, df))
  slopes.right = rowsum(upper.slopes * beta[2:(df+1)], group = rep(1, df))

  ## Now we have to multiply the derivatives with the difference
  ## in the x-values (contained in lowdiff and highdiff)
  ## lowdiff and highdiff are matrices with dimensions n x p, i.e. the
  ## dimension of the newdata object.
  ## Each column of lowdiff and highdiff is multiplied with the slope
  ## value. The result will be what we have to add beyond the boundaries.
  ## add.left and add.right will also have dimension n x p.

  ## 'as.array' is here to force a warning message if recycling would
  ## take place (see help file of sweep)
  add.left  = sweep(lowdiff, MARGIN = 2, STATS = as.array(slopes.left), FUN = "*")
  add.right = sweep(highdiff, MARGIN = 2, STATS = as.array(slopes.right), FUN = "*")

  ## Calculate the final prediction:
  ## Take the prediction of the 'cut-down' matrix and add the linear
  ## extrapolation part (add.left + add.right). We have to take the sum
  ## in each row of the linear extrapolation part (add.left + add.right)
  pred.pre = pred.pre + rowSums(add.left + add.right)
  return(pred.pre)
}


### Cross_Validation
### Function: Use cross validation to get the best knots for the semi-parametric model
### Input: DATA: First column is the n by 1 outcome, the rest is the covariate matrix
###        folds: Number of folds for cross validation
###        knots: Knots created by create_knots function
### Output: A matrix of MSE in different folds
Cross_Validation = function(DATA, folds, knots) {
  l_knots = length(knots)
  obs = length(DATA[,1])
  outcome = rep(0, folds*l_knots)
  dim(outcome) = c(folds, l_knots)
  Resample = sample(obs)
  sub_obs = floor(obs/folds)
  for (i in 1:(folds - 1))
  {
    sub_test_index = seq((i - 1) * sub_obs + 1, i * sub_obs)
    sub_test_index = Resample[sub_test_index]

    sub_training = DATA[-sub_test_index,]
    sub_testing = DATA[sub_test_index,]
    knots_inner = knots * (length(sub_training[,1])^0.8)
    for (j in 1:l_knots)
    {
      knot =as.integer(round(knots_inner[j]))
      MODEL = ESTIMATE(sub_training[,-1], sub_training[,1], knot)
      pred.resp = pred(MODEL, sub_testing[,-1])
      outcome[i, j] = sum((pred.resp - sub_testing[,1])^2) / length(sub_testing[,1])

    }
  }
  sub_train_index = Resample[1:((folds - 1)*sub_obs)]
  sub_training = DATA[sub_train_index,]
  sub_testing = DATA[-sub_train_index,]
  knots_inner = knots * length(sub_train_index^0.8)
  for (j in 1:l_knots)
  {
    knot =as.integer(round(knots_inner[j]))
    MODEL = ESTIMATE(sub_training[,-1], sub_training[,1], knot)
    pred.resp = pred(MODEL, sub_testing[,-1])
    outcome[folds, j] = sum((pred.resp - sub_testing[,1])^2) / length(sub_testing[,1])
  }
  return(outcome)
}


