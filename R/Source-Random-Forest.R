#' @title Two Stage Curvature Identification with Random Forests
#' @description This function implements Two Stage Curvature Identification with the Random Forests. It tests the IV strength and chooses the best violation form, and also constructs the confidence interval for the treatment effect with the selected violation form.
#'
#' @param Y outcome with dimension n by 1
#' @param D treatment with dimension n by 1
#' @param Z instrument variable with dimension n by 1
#' @param X baseline covariates with dimension n by p
#' @param vio.space a matrix or a list. If a matrix, then each column corresponds to a violation form of Z; If a list, then each element corresponds to a violation form of Z and must be a matrix of n rows, e.g. (Z^3,Z^2); If NULL, then default by the n by 3 matrix (Z^3, Z^2, Z). Violation form selection will be performed according to provided violation forms, for example, null violation space vs Z vs (Z^2, Z) vs (Z^3, Z^2, Z) in the default case
#' @param intercept logic, including the intercept or not in the outcome model, default by TRUE
#' @param A1.ind the indices of samples in A1, used for constructing the point estimator and the confidence interval, default by randomly selected round(2/3*n) samples from 1 to n
#' @param num.trees number of trees in Random Forests, default by 200
#' @param mtry number of covariates to possibly split at in each node of the tree in Random Forests, default by a sequence from round((p+1)/3) to round(2(p+1)/3)
#' @param max.depth maximal tree depth in Random Forests, default by 0, which refers to unlimited depth
#' @param min.node.size minimal size of each leaf node in Random Forests, default by the set {5, 10, 15}
#' @param str.thol minimal value of the threshold of IV strength test, default by 10
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
#'     \item{\code{Qmax}}{the index of largest violation space selected by IV strength test. If -1, the IV strength test fails for null violation space and run OLS. If 0, the IV Strength test fails for the null violation space and run TSRF only for null violation space. In other cases, violation space selection is performed}
#'     \item{\code{q.hat}}{the index of estimated violation space corresponding to Qmax}
#'     \item{\code{invalidity}}{invalidity of TSLS. If TRUE, the IV is invalid; Otherwise, the IV is valid}
#' @export
#'
#' @examples
#' \dontrun{
#' # dimension
#' p = 10
#' # sample size
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
#' # instrument variable
#' Z=W[,1]
#' # baseline covariates
#' X=W[,-1]
#' # generate the treatment variable D
#' D=f(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
#' # generate the outcome variable Y
#' Y=D*beta+tau*Z+X%*%gamma+Error[,2]
#'
#'
#' # Two Stage Random Forest
#' output.RF = TSRF(Y,D,Z,X)
#' # point estimates
#' output.RF$Coef.robust
#' # standard errors
#' output.RF$sd.robust
#' # confidence intervals
#' output.RF$CI.robust
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom Matrix Matrix
#' @importFrom stats coef lm poly predict qnorm quantile resid rnorm
#' @import ranger
#'
TSRF = function(Y,D,Z,X,vio.space=NULL,intercept=TRUE,A1.ind=NULL,num.trees=NULL,mtry=NULL,max.depth=NULL,min.node.size=NULL,str.thol=10,alpha=0.05) {
  if (!is.null(vio.space)) {
    if (class(vio.space)[1] != "matrix" & class(vio.space)[1] != "list") {
      stop("The violation space must be input as matrix or list")
    }
  }
  Y = as.matrix(Y); D = as.matrix(D); Z = as.matrix(Z); X = as.matrix(X);
  n = nrow(X); p = ncol(X)
  # default value for hyper-parameters tuning
  if (is.null(num.trees)) num.trees = 200
  if (is.null(mtry)) mtry = seq(round((p+1)/3),round(2*(p+1)/3),by=1)
  if (is.null(max.depth)) max.depth = 0
  if(is.null(min.node.size)) min.node.size = c(5,10,20)


  # indices of A1
  if (!is.null(A1.ind)) {
    n_A1 = length(A1.ind)
  } else {
    A1.ind = sample(1:n,round(2/3*n))
    n_A1 = length(A1.ind)
  }


  # perform Random Forest
  # Treatment model fitting
  forest = TSRF.fit(D,Z,X,num.trees=num.trees,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size,A1.ind=A1.ind)
  # Compute weight matrix
  weight = TSRF.weight(forest$nodes.A1)

  # Selection
  outputs = TSRF.Selection(Y=Y,D=D,Z=Z,X=X,vio.space=vio.space,A1.ind=A1.ind,weight=weight,intercept=intercept,str.thol=str.thol,alpha=alpha)
  return(outputs)
}


#' @title Random Forest Fitting
#' @description Implement Random forest with data splitting, with cross validation using out of bag error
#'
#' @param D treatment with dimension n by 1
#' @param Z instrument variable with dimension n by 1
#' @param X baseline covariates with dimension n by p
#' @param num.trees number of trees in Random Forest, default by 200
#' @param mtry number of covariates to possibly split at in each node of the tree in Random Forest, default by a sequence from round((p+1)/3) to round(2(p+1)/3)
#' @param max.depth maximal tree depth in Random Forest, default by 0, which refers to unlimited depth
#' @param min.node.size minimal size of each leaf node in Random Forest, default by the set {5, 10, 15}
#' @param A1.ind the indices of samples in A1, used for constructing point estimator and confidence intervals, default by randomly selected round(2/3*n) indices from 1 to n
#' @param MSE.thol a very large threshold for MSE, default by 1e12, used for the start of hyper-parameter selection
#' @param forest.save save the Random Forest output or not, default by FALSE
#'
#' @return:
#'     \item{\code{forest.A2}}{random forest built on subsample A2, available if forest.save=TRUE}
#'     \item{\code{params.A2}}{best hyper-parameters for forest.A2 selected by out-of-bag error}
#'     \item{\code{A1.ind}}{indices of subsample A1}
#'     \item{\code{nodes.A1}}{a n_A1 by num.trees matrix, rows refer to different samples, columns refer to different trees, the entrees are leaf node indices of each sample in each tree}
#'     \item{\code{MSE.oob}}{minimal out-of-bag error using the best hyper-parameters}
#' @noRd
#'
TSRF.fit = function(D,Z,X,num.trees,mtry,max.depth,min.node.size,A1.ind,MSE.thol=1e12,forest.save=FALSE) {
  W = as.matrix(cbind(Z,X)); D = as.matrix(D)
  n = nrow(W); p = ncol(X)
  Data = data.frame(cbind(D,W))
  names(Data) = c("D", paste("W", 1:(p+1), sep = ""))
  # grid search
  params.grid = expand.grid(
    num.trees = num.trees,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )
  # split the data into two parts A1 and A2
  # use A2 to build the random forest and use A1 to construct weight matrix
  n_A1 = length(A1.ind)
  Data.A1 = Data[A1.ind,]
  Data.A2 = Data[-A1.ind,]
  forest.A2 = NULL;
  MSE.oob.A2 = MSE.thol
  params.A2 = NULL
  ### use out of bag error to do hyper-parameter tuning
  for (i in 1:nrow(params.grid)) {
    temp.A2 = ranger(D~., data = Data.A2,
                      num.trees=params.grid$num.trees[i],
                      mtry=params.grid$mtry[i],
                      max.depth = params.grid$max.depth[i],
                      min.node.size = params.grid$min.node.size[i],
                      importance = "impurity"
    )
    if (temp.A2$prediction.error <= MSE.oob.A2) {
      forest.A2 = temp.A2
      params.A2 = params.grid[i,]
      MSE.oob.A2 = temp.A2$prediction.error
    }
  }
  # leaf nodes information of A1 on the Random Forest built on A2
  nodes.A1 = predict(forest.A2, data = Data.A1, type = "terminalNodes")$predictions
  returnList = list(forest.A2 = forest.A2,
                     params.A2 = params.A2,
                     A1.ind = A1.ind,
                     nodes.A1 = nodes.A1,
                     MSE.oob.A2 = MSE.oob.A2)
  if (!forest.save) returnList = returnList[-1]
  returnList
}


#' @title Weight matrix calculation for Random Forest
#' @description Calculate weight matrix for Random Forest with node information of samples in A1
#'
#' @param nodes A n_A1 by num.trees matrix, with rows referring to different samples, columns referring to different trees, the entrees referring leaf node indices of each sample in each tree.
#'
#' @return:
#'     \item{\code{out.weight}}{A n_A1 by n_A1 symmetric sparse weight matrix, with class dgCMatrix}
#' @noRd
#'
TSRF.weight = function(nodes) {
  n_A1 = nrow(nodes); num.trees = ncol(nodes)
  out.weight = matrix(0,n_A1,n_A1)
  for (j in 1:num.trees) {
    weight.mat = matrix(0,n_A1,n_A1) # weight matrix for single tree
    unique.nodes = unique(nodes[,j])
    for (i in 1:length(unique.nodes)) {
      ind = nodes[,j]==unique.nodes[i] # indices of samples in the node
      num.samples = sum(ind) # number of samples in the node
      w = 1/(num.samples-1)  # weight, to remove self-prediction
      weight.vec = ifelse(ind,yes=w,no=0)
      weight.mat[ind,] = matrix(rep(weight.vec,num.samples),num.samples,byrow=T)/num.trees
    }
    diag(weight.mat) = 0 # remove self prediction
    out.weight = out.weight + weight.mat
  }
  out.weight = Matrix(out.weight, sparse = T) # sparse matrix to save memory
  return(out.weight)
}


#' @title Statistics Calculation for violation space selection
#' @description Calculate the statistics needed for violation space selection
#'
#' @param D.rep transformed treatment of dimension n_A1 by 1 corresponding to a violation space
#' @param Cov.rep transformed augmented covariates (vio.space, X) of dimension n_A1 by (p+ncol(vio.space)) corresponding to a violation space
#' @param weight n_A1 by n_A1 weight matrix
#' @param n full sample size
#' @param eps.hat residuals in the outcome model
#' @param delta.hat residuals in Random Forest corresponding to samples in A1
#' @param str.thol the minimal value of the threshold of IV strength test
#'
#' @return:
#'     \item{\code{sd}}{estimated standard error for betaHat}
#'     \item{\code{D.resid}}{residuals for D.rep~Cov.rep}
#'     \item{\code{iv.str}}{estimated IV strength}
#'     \item{\code{iv.thol}}{threshold for IV strength test}
#'     \item{\code{diag.M}}{diagonal of matrix M_{RF}(V) corresponding to a violation space}
#' @noRd
#'
TSRF.stat = function(D.rep, Cov.rep, weight, n, eps.hat, delta.hat, str.thol) {
  n_A1 = length(D.rep)
  # compute the trace of M_{RF}(V)
  # the trace of M_{RF}(V) matrix can be computed as RSS of regressing each column of weight matrix on Cov.rep
  SigmaSqD = mean(delta.hat^2)
  diag.M = rep(NA,n_A1)
  for (j in 1:n_A1) {
    diag.M[j] = sum(resid(lm(weight[,j]~Cov.rep))^2)
  }
  trace.M = sum(diag.M)
  D.rep2 = weight%*%D.rep
  D.resid = resid(lm(D.rep~Cov.rep))
  D.RSS = sum(D.resid^2)
  iv.str = D.RSS/SigmaSqD
  sd = sqrt(sum(eps.hat^2*(weight%*%D.resid)^2))/D.RSS


  # bootstrap for the threshold of IV strength test
  boot.vec = rep(NA,300)
  delta.cent = delta.hat - mean(delta.hat)
  for (i in 1:300) {
    delta = rep(NA,n_A1)
    for (j in 1:n_A1) {
      U.j = rnorm(1)
      delta[j] = delta.cent[j]*U.j
    }

    delta.rep = weight%*%delta
    delta.resid = resid(lm(as.matrix(delta.rep)~Cov.rep))
    boot.vec[i] = sum(delta.resid^2) + 2*sum(D.rep2*delta.resid)
  }
  iv.thol = quantile(boot.vec,0.975)/SigmaSqD + max(2*trace.M, str.thol)
  returnList = list(
    sd = sd,
    D.resid = D.resid,
    iv.str = iv.str,
    iv.thol = iv.thol,
    diag.M = diag.M)
  returnList
}


#' @title Violation Space Selection of Random Forest
#' @description Select violation space for Random Forest and construct confidence intervals
#'
#' @param Y outcome with dimension n by 1
#' @param D treatment with dimension n by 1
#' @param Z instrument variable with dimension n by 1
#' @param X baseline covariates with dimension n by p
#' @param vio.space a matrix or a list. If a matrix, then each column corresponds to a violation form of Z; If a list, then each element corresponds to a violation form of Z and must be a matrix of n rows, e.g. (Z^3,Z^2); If NULL, then default by the n by 3 matrix (Z^3, Z^2, Z). Violation space selection will be performed according to provided violation space, for example, null violation space vs Z vs (Z^2, Z) vs (Z^3, Z^2, Z) in the default case
#' @param A1.ind the indices of samples in A1
#' @param weight n_A1 by n_A1 weight matrix
#' @param intercept logic, to include intercept in the outcome model or not
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
#'     \item{\code{Qmax}}{the index of largest violation space selected by IV strength test. If -1, the IV strength test fails for null violation space and run OLS. If 0, the IV Strength test fails for the first violation space and run TSRF only for null violation space. In other cases, violation space selection is performed}
#'     \item{\code{q.hat}}{the index of estimated violation space corresponding to Qmax}
#'     \item{\code{invalidity}}{invalidity of TSLS. If TRUE, the IV is invalid; Otherwise, the IV is valid}
#' @noRd
#'
TSRF.Selection = function(Y, D, Z, X, vio.space, A1.ind, weight, intercept, str.thol, alpha) {
  Y = as.matrix(Y); D = as.matrix(D); Z = as.matrix(Z); X = as.matrix(X);
  n = length(Y); n_A1 = length(A1.ind)

  if (nrow(weight) != ncol(weight)) {
    stop("Transformation matrix must be a square matrix")
  }
  if(nrow(weight) != n_A1) {
    stop("The samples to construct transformation matrix must be the same as samples in A1")
  }

  Y.A1 = Y[A1.ind]
  D.A1 = D[A1.ind]

  # define the vio.space as polynomials if not specified
  if (is.null(vio.space)) {
    Q = 4
    vio.space = matrix(NA,n,0)
    for (q in 1:(Q-1)) {
      vio.space = cbind(Z^q,vio.space)
    }
    # the indices to remove to identify violation space
    rm.ind = rep(list(NA),Q-1)
    for (i in 1:(Q-1)) {
      rm.ind[[i]] = 1:(Q-i)
    }
  }

  if (!is.null(vio.space)) {
    if (class(vio.space)[1]=="list") {
      vio.space = lapply(vio.space, as.matrix)
      Q = length(vio.space) + 1
      v.len = sapply(vio.space,dim)[2,]
      # the indices to remove to identify violation space
      rm.ind = rep(list(NA),Q-1)
      for (i in 1:(Q-1)) {
        rm.ind[[i]] = 1:sum(v.len[1:(Q-i)])
      }
      # merge the list of violation space to a matrix
      vio.space = Reduce(cbind,vio.space)
    } else if (class(vio.space)[1]=="matrix") {
      Q = ncol(vio.space) + 1
      rm.ind = rep(list(NA),Q-1)
      for (i in 1:(Q-1)) {
        rm.ind[[i]] = 1:(Q-i)
      }
    }
  }

  # define the augmentation of covariates as combination of violation space and baseline covariates
  Cov.aug = cbind(vio.space,X)
  Cov.aug.A1 = Cov.aug[A1.ind,]

  # compute the representations
  Y.rep = as.matrix(weight%*%Y.A1); D.rep = as.matrix(weight%*%D.A1)
  Cov.rep = as.matrix(weight%*%Cov.aug.A1)
  # the residuals of Random Forest
  delta.hat = D.A1 - D.rep
  SigmaSqD = mean(delta.hat^2)
  # save the estimates for selection part
  Coef.names = c(paste("TSRF-Init-q",0:(Q-1),sep=""),paste("TSRF-q",0:(Q-1),sep=""))
  Coef.all = sd.all = rep(NA,2*Q)
  names(Coef.all) = names(sd.all) = Coef.names
  # IV strength test
  iv.str = iv.thol = rep(NA,Q)
  names(iv.str) = names(iv.thol) = paste("q",0:(Q-1),sep="")
  # the residuals of outcome model
  eps.hat = rep(list(NA),Q)


  ### fixed violation space, compute necessary inputs of selection part
  D.resid = diag.M.list = rep(list(NA),Q)
  for (index in 1:Q) {
    if (index==Q) {
      if (intercept) {
        reg.rf = lm(Y.rep~D.rep+Cov.rep)
        betaHat = coef(reg.rf)[2]
      } else {
        reg.rf = lm(Y.rep~D.rep+Cov.rep-1)
        betaHat = coef(reg.rf)[1]
      }
      Coef.all[index] = betaHat
      eps.hat[[index]] = resid(lm(Y.A1-D.A1*betaHat~Cov.aug.A1))
      stat.outputs = TSRF.stat(D.rep,Cov.rep,weight,n,eps.hat[[index]],delta.hat,str.thol=str.thol)
    } else {
      if (intercept) {
        reg.rf = lm(Y.rep~D.rep+Cov.rep[,-rm.ind[[index]]])
        betaHat = coef(reg.rf)[2]
      } else {
        reg.rf = lm(Y.rep~D.rep+Cov.rep[,-rm.ind[[index]]]-1)
        betaHat = coef(reg.rf)[1]
      }
      Coef.all[index] = betaHat
      eps.hat[[index]] = resid(lm(Y.A1-D.A1*betaHat~Cov.aug.A1[,-rm.ind[[index]]]))
      stat.outputs = TSRF.stat(D.rep,Cov.rep[,-rm.ind[[index]]],weight,n,eps.hat[[index]],delta.hat,str.thol=str.thol)
    }
    # the necessary statistics
    sd.all[index] = stat.outputs$sd
    sd.all[index+Q] = stat.outputs$sd
    iv.str[index] = stat.outputs$iv.str; iv.thol[index] = stat.outputs$iv.thol;
    D.resid[[index]] = stat.outputs$D.resid
    diag.M.list[[index]] = stat.outputs$diag.M
  }
  # Residual sum of squares of D.rep~Cov.rep
  D.RSS = iv.str*SigmaSqD


  # violation space selection
  # all of the q below are from 0 to (Q-1), so use q+1 to index the columns
  # Comparison and robust estimators
  Coef.robust = sd.robust = rep(NA,2)
  names(Coef.robust) = names(sd.robust) = c("TSRF-comp","TSRF-robust")
  ivtest.vec = (iv.str>=iv.thol)
  if (sum(ivtest.vec)==0) {
    warning("Weak IV, even if the IV is assumed to be valid; run OLS")
    Qmax = -1
  } else {
    Qmax = sum(ivtest.vec)-1
    if (Qmax==0) {
      warning("Weak IV, if the IV is invalid. We still test the IV invalidity.")
    }
  }


  # Compute bias-corrected estimators
  for (i in 1:Q) {
    Coef.all[i+Q] = Coef.all[i] - sum(diag.M.list[[i]]*delta.hat*eps.hat[[i]])/D.RSS[i]
  }
  sd.all[-(1:Q)] = sd.all[1:Q]


  # If IV test fails at q0 or q1, we do not need to do selection
  if (Qmax>=1) { # selection
    eps.Qmax = eps.hat[[Qmax+1]]
    Coef.Qmax = rep(NA,Q)
    for (i in 1:Q) {
      Coef.Qmax[i] = Coef.all[i] - sum(diag.M.list[[i]]*delta.hat*eps.Qmax)/D.RSS[i]
    }

    ### Selection
    # define comparison matrix
    H = beta.diff = matrix(0,Qmax,Qmax)
    # compute H matrix
    for (q1 in 0:(Qmax-1)) {
      for (q2 in (q1+1):(Qmax)) {
        H[q1+1,q2] = as.numeric(sum((weight%*%D.resid[[q1+1]])^2*eps.Qmax^2)/(D.RSS[q1+1]^2) +
                                  sum((weight%*%D.resid[[q2+1]])^2*eps.Qmax^2)/(D.RSS[q2+1]^2) -
                                  2*sum(eps.Qmax^2*(weight%*%D.resid[[q1+1]])*(weight%*%D.resid[[q2+1]]))/(D.RSS[q1+1]*D.RSS[q2+1])
        )

      }
    }
    # compute beta difference matrix, use Qmax
    for (q in 0:(Qmax-1)) {
      beta.diff[q+1,(q+1):(Qmax)] = abs(Coef.Qmax[q+1]-Coef.Qmax[(q+2):(Qmax+1)]) # use bias-corrected estimator
    }
    # bootstrap for the quantile of the differences
    max.val = rep(NA,300)
    eps.Qmax.cent = eps.Qmax - mean(eps.Qmax)
    for (i in 1:300) {
      diff.mat = matrix(0,Qmax,Qmax)
      eps = rep(NA,n_A1)
      for (j in 1:n_A1) {
        U.j = rnorm(1)
        eps[j] = eps.Qmax.cent[j]*U.j
      }
      eps.rep = weight%*%eps
      for (q1 in 0:(Qmax-1)) {
        for (q2 in (q1+1):(Qmax)) {
          diff.mat[q1+1, q2] = sum(D.resid[[q2+1]]*eps.rep)/(D.RSS[q2+1])-sum(D.resid[[q1+1]]*eps.rep)/(D.RSS[q1+1])
        }
      }
      diff.mat = abs(diff.mat)/sqrt(H)
      max.val[i] = max(diff.mat,na.rm = TRUE)
    }
    z.alpha = 1.01*quantile(max.val,0.975)
    diff.thol = z.alpha*sqrt(H)
    # comparison matrix
    C.alpha = ifelse(beta.diff<=diff.thol,0,1)


    # vector indicating the selection of each layer
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
  Coef.all = c(Coef.all[(Q+1):(2*Q)],Coef.OLS)
  names(Coef.all) = c(Coef.names[(Q+1):(2*Q)],"OLS")
  sd.all = c(sd.all[(Q+1):(2*Q)],sd.OLS)
  names(sd.all) = c(Coef.names[(Q+1):(2*Q)],"OLS")

  # confidence intervals for all violation spaces
  CI.all = rbind(Coef.all + qnorm(alpha/2)*sd.all,Coef.all + qnorm(1-alpha/2)*sd.all)
  rownames(CI.all) = c("lower","upper")


  ### estimated violation space and corresponding estimator
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


  # confidence intervals for selected violation space
  CI.robust = rbind(Coef.robust + qnorm(alpha/2)*sd.robust,Coef.robust + qnorm(1-alpha/2)*sd.robust)
  rownames(CI.robust) = c("lower","upper")



  returnList = list(Coef.all = Coef.all,
                    sd.all = sd.all,
                    CI.all = CI.all,
                    Coef.robust = Coef.robust,
                    sd.robust = sd.robust,
                    CI.robust = CI.robust,
                    iv.str = iv.str, iv.thol = iv.thol,
                    Qmax = Qmax,
                    q.hat = q.hat,
                    invalidity = invalidity)
  returnList
}

