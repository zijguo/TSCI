library(MASS)
library(fda)
library(AER)

### This is the function to generate covariace matrix

A1gen<-function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    } 
  }
  A1
}

getDesign <- function(x, knots)
{
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


create_knots <- function(P, N, folds, begin, end, l, min_knots = 2, max_knots = 50, num = 10, exe_max = FALSE, p_max = 0.5)
{
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




ESTIMATE <- function(X, Y, knots)
{
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

pred <- function(object, newdata)
{
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




Cross_Validation <- function(DATA, folds, knots)
{
  
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








