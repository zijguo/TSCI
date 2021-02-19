library(ranger)
library(caret)
library(Matrix)
# library(xgboost)


### random forest using whole data
### tuning on parameters using out of bag error
random.forest <- function(X, y, ntree=200, mtry=NULL, max.depth=0, min.node.size=5) {
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(mtry)) mtry <- round(p/3) # can be slow if p is large
  
  Data <- data.frame(cbind(y, X))
  names(Data) <- c("y", paste("X", 1:p, sep = ""))
  params.grid <- expand.grid(
    ntree = ntree,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )
  best.model <- NULL
  best.params <- NULL
  MSE <- 1e6
  for (i in 1:nrow(params.grid)) {
    out.forest <- ranger(y~., data = Data, 
                         num.trees=params.grid$ntree[i], 
                         mtry=params.grid$mtry[i],
                         max.depth = params.grid$max.depth[i],
                         min.node.size = params.grid$min.node.size[i]
    )
    if (out.forest$prediction.error <= MSE) {
      best.model <- out.forest
      best.params <- params.grid[i,]
      MSE <- out.forest$prediction.error
    }
  }
  nodes<- predict(best.model, data = Data, type = "terminalNodes")$predictions
  # this is the in-sample prediction
  predicted.values <- predict(best.model, data = Data, type = "response")$predictions
  
  returnList <- list(ranger.obj = best.model,
                     best.params = best.params,
                     nodes = nodes,
                     predicted.values = predicted.values)
  returnList
}


### data split version of random forest
### split the data into two parts a, b and cross fit
### tuning on parameters using out of bag error
random.forest.split <- function(X, y, ntree=200, mtry=NULL, max.depth=0, min.node.size=5) {
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(mtry)) mtry <- round(p/3)
  
  Data <- data.frame(cbind(y, X))
  names(Data) <- c("y", paste("X", 1:p, sep = ""))
  params.grid <- expand.grid(
    ntree = ntree,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )
  
  # split the data into two parts a and b
  all.ind <- cut(1:n, breaks = 2, labels = FALSE)
  a.ind <- which(all.ind==1, arr.ind = TRUE)
  b.ind <- setdiff(1:n, a.ind)
  Data.a <- Data[a.ind, ]
  Data.b <- Data[b.ind, ]
  
  best.model.a <- NULL; best.model.a <- NULL
  MSE.a <- 1e6; MSE.b <- 1e6
  best.params.a <- NULL; best.params.b <- NULL
  for (i in nrow(params.grid)) {
    forest.a <- ranger(y~., data = Data.a, 
                       num.trees=params.grid$ntree[i], 
                       mtry=params.grid$mtry[i],
                       max.depth = params.grid$max.depth[i],
                       min.node.size = params.grid$min.node.size[i]
    )
    forest.b <- ranger(y~., data = Data.b, 
                       num.trees=params.grid$ntree[i], 
                       mtry=params.grid$mtry[i],
                       max.depth = params.grid$max.depth[i],
                       min.node.size = params.grid$min.node.size[i]
    )
    if (forest.a$prediction.error <= MSE.a) {
      best.model.a <- forest.a
      best.params.a <- params.grid[i,]
      MSE.a <- forest.a$prediction.error
    }
    if (forest.b$prediction.error <= MSE.b) {
      best.model.b <- forest.b
      best.params.b <- params.grid[i,]
      MSE.b <- forest.b$prediction.error
    }
  }

  pred.a <- predict(best.model.b, data = Data.a)$predictions
  pred.b <- predict(best.model.a, data = Data.b)$predictions
  predicted.values <- c(pred.a, pred.b)
  mse.crossfit <- mean((predicted.values-y)^2)
  
  ### Note: nodes.aa, nodes.bb are the nodes information of data a and b
  ### on the their own fitted models respectively
  ### nodes.ab is the nodes information of data a on model b, and nodes.ba vice versa
  nodes.aa <- predict(best.model.a, data = Data.a, type = "terminalNodes")$predictions
  nodes.bb <- predict(best.model.b, data = Data.b, type = "terminalNodes")$predictions
  nodes.ab <- predict(best.model.b, data = Data.a, type = "terminalNodes")$predictions
  nodes.ba <- predict(best.model.a, data = Data.b, type = "terminalNodes")$predictions
  
  returnList <- list(forest.a = best.model.a,
                     best.params.a = best.params.a,
                     forest.b = best.model.b,
                     best.params.b = best.params.b,
                     a.ind = a.ind,
                     b.ind = b.ind,
                     predicted.values = predicted.values,
                     mse.crossfit = mse.crossfit,
                     nodes.aa = nodes.aa,
                     nodes.bb = nodes.bb,
                     nodes.ab = nodes.ab,
                     nodes.ba = nodes.ba)
  returnList
}


### multiple split verison of random forest
### the larger the folds is, the slower the algorithm
random.forest.ksplit <- function(X,y,kfolds=5,ntree=200,mtry=NULL,max.depth=0,min.node.size=5,forest.save=F) {
  n <- nrow(X);p<-ncol(X)
  ### check mtry, it's better provide mtry
  if (is.null(mtry)) mtry <- round(p/3)
  Data <- data.frame(cbind(y, X))
  names(Data) <- c("y", paste("X", 1:p, sep = ""))
  params.grid <- expand.grid(
    ntree = ntree,
    mtry = mtry,
    max.depth = max.depth,
    min.node.size = min.node.size
  )
  
  all.ind <- cut(1:n, breaks = kfolds, labels = FALSE)
  k.ind <- rep(list(NA), kfolds)
  for (k in 1:kfolds) {
    k.ind[[k]] <- which(all.ind==k)
  }
  
  forest <- best.params <- rep(list(NA), kfolds)
  MSE <- rep(1e6, kfolds)
  for (k in 1:kfolds) {
    for (i in 1:nrow(params.grid)) {
      out.forest <- ranger(y~., data = Data[-k.ind[[k]],], 
                           num.trees=params.grid$ntree[i], 
                           mtry=params.grid$mtry[i],
                           max.depth = params.grid$max.depth[i],
                           min.node.size = params.grid$min.node.size[i]
      )
      if (out.forest$prediction.error < MSE[k]) {
        forest[[k]] <- out.forest
        best.params[[k]] <- params.grid[i,]
        MSE[k] <- out.forest$prediction.error
      }
    }
  }
  
  predicted.values <- rep(0,n)
  for (k in 1:kfolds) {
    predicted.values[k.ind[[k]]] <- predict(forest[[k]],data=Data[k.ind[[k]],])$predictions
  }
  mse.crossfit <- mean((y-predicted.values)^2)
  
  ### compute nodes information
  nodes <- nodes.ref <- rep(list(NA), kfolds)
  for (k in 1:kfolds) {
    nodes[[k]] <- predict(forest[[k]],data=Data[k.ind[[k]],],type="terminalNodes")$predictions
    nodes.ref[[k]] <- predict(forest[[k]],data=Data[-k.ind[[k]],],type="terminalNodes")$predictions
  }
  
  returnList <- list(forest = forest,
                     best.params = best.params,
                     k.ind = k.ind,
                     predicted.values = predicted.values,
                     mse.crossfit = mse.crossfit,
                     nodes = nodes,
                     nodes.ref = nodes.ref)
  if (forest.save==F) {
    returnList <- returnList[-1]
  }

  returnList
}


### weighted average transformation, fast implementation
### by not calculating the transformation weight matrix
### vec.target is the target vector to transform
### rf.obj is the model object fitted by random.forest()
weight.trans <- function(vec.target, rf.obj) {
  nodes <- rf.obj$nodes
  stopifnot(nrow(nodes)==length(vec.target))
  n <- length(vec.target)
  B <- ncol(nodes)
  
  trans.mat <- matrix(0, n, B)
  for (j in 1:B) {
    for (i in 1:n) {
      ind <- nodes[,j]==nodes[i,j]
      trans.mat[i, j] <- mean(vec.target[ind])
    }
  }
  
  trans.vec <- apply(trans.mat, 1, mean) # aggregate over all trees
  trans.vec
}


### helper function to compute weights in data split setting
### nodes.ref is the nodes information of fitted data
### nodes.target is the nodes information of target data we need to transform
### target.vec is the vector we want to transform
cross.trans <- function(nodes.ref, nodes.target, vec.ref) {
  B <- ncol(nodes.ref)
  out.n <- nrow(nodes.target)
  trans.mat <- matrix(0, out.n, B)
  for (i in 1:out.n) {
    for (j in 1:B) {
      ind <- nodes.ref[,j]==nodes.target[i, j]
      trans.mat[i, j] <- mean(vec.ref[ind])
    }
  }
  trans.vec <- apply(trans.mat, 1, mean)
  trans.vec
}


### data split version of weight.trans()
### vec.target is the target vector to transform
### rf.obj is the model object fitted by random.forest.split()
weight.trans.split <- function(vec.target, rf.obj) {
  vec.ref.a <- vec.target[rf.obj$a.ind]
  vec.ref.b <- vec.target[rf.obj$b.ind]
  
  trans.a <- cross.trans(rf.obj$nodes.bb, rf.obj$nodes.ab, vec.ref.b)
  trans.b <- cross.trans(rf.obj$nodes.aa, rf.obj$nodes.ba, vec.ref.a)
  trans <- c(trans.a, trans.b)
  trans
}


### k split version of weight.trans()
weight.trans.ksplit <- function(vec.target, rf.obj) {
  k.ind <- rf.obj$k.ind
  nodes <- rf.obj$nodes
  nodes.ref <- rf.obj$nodes.ref
  kfolds <- length(nodes)
  
  kfolds <- length(nodes)
  trans.lst <- rep(list(NA),kfolds)
  for (k in 1:kfolds) {
    trans.lst[[k]] <- cross.trans(nodes.ref[[k]],nodes[[k]],vec.target[-k.ind[[k]]])
  }
  out.trans <- Reduce(c,trans.lst)
  out.trans
}


### calculate the weight matrix without data split
weight.mat <- function(rf.obj) {
  nodes <- rf.obj$nodes
  n <- nrow(nodes)
  lst.w <- rep(list(NA), ncol(nodes))
  for (j in 1:ncol(nodes)) {
    w <- matrix(0, n, n)
    for (i in 1:n) {
      ind <- nodes[,j]==nodes[i, j]
      wt <- 1/sum(ind)
      ind[ind] <- wt
      w[i, ] <- ind/ncol(nodes)
    }
    lst.w[[j]] <- Matrix(w, sparse = TRUE)
  }
  out.w <- Reduce("+", lst.w)
  out.w
}


### the data split version of weight.mat()
weight.mat.split <- function(rf.obj) {
  nodes.aa <- rf.obj$nodes.aa
  nodes.bb <- rf.obj$nodes.bb
  nodes.ab <- rf.obj$nodes.ab
  nodes.ba <- rf.obj$nodes.ba
  n_a <- nrow(nodes.aa); n_b <- nrow(nodes.bb)
  ntree.a <- ncol(nodes.aa); ntree.b <- ncol(nodes.bb)
  
  lst.ab <- rep(list(NA),ntree.b);lst.ba <- rep(list(NA),ntree.a)
  for (j in 1:ntree.b) {
    w <- matrix(0,n_a,n_b)
    for (i in 1:n_a) {
      ind <- nodes.bb[,j] == nodes.ab[i, j]
      wt <- 1/sum(ind)
      ind[ind] <- wt
      w[i,] <- ind/ntree.b
    }
    lst.ab[[j]] <- Matrix(w, sparse=T)
  }
  out.ab <- Reduce("+", lst.ab)
  
  for (j in 1:ntree.a) {
    w <- matrix(0,n_b,n_a)
    for (i in 1:n_b) {
      ind <- nodes.aa[,j] == nodes.ba[i, j]
      wt <- 1/sum(ind)
      ind[ind] <- wt
      w[i,] <- ind/ntree.a
    }
    lst.ba[[j]] <- Matrix(w, sparse=T)
  }
  out.ba <- Reduce("+", lst.ba)
  
  mat.up <- cbind(matrix(0,n_a,n_b), as.matrix(out.ab))
  mat.lo <- cbind(as.matrix(out.ba), matrix(0,n_b,n_a))
  out.w <- Matrix(rbind(mat.up, mat.lo), sparse=T)
  out.w
}


### k split version of weight.mat()
weight.mat.ksplit <- function(rf.obj) {
  n <- length(rf.obj$predicted.values)
  k.ind <- rf.obj$k.ind
  nodes <- rf.obj$nodes
  nodes.ref <- rf.obj$nodes.ref
  kfolds <- length(nodes)
  out.lst <- rep(list(NA),kfolds)
  
  for (k in 1:kfolds) {
    n.k <- length(k.ind[[k]])
    ntree.k <- ncol(nodes[[k]])
    nodes.k <- nodes[[k]]
    nodes.ref.k <- nodes.ref[[k]]
    lst.k <- rep(list(NA),ntree.k)
    for (j in 1:ntree.k) {
      w <- matrix(0,n.k,n-n.k)
      for (i in 1:n.k) {
        ind <- nodes.ref.k[,j] == nodes.k[i, j]
        wt <- 1/sum(ind)
        ind[ind] <- wt
        w[i,] <- ind/ntree.k
      }
      lst.k[[j]] <- Matrix(w, sparse = T)
    }
    
    if (k==1) {
      out.lst[[k]] <- cbind(matrix(rep(0,n.k*n.k),n.k,n.k), Reduce("+", lst.k))
    } else if (k==kfolds) {
      out.lst[[k]] <- cbind(Reduce("+", lst.k), matrix(rep(0,n.k*n.k),n.k,n.k))
    } else {
      mat.temp <- Reduce("+", lst.k)
      persize <- floor(n/kfolds)
      mat.front <- mat.temp[, 1:(persize*(k-1))]
      mat.back <- mat.temp[, -(1:(persize*(k-1)))]
      out.lst[[k]] <- cbind(mat.front, matrix(rep(0,n.k*n.k),n.k,n.k), mat.back)
    }
    
  }
  do.call(rbind, out.lst)
}

### calculate the sum of squares of weight matrix without data split
### fast implementation without calculating the matrix
# get.sumsq <- function(rf.obj, projection=FALSE, UHat=NULL) {
#   
#   if (projection) {
#     if (is.null(UHat)) {
#       stop("UHat should be provided if projections is TRUE")
#     }
#   }
#   
#   nodes <- rf.obj$nodes
#   n <- nrow(nodes);B <- ncol(nodes)
#   sumsq <- 0
#   
#   for (i in 1:n) {
#     w.vec <- rep(0,n)
#     for (j in 1:B) {
#       k.ind <- nodes[,j]==nodes[i,j]
#       wt <- 1/sum(k.ind)
#       k.ind[k.ind] <- wt
#       w.vec <- w.vec + k.ind/B
#     }
#     sumsq <- sumsq + sum(w.vec^2)
#   }
#   sumsq
# }
# 
# ### data split version of get.sumsq()
# get.sumsq.split <- function(rf.obj, projection=FALSE, UHat=NULL) {
#   
#   if (projection) {
#     if (is.null(UHat)) {
#       stop("UHat should be provided if projections is TRUE")
#     }
#   }
#   
#   nodes.aa <- rf.obj$nodes.aa
#   nodes.bb <- rf.obj$nodes.bb
#   nodes.ab <- rf.obj$nodes.ab
#   nodes.ba <- rf.obj$nodes.ba
#   n_a <- nrow(nodes.aa); n_b <- nrow(nodes.bb)
#   ntree.a <- ncol(nodes.aa); ntree.b <- ncol(nodes.bb)
#   sumsq.a <- sumsq.b <- 0
#   
#   for (i in 1:n_a) {
#     w.vec <- rep(0, n_a)
#     for (j in 1:ntree.b) {
#       k.ind <- nodes.bb[,j] == nodes.ab[i, j]
#       wt <- 1/sum(k.ind)
#       k.ind[k.ind] <- wt
#       w.vec <- w.vec + k.ind/ntree.b
#     }
#     sumsq.a <- sumsq.a + sum(w.vec^2)
#   }
#   
#   for (i in 1:n_b) {
#     w.vec <- rep(0, n_b)
#     for (j in 1:ntree.a) {
#       k.ind <- nodes.aa[,j] == nodes.ba[i, j]
#       wt <- 1/sum(k.ind)
#       k.ind[k.ind] <- wt
#       w.vec <- w.vec + k.ind/ntree.a
#     }
#     sumsq.b <- sumsq.b + sum(w.vec^2)
#   }
#   
#   sumsq <- sumsq.a + sumsq.b
#   sumsq
# }



############## functions not in use ###############
### random forest with ranger package
### cross validation for hyper-parameter tuning
### replaced by out of bag parameter selection
# cv.forest <- function(X, y, folds, ntree=100, mtry=NULL, max.depth=0, min.node.size=5) {
#   n <- nrow(X)
#   p <- ncol(X)
#   
#   Data <- data.frame(cbind(y, X))
#   names(Data) <- c("y", paste("X", 1:p, sep = ""))
#   if (is.null(mtry)) mtry <- round(2*p/3) # use rule of thumb mtry
#   
#   params.grid <- expand.grid(
#     ntree = ntree,
#     mtry = mtry,
#     max.depth = max.depth,
#     min.node.size = min.node.size
#   )
#   MSE <- rep(0, nrow(params.grid))
#   all_folds <- cut(seq(1,n),breaks=folds,labels=FALSE)
#   for (i in 1:nrow(params.grid)) {
#     for (fold in 1:folds) {
#       testInd <- which(all_folds == fold, arr.ind = TRUE)
#       # create train and test data
#       test_data <- Data[testInd, ]
#       train_data <- Data[-testInd, ]
#       
#       out.forest <- ranger(y~., data = train_data, 
#                            num.trees=params.grid$ntree[i], 
#                            mtry=params.grid$mtry[i],
#                            max.depth = params.grid$max.depth[i],
#                            min.node.size = params.grid$min.node.size[i]
#       )
#       pred <- predict(out.forest, data = test_data[,-1])$predictions
#       MSE[i] <- MSE[i] + mean((pred-test_data$y)^2)/folds
#     }
#   }
#   best.ind <- which(MSE==min(MSE), arr.ind = T)
#   best.forest <-ranger(y~., data = Data, 
#                        num.trees=params.grid$ntree[best.ind], 
#                        mtry=params.grid$mtry[best.ind], 
#                        max.depth = params.grid$max.depth[best.ind],
#                        min.node.size = params.grid$min.node.size[best.ind]
#   )
#   
#   nodes.info <- predict(best.forest, data = Data, type = "terminalNodes")$predictions
#   
#   returnList <- list(ranger.obj = best.forest,
#                      best.params = params.grid[best.ind,],
#                      cv.mse = min(MSE),
#                      nodes.info = nodes.info)
#   returnList
# }



## boosting using xgboost, benchmark model
## weights cannot be calculated
# cv.xgb <- function(X, y, folds, params.grid=NULL) {
#   n <- nrow(X)
#   p <- ncol(X)
# 
#   Data <- data.frame(cbind(y, X))
#   names(Data) <- c("y", paste("X", 1:p, sep = ""))
# 
#   # these parameter are by rule of thumb after tuning
#   if (is.null(params.grid)) {
#     params.grid = expand.grid(
#       nrounds = 1000,
#       eta = 0.01,
#       max_depth = 2,
#       gamma = 0,
#       colsample_bytree = 1,
#       min_child_weight = 3,
#       subsample = 0.5
#     )
#   }
# 
#   xgb.trcontrol = trainControl(
#     method = "cv",
#     number = 5,
#     verboseIter = FALSE,
#     allowParallel = TRUE
#   )
# 
#   xgb.train = train(
#     y~.,
#     data = Data,
#     trControl = xgb.trcontrol,
#     tuneGrid = params.grid,
#     method = "xgbTree",
#     verbose = FALSE,
#     objective = "reg:squarederror"
#   )
# 
#   nodes.info <- predict(xgb.train$finalModel, newdata = X, predleaf = TRUE)
# 
#   returnList <- list(xgb.train = xgb.train,
#                      cv.mse = min(xgb.train$results$RMSE)^2,
#                      nodes.info = nodes.info)
#   returnList
# }

### k split version of boosting
# xgb.ksplit <- function(X,y,kfolds=5,params=NULL) {
#   n <- nrow(X); p <- ncol(X)
#   
#   Data <- cbind(y, X)
#   
#   if (is.null(params)) {
#     params = list(
#       eta = 0.01,
#       max_depth = 2,
#       gamma = 0,
#       colsample_bytree = 1,
#       min_child_weight = 3,
#       subsample = 0.5
#     )
#   }
#   nrounds = 1000
#   
#   all.ind <- cut(1:n, breaks = kfolds, labels = FALSE)
#   k.ind <- rep(list(NA), kfolds)
#   for (k in 1:kfolds) {
#     k.ind[[k]] <- which(all.ind==k)
#   }
#   
#   pred <- rep(0,n)
#   for (k in 1:kfolds) {
#     train.k <- Data[-k.ind[[k]],]
#     test.k <- Data[k.ind[[k]],]
#     dtrain.k <- xgb.DMatrix(train.k[,-1],label=train.k[,1])
#     dtest.k <- xgb.DMatrix(test.k[,-1],label=test.k[,1])
#     xgb.model <- xgb.train(params,dtrain.k,nrounds = nrounds)
#     pred[k.ind[[k]]] <- predict(xgb.model,dtest.k)
#   }
#   returnList <- list(params=params,
#                      predictions=pred)
# }

