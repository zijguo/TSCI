### TSCI: Random Forest with Binary IV and Continuous Treatment
### Setting: with interaction, compare random forest and regression model


library(MASS)
source('~/Dropbox/Projects-Supervision/Wei-Yuan/TSCI/source-forest-boosting.R', encoding = 'UTF-8')

A1gen<-function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    } 
  }
  A1
}

###### dimension change the dimension 1,2,3,4,5,8,10,15,17,20
p=20
####please change this n= 1000, 2000
n=2000
##### change the interaction 1,2
inter.val<-2
##### set up the number of folds
folds<-5
#### the number of simulation numbers
nsim<-2
##################################################
#### binary iv, TRUE(1) or FALSE(0)
binary.iv<-FALSE
### Set different ways of generating E(D|Z,X) function
if(binary.iv){
  #### change f.index across 1,2,3,4,5,6
  f.index<-4
}else{
  #### change f.index across 1,2,3,4,5,6
  f.index<-5
}
#### a denotes the IV strength, set as 1
a<-1
##### Set different violation functions
if(binary.iv){
  ##### vio.index across 0,1
  vio.index<-1
}else{
  ##### vio.index across 0,1,2,3
  vio.index<-1
}
##### tau denotes the violation strength 
##### set tau as 1 
tau<-1
##################################################
f_1 <- function(x){x+a*(x^2+0.5*x^4) -25/12}
f_2 <- function(x){exp(2*a+x+0.5*x^3) - 2/5 * sinh(5/2)}
f_3 <- function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
rho1=0.5
Cov<-(A1gen(rho1,p+1))
mu<-rep(0,p+1)
beta=1
alpha<-as.matrix(rep(-0.3,p))
gamma<-as.matrix(rep(0.2,p))
### assign different structure of interaction as p changes
if (p>5) {
  inter<-as.matrix(c(rep(inter.val,5),rep(0,p-5)))
} else {
  inter<-as.matrix(rep(inter.val,p))
}
###### number of tested basis functions, continuous iv only
Q<-5
if (binary.iv) {
  estimator.names <- c("ssreg.reg","ssreg.ls","rf.ksplit1","rf.ksplit2","ssreg.rf.ksplit","rf.full2")
  Coef.matrix.inter<-matrix(NA,nrow=nsim,ncol=6)
} else {
  estimator.names <- c("ksplit.q1","ksplit.q2","ksplit.q3","ksplit.q4","ksplit.q5",
                       "ssreg.rf.q1","ssreg.rf.q2","ssreg.rf.q3","ssreg.rf.q4","ssreg.rf.q5",
                       "full.q1","full.q2","full.q3","full.q4","full.q5")
  Coef.matrix.inter<-matrix(NA,nrow=nsim,ncol=3*Q)
}
colnames(Coef.matrix.inter) <- estimator.names

### The WD(D.rep) and W\epsilon(error.rep) for nsim simulations, W means weight matrix
### only need the average of them?
D.rep.ksplit <- D.rep.full <- matrix(NA,n,nsim)
error.rep.ksplit <- error.rep.full <- matrix(NA,n,nsim)

### correlation between D.rep and error.rep
### column means different estimators: ksplit, full
corr <- matrix(NA, nsim, 2)

### compute var(Error[,2]-error.rep) for each simulation
### for ksplit, full
error.var <- matrix(NA, nsim, 2)

### the spectral norm, LHS of equation (27)
### this is the case kfolds for data split is 5
### only ksplit is considered here
spectral.norms <- matrix(NA,nsim,5)
### the trace of RHS of equation (28), ignoring SigmaSquare here
ksplit.trace <- matrix(NA,nsim,5)

for(i in 1:nsim){
  print(i)
  #### generate the data
  mu.error<-rep(0,2)
  Cov.error<-matrix(c(1,0.5,0.5,1),2,2)
  Error<-mvrnorm(n, mu.error, Cov.error)
  W.original<-mvrnorm(n, mu, Cov)
  W<-pnorm(W.original)
  Z<-W[,1]
  if(binary.iv){
    Z<-(Z>0.6)
  }
  X<-W[,-1]
  ###### generate the data for the treatment variable D
  if(f.index==1){
    D=f_1(Z)+X%*%alpha+Error[,1]
  }
  if(f.index==2){
    D=f_2(Z)+X%*%alpha+Error[,1]
  }
  if(f.index==3){
    D=f_3(Z)+X%*%alpha+Error[,1]
  }
  if(f.index==4){
    D=f_1(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
  }
  if(f.index==5){
    D=f_2(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
  }
  if(f.index==6){
    D=f_3(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
  }
  ####### generate the outcome variable
  if(vio.index==0){
    Y=D*beta+ X%*% gamma+Error[,2]
  }
  if(vio.index==1){
    Y=D*beta+ tau*Z+ X%*% gamma+Error[,2]
  }
  if(vio.index==2){
    Y=D*beta+ tau*(Z^2-1)+ X%*% gamma+Error[,2]
  }
  if(vio.index==3){
    Y=D*beta+ tau*Z^3+ X%*% gamma+Error[,2]
  }
  if(vio.index==4){
    Y=D*beta+ tau*Z+ Z^5/40+ X%*% gamma+Error[,2]
  }
  #if(vio.index==5){
  #  Y=D*beta+ tau*(Z^2-1)+ 1/20*sin(Z/2)+ X%*% gamma+Error[,2]
  #}
  Data<-cbind(D,W)
  ###### try the regression model with interaction
  if (binary.iv) {
    D.pred<-predict(lm(D~Z+W[,-1]+Z*W[,-1]))
    ssreg.reg<-lm(Y~D.pred+W[,-1])
    ssreg.ls<-lm(Y~D.pred+Z+W[,-1])
    Coef.matrix.inter[i,1]<-coef(ssreg.reg)[2]
    Coef.matrix.inter[i,2]<-coef(ssreg.ls)[2]
  }
  ###### random forest
  forest.cov<-cbind(Z,W[,-1])
  
  ### use k split here
  ### deciding grids of mtry by p
  if (p+1<=11) {
    mtry = 1:(p+1)
  } else {
    mtry = seq(round((p+1)/3), round(2*(p+1)/3), by=1)
  }
  forest.semi <- random.forest.ksplit(forest.cov, D, mtry = mtry, max.depth = 0:5, min.node.size = c(1,5,10))
  D.rep<-forest.semi$predicted.values
  D.rep.ksplit[,i] <- D.rep
  weight.semi <- weight.mat.ksplit(forest.semi)
  error.rep.ksplit[,i] <- as.matrix(weight.semi)%*%Error[,2]
  error.var[i,1] <- var(Error[,2]-error.rep.ksplit[,i])
  corr[i,1] <- cor(D.rep,error.rep.ksplit[,i])
  if (binary.iv) {
    W.rep <- as.matrix(weight.semi)%*%forest.cov
  } else {
    Cov.aug<-cbind(Z^4,Z^3,Z^2,Z,W[,-1])
    W.rep <- as.matrix(weight.semi)%*%Cov.aug
  }
  Y.rep <- as.matrix(weight.semi)%*%Y
  if (binary.iv) {
    reg.rf1<-lm(Y.rep~D.rep+W.rep[,-1]) 
    reg.rf2<-lm(Y.rep~D.rep+W.rep)
    #summary(reg.rf1)
    #summary(reg.rf2)
    Coef.matrix.inter[i,3]<-coef(reg.rf1)[2]
    Coef.matrix.inter[i,4]<-coef(reg.rf2)[2]
  } else {
    for(q in 0:(Q-1)){
      if(q==Q-1){
        reg.rf<-lm(Y.rep~D.rep+W.rep)    
        #print(remove.set)
      }else{
        reg.rf<-lm(Y.rep~D.rep+W.rep[,-(1:(Q-1-q))])    
        #print(remove.set)
      }
      Coef.matrix.inter[i,q+1]<-coef(reg.rf)[2]
    }
  }
  ### calculate spectral norms and trace
  # k.ind <- forest.semi$k.ind
  # for (j in 1:length(k.ind)) {
  #   n.k <- length(k.ind[[j]])
  #   U.k <- W.rep[k.ind[[j]],]
  #   ### norms
  #   P.k <- (diag(1,n.k,n.k)-U.k %*% solve(t(U.k)%*%U.k) %*% t(U.k))
  #   spectral.norms[i,j] <- norm(P.k%*%D.rep[k.ind[[j]]],type = "2")
  #   ### trace
  #   Omega.k <- weight.semi[k.ind[[j]],-k.ind[[j]]]
  #   ksplit.trace[i,j] <- sum((P.k%*%Omega.k)^2)
  # }
  ###### try the random forest with D.rep defined as the predicted value
  ###### after applying k-split random forest to D
  if (binary.iv) {
    ssreg.rf<-lm(Y~D.rep+Z+W[,-1])
    Coef.matrix.inter[i,5]<-coef(ssreg.rf)[2]
  } else {
    for(q in 0:(Q-1)){
      if(q==Q-1){
        ssreg.rf<-lm(Y~D.rep+Cov.aug)    
        #print(remove.set)
      }else{
        ssreg.rf<-lm(Y~D.rep+Cov.aug[,-(1:(Q-1-q))])    
        #print(remove.set)
      }
      Coef.matrix.inter[i,q+6]<-coef(ssreg.rf)[2]
    }
  }
  ### use random forest without data split, full data
  ### use 2 to denote this
  forest.semi2 <- random.forest(forest.cov, D, mtry, max.depth = 0:5, min.node.size = c(1,5,10))
  # in-sample predictions
  D.rep2 <- forest.semi2$predicted.values
  D.rep.full[,i] <- D.rep2
  weight.semi2 <- weight.mat(forest.semi2)
  error.rep.full[,i] <- as.matrix(weight.semi2)%*%Error[,2]
  error.var[i,2] <- var(Error[,2]-error.rep.full[,i])
  corr[i,2] <- cor(D.rep2,error.rep.full[,i])
  if (binary.iv) {
    W.rep2 <- as.matrix(weight.semi2)%*%forest.cov
  } else {
    Cov.aug<-cbind(Z^4,Z^3,Z^2,Z,W[,-1])
    W.rep2 <- as.matrix(weight.semi2)%*%Cov.aug
  }
  Y.rep2 <- as.matrix(weight.semi2)%*%Y
  if (binary.iv) {
    reg2.rf2<-lm(Y.rep2~D.rep2+W.rep2)
    Coef.matrix.inter[i,6]<-coef(reg2.rf2)[2]
  } else {
    for(q in 0:(Q-1)){
      if(q==Q-1){
        reg.rf.full<-lm(Y.rep2~D.rep2+W.rep2)    
        #print(remove.set)
      }else{
        reg.rf.full<-lm(Y.rep2~D.rep2+W.rep2[,-(1:(Q-1-q))])    
        #print(remove.set)
      }
      Coef.matrix.inter[i,q+11]<-coef(reg.rf.full)[2]
    }
  }
}    
p
n
inter.val
vio.index
f.index

# lm(D.rep~Z+W[,-1])
# lm(D.pred~Z+W[,-1])
# D.resid.rf<-resid(lm(D.rep~Z+W[,-1]))
# sum(D.resid.rf*Y)/sum(D.resid.rf^2)
# D.resid<-resid(lm(D.pred~Z+W[,-1]))
# sum(D.resid*Y)/sum(D.resid^2)
# mean((D-D.rep)^2)
# mean((D-D.pred)^2)
# mean(D.rep*Y)
# mean(D.pred*Y)
# summary(lm(D.resid.rf~D.resid))

apply(Coef.matrix.inter,2,mean)
apply(Coef.matrix.inter,2,sd)

### to get the average of WD and W/epsilon
# apply(D.rep.ksplit,1,mean)
# apply(error.rep.ksplit,1,mean)
# apply(D.rep.full,1,mean)
# apply(error.rep.full,1,mean)


### compute coverage
sd <- apply(Coef.matrix.inter,2,sd)
Cov.mat <- matrix(NA,nrow=nsim,ncol=ncol(Coef.matrix.inter))
colnames(Cov.mat) <- estimator.names
for (j in 1:ncol(Coef.matrix.inter)) {
  Cov.mat[,j] <- ifelse(Coef.matrix.inter[,j]-1.96*sd[j]<=beta & beta<=Coef.matrix.inter[,j]+1.96*sd[j],yes=1,no=0)
}
apply(Cov.mat,2,mean)
