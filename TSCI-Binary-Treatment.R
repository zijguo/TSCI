### TSCI: Random Forest with Binary IV and Binary Treatment
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
p=1
####please change this n= 1000, 2000
n=1000
##### change the interaction 1,2
inter.val<-2
##### change f.index across 4,5,6
f.index<-4
##### set up the number of folds
folds<-5
#### the number of simulation numbers
nsim<-10
##################################################
#### binary iv, TRUE in this setting
binary.iv<-TRUE
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
estimator.names <- c("ssreg.reg","ssreg.ls","rf.ksplit1","rf.ksplit2","ssreg.rf.ksplit","rf.full2")
Coef.matrix.inter<-matrix(NA,nrow=nsim,ncol=6)
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
    exp.val=a*Z+X%*%alpha+Error[,1]
    prob = exp(exp.val)/(1+exp(exp.val))
    D = rbinom(n,1,prob)
  }
  if(f.index==2){
    exp.val=a*Z+X%*%alpha+Error[,1]
    prob = exp(exp.val)/(1+exp(exp.val))
    D = rbinom(n,1,prob)
  }
  if(f.index==3){
    exp.val=a*Z+X%*%alpha+Error[,1]
    prob = exp(exp.val)/(1+exp(exp.val))
    D = rbinom(n,1,prob)
  }
  if(f.index==4){
    exp.val=a*Z+X%*%alpha+Z*X%*%inter+Error[,1]
    prob = exp(exp.val)/(1+exp(exp.val))
    D = rbinom(n,1,prob)
  }
  if(f.index==5){
    exp.val=a*Z+X%*%alpha+Z*X%*%inter+Error[,1]
    prob = exp(exp.val)/(1+exp(exp.val))
    D = rbinom(n,1,prob)
  }
  if(f.index==6){
    exp.val=a*Z+X%*%alpha+Z*X%*%inter+Error[,1]
    prob = exp(exp.val)/(1+exp(exp.val))
    D = rbinom(n,1,prob)
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
  ###### try the logistic regression model with interaction
  D.pred<-predict(glm(factor(D)~Z+W[,-1]+Z*W[,-1],family = binomial),type = "response")
  ssreg.reg<-lm(Y~D.pred+W[,-1])
  ssreg.ls<-lm(Y~D.pred+Z+W[,-1])
  Coef.matrix.inter[i,1]<-coef(ssreg.reg)[2]
  Coef.matrix.inter[i,2]<-coef(ssreg.ls)[2]
  ###### random forest
  forest.cov<-cbind(Z,W[,-1])
  ### use k split here
  ### deciding grids of mtry by p
  if (p+1<=11) {
    mtry = 1:(p+1)
  } else {
    mtry = seq(round((p+1)/3), round(2*(p+1)/3), by=1)
  }
  forest.semi <- random.forest.ksplit(forest.cov, D, mtry = mtry, max.depth = 0:5, min.node.size = c(5,10))
  D.rep<-forest.semi$predicted.values
  D.rep.ksplit[,i] <- D.rep
  weight.semi <- weight.mat.ksplit(forest.semi)
  error.rep.ksplit[,i] <- as.matrix(weight.semi)%*%Error[,2]
  error.var[i,1] <- var(Error[,2]-error.rep.ksplit[,i])
  corr[i,1] <- cor(D.rep,error.rep.ksplit[,i])
  W.rep <- as.matrix(weight.semi)%*%forest.cov
  Y.rep <- as.matrix(weight.semi)%*%Y
  reg.rf1<-lm(Y.rep~D.rep+W.rep[,-1]) 
  reg.rf2<-lm(Y.rep~D.rep+W.rep) 
  Coef.matrix.inter[i,3]<-coef(reg.rf1)[2]
  Coef.matrix.inter[i,4]<-coef(reg.rf2)[2]
  
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
  
  #summary(reg.rf1)
  #summary(reg.rf2)
  ###### try the random forest with D.rep defined as the predicted value
  ###### after applying k-split random forest to D
  ssreg.rf<-lm(Y~D.rep+Z+W[,-1])
  Coef.matrix.inter[i,5]<-coef(ssreg.rf)[2]
  ###### try the random forest with D.rep defined as the predicted value
  ###### after applying k-split random forest to D
  
  ### use random forest without data split
  ### use 2 to denote this
  forest.semi2 <- random.forest(forest.cov, D, mtry, max.depth = 0:5, min.node.size = c(5,10))
  # in-sample predictions
  D.rep2 <- forest.semi2$predicted.values
  D.rep.full[,i] <- D.rep2
  weight.semi2 <- weight.mat(forest.semi2)
  error.rep.full[,i] <- as.matrix(weight.semi2)%*%Error[,2]
  error.var[i,2] <- var(Error[,2]-error.rep.full[,i])
  corr[i,2] <- cor(D.rep2,error.rep.full[,i])
  W.rep2 <- as.matrix(weight.semi2)%*%forest.cov
  Y.rep2 <- as.matrix(weight.semi2)%*%Y
  reg2.rf2<-lm(Y.rep2~D.rep2+W.rep2)
  Coef.matrix.inter[i,6]<-coef(reg2.rf2)[2]
  
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
Cov.mat <- matrix(NA,nsim,ncol(Coef.matrix.inter))
colnames(Cov.mat) <- estimator.names
for (j in 1:ncol(Coef.matrix.inter)) {
  Cov.mat[,j] <- ifelse(Coef.matrix.inter[,j]-1.96*sd[j]<=beta & beta<=Coef.matrix.inter[,j]+1.96*sd[j],yes=1,no=0)
}
apply(Cov.mat,2,mean)