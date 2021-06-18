### Start Date: Mar 6, 2021; Recently updated Date: June 18, 2021
### TSCI: TSCI using Basis Approach and Random Forest


library(MASS)
library(AER)
source('Source-Random-Forest.R', encoding = 'UTF-8')
source('Source-Basis-Approach.R', encoding = "UTF-8")


# function to generate covariance matrix
A1gen<-function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    }
  }
  A1
}


###### dimension change the dimension 10,20
p = 10
####please change this n = 1000, 3000, 5000
n = 1000
### setting, change across 1(for previous 1 and 4), 2(for previous 3 and 6)
f.index = 2
##### change the interaction 0, 0.5, 1, 1.5
inter.val = 0.5
#### a denotes the IV strength, set as 1
a = 1
#### violation index, set as 1 or 2
vio.index = 2
##### tau denotes the violation strength
##### set tau as 1
tau = 1
#### the number of simulation numbers
nsim = 10


##############################
f_1 <- function(x){x+a*(x^2+0.5*x^4) -25/12}
# f_2 <- function(x){exp(2*a+x+0.5*x^3)+a*(x^2+0.5*x^4)}
f_2 <- function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
rho1=0.5
Cov<-(A1gen(rho1,p+1))
mu<-rep(0,p+1)
beta=1
alpha<-as.matrix(rep(-0.3,p))
gamma<-as.matrix(rep(0.2,p))
inter<-as.matrix(c(rep(inter.val,5),rep(0,p-5)))


### number of violation space to consider
Q = 5
estimator.names <- c(paste("RF-q",0:(Q-1),sep=""),paste("Cor-q",0:(Q-1),sep=""))
### the fixed violation space for RF
Coef.matrix.rf<-matrix(NA,nrow=nsim,ncol=2*Q)
sd.matrix.rf<-matrix(NA,nrow=nsim,ncol=2*Q)
colnames(Coef.matrix.rf) <- colnames(sd.matrix.rf) <- estimator.names
### the fixed violation space for basis
Coef.matrix.basis<-matrix(NA,nrow=nsim,ncol=Q)
sd.matrix.basis<-matrix(NA,nrow=nsim,ncol=Q)
colnames(Coef.matrix.basis) <- colnames(sd.matrix.basis) <- paste("Basis-q",0:(Q-1),sep="")


# IV strength test
iv.str <- iv.thol <- matrix(NA,nrow=nsim,ncol=Q)
colnames(iv.str) <- colnames(iv.thol) <- paste("q",0:(Q-1),sep="")


# variance of epsilon
SigmaSqY <- matrix(NA, nsim, Q)
colnames(SigmaSqY) <- estimator.names[1:Q]
# variance of delta
SigmaSqD <- matrix(NA, nsim, 1)
# covariance of epsilon and delta
SigmaYD <- matrix(NA, nsim, Q)
colnames(SigmaYD) <- estimator.names[1:Q]


### selection
Q.max  <- matrix(NA,nsim,2)
colnames(Q.max) <- c("Basis","RF")
q.comp <- q.robust <- SigmaSqY.Qmax <- Q.max
# H <- matrix(NA,nsim,Q-1)
# C.alpha <- matrix(NA,nsim,Q-1)
Coef.robust <- matrix(NA,nsim,6)
sd.robust <- matrix(NA,nsim,6)
colnames(Coef.robust) <- colnames(sd.robust) <-
  c("Basis-comp","Basis-robust","RF-comp","RF-Cor-comp","RF-robust","RF-Cor-robust")

### oracle estimator
Coef.oracle <- sd.oracle <- matrix(NA,nsim,3)
colnames(Coef.oracle) <- colnames(sd.oracle) <- c("Basis","RF","RF-Cor")


### naive/plain estimators
Coef.naive <- sd.naive <- matrix(NA,nsim,2)
colnames(Coef.naive) <- colnames(sd.naive) <- c("naiveRF","RF-full")


### TSLS estimator
Coef.TSLS <- sd.TSLS <- matrix(NA,nsim,3)
colnames(Coef.TSLS) <- colnames(sd.TSLS) <- c("TSLS","TSLS-Inter","TSLS-Inter-NoAdj")


### weak IV problem
run.OLS <- weak.iv <- matrix(NA,nsim,2)
colnames(run.OLS) <- colnames(weak.iv) <- c("Basis","RF")


for(i in 1:nsim) {
  print(i)
  #### generate the data
  mu.error<-rep(0,2)
  Cov.error<-matrix(c(1,0.5,0.5,1),2,2)
  Error<-mvrnorm(n, mu.error, Cov.error)
  W.original<-mvrnorm(n, mu, Cov)
  W<-pnorm(W.original)
  Z<-W[,1]
  X<-W[,-1]
  ###### generate the data for the treatment variable D
  if(f.index==1){
    D=f_1(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
  }
  if(f.index==2){
    D=f_2(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
  }
  ####### generate the outcome variable
  # if(vio.index==0){
  #   Y=D*beta+ X%*% gamma+Error[,2]
  # }
  if(vio.index==1){
    Y=D*beta+ tau*Z+ X%*% gamma+Error[,2]
  }
  if(vio.index==2){
    Y=D*beta+ tau*(Z^2+Z-1)+ X%*% gamma+Error[,2] # difficult if no Z
  }
  # if(vio.index==3){
  #   Y=D*beta+ tau*Z^3+ X%*% gamma+Error[,2]
  # }
  # if(vio.index==4){
  #   Y=D*beta+ tau*Z+ Z^5/40+ X%*% gamma+Error[,2]
  # }
  #if(vio.index==5){
  #  Y=D*beta+ tau*(Z^2-1)+ 1/20*sin(Z/2)+ X%*% gamma+Error[,2]
  #}


  ### TSLS Estimators
  TSLS <- ivreg(Y~D+X|Z+X)
  TSLS.Inter <- ivreg(Y~D+Z+X|Z+X+Z*X)
  TSLS.Inter.noadj <- ivreg(Y~D+X|Z+X+Z*X)

  Coef.TSLS[i,1] <- coef(TSLS)[2]
  Coef.TSLS[i,2] <- coef(TSLS.Inter)[2]
  Coef.TSLS[i,3] <- coef(TSLS.Inter.noadj)[2]
  sd.TSLS[i,1] <- summary(TSLS)$coefficients[2,2]
  sd.TSLS[i,2] <- summary(TSLS.Inter)$coefficients[2,2]
  sd.TSLS[i,3] <- summary(TSLS.Inter.noadj)$coefficients[2,2]


  ## Basis Approach
  basis.fit <- TSCI.basis.fit(W,D)

  D.rep.basis <- basis.fit$D.rep
  M <- basis.fit$M
  knot <- basis.fit$knot

  outputs.basis <- TSCI.basis.selection(Y, D, W, D.rep.basis, knot, M, Q=Q)
  if (length(outputs.basis$Coef.vec)>=Q) {
    Coef.matrix.basis[i,] <- outputs.basis$Coef.vec
    sd.matrix.basis[i,] <- outputs.basis$sd.vec
  } else {
    Coef.matrix.basis[i,] <- c(outputs.basis$Coef.vec, rep(NA,Q-length(outputs.basis$Coef.vec)))
    sd.matrix.basis[i,] <- c(outputs.basis$sd.vec, rep(NA,Q-length(outputs.basis$sd.vec)))
  }
  Coef.robust[i,1:2] <- outputs.basis$Coef.robust
  sd.robust[i,1:2] <- outputs.basis$sd.robust
  SigmaSqY.Qmax[i,1] <- outputs.basis$SigmaSqY.Qmax
  Q.max[i,1] <- outputs.basis$Q.max; q.comp[i,1] <- outputs.basis$q.comp; q.robust[i,1] <- outputs.basis$q.robust;
  Coef.oracle[i,1] <- outputs.basis$Coef.vec[vio.index+1]
  sd.oracle[i,1] <- outputs.basis$sd.vec[vio.index+1]
  run.OLS[i,1] <- outputs.basis$run.OLS; weak.iv[i,1] <- outputs.basis$weak.iv


  ### random forest based methods
  forest.cov<-cbind(Z,W[,-1])
  # mtry from p/3 to 2*p/3
  mtry = seq(round((p+1)/3), round(2*(p+1)/3), by=1)
  ### set max.depth and min.node.size for tuning
  ### larger max.depth and smaller min.node.size means more complex trees
  max.depth <- 0; min.node.size <- c(5,10,20)


  ### Data splitting random forest
  # use 2 to denote 2split
  forest.2 <- TSRF.fit(forest.cov,D,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size,split.prop=2/3)
  A1.ind <- forest.2$A1.ind
  # weight matrix
  weight.2 <- TSRF.weight(forest.2$nodes.A1)
  rm(forest.2)

  Cov.aug <- W[,-1]
  for (q in 1:(Q-1)) {
    Cov.aug<-cbind(Z^q,Cov.aug)
  }


  ### selection
  outputs.2 <- TSRF.Selection(Y, D, Cov.aug, A1.ind, weight.2, Q=Q)
  ### outputs
  Coef.matrix.rf[i,] <- outputs.2$Coef.vec
  sd.matrix.rf[i,] <- outputs.2$sd.vec
  SigmaSqY[i,] <- outputs.2$SigmaSqY
  SigmaSqD[i] <- outputs.2$SigmaSqD
  SigmaYD[i,] <- outputs.2$SigmaYD
  iv.str[i,] <- outputs.2$iv.str; iv.thol[i,] <- outputs.2$iv.thol;

  Coef.robust[i,3:6] <- outputs.2$Coef.robust; sd.robust[i,3:6] <- outputs.2$sd.robust;
  SigmaSqY.Qmax[i,2] <- outputs.2$SigmaSqY.Qmax
  Q.max[i,2] <- outputs.2$Q.max; q.comp[i,2] <- outputs.2$q.comp; q.robust[i,2] <- outputs.2$q.robust;
  Coef.oracle[i,2] <- Coef.matrix.rf[i,vio.index+1]
  Coef.oracle[i,3] <- Coef.matrix.rf[i,vio.index+Q+1]
  sd.oracle[i,2] <- sd.matrix.rf[i,vio.index+1]
  sd.oracle[i,3] <- sd.matrix.rf[i,vio.index+Q+1]
  run.OLS[i,2] <- outputs.2$run.OLS; weak.iv[i,2] <- outputs.2$weak.iv


  ## random forest using plug-in method with data in A1
  naiveRF.outputs <- naiveRF.stat(Y, D, Cov.aug[,-(1:(Q-1-vio.index))], A1.ind, weight.2)
  Coef.naive[i,1] <- naiveRF.outputs$betaHat
  sd.naive[i,1] <- naiveRF.outputs$sd
  rm(weight.2)


  ### random forest with full data
  forest.full <- TSRF.full(Cov.aug[,-(1:(Q-1-vio.index))],D,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size)
  weight.full <- TSRF.weight(forest.full$nodes)
  full.outputs <- TSRF.stat.full(Y,D,forest.cov,weight.full)
  Coef.naive[i,2] <- full.outputs$betaHat
  sd.naive[i,2] <- full.outputs$sd
  rm(weight.full)
  rm(forest.full)

}


apply(Coef.matrix.rf,2,mean)
apply(sd.matrix.rf,2,mean)
### coverage of fixed space for ramdom forest
Cov.mat.rf <- matrix(NA,nsim,ncol(sd.matrix.rf))
colnames(Cov.mat.rf) <- colnames(sd.matrix.rf)
for (j in 1:ncol(sd.matrix.rf)) {
  Cov.mat.rf[,j] <- ifelse(Coef.matrix.rf[,j]-1.96*sd.matrix.rf[,j]<=beta &
                          beta<=Coef.matrix.rf[,j]+1.96*sd.matrix.rf[,j],1,0)
}
apply(Cov.mat.rf,2,mean)


Coef.matrix.basis <- na.omit(Coef.matrix.basis)
sd.matrix.basis <- na.omit(sd.matrix.basis)
apply(Coef.matrix.basis,2,mean)
apply(sd.matrix.basis,2,mean)
### coverage of fixed space for basis
Cov.mat.basis <- matrix(NA,nrow(sd.matrix.basis),ncol(sd.matrix.basis))
colnames(Cov.mat.basis) <- colnames(sd.matrix.basis)
for (j in 1:ncol(sd.matrix.basis)) {
  Cov.mat.basis[,j] <- ifelse(Coef.matrix.basis[,j]-1.96*sd.matrix.basis[,j]<=beta &
                          beta<=Coef.matrix.basis[,j]+1.96*sd.matrix.basis[,j],1,0)
}
apply(Cov.mat.basis,2,mean)


### coverage of oracle violation space
Cov.oracle <- matrix(NA,nsim,ncol(sd.oracle))
colnames(Cov.oracle) <- colnames(sd.oracle)
for (j in 1:ncol(sd.oracle)) {
  Cov.oracle[,j] <- ifelse(Coef.oracle[,j]-1.96*sd.oracle[,j]<=beta &
                          beta<=Coef.oracle[,j]+1.96*sd.oracle[,j],1,0)
}
apply(Cov.oracle,2,mean)


### coverage of selection method
apply(Coef.robust,2,mean)
apply(sd.robust,2,mean)

Cov.robust <- matrix(NA,nsim,ncol(sd.robust))
colnames(Cov.robust) <- colnames(sd.robust)
for (j in 1:ncol(sd.robust)) {
  Cov.robust[,j] <- ifelse(Coef.robust[,j]-1.96*sd.robust[,j]<=beta &
                             beta<=Coef.robust[,j]+1.96*sd.robust[,j],1,0)
}
apply(Cov.robust,2,mean)


## coverage of naive methods
Cov.naive <- matrix(NA,nsim,ncol(sd.naive))
colnames(Cov.naive) <- colnames(sd.naive)
for (j in 1:ncol(sd.naive)) {
  Cov.naive[,j] <- ifelse(Coef.naive[,j]-1.96*sd.naive[,j]<=beta &
                             beta<=Coef.naive[,j]+1.96*sd.naive[,j],1,0)
}
apply(Cov.naive,2,mean)


### coverage of TSLS
Cov.TSLS <- matrix(NA,nsim,ncol(sd.TSLS))
colnames(Cov.TSLS) <- colnames(sd.TSLS)
for (j in 1:ncol(sd.TSLS)) {
  Cov.TSLS[,j] <- ifelse(Coef.TSLS[,j]-1.96*sd.TSLS[,j]<=beta &
                           beta<=Coef.TSLS[,j]+1.96*sd.TSLS[,j],1,0)
}
apply(Cov.TSLS,2,mean)


filename <- paste("TSCI-ConIV-Setting",f.index,"-Interaction",inter.val,"-Violation",vio.index,"-p",p,"-n",n,".RData",sep="")
save.image(filename)

