### Start Date: Mar 6, 2021; Recently updated Date: Mar 6, 2021
### TSCI: Random Forest with Binary IV and Continuous Treatment
### Setting: with interaction, compare random forest and regression model


library(MASS)
library(AER)
source('Source-Random-Forest.R', encoding = 'UTF-8')


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

# function to easily compute the SigmaSqY
get.sigma <- function(beta, Y, D, VW) {
  n <- length(Y)
  Portho.VW <- diag(1,n,n) - VW%*%solve(t(VW)%*%VW)%*%t(VW)
  sum((Portho.VW%*%(Y-beta*D))^2)/(n-ncol(VW))
}

###### dimension change the dimension 5,10,20
p = 10
####please change this n = 1000, 2000
n = 1000
##### change the interaction 0.5, 1, 1.5
inter.val = 1
#### a denotes the IV strength, set as 1
a = 1
#### violation index, set as 1 or 2
vio.index = 1
##### tau denotes the violation strength 
##### set tau as 1 
tau = 1
#### the number of simulation numbers
nsim = 20


##############################
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


estimator.names <- c("RF-2split-1","RF-2split-2","Cor-2split-1","Cor-2split-2","RF-kcross-1","RF-kcross-2",
                     "RF-TSLS-2","RF-full-2","TSLS-1","TSLS-2")
Coef.matrix.inter<-matrix(NA,nrow=nsim,ncol=10)
sd.matrix.inter<-matrix(NA,nrow=nsim,ncol=10)
colnames(Coef.matrix.inter) <- colnames(sd.matrix.inter) <- estimator.names


ivtest.2split <- matrix(NA,nrow=nsim,ncol=4)
colnames(ivtest.2split) <- c("2split-1-LHS","2split-1-RHS", "2split-2-LHS","2split-2-RHS")
SignalTest.2split <- matrix(NA,nrow=nsim,ncol=4)
colnames(SignalTest.2split) <- c("RF-2split-1-LHS","RF-2split-1-RHS", "RF-2split-2-LHS","RF-2split-2-RHS")
SignalTest.cor.2split <- matrix(NA,nrow=nsim,ncol=4)
colnames(SignalTest.cor.2split) <- c("Cor-2split-1-LHS","Cor-2split-1-RHS", "Cor-2split-2-LHS","Cor-2split-2-RHS")


SigmaSqY <- matrix(NA, nsim, 2)
colnames(SigmaSqY) <- c("2split-1","2split-2")
var.delta <- matrix(NA, nsim, 5)
colnames(var.delta) <- c("RF-2split","RF-kcross","RF-TSLS","RF-full","TSLS")
SigmaYD <- matrix(NA, nsim, 2)
colnames(SigmaYD) <- c("2split-1","2split-2")


decomp <-array(NA,dim = c(nsim,3,2))


tr.TRF2 <- tr.TRF <- matrix(NA,nsim,2)
colnames(tr.TRF) <- colnames(tr.TRF2) <- c("2split-1","2split-2")


### simulations
for(i in 1:nsim){
  print(i)
  #### generate the data
  mu.error<-rep(0,2)
  Cov.error<-matrix(c(1,0.5,0.5,1),2,2)
  Error<-mvrnorm(n, mu.error, Cov.error)
  W.original<-mvrnorm(n, mu, Cov)
  W<-pnorm(W.original)
  Z<-W[,1]
  Z<-(Z>0.6)
  X<-W[,-1]
  ###### generate the data for the treatment variable D
  D=Z*a+X%*%alpha+Z*X%*%inter+Error[,1]
  
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
  
  
  ### TSLS estimator
  TSLS1 <- ivreg(Y~D+W[, -1]|Z+W[, -1]+Z*W[, -1])
  TSLS2 <- ivreg(Y~D+Z+W[, -1]|Z+W[, -1]+Z*W[, -1])
  Coef.matrix.inter[i,9] <- coef(TSLS1)[2]
  Coef.matrix.inter[i,10] <- coef(TSLS2)[2]
  sd.matrix.inter[i,9] <- summary(TSLS1)$coefficients[2,2]
  sd.matrix.inter[i,10] <- summary(TSLS2)$coefficients[2,2]
  
  ### random forest based methods
  forest.cov<-cbind(Z,W[,-1])
  ### set mtry according to p to save time
  if (p+1<=11) {
    mtry = 1:(p+1)
  } else {
    mtry = seq(round((p+1)/3), round(2*(p+1)/3), by=1)
  }
  ### set max.depth and min.node.size for tuning
  max.depth <- 0; min.node.size <- c(5,10,20)
  
  
  ### Data split random forest
  forest.2 <- rf.2split(forest.cov,D,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size)
  weight.2 <- weight.2split(forest.2$nodes.A1)
  D.rep.2 <- as.matrix(weight.2)%*%D[forest.2$A1.ind]
  var.delta[i,1] <- mean((D[forest.2$A1.ind]-D.rep.2)^2)
  print(var.delta[i,1])
  W.rep.2 <- as.matrix(weight.2)%*%forest.cov[forest.2$A1.ind,]
  Y.rep.2 <- as.matrix(weight.2)%*%Y[forest.2$A1.ind]
  ### reg2 means using two data splits
  reg2.rf1 <- lm(Y.rep.2~D.rep.2+W.rep.2[,-1])
  reg2.rf2 <- lm(Y.rep.2~D.rep.2+W.rep.2)
  Coef.matrix.inter[i,1] <- coef(reg2.rf1)[2]
  Coef.matrix.inter[i,2] <- coef(reg2.rf2)[2]
  
  SigmaSqY[i,1] <- get.sigma(coef(reg2.rf1)[2], Y, D, cbind(1,forest.cov[,-1]))
  SigmaSqY[i,2] <- get.sigma(coef(reg2.rf2)[2], Y, D, cbind(1,forest.cov))
  
  stat.result21 <- statRF.2split(Y.rep.2, D.rep.2, Y[forest.2$A1.ind], D[forest.2$A1.ind], 
                                 cbind(1,W.rep.2[,-1]), cbind(1,forest.cov[forest.2$A1.ind,-1]), 
                                 coef(reg2.rf1)[2], weight.2, n, 
                                 SigmaSqY[i,1], var.delta[i,1])
  stat.result22 <- statRF.2split(Y.rep.2, D.rep.2, Y[forest.2$A1.ind], D[forest.2$A1.ind], 
                                 cbind(1,W.rep.2), cbind(1,forest.cov[forest.2$A1.ind,]), 
                                 coef(reg2.rf2)[2], weight.2, n, 
                                 SigmaSqY[i,2], var.delta[i,1])
  Coef.matrix.inter[i,3] <- stat.result21$beta.cor; Coef.matrix.inter[i,4] <- stat.result22$beta.cor;
  SigmaYD[i,1] <- stat.result21$SigmaYD; SigmaYD[i,2] <- stat.result22$SigmaYD;
  sd.matrix.inter[i,1] <- stat.result21$Sd; sd.matrix.inter[i,2] <- stat.result22$Sd;
  sd.matrix.inter[i,3] <- stat.result21$Sd.cor; sd.matrix.inter[i,4] <- stat.result22$Sd.cor;
  ivtest.2split[i,1] <- stat.result21$LHS; ivtest.2split[i,2] <- stat.result21$RHS;
  ivtest.2split[i,3] <- stat.result22$LHS; ivtest.2split[i,4] <- stat.result22$RHS;
  SignalTest.2split[i,1] <- stat.result21$SignalTest.LHS; SignalTest.2split[i,2] <- stat.result21$SignalTest.RHS;
  SignalTest.2split[i,3] <- stat.result22$SignalTest.LHS; SignalTest.2split[i,4] <- stat.result22$SignalTest.RHS;
  SignalTest.cor.2split[i,1] <- stat.result21$SignalTest.cor.LHS; SignalTest.cor.2split[i,2] <- stat.result21$SignalTest.cor.RHS;
  SignalTest.cor.2split[i,3] <- stat.result22$SignalTest.cor.LHS; SignalTest.cor.2split[i,4] <- stat.result22$SignalTest.cor.RHS;
  decomp[i,1,1] <- stat.result21$term1;decomp[i,2,1] <- stat.result21$term2;decomp[i,3,1] <- stat.result21$term3
  decomp[i,1,2] <- stat.result22$term1;decomp[i,2,2] <- stat.result22$term2;decomp[i,3,2] <- stat.result22$term3
  tr.TRF[i,1] <- stat.result21$tr.TRF;tr.TRF[i,2] <- stat.result22$tr.TRF;
  tr.TRF2[i,1] <- stat.result21$tr.TRF2;tr.TRF2[i,2] <- stat.result22$tr.TRF2;
  
  
  ### K split cross-fitting random forest
  # forest.k <- rf.kcross(forest.cov,D,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size)
  # weight.k <- weight.kcross(forest.k$nodes)
  # D.rep.k <- as.matrix(weight.k)%*%D
  # var.delta[i,2] <- mean((D-D.rep.k)^2)
  # print(var.delta[i,2])
  # W.rep.k <- as.matrix(weight.k)%*%forest.cov
  # Y.rep.k <- as.matrix(weight.k)%*%Y
  # ### reg2 means using k splits cross-fitting
  # regk.rf1 <- lm(Y.rep.k~D.rep.k+W.rep.k[,-1])
  # regk.rf2 <- lm(Y.rep.k~D.rep.k+W.rep.k)
  # Coef.matrix.inter[i,3] <- coef(regk.rf1)[2]
  # Coef.matrix.inter[i,4] <- coef(regk.rf2)[2]
  # 
  # SigmaSqY[i,3] <- get.sigma(coef(regk.rf1)[2], Y, D, cbind(1,forest.cov[,-1]))
  # SigmaSqY[i,4] <- get.sigma(coef(regk.rf2)[2], Y, D, cbind(1,forest.cov))
  # 
  # stat.resultk1 <- statRF.kcross(D.rep.k, D, cbind(1,W.rep.k[,-1]), weight.k, SigmaSqY[i,3], var.delta[i,2])
  # stat.resultk2 <- statRF.kcross(D.rep.k, D, cbind(1,W.rep.k), weight.k, SigmaSqY[i,4], var.delta[i,2])
  # sd.matrix.inter[i,3] <- stat.resultk1$Sd; sd.matrix.inter[i,4] <- stat.resultk2$Sd;
  # ivtest.kcross[i,1] <- stat.resultk1$LHS; ivtest.kcross[i,2] <- stat.resultk1$RHS;
  # ivtest.kcross[i,3] <- stat.resultk2$LHS; ivtest.kcross[i,4] <- stat.resultk2$RHS;
  
  
  # Two Stage Least Squares Random Forest
  # regTSLS.rf2 <- lm(Y[forest.2$A1.ind]~D.rep.2+forest.cov[forest.2$A1.ind,])
  # Coef.matrix.inter[i,5] <- coef(regTSLS.rf2)[2]
  # SigmaSqY.RF.TSLS <- get.sigma(coef(regTSLS.rf2)[2], Y, D, cbind(1, forest.cov))
  # sd.matrix.inter[i,5] <- statRF.TSLS(D.rep.2, cbind(1,forest.cov[forest.2$A1.ind,]), SigmaSqY.RF.TSLS)$Sd
  # var.delta[i,3] <- var.delta[i,1]
  
  
  # random forest with full data
  # forest.full <- rf.full(forest.cov,D,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size)
  # weight.full <- weight.2split(forest.full$nodes)
  # D.rep.full <- as.matrix(weight.full)%*%D
  # var.delta[i,4] <- mean((D.rep.full-D)^2)
  # W.rep.full <- as.matrix(weight.full)%*%forest.cov
  # Y.rep.full <- as.matrix(weight.full)%*%Y
  # regfull.rf2 <- lm(Y.rep.full~D.rep.full+W.rep.full)
  # Coef.matrix.inter[i,6] <- coef(regfull.rf2)[2]
  # SigmaSqY.full <- get.sigma(coef(regfull.rf2)[2], Y, D, cbind(1, forest.cov))
  # sd.matrix.inter[i,6] <- statRF.full(D.rep.full, D, W.rep.full, weight.full, SigmaSqY.full)$Sd
}


apply(Coef.matrix.inter,2,mean)
apply(Coef.matrix.inter,2,sd)
apply(sd.matrix.inter,2,mean)


### compute coverage
Cov.mat <- matrix(NA,nsim,ncol(sd.matrix.inter))
colnames(Cov.mat) <- colnames(sd.matrix.inter)
for (j in 1:ncol(sd.matrix.inter)) {
  Cov.mat[,j] <- ifelse(Coef.matrix.inter[,j]-1.96*sd.matrix.inter[,j]<=beta & 
                        beta<=Coef.matrix.inter[,j]+1.96*sd.matrix.inter[,j],1,0)

}
apply(Cov.mat,2,mean)


### empirical coverage using oracle sd
Cov.orac <- matrix(NA,nsim,ncol(sd.matrix.inter))
colnames(Cov.orac) <- colnames(sd.matrix.inter)
for (j in 1:ncol(sd.matrix.inter)) {
  Cov.orac[,j] <- ifelse(Coef.matrix.inter[,j]-1.96*apply(Coef.matrix.inter,2,sd)[j]<=beta & 
                         beta<=Coef.matrix.inter[,j]+1.96*apply(Coef.matrix.inter,2,sd)[j],1,0)
}
apply(Cov.orac,2,mean)


