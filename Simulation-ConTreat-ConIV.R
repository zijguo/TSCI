### Start Date: Mar 6, 2021; Recently updated Date: Mar 6, 2021
### TSCI: Random Forest with Binary IV and Continuous Treatment
### Setting: with interaction, compare random forest and regression model


library(MASS)
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
n = 2000
### setting, change across 1, 2, 3, 4, 5, 6
f.index = 3
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


Q = 5
estimator.names <- c(paste("RF-2split-q",1:Q,sep=""),paste("Cor-2split-q",1:Q,sep=""),paste("RF-kcross-q",1:Q,sep=""),"RF-TSLS-2","RF-full-2")
Coef.matrix.inter<-matrix(NA,nrow=nsim,ncol=3*Q+2)
sd.matrix.inter<-matrix(NA,nrow=nsim,ncol=3*Q+2)
colnames(Coef.matrix.inter) <- colnames(sd.matrix.inter) <- estimator.names


# sd.matrix.boot <- matrix(NA,nsim,5)
# colnames(sd.matrix.boot) <- estimator.names[(Q+1):(2*Q)]


ivtest.2split <- matrix(NA,nrow=nsim,ncol=2*Q)
colnames(ivtest.2split) <- paste(rep(paste("q",1:Q,sep=""),each=2),c("LHS","RHS"), sep = "-")
SignalTest.2split <- matrix(NA,nrow=nsim,ncol=2*Q)
colnames(SignalTest.2split) <- paste(rep(estimator.names[1:Q], each = 2), c("LHS","RHS"), sep = "-")
SignalTest.cor.2split <- matrix(NA,nrow=nsim,ncol=2*Q)
colnames(SignalTest.cor.2split) <- paste(rep(estimator.names[(Q+1):(2*Q)], each = 2), c("LHS","RHS"), sep = "-")


SigmaSqY <- matrix(NA, nsim, 2*Q)
colnames(SigmaSqY) <- estimator.names[1:(2*Q)]
var.delta <- matrix(NA, nsim, 4)
colnames(var.delta) <- c("RF-2split","RF-kcross","RF-TSLS","RF-full")
SigmaYD <- matrix(NA, nsim, 2*Q)
colnames(SigmaYD) <- estimator.names[1:(2*Q)]


decomp <-array(NA,dim = c(nsim,3,5))

tr.TRF2 <- tr.TRF <- matrix(NA,nsim,5)
colnames(tr.TRF) <- colnames(tr.TRF2) <- paste("q",1:5,sep = "")


for(i in 1:nsim){
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
  Cov.aug<-cbind(Z^4,Z^3,Z^2,Z,W[,-1])
  W.rep.2 <- as.matrix(weight.2)%*%Cov.aug[forest.2$A1.ind,]
  Y.rep.2 <- as.matrix(weight.2)%*%Y[forest.2$A1.ind]

  for (q in 0:(Q-1)) {
    if (q==Q-1) {
      reg2.rf<-lm(Y.rep.2~D.rep.2+W.rep.2)
      Coef.matrix.inter[i,q+1]<-coef(reg2.rf)[2]
      SigmaSqY[i,q+1] <- get.sigma(coef(reg2.rf)[2], Y, D, cbind(1,Cov.aug))
      stat.result2 <- statRF.2split(Y.rep.2, D.rep.2, Y[forest.2$A1.ind], D[forest.2$A1.ind],
                                    cbind(1,W.rep.2), cbind(1,Cov.aug[forest.2$A1.ind,]),
                                    coef(reg2.rf)[2], weight.2, n,
                                    SigmaSqY[i,q+1], var.delta[i,1])
      Coef.matrix.inter[i,q+1+Q] <- stat.result2$beta.cor
      SigmaYD[i, q+1] <- stat.result2$SigmaYD
      sd.matrix.inter[i,q+1] <- stat.result2$Sd
      sd.matrix.inter[i,q+1+Q] <- stat.result2$Sd.cor
      # sd.matrix.boot[i,q+1] <- stat.result2$Sd.cor.boot
      ivtest.2split[i,2*q+1] <- stat.result2$LHS; ivtest.2split[i,2*(q+1)] <- stat.result2$RHS;
      SignalTest.2split[i,2*q+1] <- stat.result2$SignalTest.LHS; SignalTest.2split[i,2*(q+1)] <- stat.result2$SignalTest.RHS;
      SignalTest.cor.2split[i,2*q+1] <- stat.result2$SignalTest.cor.LHS; SignalTest.cor.2split[i,2*(q+1)] <- stat.result2$SignalTest.cor.RHS;
      decomp[i,1,q+1] <- stat.result2$term1
      decomp[i,2,q+1] <- stat.result2$term2
      decomp[i,3,q+1] <- stat.result2$term3
      tr.TRF[i,q+1] <- stat.result2$tr.TRF
      tr.TRF2[i,q+1] <- stat.result2$tr.TRF2
    } else {
      reg2.rf<-lm(Y.rep.2~D.rep.2+W.rep.2[,-(1:(Q-1-q))])
      Coef.matrix.inter[i,q+1]<-coef(reg2.rf)[2]
      SigmaSqY[i,q+1] <- get.sigma(coef(reg2.rf)[2], Y, D, cbind(1,Cov.aug[,-(1:(Q-1-q))]))
      stat.result2 <- statRF.2split(Y.rep.2, D.rep.2, Y[forest.2$A1.ind], D[forest.2$A1.ind],
                                    cbind(1,W.rep.2[,-(1:(Q-1-q))]), cbind(1,Cov.aug[forest.2$A1.ind,-(1:(Q-1-q))]),
                                    coef(reg2.rf)[2], weight.2, n,
                                    SigmaSqY[i,q+1],var.delta[i,1])
      Coef.matrix.inter[i,q+1+Q] <- stat.result2$beta.cor
      SigmaYD[i, q+1] <- stat.result2$SigmaYD
      sd.matrix.inter[i,q+1] <- stat.result2$Sd
      sd.matrix.inter[i,q+1+Q] <- stat.result2$Sd.cor
      # sd.matrix.boot[i,q+1] <- stat.result2$Sd.cor.boot
      ivtest.2split[i,2*q+1] <- stat.result2$LHS; ivtest.2split[i,2*(q+1)] <- stat.result2$RHS;
      SignalTest.2split[i,2*q+1] <- stat.result2$SignalTest.LHS; SignalTest.2split[i,2*(q+1)] <- stat.result2$SignalTest.RHS;
      SignalTest.cor.2split[i,2*q+1] <- stat.result2$SignalTest.cor.LHS; SignalTest.cor.2split[i,2*(q+1)] <- stat.result2$SignalTest.cor.RHS;
      decomp[i,1,q+1] <- stat.result2$term1
      decomp[i,2,q+1] <- stat.result2$term2
      decomp[i,3,q+1] <- stat.result2$term3
      tr.TRF[i,q+1] <- stat.result2$tr.TRF
      tr.TRF2[i,q+1] <- stat.result2$tr.TRF2
    }
  }


  ### k split cross-fitting random forest
  # forest.k <- rf.kcross(forest.cov,D,k=2,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size)
  # weight.k <- weight.kcross(forest.k$nodes)
  # D.rep.k <- as.matrix(weight.k)%*%D
  # var.delta[i,2] <- mean((D-D.rep.k)^2)
  # print(var.delta[i,2])
  # Cov.aug<-cbind(Z^4,Z^3,Z^2,Z,W[,-1])
  # W.rep.k <- as.matrix(weight.k)%*%Cov.aug
  # Y.rep.k <- as.matrix(weight.k)%*%Y
  #
  # for (q in 0:(Q-1)) {
  #   if (q==Q-1) {
  #     regk.rf<-lm(Y.rep.k~D.rep.k+W.rep.k)
  #     Coef.matrix.inter[i,q+1+2*Q]<-coef(regk.rf)[2]
  #     SigmaSqY[i,q+1+2*Q] <- get.sigma(coef(regk.rf)[2], Y, D, cbind(1,Cov.aug))
  #     stat.resultk <- statRF.kcross(D.rep.k, D, cbind(1,W.rep.k), weight.k, SigmaSqY[i,q+1+2*Q],var.delta[i,2])
  #     sd.matrix.inter[i,q+1+2*Q] <- stat.resultk$Sd
  #     # ivtest.kcross[i,2*q+1] <- stat.resultk$LHS; ivtest.kcross[i,2*(q+1)] <- stat.resultk$RHS;
  #   } else {
  #     regk.rf<-lm(Y.rep.k~D.rep.k+W.rep.k[,-(1:(Q-1-q))])
  #     Coef.matrix.inter[i,q+1+2*Q]<-coef(regk.rf)[2]
  #     SigmaSqY[i,q+1+2*Q] <- get.sigma(coef(regk.rf)[2], Y, D, cbind(1,Cov.aug[,-(1:(Q-1-q))]))
  #     stat.resultk <- statRF.kcross(D.rep.k, D, cbind(1,W.rep.k[,-(1:(Q-1-q))]), weight.k, SigmaSqY[i,q+1+2*Q],var.delta[i,2])
  #     sd.matrix.inter[i,q+1+2*Q] <- stat.resultk$Sd
  #     # ivtest.kcross[i,2*q+1] <- stat.resultk$LHS; ivtest.kcross[i,2*(q+1)] <- stat.resultk$RHS;
  #   }
  # }


  # Two Stage Least Squares Random Forest
  # regTSLS.rf2 <- lm(Y[forest.2$A1.ind]~D.rep.2+forest.cov[forest.2$A1.ind,])
  # Coef.matrix.inter[i,3*Q+1] <- coef(regTSLS.rf2)[2]
  # SigmaSqY.RF.TSLS <- get.sigma(coef(regTSLS.rf2)[2], Y, D, cbind(1, forest.cov))
  # sd.matrix.inter[i,3*Q+1] <- statRF.TSLS(D.rep.2, cbind(1,forest.cov[forest.2$A1.ind,]), SigmaSqY.RF.TSLS)$Sd
  # var.delta[i,3] <- var.delta[i,1]


  # random forest with full data
  # forest.full <- rf.full(forest.cov,D,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size)
  # weight.full <- weight.2split(forest.full$nodes)
  # D.rep.full <- as.matrix(weight.full)%*%D
  # var.delta[i,4] <- mean((D.rep.full-D)^2)
  # W.rep.full <- as.matrix(weight.full)%*%forest.cov
  # Y.rep.full <- as.matrix(weight.full)%*%Y
  # regfull.rf2 <- lm(Y.rep.full~D.rep.full+W.rep.full)
  # Coef.matrix.inter[i,3*Q+2] <- coef(regfull.rf2)[2]
  # SigmaSqY.full <- get.sigma(coef(regfull.rf2)[2], Y, D, cbind(1, forest.cov))
  # sd.matrix.inter[i,3*Q+2] <- statRF.full(D.rep.full, D, W.rep.full, weight.full, SigmaSqY.full)$Sd


}


apply(Coef.matrix.inter,2,mean)
apply(Coef.matrix.inter,2,sd)
apply(sd.matrix.inter,2,mean)
# apply(sd.matrix.boot,2,mean)


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


setwd("D:/Desktop/test")
filename <- paste("RF-Continuous-IV-Setting",f.index,"-Interaction",inter.val,"-p",p,"-n",n,".RData",sep="")
save.image(filename)

