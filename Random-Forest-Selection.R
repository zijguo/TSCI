### Start Date: Mar 6, 2021; Recently updated Date: April 3, 2021
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
n = 1000
### setting, change across 1, 2, 3, 4, 5, 6, 7
f.index = 5
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
f_4 <- function(x){exp(a*x)+a*(x^2+0.5*x^4)}
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
estimator.names <- c(paste("RF-2split-q",0:(Q-1),sep=""),paste("Cor-2split-q",0:(Q-1),sep=""))
Coef.matrix.inter<-matrix(NA,nrow=nsim,ncol=2*Q)
sd.matrix.inter<-matrix(NA,nrow=nsim,ncol=2*Q)
colnames(Coef.matrix.inter) <- colnames(sd.matrix.inter) <- estimator.names


# IV strength test
iv.str <- iv.thol <- matrix(NA,nrow=nsim,ncol=Q)
colnames(iv.str) <- colnames(iv.thol) <- paste("q",0:(Q-1),sep="")
# signal strength test for data split estimator
signal.str <- signal.thol <- matrix(NA,nrow=nsim,ncol=Q)
colnames(signal.str) <- colnames(signal.thol) <- paste("q",0:(Q-1),sep="")
# signal strength test for bias corrected estimator
signal.str.cor <- signal.thol.cor <- matrix(NA,nrow=nsim,ncol=Q)
colnames(signal.str.cor) <- colnames(signal.thol.cor) <- paste("q",0:(Q-1),sep="")


# variance of epsilon
SigmaSqY <- matrix(NA, nsim, Q)
colnames(SigmaSqY) <- estimator.names[1:Q]
# variance of delta
SigmaSqD <- matrix(NA, nsim, 1)
colnames(SigmaSqD) <- c("RF-2split")
# covariance of epsilon and delta
SigmaYD <- matrix(NA, nsim, Q)
colnames(SigmaYD) <- estimator.names[1:Q]


# the numerator of variance estimate, saved to speed up computation
numerator <- matrix(NA, nsim, Q)
colnames(numerator) <- paste("q",0:(Q-1),sep = "")
# trace of TRF.V and (TRF.V)^2
tr.TRF2 <- tr.TRF <- matrix(NA,nsim,Q)
colnames(tr.TRF) <- colnames(tr.TRF2) <- paste("q",0:(Q-1),sep = "")


### selection
Q.max <- rep(NA,nsim)
qhat.c <- rep(NA,nsim)
qhat.r <- rep(NA,nsim)
H <- matrix(NA,nsim,Q-1)
SigmaSqY.max <- rep(NA,nsim)
C.alpha <- matrix(NA,nsim,Q-1)
Coef.robust <- matrix(NA,nsim,2)
sd.robust <- matrix(NA,nsim,2)
colnames(Coef.robust) <- colnames(sd.robust) <- c("RF-2split","Cor-2split")
# TSLS validity 
validity <- rep(0,nsim)


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
  if(f.index==7){
    D=f_4(Z)+X%*%alpha+Error[,1]
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
  # mtry from 1 to p if p<10, from p/3 to 2*p/3 if p>10
  if (p+1<=11) {
    mtry = 1:(p+1)
  } else {
    mtry = seq(round((p+1)/3), round(2*(p+1)/3), by=1)
  }
  ### set max.depth and min.node.size for tuning
  ### larger max.depth and smaller min.node.size means more complex trees
  max.depth <- 0; min.node.size <- c(5,10,20)
  
  
  ### Data split random forest
  # use 2 to denote 2split
  forest.2 <- rf.2split(forest.cov,D,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size)
  # weight matrix
  weight.2 <- weight.2split(forest.2$nodes.A1)
  D.rep.2 <- as.matrix(weight.2)%*%D[forest.2$A1.ind]
  SigmaSqD[i,1] <- mean((D[forest.2$A1.ind]-D.rep.2)^2)
  print(SigmaSqD[i,1])
  Cov.aug<-cbind(Z^4,Z^3,Z^2,Z,W[,-1])
  W.rep.2 <- as.matrix(weight.2)%*%Cov.aug[forest.2$A1.ind,]
  Y.rep.2 <- as.matrix(weight.2)%*%Y[forest.2$A1.ind]
  
  
  
  TRF.V <- rep(list(NA),Q)
  # conduct the computation for each violation space
  for (q in 0:(Q-1)) {
    if (q==Q-1) {
      reg2.rf<-lm(Y.rep.2~D.rep.2+W.rep.2)
      Coef.matrix.inter[i,q+1]<-coef(reg2.rf)[2]
      SigmaSqY[i,q+1] <- get.sigma(coef(reg2.rf)[2], Y, D, cbind(1,Cov.aug))
      stat.result2 <- statRF.2split(Y.rep.2, D.rep.2, Y[forest.2$A1.ind], D[forest.2$A1.ind],
                                    cbind(1,W.rep.2), cbind(1,Cov.aug[forest.2$A1.ind,]),
                                    coef(reg2.rf)[2], weight.2, n,
                                    SigmaSqY[i,q+1], SigmaSqD[i,1])
      Coef.matrix.inter[i,q+1+Q] <- stat.result2$beta.cor
      SigmaYD[i, q+1] <- stat.result2$SigmaYD
      sd.matrix.inter[i,q+1] <- stat.result2$Sd
      iv.str[i,q+1] <- stat.result2$iv.str; iv.thol[i,q+1] <- stat.result2$iv.thol;
      signal.str[i,q+1] <- stat.result2$signal.str; signal.thol[i,q+1] <- stat.result2$signal.thol;
      signal.str.cor[i,q+1] <- stat.result2$signal.str.cor; signal.thol.cor[i,q+1] <- stat.result2$signal.thol.cor;
      tr.TRF[i,q+1] <- stat.result2$tr.TRF
      tr.TRF2[i,q+1] <- stat.result2$tr.TRF2
      numerator[i,q+1] <- stat.result2$numerator
      TRF.V[[q+1]] <- stat.result2$TRF.V
    } else {
      reg2.rf<-lm(Y.rep.2~D.rep.2+W.rep.2[,-(1:(Q-1-q))])
      Coef.matrix.inter[i,q+1]<-coef(reg2.rf)[2]
      SigmaSqY[i,q+1] <- get.sigma(coef(reg2.rf)[2], Y, D, cbind(1,Cov.aug[,-(1:(Q-1-q))]))
      stat.result2 <- statRF.2split(Y.rep.2, D.rep.2, Y[forest.2$A1.ind], D[forest.2$A1.ind],
                                    cbind(1,W.rep.2[,-(1:(Q-1-q))]), cbind(1,Cov.aug[forest.2$A1.ind,-(1:(Q-1-q))]),
                                    coef(reg2.rf)[2], weight.2, n,
                                    SigmaSqY[i,q+1],SigmaSqD[i,1])
      Coef.matrix.inter[i,q+1+Q] <- stat.result2$beta.cor
      SigmaYD[i, q+1] <- stat.result2$SigmaYD
      sd.matrix.inter[i,q+1] <- stat.result2$Sd
      sd.matrix.inter[i,q+1+Q] <- stat.result2$Sd.cor
      iv.str[i,q+1] <- stat.result2$iv.str; iv.thol[i,q+1] <- stat.result2$iv.thol;
      signal.str[i,q+1] <- stat.result2$signal.str; signal.thol[i,q+1] <- stat.result2$signal.thol;
      signal.str.cor[i,q+1] <- stat.result2$signal.str.cor; signal.thol.cor[i,q+1] <- stat.result2$signal.thol.cor;
      tr.TRF[i,q+1] <- stat.result2$tr.TRF
      tr.TRF2[i,q+1] <- stat.result2$tr.TRF2
      numerator[i,q+1] <- stat.result2$numerator
      TRF.V[[q+1]] <- stat.result2$TRF.V
    }
  }
  
  
  ### selection
  ### all of q here are from 0 to 4, so use q+1 to index the columns
  Q.max[i] <- sum(iv.str[i,]>=iv.thol[i,])-1
  ### noise level defined by V_{Q.max}
  SigmaSqY.max[i] <- SigmaSqY[i,Q.max[i]+1]
  ### compute CRF.alpha
  # difference between beta
  beta.diff <- rep(NA,Q-1)
  # the threshold of difference
  thol <- rep(NA,Q-1)
  for (q in 0:(Q-2)) {
    ### H can be negative in some settings...
    H[i,q+1] <- numerator[i,q+1]/(iv.str[i,q+1]^2)+numerator[i,q+2]/(iv.str[i,q+2]^2)-
      2*t(D[forest.2$A1.ind])%*%TRF.V[[q+2]]%*%TRF.V[[q+1]]%*%D[forest.2$A1.ind]/(iv.str[i,q+1]*iv.str[i,q+2])
    beta.diff[q+1] <- abs(Coef.matrix.inter[i,q+2]-Coef.matrix.inter[i,q+1])
    thol[q+1] <- qnorm(1-0.05/log(n))*SigmaSqY.max[i]*sqrt(H[i,q+1])
  }
  C.alpha[i,] <- ifelse(beta.diff<=thol,0,1)
  if (sum(C.alpha[i,])==Q-1) {
    qhat.c[i] <- Q.max[i]-1
  } else {
    qhat.c[i] <- min(which(C.alpha[i,]==0))-1
  }
  qhat.r[i] <- min(qhat.c[i]+1, Q.max[i])
  Coef.robust[i,1] <- Coef.matrix.inter[i,qhat.r[i]+1]
  Coef.robust[i,2] <- Coef.matrix.inter[i,qhat.r[i]+Q+1]
  sd.robust[i,1] <- sd.matrix.inter[i,qhat.r[i]+1]
  sd.robust[i,2] <- sd.matrix.inter[i,qhat.r[i]+Q+1]
  
  if(qhat.c[i]>=1) {
    validity[i] <- 1
  }
  
  
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


### coverage of robust selection
Cov.robust <- matrix(NA,nsim,ncol(sd.robust))
colnames(Cov.robust) <- colnames(sd.robust)
for (j in 1:ncol(sd.robust)) {
  Cov.robust[,j] <- ifelse(Coef.robust[,j]-1.96*sd.robust[,j]<=beta &
                           beta<=Coef.robust[,j]+1.96*sd.robust[,j],1,0)
}
apply(Cov.orac,2,mean)



filename <- paste("RF-Continuous-IV-Setting",f.index,"-Interaction",inter.val,"-p",p,"-n",n,".RData",sep="")
save.image(filename)

