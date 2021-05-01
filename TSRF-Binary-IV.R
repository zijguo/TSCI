### Start Date: April 30, 2021; Recently updated Date: April 30, 2021
### TSCI: Random Forest with Binary IV

### Start Date: Mar 6, 2021; Recently updated Date: April 30, 2021
### TSCI: TSCI using Basis Approach and Random Forest


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


###### dimension change the dimension 5,10,20
p = 10
####please change this n = 1000, 3000, 5000
n = 1000
##### change the interaction 0.5, 1, 1.5
inter.val = 0.5
#### a denotes the IV strength, change over 0.5, 1, 1.5
a = 1
#### violation index, set as 1
vio.index = 1
##### tau denotes the violation strength
##### set tau as 1
tau = 1
#### the number of simulation numbers
nsim = 10

### save difference testing
H <- beta.diff <- rep(list(NA),nsim)
z.alpha <- rep(NA,nsim)

##############################
f_1 <- function(x){x+a*(x^2+0.5*x^4) -25/12}
f_2 <- function(x){exp(2*a+x+0.5*x^3) - 2/5 * sinh(5/2)}
f_3 <- function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
f_4 <- function(x){exp(2*a+x+0.5*x^3)+a*(x^2+0.5*x^4)}
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


Q = 2
estimator.names <- c(paste("RF-q",0:(Q-1),sep=""),paste("Cor-q",0:(Q-1),sep=""))
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
# covariance of epsilon and delta
SigmaYD <- matrix(NA, nsim, Q)
colnames(SigmaYD) <- estimator.names[1:Q]


# the numerator of variance estimate, saved to speed up computation
DT.Sq <- matrix(NA, nsim, Q)
colnames(DT.Sq) <- paste("q",0:(Q-1),sep = "")
# trace of TRF.V and (TRF.V)^2
trace.T2 <- trace.T <- matrix(NA,nsim,Q)
colnames(trace.T) <- colnames(trace.T2) <- paste("q",0:(Q-1),sep = "")


### selection
Q.max  <- matrix(NA,nsim,1)
colnames(Q.max) <- "RF"
qhat.c <- qhat.r <- validity <- SigmaSqY.Qmax <- Q.max
# H <- matrix(NA,nsim,Q-1)
# C.alpha <- matrix(NA,nsim,Q-1)
Coef.robust <- matrix(NA,nsim,4)
sd.robust <- matrix(NA,nsim,4)
colnames(Coef.robust) <- colnames(sd.robust) <- 
  c("RF-comp","RF-Cor-comp","RF-robust","RF-Cor-robust")

### oracle estimator
oracle.est <- oracle.sd <- matrix(NA,nsim,2)
colnames(oracle.est) <- colnames(oracle.sd) <- c("RF","RF-Cor")


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
  
  
  ### Data splitting random forest
  # use 2 to denote 2split
  forest.2 <- TSRF.fit(forest.cov,D,mtry=mtry,max.depth=max.depth,min.node.size=min.node.size)
  A1.ind <- forest.2$A1.ind
  # weight matrix
  weight.2 <- TSRF.weight(forest.2$nodes.A1)
  
  Cov.aug <- W[,-1]
  for (q in 1:(Q-1)) {
    Cov.aug<-cbind(Z^q,Cov.aug)
  }
  
  
  ### selection
  outputs <- TSRF.Selection(Y, D, Cov.aug, A1.ind, weight.2, Q=Q)
  ### outputs
  Coef.matrix.inter[i,] <- outputs$Coef.vec
  sd.matrix.inter[i,] <- outputs$sd.vec
  SigmaSqY[i,] <- outputs$SigmaSqY
  SigmaSqD[i] <- outputs$SigmaSqD
  SigmaYD[i,] <- outputs$SigmaYD
  iv.str[i,] <- outputs$iv.str; iv.thol[i,] <- outputs$iv.thol; DT.Sq[i,] <- outputs$DT.Sq
  signal.str[i,] <- outputs$signal.str; signal.thol[i,] <- outputs$signal.thol;
  signal.str.cor[i,] <- outputs$signal.str.cor; signal.thol.cor[i,] <- outputs$signal.thol.cor;
  trace.T[i,] <- outputs$trace.T;   trace.T2[i,] <- outputs$trace.T2
  
  Coef.robust[i,] <- outputs$Coef.robust; sd.robust[i,] <- outputs$sd.robust;
  # H[i,] <- outputs$H; C.alpha[i,] <- outputs$C.alpha;
  SigmaSqY.Qmax[i,1] <- outputs$SigmaSqY.Qmax
  Q.max[i,1] <- outputs$Q.max; qhat.c[i,1] <- outputs$qhat.c; qhat.r[i,1] <- outputs$qhat.r;
  validity[i,1] <- outputs$validity
  oracle.est[i,1] <- Coef.matrix.inter[i,vio.index+1]
  oracle.est[i,2] <- Coef.matrix.inter[i,vio.index+Q+1]
  oracle.sd[i,1] <- sd.matrix.inter[i,vio.index+1]
  oracle.sd[i,2] <- sd.matrix.inter[i,vio.index+Q+1]
  H[[i]] <- outputs$H
  beta.diff[[i]] <- outputs$beta.diff
  z.alpha[i] <- outputs$z.alpha
  
}


apply(Coef.matrix.inter,2,mean)
# apply(Coef.matrix.inter,2,sd)
apply(sd.matrix.inter,2,mean)


### coverage of fixed space for ramdom forest
Cov.mat <- matrix(NA,nsim,ncol(sd.matrix.inter))
colnames(Cov.mat) <- colnames(sd.matrix.inter)
for (j in 1:ncol(sd.matrix.inter)) {
  Cov.mat[,j] <- ifelse(Coef.matrix.inter[,j]-1.96*sd.matrix.inter[,j]<=beta &
                          beta<=Coef.matrix.inter[,j]+1.96*sd.matrix.inter[,j],1,0)
}
apply(Cov.mat,2,mean)


### coverage of oracle method
Cov.oracle <- matrix(NA,nsim,ncol(oracle.sd))
colnames(Cov.oracle) <- colnames(oracle.sd)
for (j in 1:ncol(oracle.sd)) {
  Cov.oracle[,j] <- ifelse(oracle.est[,j]-1.96*oracle.sd[,j]<=beta &
                             beta<=oracle.est[,j]+1.96*oracle.sd[,j],1,0)
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


filename <- paste("RF-BiIV-Strength",a,"-Violation",vio.index,"-Interaction",inter.val,"-p",p,"-n",n,".RData",sep="")
save.image(filename)