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
###### correlation paramaeter 
rho1=0.5
###### dimension
p=20
folds<-5
##################################################
####please change this n=1000, 5000, 10000 and 20000
n=5000
#### we generate 20 candidate notes
knots <- unlist(create_knots(1, n, folds, 0.01, 0.1, 20, min_knots = 2, max_knots = 100, num = 10, exe_max = FALSE)) 
#### to speed up, you can use nsim=200
nsim<-6
#### change f.index across 1,2,3
f.index<-2
##################################################
a<-0.25
Cov<-(A1gen(rho1,p+1))
mu<-rep(0,p+1)
W.original<-mvrnorm(n, mu, Cov)
W<-pnorm(W.original)
Z<-W[,1]
X<-W[,-1]
f_1 <- function(x){x+a*(2*x^2+2*x^3) -25/12}
f_2 <- function(x){exp(2*a+x) - 2/5 * sinh(5/2)}
f_3 <- function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
hist(X[,1])
beta=1

#alpha <- runif(p-1, -0.5, 0.5)
#gamma <-runif(p,-1, 1)
alpha<-rep(-0.3,p)
gamma<-rep(0.2,p)
tau<-0.5
Coef.matrix<-matrix(NA,nrow=nsim,ncol=5)
for(i in 1: nsim){
  print(i)
  mu.error<-rep(0,2)
  Cov.error<-matrix(c(1,0.5,0.5,1),2,2)
  Error<-mvrnorm(n, mu, Cov)
  if(f.index==1){
    D=f_1(Z)+X%*%alpha+Error[,1]
  }
  if(f.index==2){
    D=f_2(Z)+X%*%alpha+Error[,1]
  }
  if(f.index==3){
    D=f_3(Z)+X%*%alpha+Error[,1]
  }
  Y=D*beta+ tau*Z+ X%*% gamma+Error[,2]
  Data<-cbind(D,W)
  knot.values<-apply(Cross_Validation(Data, folds, knots),2,mean)
  opt.index<-which.min(knot.values)
  knot<-round(knots[opt.index] * n^0.8)
  opt.index
  knot
  MODEL.D <- ESTIMATE(W, D, knot)
  D.rep<- pred(MODEL.D, W)
  MODEL.Y <- ESTIMATE(W, Y, knot)
  Y.rep<- pred(MODEL.Y, W)
  V<-Z
  MODEL.V <- ESTIMATE(W,V,knot)
  V.rep <- pred(MODEL.V, W)
  Cov.total<-cbind(V.rep, X)
  D.resid<-resid(lm(D.rep~Cov.total,-1))
  Y.resid<-resid(lm(Y.rep~Cov.total,-1))
  prop.est<-coef(lm(Y.resid~D.resid))[2]
  D.res<-D-pred(MODEL.D, W)
  Y.res<-Y-pred(MODEL.Y, W)
  sigsq.hat<-mean((Y.res-prop.est*D.res)^2)
  prop.sd<-sqrt(sigsq.hat/sum(D.resid^2))
  prop.sd
  ols.est<-coef(lm(Y~D+X))[2]
  tsls.est<-coef(lm(Y~D.rep+X))[2]
  #tsls<-ivreg(Y ~ D + X | X + Z)
  #coef(tsls)[2]
  Coef.matrix[i,1]<-prop.est
  Coef.matrix[i,2]<-prop.sd
  Coef.matrix[i,3]<-ols.est
  Coef.matrix[i,4]<-tsls.est
  Coef.matrix[i,5]<-(prop.est-1.96*prop.sd<beta)*(prop.est+1.96*prop.sd>beta)
}
#apply(Coef.matrix,2,mean)
#apply(Coef.matrix,2,sd)
sel<-c(1,3,4)
report.vec<-rep(NA,10)
report.vec[1:3]<-abs(apply(Coef.matrix[,sel],2,mean)-beta)
report.vec[4:6]<-apply(Coef.matrix[,sel],2,sd)
report.vec[7:9]<-sqrt(apply((Coef.matrix[,sel]-beta)^2,2,mean))
report.vec[10]<-apply(Coef.matrix[,-sel],2,mean)[2]
library(xtable)
round(report.vec,4)
f.index
a
n
tau