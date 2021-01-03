### Start Date: Aug 11, 2020; Recently updated Date: Aug 11, 2020
### For this code, we conduct the violation space selection
### We set different ways of generating E(D|Z,X) function and set a linear violation
### We use a data-dependent way of choosing among (1)no violation, (2)linear violation, (3) quadratic violation and (4) cubic violation. 
source('~/Dropbox/Projects-Collaboration/Peter-Zijian-Causal/Invalid-IV-Invariance/Simulation/Source-Semiparametric.R', encoding = 'UTF-8')
###### dimension
p=20
folds<-5
##################################################
####please change this n=1000, 5000, 10000 and 20000
n=5000
#### change f.index across 1,2,3
f.index<-3
##### setting.vio across 0,1,2,3,4,5
setting.vio<-4
#### to speed up, you can use nsim=200
nsim<-1000
#### we generate 20 candidate notes
knots <- unlist(create_knots(1, n, folds, 0.01, 0.1, 20, min_knots = 2, max_knots = 100, num = 10, exe_max = FALSE)) 
##################################################
a<-0.25
tau<-2
rho1=0.5
Cov<-(A1gen(rho1,p+1))
mu<-rep(0,p+1)
W.original<-mvrnorm(n, mu, Cov)
W<-pnorm(W.original)
Z<-W[,1]+0.5
X<-W[,-1]
f_1 <- function(x){x+a*(2*x^2+4*x^3) -25/12}
f_2 <- function(x){exp(2*a+x+0.5*x^3) - 2/5 * sinh(5/2)}
f_3 <- function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
beta=1
alpha<-rep(-0.3,p)
gamma<-rep(0.2,p)
###### number of tested basis functions
Q<-5
order<-seq(1,Q,by=1)-1
#max.matrix<-matrix(0,nrow=nsim,ncol=2)
Coef.matrix<-matrix(NA,nrow=nsim,ncol=2+Q)
Prop.sel.matrix<-matrix(0,nrow=nsim,ncol=14)
SE.matrix<-matrix(NA,nrow=nsim,ncol=2+Q)
Cov.matrix<-matrix(NA,nrow=nsim,ncol=2+Q)
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
  if(setting.vio==0){
    Y=D*beta+ X%*% gamma+Error[,2]
  }
  if(setting.vio==1){
  Y=D*beta+ tau*Z+ X%*% gamma+Error[,2]
  }
  if(setting.vio==2){
    Y=D*beta+ tau*(Z^2-1)+ X%*% gamma+Error[,2]
  }
  if(setting.vio==3){
    Y=D*beta+ tau*Z^3+ X%*% gamma+Error[,2]
  }
  if(setting.vio==4){
    Y=D*beta+ tau*(Z^2-1)+ Z^5/40+ X%*% gamma+Error[,2]
  }
  if(setting.vio==5){
    Y=D*beta+ tau*(Z^2-1)+ 1/20*sin(Z/2)+ X%*% gamma+Error[,2]
  }
  ##### ols and TSLS
  ols<-lm(Y~D+X)
  tsls<-ivreg(Y~D+X|Z+X)
  Coef.matrix[i,1]<-coef(ols)[2]
  Coef.matrix[i,2]<-coef(tsls)[2]
  SE.matrix[i,1]<-sqrt(vcov(ols)[2,2])
  SE.matrix[i,2]<-sqrt(vcov(tsls)[2,2])
  ##### choose the best knots for D~Z+ X
  Data<-cbind(D,W)
  knot.values<-apply(Cross_Validation(Data, folds, knots),2,mean)
  opt.index<-which.min(knot.values)
  knot<-round(knots[opt.index] * n^0.8)
  ##### Fit the semiparametric D model
  MODEL.D <- ESTIMATE(W, D, knot)
  D.rep<- pred(MODEL.D, W)
  resid<-D-D.rep
  sd.D=sqrt(sum(resid^2)/length(resid))
  ###### test signal strength in TSCI
  M<-knot+2
  Q.max<-min(Q, M-1)
  str.vec<-rep(NA,Q)
  thre.vec<-rep(NA,Q)
  prop.vec<-rep(NA,Q)
  se.vec<-rep(NA,Q)
  inver.design<-rep(NA,Q)
  noise.vec<-rep(NA,Q)
  pred.error<-rep(NA,Q)
  for(index in 1:Q.max){
    q=order[index]
    #print(q)
    if(q==0){
      V=NULL
    }else{
      V=poly(Z,q)
    }
    if(is.null(V)){
      Cov.total<-W[,-1]
    }else{
      V.rep=V
      q=dim(V)[2]
      for(l in 1:q){
        MODEL.V<-ESTIMATE(W,V[,l],knot)
        V.rep[,l] <- pred(MODEL.V, W)
      }
      Cov.total<-cbind(V.rep, W[,-1])
    }
    D.resid<-resid(lm(D.rep~Cov.total,-1))
    str.vec[index]<-sum(D.resid^2)/(sd.D^2)
    #print(str.vec)
    taun<-1/log(M-q)
    thre.vec[index]<-(M-q)^2+sqrt((M-q)*2*log(M-q))*(1+2*max(taun,sqrt(str.vec[index]/(M-q))))
    if(str.vec[index]>min(thre.vec[index],10)){
      MODEL.Y <- ESTIMATE(W, Y, knot)
      Y.rep<- pred(MODEL.Y, W)
      D.resid<-resid(lm(D.rep~Cov.total,-1))
      Y.resid<-resid(lm(Y.rep~Cov.total,-1))
      ### estimate the point estimator
      prop.vec[index]<-sum(Y.resid*D.resid)/sum(D.resid^2)
      ### estimate the standard error
      D.res<-D-pred(MODEL.D, W)
      Y.res<-Y-pred(MODEL.Y, W)
      sigsq.hat<-mean((Y.res-prop.vec[index]*D.res)^2)
      noise.vec[index]<-sigsq.hat
      inver.design[index]<-1/sum(D.resid^2)
      scale<-1.05
      se.vec[index]<-scale*sqrt(sigsq.hat/sum(D.resid^2))
      pred.error[index]<-sum((Y.resid-prop.vec[index]*D.resid)^2)
    }
  }  
  Coef.matrix[i,3:(2+Q)]<-prop.vec
  SE.matrix[i,3:(2+Q)]<-se.vec
  Cov.matrix[i,]<-(Coef.matrix[i,]-1.96*SE.matrix[i,]<beta)*(Coef.matrix[i,]+1.96*SE.matrix[i,]>beta)
  test.diff<-abs(prop.vec[-1]-prop.vec[-Q])
  test.number<-length(na.omit(test.diff))
  test.threshold<-qnorm(1-0.025)*na.omit(noise.vec)[length(na.omit(noise.vec))]*sqrt(inver.design[-1]-inver.design[-Q])
  if(length(which(test.diff<test.threshold))==0){
  q.comp=length(na.omit(prop.vec))
  Prop.sel.matrix[i,5]<-1
  }else{
  q.comp<-min(which(test.diff<test.threshold))
  }
  Prop.sel.matrix[i,1]<-prop.vec[q.comp]
  Prop.sel.matrix[i,2]<-se.vec[q.comp]
  Prop.sel.matrix[i,3]<-(Prop.sel.matrix[i,1]-1.96*Prop.sel.matrix[i,2]<beta)*(Prop.sel.matrix[i,1]+1.96*Prop.sel.matrix[i,2]>beta)
  Prop.sel.matrix[i,4]<-q.comp
  test.pred<-pred.error[1:length(na.omit(noise.vec))]/na.omit(noise.vec)[length(na.omit(noise.vec))]
  test.pred.thre<-qchisq(1-0.05, M-order[1:length(na.omit(noise.vec))])
  if(length(which(test.pred<test.pred.thre))==0){
    q.pred=length(na.omit(prop.vec))
    Prop.sel.matrix[i,10]<-1
  }else{
    q.pred<-min(which(test.pred<test.pred.thre))
  }
  Prop.sel.matrix[i,6]<-prop.vec[q.pred]
  Prop.sel.matrix[i,7]<-se.vec[q.pred]
  Prop.sel.matrix[i,8]<-(Prop.sel.matrix[i,6]-1.96*Prop.sel.matrix[i,7]<beta)*(Prop.sel.matrix[i,6]+1.96*Prop.sel.matrix[i,7]>beta)
  Prop.sel.matrix[i,9]<-q.pred
  q.conser<-min(q.comp+1,length(na.omit(prop.vec)))
  Prop.sel.matrix[i,11]<-prop.vec[q.conser]
  Prop.sel.matrix[i,12]<-se.vec[q.conser]
  Prop.sel.matrix[i,13]<-(Prop.sel.matrix[i,11]-1.96*Prop.sel.matrix[i,12]<beta)*(Prop.sel.matrix[i,11]+1.96*Prop.sel.matrix[i,12]>beta)
  Prop.sel.matrix[i,14]<-q.conser
}
abs(apply(Coef.matrix,2,mean)-beta)
apply(Coef.matrix,2,sd)
apply(SE.matrix,2,mean)
sqrt(apply((Coef.matrix-beta)^2,2,mean))
apply(Cov.matrix,2,mean)
apply(Prop.sel.matrix,2,mean)
f.index
setting.vio
a
n
tau


###### save the R.data
filename <- paste("TSCI-Selection-Setting", f.index,"-Violation",setting.vio,"-Level",tau, "-n", n, ".RData", sep="")
save.image(filename)