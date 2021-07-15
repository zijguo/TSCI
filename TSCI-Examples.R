library(MASS)
library(AER)
source('Source-RF.R', encoding = 'UTF-8')
source('Source-Basis.R', encoding = "UTF-8")


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


### dimension change the dimension 10,20
p = 10
### please change this n = 1000, 3000, 5000
n = 3000
### setting, change across 1(for previous 1 and 4), 2(for previous 3 and 6)
f.index = 2
### change the interaction 0, 0.5, 1, 1.5
inter.val = 0
### a denotes the IV strength, set as 1
a = 1
### violation index, set as 1 or 2
vio.index = 1
### tau denotes the violation strength
### set tau as 1
tau = 1


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


### random forest
output.RF <- TSCI.RF(Y,D,Z,X)
# point estimate
output.RF$Coef.robust["RF-Cor-robust"]
# standard error
output.RF$sd.robust["RF-Cor-robust"]

### basis
output.basis <- TSCI.basis(Y,D,Z,X)
# point estimate
output.basis$Coef.robust["Basis-robust"]
# standard error
output.basis$sd.robust["Basis-robust"]
