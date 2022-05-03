# TSCI
This package implements the Two Stage Curvature Idetification (TSCI) method proposed in the paper https://arxiv.org/abs/2203.12808 using Ramdom Forests and Basis Approach. It constructs a robust point estimator of the treatment effect even the instrumental variabls are possibly invalid. Confidence intervals are futher constructed.


## Installation
The package can be installed from Github using the following code.
```
# install.packages("devtools")
library(devtools)
devtools::install_github("https://github.com/zijguo/TSCI")
```

## Data generation 
#We generate the data where the instrumental variable has a linear violation form.
``` r
library(TSCI)
library(MASS)

# dimension
p = 20
# sample size
n = 1000
# interaction value of Z and X
inter.val = 0.5
# the IV strength
a = 1
# the conditional mean function for the treatment model
f = function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
# function to generate covariance matrix
A1gen=function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]=rho^(abs(i-j))
    }
  }
  A1
}
rho1=0.5
Cov=(A1gen(rho1,p+1))
mu=rep(0,p+1)
# true effect
beta=1
# violation strength
tau = 1
# other model parameters
alpha=as.matrix(rep(-0.3,p))
gamma=as.matrix(rep(0.2,p))
inter=as.matrix(c(rep(inter.val,5),rep(0,p-5)))


# generate the data
mu.error=rep(0,2)
Cov.error=matrix(c(1,0.5,0.5,1),2,2)
Error=mvrnorm(n, mu.error, Cov.error)
W.original=mvrnorm(n, mu, Cov)
W=pnorm(W.original)
# instrument variable
Z=W[,1]
# baseline covariates
X=W[,-1]
# generate the treatment variable D
D=f(Z)+X%*%alpha+Z*X%*%inter+Error[,1]
# generate the outcome variable Y, linear violation of Z
Y=D*beta+tau*Z+X%*%gamma+Error[,2]
```



# TSCI with Random Forests
``` r
# implement TSCI with Random Forests 
output.RF = TSRF(Y,D,Z,X)
# point estimates
output.RF$Coef.robust
# standard errors
output.RF$sd.robust
# confidence intervals
output.RF$CI.robust
# estimated violation space
output.RF$q.hat
```

# TSCI with Basis Approaches
``` r
# implement TSCI with Basis Approach
output.BA = TSBA(Y,D,Z,X)
# point estimates
output.BA$Coef.robust
# standard errors
output.BA$sd.robust
# confidence intervals
output.BA$CI.robust
# estimated violation space
output.BA$q.hat
```


# Analysis of Card Education Data by TSCI 

``` r
library(TSCI)
# card data
library(ivmodel)
data(card.data)

Xname2=c("exper", "expersq", "black", "south", "smsa", "reg661","reg662", "reg663", "reg664", "reg665", "reg666", "reg667","reg668", "smsa66")

Y=card.data[,"lwage"]
D=card.data[,"educ"]
Z=card.data[,"nearc4"]
X2=card.data[,Xname2]

# results for TSLS
TSLS = ivmodel(Y=Y,D=D,Z=Z,X=X2)
summary(TSLS)

# results for TSCI with Random Forest
card.RF = TSRF(Y,D,Z,X2)
# point estimates
card.RF$Coef.robust
# standard errors
card.RF$sd.robust
# confidence intervals
card.RF$CI.robust
```
