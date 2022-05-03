# TSCI
This package implements the Two-Stage Curvature Idetification method proposed in the paper [https://arxiv.org/abs/2203.12808] using Ramdom Forest and Bais Approach. It constructs the estimator for treatment effect in the presence of invalid IVs by doing a IV strength test and the violation space selection. Confidence intervals are futher constructed for this estimator.


## Installation
The package can be installed from Github using the following code.
```
# install.packages("devtools")
library(devtools)
devtools::install_github("https://github.com/zijguo/TSCI")
```

## Examples
This example shows the point estimators and confidence intervals for TSCI with random forest and basis approach, The true outcome model has a linear violation of the instrument variable.

``` r
library(MASS)

# dimension
p = 20
# sample size
n = 1000
# interaction value
inter.val = 0.5
# the IV strength
a = 1
# violation strength
tau = 1
f = function(x){a*(1*sin(2*pi*x) + 1.5*cos(2*pi*x))}
rho1=0.5
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
Cov=(A1gen(rho1,p+1))
mu=rep(0,p+1)
# true effect
beta=1
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
# generate the outcome variable Y
Y=D*beta+tau*Z+X%*%gamma+Error[,2]

# Two Stage Random Forest
output.RF = TSRF(Y,D,Z,X)
# point estimates
output.RF$Coef.robust
# standard errors
output.RF$sd.robust
# confidence intervals
output.RF$CI.robust


# Two Stage Basis Approach
output.BA = TSBA(Y,D,Z,X)
# point estimates
output.BA$Coef.robust
# standard errors
output.BA$sd.robust
# confidence intervals
output.BA$CI.robust
```
