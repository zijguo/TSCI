# TSCI
This package implements the Two-Stage Curvature Idetification method proposed in the paper [https://arxiv.org/abs/2203.12808] using Ramdom Forest and Bais Approach. It constructs the estimator for treatment effect in the presence of invalid IVs by doing a IV strength test and the violation space selection. Confidence intervals are futher constructed for this estimator.


## Installation
The package can be installed from Github using the following code:(cannot be installed for now since this is a private repo.)
```
# install.packages("devtools")
library(devtools)
devtools::install_github("https://github.com/zijguo/TSCI")
```

## Examples
