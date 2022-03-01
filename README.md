# Bayes-Correlation-Test Tutorial


## Load the following libraries
```
library(Rcpp)
library(RcppDist)
```

## Load the cpp function BayesDMC

The function is located in the code folder in this repository.

```
sourceCpp("BayesDMC.cpp")
```


## C++ implementation of the Bayes Correlation Test

Supply the counts and target from the same group. The function can only run one group at a time so be sure to subset your data by group.  This will be updated later to handle multiple groups.  

Counts: a matrix of count data.  The rows must be samples; columns are taxa or other feature.

Target: a vector of real-valued data that will be tested for correlation with each taxa.

```
res <- BayesDMC_cpp(counts, target, store = T)
```

## Output of BayesDMC

The function will return two 3-dimensional arrays of matrices: the estimated values of the DM parameter, alpha, for each taxa or feature at each thinned after burn-in MCMC iteration; the estimated correlation between each taxa and the target at each thinned after burn-in MCMC iteration.

The dimensions of each matrix are (sample i, taxa j, thinned after burn-in MCMC sample t) where each matrix pertains to one sample. You can check the dimensions of each using:
```
dim(res$A_store)
dim(res$cor_store)
```

# Summary of the Correlations
The posterior mean and 95% credible interval (lower bound = LL, upper bound = UL) between target and taxa are displayed by running the command below.  Each row is a taxon or feature.   The columns are the posterior mean, LL, and UL for the correlation between the taxaon and the target.

```
cbind.data.frame(Mean = apply(res$cor_store, 1, mean),
                 LL = apply(res$cor_store, 1, function(i) quantile(i, c(0.025))),
                 UL = apply(res$cor_store, 1, function(i) quantile(i, c(0.975))))
```

