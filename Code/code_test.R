# Load libraries
library(tidyverse)
library(dirmult)
library(Rcpp)
library(RcppDist)

sourceCpp("BayesDMC.cpp")
source("DM_Cor_Test.R")

# Simulating data

# Using dirichlet and multinomial separately
DM_simulator <- function(N, p,
                           n_reads_min = 1000, n_reads_max = 2000, seed = 1){
  #browser()
  set.seed(seed)
  
  #browser()
  p_mat <- Y_mat <- alpha_mat <- matrix(0, nrow = N, ncol = p)
  
  for(i in 1:N){
    #browser()
    alpha_mat[i, ] <- runif(p)
    p_mat[i, ] <- rdirichlet(n=1, alpha=alpha_mat[i, ])
    while(any(is.na(p_mat[i, ]))){
      p_mat[i, ] <- rdirichlet(n=1, alpha=alpha_mat[i, ])
    }
    Y_mat[i, ] <- rmultinom(1, round(runif(p, n_reads_min, n_reads_max)), p_mat[i, ])
  }
  return(list(Y = Y_mat, P = p_mat, alpha = alpha_mat))
  
}

N <- 100 # Samples
p <- 5 # Taxons
sim_res <- DM_simulator(N, p)
sim_res$alpha
sim_res$P
sim_res$Y

# C++ implementation
# Assuming X to be a linear function of Taxon 1 
res <- BayesDMC_cpp(sim_res$Y, 2*sim_res$Y[,1] + runif(N), store = T)
system.time(BayesDMC_cpp(sim_res$Y, 2*sim_res$Y[,1] + runif(N), store = T)) # 15.132 seconds elapsed

dim(res$A_store)
dim(res$cor_store)

# observed high correlation b/w alpha[,1] and X as expected
cbind.data.frame(Mean = apply(res$cor_store, 1, mean),
                 LL = apply(res$cor_store, 1, function(i) quantile(i, c(0.025))),
                 UL = apply(res$cor_store, 1, function(i) quantile(i, c(0.975))))

# R implementation
BayesDMC(sim_res$Y, rep(1:2, each = 50), rnorm(100))
system.time(BayesDMC(sim_res$Y, rep(1:2, each = 50), rnorm(100))) # 67.710 seconds elapsed

# Checking correlation between outputs from R and C++ codes
## A_means_1: mean alpha matrix from cpp code
## A_means_2: mean alpha matrix from R code

A_means_1 <- A_means_2 <- matrix(nrow = N, ncol = p)
cor.mat <- NULL
for(i in 1:N){
  for(j in 1:p){
    #browser()
    A_means_1[i, j] <- mean(res$A_store[i, j, ])
    A_means_2[i, j] <- mean(alpha[, i, j])
  }
  #cor1 <- cor(A_means_1[i, ], sim_res$alpha[i, ])
  #cor2 <- cor(A_means_2[i, ], sim_res$alpha[i, ])
  cor3 <- cor(A_means_1[i, ], A_means_2[i, ])
  #cor.mat <- rbind(cor.mat, cbind(cor1, cor2, cor3))
  cor.mat <- rbind(cor.mat, cor3)
}

round(A_means_1)
round(A_means_2)
cor.mat # close to 1 (results match!)
