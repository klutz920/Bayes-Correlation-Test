// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double cor_cpp(Col<double> Y, colvec rank_X){
  
  // Function to compute spearman correlation in cpp
  
  int N = Y.n_rows;
  // Removing ties in Y at random
  for(int i = 0; i < N; i++){
    if(Y(i) == 0.0){
      Y(i) = Y(i) + runif(1, 0, 0.0001)(0);
    }
  }
  
  uvec rank_Y_temp = sort_index(sort_index(Y)) + 1;
  colvec rank_Y = conv_to<colvec>::from(rank_Y_temp);
  
  colvec D(N);
  D = rank_Y - rank_X;
  
  double cor_res = 0.0, D_sum = 0.0;
  
  for(int i = 0; i < N; i++){
    D_sum = D_sum + pow(D(i), 2.0);
  }
  
  cor_res = 1.0 - (6*D_sum / (N*(pow(N, 2) - 1)));
  
  return cor_res;
}

