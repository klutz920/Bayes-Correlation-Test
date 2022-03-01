#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;
using namespace arma;

//static Mat<double> logDM(Row<double> Y_i, Row<double> A_i, int tax_id);
static Mat<double> logDM(Row<double> Y_i, Row<double> A_i);
static double cor_cpp(Col<double> Y, colvec rank_X);


// DM model with linear regression model

// [[Rcpp::export]]
Rcpp::List BayesDMC_cpp(Mat<double> Y, colvec X, bool store) {
  
  // Input: Count Matrix Y, Vector X (with removed ties) for correlation analysis 
  
  // Output: A_store - Normalized abundance cube with N (samples) rows, p (taxons) columns and [(iter - burn)/thin] slices
  // cor_store - Spearman correlation matrix with p (taxons) rows and [(iter - burn)/thin] columns
  
  // Read data info
  int N = Y.n_rows;
  int p = Y.n_cols;
  
  // Initializing model parameters & hyperparameters
  
  Mat<double> A(N, p, fill::none);
  A.fill(1);
  //A = Y;
  
  // MCMC settings
  int iter = 10000;
  int burn = 0.5*iter;
  int thin = 10;
  
  double tau = 0.1;
  
  // Store & Temp variables
  double hastings;
  int count = 0;
  Col<int> count_2(p, fill::zeros);
  //Mat<double> lnA_temp = lnA;
  Mat<double> A_temp = A;
  Mat<double> A_map = A;
  double A_sum = 0, A_sum_temp = 0;
  //Cube<double> A_store(N, p, iter - burn, fill::zeros);
  Cube<double> A_store(N, p, (iter - burn)/thin, fill::zeros);
  Mat<double> cor_store(p, (iter - burn)/thin, fill::zeros);
  
  // Performing MCMC
  
  for(int i = 0; i < N; i++){
    
    // Update alpha
    // Fix conflict with A_sum between iterations
    // Try calculating entire log posterior and then subtracting it with old one (slow method)
    
    if(i*100/N == count){
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    count_2.fill(0);
    for(int it = 0; it < iter; it++){
      
      for(int j = 0; j < p; j++){
        A_temp(i, j) = exp(rnorm(1, log(A(i, j)), tau)(0));
        
        hastings = 0.0;
        
        //hastings = hastings + logDM(Y.row(i), A_temp.row(i), j)(0, 0);
        //hastings = hastings - logDM(Y.row(i), A.row(i), j)(0, 0);
        hastings = hastings + logDM(Y.row(i), A_temp.row(i))(0, 0);
        hastings = hastings - logDM(Y.row(i), A.row(i))(0, 0);
        
        if(hastings >= log(double(rand()%10001)/10000)){
          A(i, j) = A_temp(i, j);
        }
        
        
        if(store){
          if(it > burn){
            if((it - burn)%thin == 1){
              A_store(i, j, count_2(j)) = A(i, j);
              count_2(j)++;
            }
          }
        }

        
      }
      
    }

    
  }
  
  // Correlation analysis
  
  uvec rank_X_temp = sort_index(sort_index(X)) + 1;
  colvec rank_X = conv_to<colvec>::from(rank_X_temp);
  
  for(int j = 0; j < p; j++){
    for(int k = 0; k < A_store.n_slices; k++){
      cor_store(j, k) = cor_cpp(A_store.slice(k).col(j), rank_X);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("A_store") = A_store,
                            Rcpp::Named("cor_store") = cor_store);
  
}


// [[Rcpp::export]]
Mat<double> logDM(Row<double> Y_i, Row<double> A_i) { 
  Row<double> res(1, 1, fill::zeros);
  res = res + lgamma(sum(A_i)) - lgamma(sum(Y_i + A_i));
  res = res + sum(lgamma(Y_i + A_i) - lgamma(A_i));
  return res;
}

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

/*
// [[Rcpp::export]]
Mat<double> logDM(Row<double> Y_i, Row<double> A_i, int tax_id) { 
  Row<double> res(1, 1, fill::zeros);
  res = res + lgamma(sum(A_i)) - lgamma(sum(Y_i) + sum(A_i));
  res = res + lgamma(Y_i(tax_id) + A_i(tax_id)) - lgamma(A_i(tax_id));
  return res;
}*/







