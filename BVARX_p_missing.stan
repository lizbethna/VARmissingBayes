/*
BVARX(p)-missing
Bayesian vector autoregressive approach with exogenous variables
Missing data 
Parallel: Reduce-Sum
*/ 

functions {
// function used to apply parallelization by "reduce_sum" function
real partial_sum_lpdf(vector[] y_slice, 
                      int start, int end, 
                      vector[] eta, matrix L_sigma ){  
    
   return multi_normal_cholesky_lpdf(y_slice | eta[start:end], L_sigma);    
}

// function to compute the linear predictor in the BVARX(p)
vector[] fn_eta(int N, int P, int D, 
                vector MU,             
                matrix Y_imp, matrix[] PHI, 
                matrix GAMMA, matrix Ct){
   vector[D] eta[N-P];
   matrix[N-P,D] eta_aux;
//   eta_aux = rep_matrix(MU,N-P)' + Ct[(P+1):N,] * GAMMA';
   eta_aux = Ct[(P+1):N,] * GAMMA';
   for (p in 1:P){ 
      eta_aux =  eta_aux + Y_imp[p:(N-P+p-1),] * PHI[p]';
   }
   for(t in 1:(N-P)){
      eta[t] = to_vector(eta_aux[t,]);
   }
   return eta;
}
} 

// Data block
data {
   int<lower=0> N;          // No. of time periods observed 
   int<lower=0> D;          // No. of variables in y_t 
   int<lower=1> P;          // Lag order 
   int<lower=1> M;          // No. of covariates 
   int<lower=1> Nobs;       // No. of observed data
   int<lower=1> Nmis;       // No. of missing data 
   matrix[Nobs,D] Y_obs;    // time series vector
   matrix[N,M] Ct;          // Covariates 
   int<lower=0> Nnew;       // No. of time periods for prediction 
   vector[M] Ct_test[Nnew]; // Covariates for prediction 
   int<lower=1, upper=Nobs+Nmis> ii_obs[Nobs]; // index for observed data
   int<lower=1, upper=Nobs+Nmis> ii_mis[Nmis]; // index for missing data
}

// Parameters block
parameters {
   matrix[D, D] PHI[P];   // autoregressive coefficient
   vector[D] MU;          // intercept
   matrix[D, M] GAMMA;    // coefficients for the exogenous inputs
   cholesky_factor_corr[D] L_corr_noise; // Cholesky factor: lower triangular
   vector<lower=0>[D] sd_noise;          // Cholesky factor: scale parameters
   matrix[Nmis,D] Y_mis;  // Imputed data
}

// Transformed parameters block
transformed parameters {
   matrix[D,D] L_sigma;
   matrix[N,D] Y_imp;
   vector[D] Y[N-P];
   L_sigma = diag_pre_multiply(sd_noise, L_corr_noise);
   Y_imp[ii_obs ,] = Y_obs ;
   Y_imp[ii_mis ,] = Y_mis ;
   for (t in 1:(N-P)){
      Y[t] = to_vector(Y_imp[t+P,]) ; 
   }
}

// Model block
model {
   int grainsize=1;
   matrix[N-P,D] eta_aux;
   vector[D] eta[N-P]; 
   
// Prior
   L_corr_noise ~ lkj_corr_cholesky(2.0);
   sd_noise ~ normal(0,1);
   MU ~ normal(0,1);
   for(p in 1:P){
      to_vector(PHI[p]) ~ normal(0, 1);
   }

// Linear predictor for the BVARX(p)
//   eta_aux =  rep_matrix(MU,N-P)' + Ct[(P+1):N,] * GAMMA';
   eta_aux =  Ct[(P+1):N,] * GAMMA';
   for (p in 1:P){ 
      eta_aux =  eta_aux + Y_imp[p:(N-P+p-1),] * PHI[p]';
   }
   for(t in 1:(N-P)){
      eta[t] = to_vector(eta_aux[t,]);
   }

// Likelihood, parallelization by "reduce_sum" function  
   target += reduce_sum(partial_sum_lpdf, Y, 
                        grainsize, 
                        eta, L_sigma);

}

// Generated quantities block
generated quantities {
   matrix[D,D] Sigma;          // Covariance 
   vector[D] Y_fit[N-P];       // Fitted time series
   vector[D] eta_fit[N-P];     // Linear predictor for training dataset
   vector[D] Y_test[Nnew+P];   // Predicted time series
   vector[D] eta_test[Nnew];   // Linear predictor for testing dataset 

   Sigma = L_sigma * L_sigma';

// Fitted values for Training dataset
   eta_fit = fn_eta(N,P,D, MU, Y_imp, PHI, GAMMA,Ct);
   Y_fit = multi_normal_cholesky_rng(eta_fit,L_sigma);

// New Predictions for Testing dataset
   for (t in 1:P){ 
      Y_test[t] = to_vector(Y_imp[N-P+t,]);
   }
   for (t in 1:Nnew) {
//      eta_test[t] = MU + GAMMA * Ct_test[t];
      eta_test[t] = GAMMA * Ct_test[t];
      for (p in 1:P){ 
         eta_test[t] = eta_test[t] + PHI[p] * Y_test[t+p-1];
      }
      Y_test[P+t] = multi_normal_cholesky_rng(eta_test[t],L_sigma);
   }
}
