data {
  int<lower=1> N;             // number of observations
  int<lower=1> P;             // number of variables
  matrix[N,P] Y;              // data matrix of order [N,P]
  int<lower=1> D;             // number of factors
}
transformed data {
  int<lower=1> M;
  vector[P] mu;
  M  = D*P;  // number of non-zero loadings
  mu = rep_vector(0.0,P);
}
parameters {    
  matrix[P, D] L; 	      // factor loadings matrix
  vector<lower=0>[P] psi;     // vector of variances
}
transformed parameters{
  cov_matrix[P] Q;   //Covariance mat
  Q=L*L'+diag_matrix(psi); 
}
model {
// the priors 
psi ~ inv_gamma(0.0005, 0.0005);
//The likelihood
for( j in 1:N)
    Y[j] ~ multi_normal(mu,Q); 
}
