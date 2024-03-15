data {
  int<lower=1> N; // Number of observations
  int<lower=1> D; // Number of variables (dimension of the data)
  int<lower=1> K; // Truncation level (max number of clusters)
  matrix[N, D] y; // Data matrix
}

parameters {
  simplex[K] weights; // Mixture weights
  vector[D] mu[K]; // Array of means for each cluster
  cholesky_factor_corr[D] L_Omega[K]; // Cholesky factor of correlation matrix for each cluster
  vector<lower=0>[D] sigma[K]; // Array of standard deviations for each cluster
}

model {
  // Priors
  for (k in 1:K) {
    mu[k] ~ normal(0, 10);
    L_Omega[k] ~ lkj_corr_cholesky(2.0);
    sigma[k] ~ cauchy(0, 2.5);
  }
  
  // Likelihood
  for (n in 1:N) {
    vector[K] logp;
    for (k in 1:K) {
      matrix[D, D] cov_k = diag_pre_multiply(sigma[k], L_Omega[k]); // Construct covariance matrix
      logp[k] = log(weights[k]) + multi_normal_cholesky_lpdf(y[n] | mu[k], L_Omega[k]);
    }
    target += log_sum_exp(logp);
  }
}

generated quantities {
  int<lower=1, upper=K> cluster_assignments[N];
  for (n in 1:N) {
    vector[K] log_prob; // Log probabilities for each cluster
    for (k in 1:K) {
      matrix[D, D] cov_k = diag_pre_multiply(sigma[k], L_Omega[k]); // Consistent with model block
      log_prob[k] = multi_normal_cholesky_lpdf(y[n] | mu[k], L_Omega[k]); // Corrected for multivariate
    }
    cluster_assignments[n] = categorical_logit_rng(log_prob);
  }
}
