data {
  int<lower=1> K;  // number of clusters, known at the time of post-processing
  int<lower=1> p;  // number of food items
  int<lower=1> d;  // number of consumption levels
  int<lower=1> n;  // number of subjects
  int<lower=1> q;  // number of covariates in probit regression
  
  int X[n, p];                // categorical food data
  int<lower=0, upper=1> y[n]; // binary outcome data
  real V[n, q];               // covariate matrix excluding class assignment
  vector<lower=0>[n] weights;       // individual-level survey weights
  
  vector[K] alpha;         // hyperparameter for pi prior
  vector[d] eta;           // hyperparameter for theta prior
  vector[q] mu0[K];           // hyperparameter for mean of xi prior
  cov_matrix[q] Sig0[K];      // hyperparameter for covariance of xi prior
}
parameters {
  simplex[K] pi;                 // cluster probabilities
  simplex[d] theta[p, K];  // cluster-specific item consumption probabilities
  vector[q] xi[K];         // regression coefficients
}
transformed parameters {
  matrix[p, d] theta_prod;  // to check convergence up to permutation of labels
  vector[q] xi_prod;  // to check convergence up to permutation of labels
  for (j in 1:p) {
    for (r in 1:d) {
      theta_prod[j, r] = prod(theta[j, ,r]);
    }
  }
  for (v in 1:q) {
    xi_prod[v] = prod(xi[, v]);
  }
}
model {
  vector[K] log_cond_c[n]; // log p(c_i=k| -)
  
  pi ~ dirichlet(alpha);  // prior for pi
  for (j in 1:p) {        // prior for theta
    for (k in 1:K) {
      theta[j, k] ~ dirichlet(eta);
    }
  }
  for (k in 1:K) {
    xi[k] ~ multi_normal(mu0[k], Sig0[k]);
  }
  
  for (i in 1:n) {
    for (k in 1:K) {
      log_cond_c[i, k] = log(pi[k]) + bernoulli_lpmf(y[i] | Phi(to_row_vector(V[i, ]) * xi[k]));
      for (j in 1:p) {
        log_cond_c[i, k] = log_cond_c[i, k] + categorical_lpmf(X[i,j] | theta[j, k, ]);
      }
    }
    // increment log-probability
    target += weights[i] * log_sum_exp(log_cond_c[i]);
  }
}
generated quantities {
  simplex[K] pred_class_probs[n];  // posterior class probs for each individual
  int<lower=1> pred_class[n];  // posterior predicted class for individuals
  vector[K] log_cond_c[n]; // log p(c_i=k| -)
  
  for (i in 1:n) {
    for (k in 1:K) {
      // need to recalculate log-density to avoid scope issues for 'log_cond_c'
      log_cond_c[i, k] = log(pi[k]) + bernoulli_lpmf(y[i] | Phi(to_row_vector(V[i, ]) * xi[k]));
      for (j in 1:p) {
        log_cond_c[i, k] = log_cond_c[i, k] + categorical_lpmf(X[i,j] | theta[j, k, ]);
      }
      // print("log_cond_c[", i,",", k,"]:", log_cond_c[i,k]);
    }
    pred_class_probs[i, ] = exp(log_cond_c[i, ] - log_sum_exp(log_cond_c[i]));
    // print("log_sum_exp[", i,"]:", log_sum_exp(log_cond_c[i]));
    // print(i, " indiv:", pred_class_probs[i, ]);
    pred_class[i] = categorical_rng(pred_class_probs[i, ]);
  }
}


