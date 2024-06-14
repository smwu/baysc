data {
  int<lower=1> K;  // number of clusters, known at the time of post-processing
  int<lower=1> J;  // number of food items
  int<lower=1> R;  // number of consumption levels
  int<lower=1> n;  // number of subjects

  int X[n, J];                // categorical food data
  vector<lower=0>[n] weights;       // individual-level survey weights
  
  vector[K] alpha;         // hyperparameter for pi prior
  matrix[J, R] eta;           // hyperparameter for theta prior
}
parameters {
  simplex[K] pi;                 // cluster probabilities
  simplex[R] theta[J, K];  // cluster-specific item consumption probabilities
}
transformed parameters {
  matrix[J, R] theta_prod;  // to check convergence up to permutation of labels
  for (j in 1:J) {
    for (r in 1:R) {
      theta_prod[j, r] = prod(theta[j, ,r]);
    }
  }
}
model {
  vector[K] log_cond_c[n]; // log p(c_i=k| -)
  
  pi ~ dirichlet(alpha);  // prior for pi
  for (j in 1:J) {        // prior for theta
    for (k in 1:K) {
      theta[j, k] ~ dirichlet(to_vector(eta[j, ]));
    }
  }
  
  for (i in 1:n) {
    for (k in 1:K) {
      log_cond_c[i, k] = log(pi[k]);
      for (j in 1:J) {
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
      log_cond_c[i, k] = log(pi[k]);
      for (j in 1:J) {
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


