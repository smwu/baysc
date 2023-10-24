#include <RcppArmadillo.h>
#include <RcppTN.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppTN)]]

// Single draw from Dirichlet
// [[Rcpp::export]]
vec rdirichlet_cpp(vec alpha) {
  int distribution_size = alpha.size();
  // Draw from a Dirichlet
  vec distribution(distribution_size);
  
  double sum_term = 0;
  // Loop through the distribution and draw Gamma variables
  for (int i = 0; i < distribution_size; i++) {
    double cur = R::rgamma(alpha[i],1.0);
    distribution(i) = cur;
    sum_term += cur;
  }
  // Now normalize
  for (int i = 0; i < distribution_size; i++) {
    distribution(i) = distribution(i)/sum_term;
  }
  return distribution;
}

// Single draw from Categorical distribution
// [[Rcpp::export]]
int rcat_cpp(rowvec probs) {
  int num_categs = probs.size();
  IntegerVector draw(num_categs);
  rmultinom(1, probs.begin(), num_categs, draw.begin());
  int categ = which_max(draw);
  return categ; 
}

// Draw from multivariate Normal distribution
// [[Rcpp::export]]
mat mvrnorm_cpp(const int& n, const vec& mu, const mat& sigma) {
  int ncols = sigma.n_cols;
  mat Y = randn(n, ncols);
  return repmat(mu, 1, n).t() + Y * chol(sigma);
}

// Draw from multivariate Normal distribution
// [[Rcpp::export]]
mat mvrnorm_cpp2(const int& n, const vec& mu, const mat& sigma) {
  int ncols = sigma.n_cols;
  mat z = randn(n, ncols) * chol(sigma);
  return mu.t() + z;
}

// Draw from multivariate Normal distribution
// [[Rcpp::export]]
mat mvrnorm_cpp3(const int& n, const rowvec& mu, const mat& sigma) {
  Environment pkg = Environment::namespace_env("LaplacesDemon");
  Function f("rmvn");
  int ncols = sigma.n_cols;
  // NumericMatrix temp = f(_["n"] = 1, _["mu"] = mu, _["Sigma"] = sigma);
  mat temp = as<arma::mat>(f(1, _["mu"] = mu, _["Sigma"] = sigma));
  return temp;
}

// Draw from truncated Normal distribution
// [[Rcpp::export]]
double rtruncnorm_cpp(const int& n, const double& a, const double& b,
                      const double& mean, const double& sd) {
  Environment pkg = Environment::namespace_env("truncnorm");
  Function f("rtruncnorm");
  double rtn = as<double>(f(_["n"] = n, _["a"] = a, _["b"] = b, _["mean"] = mean, 
                            _["sd"] = sd));
  return rtn;
}

// log-sum-exp trick from https://github.com/helske/seqHMM/blob/main/src/logSumExp.cpp
// [[Rcpp::export]]
double logSumExp_cpp(const arma::rowvec& x) {
  int maxi = x.index_max();
  double maxv = x(maxi);
  if (!(maxv > -arma::datum::inf)) {
    return -arma::datum::inf;
  }
  double cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) && (x(i) > -arma::datum::inf)) {
      cumsum += exp(x(i) - maxv);
    }
  }
  return maxv + log1p(cumsum);
}

// Update pi
// [[Rcpp::export]]
void update_pi(vec& pi, const vec& w_all, const vec& c_all, const int& K, 
               const vec& alpha) {
  NumericVector w_all_copy = as<NumericVector>(wrap(w_all));
  // Posterior parameters for pi
  vec alpha_post(K);  
  for (int k = 0; k < K; k++) {
    // Add sum of normalized weights for those assigned to class k, equiv. to
    // weighted number of individuals assigned to each class
    // Careful with 0-based indexing
    LogicalVector indiv_k = (as<IntegerVector>(wrap(c_all)) == (k + 1));
    NumericVector weights_k = w_all_copy[(indiv_k)];
    alpha_post[k] = alpha[k] + sum(weights_k);
  }
  // Rcout <<"alpha_post: " << alpha_post;  // Print out posterior alpha
  // Draw pi from posterior
  vec out = rdirichlet_cpp(alpha_post);
  // Rcout << "rdir " << out;
  pi = out;
  // return pi;
}

// Update c
// [[Rcpp::export]]
void update_c(vec& c_all, const int& n, const int& K, const int& p,
              const cube& theta, const mat& x_mat, const vec& pi,
              const vec& z_all, const mat& V, const mat& xi, const vec& y_all) {
  mat log_cond_c(n, K);        // Individual log-likelihood for each class
  mat pred_class_probs(n, K);  // Posterior class membership probabilities

  // Calculate posterior class membership, p(c_i=k|-), for each class k and
  // update class assignments
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < K; k++) {
      // Calculate theta component of individual log-likelihood for class k
      double log_theta_comp_k = 0.0;
      double log_probit_comp_k = 0.0;
      for (int j = 0; j < p; j++) {
        // Subtract 1 from exposure value due to 0-based indexing
        log_theta_comp_k += log(theta(j, k, x_mat(i, j) - 1));
      }
      // Calculate and control extremes for probit component
      log_probit_comp_k = log(R::dnorm(z_all(i), (V.row(i) * xi.row(k).t()).eval()(0,0), 1.0, false));
      if (log_probit_comp_k == R_NegInf) {
        log_probit_comp_k = log(1e-16);
      }
      log_probit_comp_k += log((y_all(i) * (z_all(i) > 0)) + ((1 - y_all(i)) * (z_all(i) <= 0)));
      // Individual log-likelihood for class k
      log_cond_c(i, k) = log(pi(k)) + log_theta_comp_k + log_probit_comp_k;
      // log_cond_c(i, k) = log(pi(k)) + log_theta_comp_k + 
      //   log(R::dnorm(z_all(i), (V.row(i) * xi.row(k).t()).eval()(0,0), 1.0, false)) + 
      //   log((y_all(i) * (z_all(i) > 0)) + ((1 - y_all(i)) * (z_all(i) <= 0)));
    }
    // Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
    pred_class_probs.row(i) = exp(log_cond_c.row(i) - logSumExp_cpp(log_cond_c.row(i)));
    // Update class assignment using the posterior probabilities
    // Be careful of 0-based indexing
    c_all(i) = rcat_cpp(pred_class_probs.row(i)) + 1;
  }
  // return c_all;
}


// Update c for WOLCA
// [[Rcpp::export]]
void update_c_WOLCA(vec& c_all, const int& n, const int& K, const int& p, 
                  const cube& theta, const mat& x_mat, const vec& pi) {
  mat log_cond_c(n, K);        // Individual log-likelihood for each class
  mat pred_class_probs(n, K);  // Posterior class membership probabilities
  
  // Calculate posterior class membership, p(c_i=k|-), for each class k and
  // update class assignments
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < K; k++) {
      // Calculate theta component of individual log-likelihood for class k
      double log_theta_comp_k = 0.0;
      for (int j = 0; j < p; j++) {
        // Subtract 1 from exposure value due to 0-based indexing
        log_theta_comp_k += log(theta(j, k, x_mat(i, j) - 1));
      }
      // Individual log-likelihood for class k
      log_cond_c(i, k) = log(pi(k)) + log_theta_comp_k;
    }
    // Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
    pred_class_probs.row(i) = exp(log_cond_c.row(i) - logSumExp_cpp(log_cond_c.row(i)));
    // Update class assignment using the posterior probabilities
    // Be careful of 0-based indexing
    c_all(i) = rcat_cpp(pred_class_probs.row(i)) + 1;
  }
  // return c_all;
}

// Update theta
// [[Rcpp::export]]
void update_theta(cube& theta, const int& p, const int& K, const int& d, 
                  const vec& eta, const vec& w_all, const vec& c_all, mat x_mat) {
  NumericVector w_all_copy = as<NumericVector>(wrap(w_all));
  // Posterior parameters for theta
  NumericVector eta_post(d);  
  for (int j = 0; j < p; j++) {
    for (int k = 0; k < K; k++) {
      for (int r = 0; r < d; r++) {
        // Add sum of normalized weights for those assigned to class k with x_ij = r
        // Careful with 0-based indexing
        LogicalVector indiv_r = (as<NumericVector>(wrap(x_mat.col(j))) == (r + 1));
        LogicalVector indiv_k = (as<IntegerVector>(wrap(c_all)) == (k + 1));
        NumericVector weights_r = w_all_copy[(indiv_r & indiv_k)];
        eta_post[r] = eta[r] + sum(weights_r);
      }
      // Rcout <<"eta_post: " << eta_post << "\n";  // Print out posterior alpha
      // Draw theta from posterior
      theta.tube(j, k) = rdirichlet_cpp(eta_post);
    }
  }
  // return theta;
}

// Update xi
// [[Rcpp::export]]
mat update_xi(mat& xi, const int& n, const int& K, const vec& w_all, 
               const vec& c_all, const vec& z_all, const mat& V, 
               const List& mu0, const List& Sig0) {
  // Sparse diagonal normalized weight matrix
  sp_mat W_tilde(n, n);
  W_tilde.diag() = w_all;
  // Sparse design matrix without class assignments
  sp_mat V_sparse = sp_mat(V);
  
  for (int k = 0; k < K; k++) {
    // Sparse diagonal matrix subsetting to obs in class k
    sp_mat C_k(n, n);
    LogicalVector indiv_k = (as<IntegerVector>(wrap(c_all)) == (k + 1));
    C_k.diag() = as<arma::vec>(wrap(indiv_k));
    
    // Draw xi from conditional posterior distribution
    mat Sig0_k = as<arma::mat>(Sig0[k]);
    mat Sig_post = inv(Sig0_k) + mat(V.t() * C_k * W_tilde * V);
    vec mu_right = inv(Sig0_k) * as<arma::vec>(wrap(mu0[k])) + 
      mat(V.t() * C_k * W_tilde * z_all);
    vec mu_post = inv(Sig_post) * mu_right;
    // Rcout <<"mu_post: " << mu_post << "\n";
    
    // Update xi
    xi.row(k) = mvrnorm_cpp(1, mu_post, inv(Sig_post));
    // // Change to R function rmvn
    // xi.row(k) = mvrnorm_cpp3(1, mu_post.t(), inv(Sig_post));
  }
  return xi;
}

// Update z
// [[Rcpp::export]]
vec update_z(vec& z_all, const int& n, const mat& V, const mat& xi, 
              const vec& c_all, const vec& y_all) {
  // Linear predictor using covariate values and class assignment for each individual
  vec lin_pred(n);
  for (int i = 0; i < n; i++) {
    // Be careful of 0-based indexing
    lin_pred(i) = (V.row(i) * xi.row(c_all(i) - 1).t()).eval()(0,0);
    if (y_all(i) == 1) {
      // Probit model latent variable z, following Albert and Chib (1993)
      // For cases, z_i ~ TruncNormal(low=0, high=Inf, mean=V*xi[c_i], var=1)
      z_all(i) = RcppTN::rtn1(lin_pred(i), 1.0, 0.0, R_PosInf);
      // z_all(i) = rtruncnorm_cpp(1, 0.0, R_PosInf, lin_pred(i), 1.0);
    } else {
      // For controls, z_i ~ TruncNormal(low=-Inf, high=0, mean=V*xi[c_i], var=1)
      z_all(i) = RcppTN::rtn1(lin_pred(i), 1.0, R_NegInf, 0.0);
      // z_all(i) = rtruncnorm_cpp(1, R_NegInf, 0.0, lin_pred(i), 1.0);
    }
    if (z_all(i) == R_PosInf) {
      z_all(i) = R::qnorm(1.0 - 1e-10, 0.0, 1.0, true, false);
    } else if (z_all(i) == R_NegInf) {
      z_all(i) = R::qnorm(1e-10, 0.0, 1.0, true, false);
    }
  }
  // Rcout << "lin_pred: " << lin_pred << "\n";
  
  return z_all;
}

// Update individual log-likelihood
// [[Rcpp::export]]
void update_loglik(vec& loglik, const int& n, const int& p, const vec& c_all, 
                   const cube& theta, const mat& x_mat, const vec& pi, 
                   const vec& z_all, const mat& V, const mat& xi, const vec& y_all) {
  for (int i = 0; i < n; i++) {
    int c_i = c_all(i);
    // Calculate theta component of individual log-likelihood
    // Calculate theta component of individual log-likelihood for class k
    double log_theta_comp = 0.0;
    for (int j = 0; j < p; j++) {
      // Subtract 1 due to 0-based indexing
      log_theta_comp += log(theta(j, c_i - 1, x_mat(i, j) - 1));
    }
    // Update individual log-likelihood for class k
    loglik(i) = log(pi(c_i - 1)) + log_theta_comp +
      log(R::dnorm(z_all(i), (V.row(i) * xi.row(c_i - 1).t()).eval()(0,0), 1.0, false)) + 
      log((y_all(i) * (z_all(i) > 0)) + ((1 - y_all(i)) * (z_all(i) <= 0)));
  }
  // return loglik;
}

// Update c test
// [[Rcpp::export]]
void update_c_test(vec& c_all, const int& n, const int& K, const int& p,
                   const cube& theta, const mat& x_mat, const vec& pi,
                   const vec& z_all, const mat& V, const mat& xi, const vec& y_all) {
  mat log_cond_c(n, K);        // Individual log-likelihood for each class
  mat pred_class_probs(n, K);  // Posterior class membership probabilities
  
  // Calculate posterior class membership, p(c_i=k|-), for each class k and
  // update class assignments
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < K; k++) {
      // Calculate theta component of individual log-likelihood for class k
      double log_theta_comp_k = 0.0;
      for (int j = 0; j < p; j++) {
        // Subtract 1 from exposure value due to 0-based indexing
        log_theta_comp_k += log(theta(j, k, x_mat(i, j) - 1));
      }
      // Individual log-likelihood for class k
      log_cond_c(i, k) = log(pi(k)) + log_theta_comp_k +
        log(R::dnorm(z_all(i), (V.row(i) * xi.row(k).t()).eval()(0,0), 1.0, false)) +
        log((y_all(i) * (z_all(i) > 0)) + ((1 - y_all(i)) * (z_all(i) <= 0)));
    }
    // Calculate p(c_i=k|-) = p(x,y,c_i=k) / p(x,y)
    double denom = sum(exp(log_cond_c.row(i)));
    pred_class_probs.row(i) = exp(log_cond_c.row(i)) / denom;
    // pred_class_probs.row(i) = exp(log_cond_c.row(i) - logSumExp_cpp(log_cond_c.row(i)));
    // Update class assignment using the posterior probabilities
    // Be careful of 0-based indexing
    c_all(i) = rcat_cpp(pred_class_probs.row(i)) + 1;
  }
  // return c_all;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# library(R.matlab)
# library(stringr)
# library(fastDummies)
# set.seed(1)
# wd = "/n/holyscratch01/stephenson_lab/Users/stephwu18/wsOFMM/"
# data_dir = "Data/"
# scen_samp = 111211
# iter_pop = 1
# samp_n = 3
# data_path = paste0(wd, data_dir, "simdata_scen", scen_samp, "_iter", iter_pop,
#                     "_samp", samp_n, ".RData")   # Input dataset
# load(data_path)
# data_vars <- sim_data
# # data_vars = readMat(data_path)$sim.data
# # names(data_vars) = str_replace_all(dimnames(data_vars)[[1]], "[.]", "_")
# 
# K = 3
# n = 10
# p = 5
# d = 4
# S = 2
# q = 2
# x_mat = data_vars$X_data[1:n, 1:p]
# y_all = data_vars$Y_data[1:n]
# z_all = rnorm(n)
# z_all = ifelse(y_all == 1, abs(z_all), -abs(z_all))
# s_all = data_vars$true_Si[1:n]
# V = as.matrix(dummy_cols(data.frame(x = factor(s_all, levels = 1:S)),
#                           remove_selected_columns = TRUE))
# kappa = sum(data_vars$sample_wt[1:n]) / n
# w_all = c(data_vars$sample_wt[1:n] / kappa)
# c_all = data_vars$true_Ci[1:n]
# 
# alpha = rep(1, 3) / 3
# eta = rep(1, d)
# mu0 = Sig0 = vector("list", K)
# for (k in 1:K) {
#   mu0[[k]] = rnorm(n = q)
#   Sig0[[k]] = diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
# }
# 
# pi = rdirichlet(1, alpha)
# theta = data_vars$true_global_thetas[1:p, , ]
# xi = matrix(data_vars$true_xi, nrow = K, ncol = q, byrow = FALSE)
# loglik = numeric(n)
*/

/*** R
# print(pi)
# print(theta)
# print(c_all)
# print(z_all)
# print(x_mat)
# print(V)
# print(y_all)
# print(xi)
# update_pi(pi, w_all, c_all, K, alpha)
# pi
# update_c(c_all, n, K, p, theta, x_mat, pi, z_all, V, xi, y_all)
# c_all
# update_theta(theta, p, K, d, eta, w_all, c_all, x_mat)
# theta
# update_xi(xi, n, K, w_all, c_all, z_all, V, y_all, mu0, Sig0)
# xi
# xi_test = array(NA, dim=c(100, 3, 2))
# for (i in 1:100) {
#   xi_test[i,,] = update_xi(n, K, w_all, c_all, z_all, V, xi, y_all, mu0, Sig0)
# }
# par(mfrow=c(3,2))
# plot(density(xi_test[,1,1]), xlab = "1,1")
# plot(density(xi_test[,1,2]), xlab = "1,2")
# plot(density(xi_test[,2,1]), xlab = "2,1")
# plot(density(xi_test[,2,2]), xlab = "2,2")
# plot(density(xi_test[,3,1]), xlab = "3,1")
# plot(density(xi_test[,3,2]), xlab = "3,2")
# update_z(z_all, n, V, xi, c_all, y_all)
# z_all
# par(mfrow=c(1,1))
# plot(density(z_all))
# update_loglik(loglik, n, p, c_all, theta, x_mat, pi, z_all, V, xi, y_all)
# loglik
*/
