# Load data and obtain relevant variables
data("sim_data")
data_vars <- sim_data
x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
cluster_id <- data_vars$cluster_id   # Cluster indicators, nx1
stratum_id <- data_vars$true_Si      # Stratum indicators, nx1
sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
n <- dim(x_mat)[1]                   # Number of individuals

# Probit model only includes latent class
V <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
glm_form <- "~ 1"

# Run swolca
res_adapt <- swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt,
                   cluster_id = cluster_id, stratum_id = stratum_id, V = V,
                   run_sampler = "adapt", glm_form = glm_form, adapt_seed = 1,
                   n_runs = 5, burn = 1, thin = 1, save_res = FALSE)
res_fixed <- swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt,
                    cluster_id = cluster_id, stratum_id = stratum_id, V = V,
                    run_sampler = "fixed", glm_form = glm_form, fixed_seed = 1, 
                    K_fixed = 3, n_runs = 5, burn = 1, thin = 1, save_res = FALSE)

test_that("adaptive sampler works", {
  expect_equal(res_adapt$K_fixed, 7)
  expect_equal(max(table(res_adapt$MCMC_out$c_all_MCMC[4, ])), 1002) 
  expect_equal(min(table(res_adapt$MCMC_out$c_all_MCMC[4, ])), 258)  
})

test_that("fixed sampler works", {
  expect_equal(round(res_fixed$estimates$pi_med, 2), c(0.48, 0.31, 0.21))
  expect_equal(max(table(res_fixed$estimates$c_all)), 2156) 
  expect_equal(min(table(res_fixed$estimates$c_all)), 842)  
})


# Run swolca with stratum covariate in probit model
V <- data.frame(stratum_id = as.factor(stratum_id))
glm_form <- "~ stratum_id"
res_fixed_strat <- swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt,
                          cluster_id = cluster_id, stratum_id = stratum_id, V = V,
                          run_sampler = "fixed", glm_form = glm_form, 
                          fixed_seed = 1, K_fixed = 3, n_runs = 5, burn = 1, 
                          thin = 1, save_res = FALSE)

test_that("stratum covariate works", {
  expect_equal(round(res_fixed_strat$estimates$pi_med, 2), c(0.51, 0.26, 0.23))
  expect_equal(max(table(res_fixed_strat$estimates$c_all)), 2156) 
  expect_equal(min(table(res_fixed_strat$estimates$c_all)), 842) 
})


# Run swolca with different categories for different exposure variables
# Convert first food items to be binary: 1 or 2
x_mat[, 1] <- ifelse(x_mat[, 1] >= 3, 2, 1)
# Convert second food item to have 3 levels: 1, 2, or 3
x_mat[, 2] <- ifelse(x_mat[, 2] >= 3, 3, x_mat[, 2])

res_R_j <- swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt,
                  cluster_id = cluster_id, stratum_id = stratum_id, V = V,
                  run_sampler = "fixed", glm_form = glm_form, 
                  fixed_seed = 1, K_fixed = 3, n_runs = 5, burn = 1, 
                  thin = 1, save_res = FALSE)
res_R_j <- swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt,
                  cluster_id = cluster_id, stratum_id = stratum_id, V = V,
                  run_sampler = "fixed", glm_form = glm_form, 
                  fixed_seed = 1, K_fixed = 3, n_runs = 50, burn = 25, 
                  thin = 1, save_res = FALSE)
# Apply variance adjustment to posterior estimates
res_R_j_adjust <- swolca_var_adjust(res = res_R_j, num_reps = 100, 
                                    save_res = FALSE, adjust_seed = 1)

test_that("R_j works", {
  expect_equal(round(res_R_j$estimates$pi_med, 2), c(0.52, 0.25, 0.23))
  expect_equal(max(table(res_R_j$estimates$c_all)), 2156) 
  expect_equal(min(table(res_R_j$estimates$c_all)), 842) 
})
