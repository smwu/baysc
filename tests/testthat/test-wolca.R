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
V_data <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
glm_form <- "~ 1"

# Run wolca
res_adapt <- wolca(x_mat = x_mat, sampling_wt = sampling_wt,
                    cluster_id = cluster_id, stratum_id = stratum_id, 
                    run_sampler = "adapt", adapt_seed = 1,
                    n_runs = 5, burn = 1, thin = 1, save_res = FALSE)
res_fixed <- wolca(x_mat = x_mat, sampling_wt = sampling_wt,
                    cluster_id = cluster_id, stratum_id = stratum_id,
                    run_sampler = "fixed", fixed_seed = 1, K_fixed = 3, 
                   n_runs = 50, burn = 25, thin = 1, save_res = FALSE)

# Apply variance adjustment to posterior estimates
res_adjust <- wolca_var_adjust(res = res_fixed, num_reps = 100, save_res = FALSE,
                               adjust_seed = 1)

# Run weighted outcome regression model
res_svyglm <- wolca_svyglm(res = res_adjust, y_all = y_all,
                           glm_form = glm_form, ci_level = 0.95,
                           V_data = V_data, save_res = FALSE)

test_that("adaptive sampler works", {
  expect_equal(res_adapt$K_fixed, 2)
  expect_equal(max(table(res_adapt$MCMC_out$c_all_MCMC[4, ])), 658) 
  expect_equal(min(table(res_adapt$MCMC_out$c_all_MCMC[4, ])), 1)  
})

test_that("fixed sampler works", {
  expect_equal(round(res_fixed$estimates$pi_med, 2), c(0.27, 0.52, 0.21))
  expect_equal(max(table(res_fixed$estimates$c_all)), 442) 
  expect_equal(min(table(res_fixed$estimates$c_all)), 173)  
})

test_that("variance adjustment works", {
  expect_equal(round(res_adjust$estimates_adjust$pi_med, 2), c(0.27, 0.52, 0.21))
  expect_equal(max(table(res_adjust$estimates_adjust$c_all)), 442) 
  expect_equal(min(table(res_adjust$estimates_adjust$c_all)), 173)  
})

test_that("wolca svyglm works", {
  expect_equal(round(c(res_svyglm$estimates_svyglm$xi_est), 2), c(0.43, -0.46, -1.39))
  expect_equal(res_svyglm$data_vars$Q, 1)
})

# Run wolca with stratum covariate in probit model
V_data <- data.frame(stratum_id = as.factor(stratum_id))
glm_form <- "~ stratum_id"
res_svyglm_strat <- wolca_svyglm(res = res_adjust, y_all = y_all,
                                 glm_form = glm_form, ci_level = 0.95,
                                 V_data = V_data, save_res = FALSE)

test_that("wolca svyglm with stratum covariate works", {
  expect_equal(round(c(res_svyglm_strat$estimates_svyglm$xi_est[, 1]), 2), 
               c(0.41, -0.36, -1.38))
  expect_equal(res_svyglm_strat$data_vars$Q, 2)
})

# Run wolca with different categories for different exposure variables
# Convert first food items to be binary: 1 or 2
x_mat[, 1] <- ifelse(x_mat[, 1] >= 3, 2, 1)
# Convert second food item to have 3 levels: 1, 2, or 3
x_mat[, 2] <- ifelse(x_mat[, 2] >= 3, 3, x_mat[, 2])
res_R_j <- wolca(x_mat = x_mat, sampling_wt = sampling_wt,
                  cluster_id = cluster_id, stratum_id = stratum_id, 
                  run_sampler = "fixed", 
                  fixed_seed = 888, K_fixed = 3, n_runs = 100, burn = 50, 
                  thin = 5, update = 20, save_res = FALSE)
res_R_j_adjust <- wolca_var_adjust(res = res_R_j, num_reps = 100,
                                    save_res = FALSE, adjust_seed = 1)

test_that("R_j works", {
  expect_equal(round(res_R_j$estimates$pi_med, 2), c(0.27, 0.53, 0.20))
  expect_equal(max(table(res_R_j$estimates$c_all)), 442) 
  expect_equal(min(table(res_R_j$estimates$c_all)), 173) 
})

