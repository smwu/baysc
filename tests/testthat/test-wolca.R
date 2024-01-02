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
  expect_equal(max(table(res_adapt$MCMC_out$c_all_MCMC[4, ])), 3290) 
  expect_equal(min(table(res_adapt$MCMC_out$c_all_MCMC[4, ])), 2)  
})

test_that("fixed sampler works", {
  expect_equal(round(res_fixed$estimates_unadj$pi_med, 2), c(0.51, 0.26, 0.23))
  expect_equal(max(table(res_fixed$estimates_unadj$c_all)), 2156) 
  expect_equal(min(table(res_fixed$estimates_unadj$c_all)), 842)  
})

test_that("variance adjustment works", {
  expect_equal(round(res_adjust$estimates$pi_med, 2), c(0.51, 0.26, 0.22))
  expect_equal(max(table(res_fixed$estimates_unadj$c_all)), 2156) 
  expect_equal(min(table(res_fixed$estimates_unadj$c_all)), 842)  
})

test_that("wolca svyglm works", {
  expect_equal(round(c(res_svyglm$estimates$xi_est), 2), c(-0.02, 0.74, -0.72))
  expect_equal(res_svyglm$data_vars$q, 1)
})

# Run wolca with stratum covariate in probit model
V_data <- data.frame(stratum_id = as.factor(stratum_id))
glm_form <- "~ stratum_id"
res_svyglm_strat <- wolca_svyglm(res = res_adjust, y_all = y_all,
                                 glm_form = glm_form, ci_level = 0.95,
                                 V_data = V_data, save_res = FALSE)

test_that("wolca svyglm with stratum covariate works", {
  expect_equal(round(c(res_svyglm_strat$estimates$xi_est[, 1]), 2), 
               c(0.28, 0.79, -0.49))
  expect_equal(res_svyglm_strat$data_vars$q, 2)
})

