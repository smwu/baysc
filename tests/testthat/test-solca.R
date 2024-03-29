# # Load data and obtain relevant variables
# data("sim_data")
# data_vars <- sim_data
# x_mat <- data_vars$X_data            # Categorical exposure matrix, nxp
# y_all <- c(data_vars$Y_data)         # Binary outcome vector, nx1
# cluster_id <- data_vars$cluster_id   # Cluster indicators, nx1
# stratum_id <- data_vars$true_Si      # Stratum indicators, nx1
# sampling_wt <- data_vars$sample_wt   # Survey sampling weights, nx1
# n <- dim(x_mat)[1]                   # Number of individuals
# 
# # Probit model only includes latent class
# V <- as.data.frame(matrix(1, nrow = n)) # Additional regression covariates
# glm_form <- "~ 1"
# 
# # Run swolca
# res_adapt <- solca(x_mat = x_mat, y_all = y_all, V = V,
#                     run_sampler = "adapt", glm_form = glm_form, adapt_seed = 1,
#                     n_runs = 5, burn = 1, thin = 1, save_res = FALSE)
# res_fixed <- solca(x_mat = x_mat, y_all = y_all, V = V,
#                     run_sampler = "fixed", glm_form = glm_form, fixed_seed = 1, 
#                     K_fixed = 3, n_runs = 5, burn = 1, thin = 1, save_res = FALSE)
# 
# test_that("adaptive sampler works", {
#   expect_equal(res_adapt$K_fixed, 7)
#   expect_equal(max(table(res_adapt$MCMC_out$c_all_MCMC[4, ])), 1002) 
#   expect_equal(min(table(res_adapt$MCMC_out$c_all_MCMC[4, ])), 242)  
# })
# 
# test_that("fixed sampler works", {
#   expect_equal(round(res_fixed$estimates$pi_med, 2), c(0.69, 0.31))
#   expect_equal(max(table(res_fixed$estimates$c_all)), 2600) 
#   expect_equal(min(table(res_fixed$estimates$c_all)), 1400)  
# })
# 
# 
# # Run swolca with stratum covariate in probit model
# V <- data.frame(stratum_id = as.factor(stratum_id))
# glm_form <- "~ stratum_id"
# res_fixed_strat <- solca(x_mat = x_mat, y_all = y_all, V = V,
#                           run_sampler = "fixed", glm_form = glm_form, 
#                           fixed_seed = 1, K_fixed = 3, n_runs = 5, burn = 1, 
#                           thin = 1, save_res = FALSE)
# 
# test_that("stratum covariate works", {
#   expect_equal(round(res_fixed_strat$estimates$pi_med, 2), c(0.54, 0.21, 0.25))
#   expect_equal(max(table(res_fixed_strat$estimates$c_all)), 2156) 
#   expect_equal(min(table(res_fixed_strat$estimates$c_all)), 842) 
# })
# 
