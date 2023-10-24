#===================================================
## Weighted Supervised Overfitted Latent Class Model
## Programmer: SM Wu   
## Data: NHANES Application with Covariates  
## Date Updated: 2023/07/15
#===================================================

# Load libraries
library(plyr)
library(dplyr)
library(tidyr)
library(forcats)
library(LaplacesDemon)
library(truncnorm)
library(stringr)
library(R.matlab)
library(fastDummies)
library(matrixStats)
library(Matrix)
library(gtools)
library(e1071)
library(rstan)
library(survey)
library(Rcpp)
library(RcppArmadillo)
library(RcppTN)
library(ggpubr)
library(tableone)
library(ggdendro)

#========================= MAIN FUNCTION =======================================

# 'WSOLCA_app_covs_Rcpp' runs the WSOLCA model with covariates and saves and 
# returns results
# Inputs:
#   data_vars: Output list from "process_data" function 
#   adapt_path: String path for adaptive sampler output file
#   adj_path: String path for adjusted output file
#   stan_path: String path for Stan file
#   save_res: Boolean specifying if results should be saved. Default = TRUE
#   n_runs: Number of MCMC iterations
#   burn: Burn-in period
#   thin: Thinning factor
#   K_known: Number of latent classes. If NULL (default), adaptive sampler runs
# Outputs: Saves and returns list `res` containing:
#   analysis_adj: List of posterior model results
#   runtime: Total runtime for model
#   data_vars: Input dataset
#   MCMC_out: List of full MCMC output
#   post_MCMC_out: List of MCMC output after relabeling
#   K_MCMC: Adaptive sampler MCMC output for K
# Also saves list `adapt_MCMC` containing:
#   MCMC_out: List of full MCMC output
#   K_fixed: Number of classes to use for fixed sampler; output from adaptive sampler
#   K_MCMC: Adaptive sampler MCMC output for K
WSOLCA_app_covs_Rcpp <- function(data_vars, adapt_path, adj_path, stan_path, 
                                 save_res = TRUE, n_runs, burn, thin, K_known = NULL) {
  start_time <- Sys.time()
  
  #================= Read in data ==============================================
  print("Read in data")
  x_mat <- data_vars$x_mat
  y_all <- data_vars$y_all
  s_all <- data_vars$s_all
  clus_id_all <- data_vars$clus_id_all
  sample_wt <- data_vars$sample_wt
  V <- data_vars$V
  
  # Obtain dimensions
  n <- dim(x_mat)[1]        # Number of individuals
  p <- dim(x_mat)[2]        # Number of exposure items
  d <- max(apply(x_mat, 2,  # Number of exposure categories 
                 function(x) length(unique(x))))  # CHANGE TO ADAPT TO ITEM
  q <- ncol(V)                         # Number of regression covariates excluding class assignment
  
  # Obtain normalized weights
  kappa <- sum(sample_wt) / n   # Weights norm. constant. If sum(weights)=N, this is 1/(sampl_frac)
  w_all <- c(sample_wt / kappa) # Weights normalized to sum to n, nx1
  
  #================= ADAPTIVE SAMPLER ==========================================
  if (!is.null(K_known)) {
    K_fixed <- K_known
  } else {
    print("Adaptive sampler")
    #================= Initialize priors and variables for OLCA model ============
    K_max <- 30                      # Upper limit for number of classes
    alpha <- rep(1, K_max) / K_max   # Hyperparameter for prior for pi
    eta <- rep(1, d)                 # Hyperparameter for prior for theta
    # Obtain pi, theta, c_all
    OLCA_params <- init_OLCA(alpha = alpha, eta = eta, n = n, K = K_max, p = p, 
                             d = d)
    
    #================= Initialize priors and variables for probit model ==========
    # Initialize hyperparameters for xi
    mu0 <- Sig0 <- vector("list", K_max)
    for (k in 1:K_max) {
      # MVN(0,1) hyperprior for prior mean of xi
      mu0[[k]] <- rnorm(n = q)
      # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated 
      # components and mean variance 2.5 for a weakly informative prior on xi
      Sig0[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
    }
    # Obtain xi, z_all                        
    probit_params <- init_probit(mu0 = mu0, Sig0 = Sig0, K = K_max, q = q, n = n, 
                                 V = V, y_all = y_all, c_all = OLCA_params$c_all)
    
    #================= Run adaptive sampler to obtain number of classes ==========
    # Obtain pi_MCMC, theta_MCMC, xi_MCMC, c_all_MCMC, z_all_MCMC, loglik_MCMC
    MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, probit_params = probit_params, 
                              n_runs = n_runs, burn = burn, thin = thin, K = K_max, 
                              p = p, d = d, n = n, q = q, w_all = w_all, x_mat = x_mat, 
                              y_all = y_all, V = V, alpha = alpha, eta = eta, 
                              Sig0 = Sig0, mu0 = mu0)
    
    #================= Post-processing for adaptive sampler ======================
    # Get median number of classes with >= 5% of individuals, over all iterations
    M <- dim(MCMC_out$pi_MCMC)[1]  # Number of stored MCMC iterations
    K_MCMC <- rowSums(MCMC_out$pi_MCMC >= 0.05)
    K_med <- round(median(K_MCMC))
    # Get number of unique classes for fixed sampler
    K_fixed <- K_med
    print(paste0("K_fixed: ", K_fixed))
    # Save adaptive output
    adapt_MCMC <- list(MCMC_out = MCMC_out, K_fixed = K_fixed, K_MCMC = K_MCMC)
    if (save_res) {
      save(adapt_MCMC, file = adapt_path)
    }
    # Reduce memory burden
    rm(OLCA_params, probit_params, MCMC_out)
  }
  
  #================= FIXED SAMPLER =============================================
  print("Fixed sampler")
  #================= Run fixed sampler to obtain posteriors ====================
  # Initialize OLCA model using fixed number of classes
  alpha <- rep(1, K_fixed) / K_fixed  # Hyperparameter for prior for pi
  eta <- rep(1, d)                 # Hyperparameter for prior for theta
  # Obtain pi, theta, c_all
  OLCA_params <- init_OLCA(alpha = alpha, eta = eta, n = n, K = K_fixed, p = p, 
                           d = d)
  
  # Initialize probit model using fixed number of classes
  # Initialize hyperparameters for xi
  mu0 <- Sig0 <- vector("list", K_fixed)
  for (k in 1:K_fixed) {
    # MVN(0,1) hyperprior for prior mean of xi
    mu0[[k]] <- rnorm(n = q)
    # InvGamma(3.5, 6.25) hyperprior for prior variance of xi. Assume uncorrelated 
    # components and mean variance 2.5 for a weakly informative prior on xi
    Sig0[[k]] <- diag(rinvgamma(n = q, shape = 3.5, scale = 6.25), nrow = q, ncol = q)
  }
  # Obtain xi, z_all                        
  probit_params <- init_probit(mu0 = mu0, Sig0 = Sig0, K = K_fixed, q = q, n = n, 
                               V = V, y_all = y_all, c_all = OLCA_params$c_all)
  
  # Run MCMC algorithm using fixed number of classes
  # Obtain pi_MCMC, theta_MCMC, xi_MCMC, c_all_MCMC, z_all_MCMC, loglik_MCMC
  MCMC_out <- run_MCMC_Rcpp(OLCA_params = OLCA_params, probit_params = probit_params, 
                            n_runs = n_runs, burn = burn, thin = thin, K = K_fixed, 
                            p = p, d = d, n = n, q = q, w_all = w_all, x_mat = x_mat, 
                            y_all = y_all, V = V, alpha = alpha, eta = eta, 
                            Sig0 = Sig0, mu0 = mu0)
  
  # Post-processing to recalibrate labels and remove extraneous empty classes
  # Obtain K_med, pi, theta, xi, loglik_MCMC
  post_MCMC_out <- post_process(MCMC_out = MCMC_out, p = p, d = d, q = q)
  
  # Obtain posterior estimates, reduce number of classes, analyze results
  # Obtain K_red, pi_red, theta_red, xi_red, pi_med, theta_med, xi_med, Phi_med, 
  # c_all, pred_class_probs, loglik_med
  analysis <- analyze_results(MCMC_out = MCMC_out, post_MCMC_out = post_MCMC_out, 
                              n = n, p = p, V = V, y_all = y_all, x_mat = x_mat)
  
  #================= VARIANCE ADJUSTMENT =======================================
  print("Variance adjustment")
  # Create Stan model
  mod_stan <- stan_model(stan_path)
  # Apply variance adjustment for correct coverage
  # Obtain pi_red_adj, theta_red_adj, xi_red_adj, pi_med_adj, theta_med_adj, 
  # xi_med_adj, Phi_med_adj, c_all, pred_class_probs, log_lik_med
  analysis_adj <- var_adjust(mod_stan = mod_stan, analysis = analysis, 
                             K = analysis$K_red, p = p, d = d, n = n, q = q, 
                             x_mat = x_mat, y_all = y_all, V = V, w_all = w_all, 
                             s_all = s_all, clus_id_all = clus_id_all)
  
  runtime <- Sys.time() - start_time
  
  #================= Save and return output ====================================
  res <- list(analysis_adj = analysis_adj, runtime = runtime, 
              data_vars = data_vars, MCMC_out = MCMC_out, 
              post_MCMC_out = post_MCMC_out)
  if (is.null(K_known)) {
    res$K_MCMC <- K_MCMC
  }
  if (save_res) {
    save(res, file = adj_path)
  }
  return(res)
}



#===================== RUN MAIN WSOLCA FUNCTION =================================

# Define directories
wd <- "/n/holyscratch01/stephenson_lab/Users/stephwu18/SWOLCA/"
# wd <- "~/Documents/Github/SWOLCA/"
data_dir <- "Data/"
res_dir <- "Results/July6/"
model_dir <- "Model_Code/"
model <- "wsOFMM"

# Define paths
data_path <- paste0(wd, data_dir, "nhanes1518_adult_low_f_12jul2023.csv")   # Input dataset
adapt_path <- paste0(wd, res_dir, model, "_adapt_nhanesNOEDUC_NOLEG", ".RData")  # Output file
adj_path <- paste0(wd, res_dir, model, "_results_adj_nhanesNOEDUC_NOLEG", ".RData")  # Adjusted output file
stan_path <- paste0(wd, model_dir, "WSOLCA_main.stan")  # Stan file

# Check if results already exist
already_done <- file.exists(adj_path)
if (already_done) {
  print(paste0('NHANES results already exist.'))
} else {
  # Source application preparation helper functions
  source(paste0(wd, model_dir, "app_helper_functions.R"))
  # Read and process data
  data_vars <- process_data(data_path = data_path,
                            covs = c("age_cat", "racethnic", "smoker", "physactive"),
                            formula = "~ age_cat + racethnic + smoker + physactive")
  
  # Source R helper functions
  source(paste0(wd, model_dir, "helper_functions.R"))
  # Source Rcpp functions
  Rcpp::sourceCpp(paste0(wd, model_dir, "main_Rcpp_functions.cpp"))
  
  # Set seed and run model
  set.seed(20230225)
  print(paste0("Running WSOLCA_application..."))
  results_adj <- WSOLCA_app_covs_Rcpp(data_vars = data_vars, 
                                      adapt_path = adapt_path,
                                      adj_path = adj_path, 
                                      stan_path = stan_path, save_res = TRUE, 
                                      n_runs = 20000, burn = 10000, thin = 5,
                                      K_known = NULL)
  print(paste0("Runtime: ", results_adj$runtime))
}


#===================== PLOT OUTPUT =============================================

load(adj_path)
age_categs <- c("[20,40)", "[40,60)", ">=60")
educ_categs <- c("Some College", "HS/GED", "<HS")
racethnic_categs <- c("NH White", "NH Black", "NH Asian", "Hispanic/Latino", "Other/Mixed")
smoker_categs <- c("Non-Smoker", "Smoker")
physactive_categs <- c("Inactive", "Active")
K <- length(res$analysis_adj$pi_med_adj)

# Reorder classes
new_order <- c(3, 2, 5, 4, 1)
# new_order <- c(5, 2, 4, 1, 3)
res <- reorder_classes(res = res, model = "wsOFMM", new_order = new_order)

# Plot pi boxplots for each latent class 
plot_pi_boxplots(res, model = "wsOFMM")

# Plot xi boxplots 
plot_xi_boxplots(res, model = "wsOFMM", age_categs = age_categs, 
                 racethnic_categs = racethnic_categs,
                 educ_categs = educ_categs, smoker_categs = smoker_categs, 
                 physactive_categs = physactive_categs)

# Plot Phi line plots, separately for each covariate 
plot_Phi_line(res, model = "wsOFMM", age_categs = age_categs, 
              racethnic_categs = racethnic_categs,
              educ_categs = educ_categs, smoker_categs = smoker_categs, 
              physactive_categs = physactive_categs)
plot_Phi_line_cis(res, model = "wsOFMM", age_categs = age_categs, 
                  racethnic_categs = racethnic_categs,
                  educ_categs = educ_categs, smoker_categs = smoker_categs, 
                  physactive_categs = physactive_categs)

# Output reference cell coefficients table for xi 
convert_to_ref(xi_med = res$analysis_adj$xi_med_adj, 
               xi_red = res$analysis_adj$xi_red_adj,
               age_categs = age_categs, 
               racethnic_categs = racethnic_categs,
               smoker_categs = smoker_categs, 
               physactive_categs = physactive_categs, format = "latex")

# Plot theta modal consumption levels
p1 <- plot_theta_modes(res = res, model = "wsOFMM")
# Plot theta probabilities
p2 <- plot_theta_probs(res = res, model = "wsOFMM")
ggarrange(p1, p2, common.legend = TRUE, widths = c(0.6, 1), legend = "top")
p1 + theme(legend.position = "top")
p2 + theme(legend.position = "top")

# Table of demographics by hypertension 
res_demog <- res$data_vars$V_data %>%
  mutate(racethnic = factor(racethnic),
         smoker = factor(smoker),
         # educ = factor(educ),
         age_cat = factor(age_cat))
res_demog$Class <- factor(res$analysis_adj$c_all)
res_demog$y <- res$data_vars$y_all
res_demog$w <- res$data_vars$sample_wt
res_demog$clus_id_all <- res$data_vars$clus_id_all
res_demog$sample_wt <- res$data_vars$sample_wt
res_demog$s_all <- res$data_vars$s_all
# Read in HEI scores
hei_data <- read.csv(paste0(wd, data_dir, "nhanes1518_hei_tert_lowF.csv")) %>% 
  filter(!DMDEDUC2 %in% c(7, 9))
res_demog$hei <- hei_data$HEI2015_TOTAL_SCORE
# Rename variables and variable names
res_demog_rename <- res_demog %>%
  mutate(age_cat = factor(age_cat, labels = c("[20, 40)", "[40, 60)", ">=60")),
         racethnic = factor(racethnic, 
                            labels = c("NH White", "NH Black", "NH Asian", 
                                       "Hispanic/Latino", "Other/Mixed")),
         smoker = factor(smoker, labels = c("No", "Yes")),
         y = factor(y, labels = c("No", "Yes"))) %>%
  rename("Age Group" = "age_cat", "Race/Ethnicity" = "racethnic", 
         "Current Smoker" = "smoker", "Physical Activity" = "physactive",
         "Hypertension" = "y")
demogs <- CreateTableOne(vars = c("Age Group", "Race/Ethnicity", "Current Smoker", 
                                  "Physical Activity"),
                         factorVars = c("Age Group", "Race/Ethnicity", 
                                        "Current Smoker", "Physical Activity"),
                         strata = "Hypertension", addOverall = TRUE, 
                         data = res_demog_rename)
demogs_df <- print(demogs, noSpaces = TRUE, showAllLevels = TRUE)
kable(demogs_df, format = "latex")


# Table of demographics by latent class
create_demog_table_pop(res_demog, age_categs, racethnic_categs,
                       smoker_categs, physactive_categs, res = res)
# create_demog_table(res_demog, age_categs, racethnic_categs,
#                    educ_categs, smoker_categs, physactive_categs)


# # Troubleshooting PEP and CrI incongruency
# library(HDInterval)
# est_mcmc <- res$analysis_adj$xi_red_adj[,2,9]
# est_mcmc <- res$analysis_adj$xi_red_adj[, 2, 9] - res$analysis_adj$xi_med_adj[1, 9]
# hdi <- hdi(est_mcmc)
# hist(est_mcmc, breaks = 50)
# hist(est_mcmc[est_mcmc > 0], add = TRUE, col = "orange")
# abline(v=quantile(est_mcmc, 0.025), col = "blue")
# abline(v=quantile(est_mcmc, 0.975), col = "blue")
# abline(v = hdi[1], col = "red")
# abline(v = hdi[2], col = "red")
# quantile(est_mcmc, c(0.025, 0.975))
# hdi(est_mcmc)
# mean(est_mcmc > 0)

#===================== Create colored dendrogram ===============================

dendrogram <- res$post_MCMC_out$dendrogram
dendro_k <- dendro_data_k(dendrogram, K_med)
plot_ggdendro(dendro_k, direction = "tb", branch.size = 0.5) + 
  theme_bw() + scale_color_brewer(palette="Set2") + 
  geom_hline(yintercept=1750, linetype = "dashed") +
  xlab("") + ylab("")

#================= OLD MISCELLANEOUS CODE ======================================
# #================= Plot theta modes ============================================
# 
# # Old Code grouping food items together
# lcmodel %>%
#   ggplot(aes(x = Probability, y = factor(Item, levels = rev(food_items)), 
#              fill = factor(Level))) + 
#   geom_bar(stat = "identity", position = "stack") + 
#   facet_grid(. ~ Class) + 
#   scale_fill_brewer(type="seq", palette="Greys") + 
#   theme_bw() + 
#   labs(x="Item consumption probabilities", y = "Food items",
#        fill ="Item \nconsumption \nlevels") + 
#   theme( axis.text.x=element_text(size=7),
#          panel.grid.major.x=element_blank())
# # axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# # Plot MEAN consumption levels
# est_item_probs <- res$analysis_adj$theta_med_adj
# d <- dim(est_item_probs)[3]
# mean_item_probs <- as.data.frame(apply(est_item_probs, c(1, 2), 
#                                        function(x) round(sum(x * c(1:d)), 1)))
# class_names <- paste0("Class ", 1:(dim(est_item_probs)[2]))
# rownames(mean_item_probs) <- food_items
# colnames(mean_item_probs) <- class_names
# mean_item_probs$Item <- rownames(mean_item_probs)
# mode_plot <- mean_item_probs %>% gather("Class", "Level", -Item) 
# mode_plot %>% ggplot(aes(x=Class, y=factor(Item, levels = rev(food_items)), 
#                          fill=Level)) + 
#   geom_tile(color="black") + 
#   geom_text(aes(label = Level), col="white", cex=3) +
#   scale_fill_gradient(trans="reverse") + 
#   theme(legend.position="none") +
#   ylab("Item") + xlab("Latent Class") +
#   ggtitle("Estimated modal item \nconsumption levels") + 
#   theme(text = element_text(size = 15))
# 
# 
# #================= Plot marginal model xi ======================================
# xi_dims <- dim(res$analysis_adj$xi_red_adj)
# K <- xi_dims[2]
# xi_red <- as.data.frame(t(matrix(res$analysis_adj$xi_red_adj, 
#                                  nrow = xi_dims[1], 
#                                  ncol = K*xi_dims[3], byrow = FALSE)))
# # Class variable: 1-K, 1-K, 1-K...
# xi_red$Class <- as.character(factor(c(rep(1:K, times = xi_dims[3]))))
# # Covariate level: RefxK, [40,60)xK, >=60xK... 
# xi_red$Covariate_Level <- rep("Ref", K)
# xi_red$Covariate <- c(rep("Ref", K))
# xi_red_plot <- xi_red %>% 
#   pivot_longer(cols = -c("Class", "Covariate", "Covariate_Level"), 
#                names_to = "iter", values_to = "value")
# xi_red_plot %>% ggplot(aes(x = Covariate_Level, y = value, group = Class, fill = Class)) +
#   theme_bw() +
#   geom_boxplot() +
#   facet_grid(.~Covariate, labeller = label_both) +
#   ggtitle(as.expression(bquote("Parameter estimation for "~xi~" ")))
# 
# 
# #================= Plot Phi using Phi_med ======================================
# df_Phi <- data.frame(Stratum = factor(res$data_vars$s_all), 
#                      Class = factor(res$analysis_adj$c_all), 
#                      Phi = res$analysis_adj$Phi_med_adj)
# ### Plot marginal Phi
# df_Phi %>% ggplot(aes(x = Class, y = Phi, col = Class)) +
#   theme_bw() +
#   geom_point() +
#   ylim(0, 1) +
#   ggtitle("Parameter estimation for conditional outcome probabilities") +
#   ylab("P(Y=1|Class)")
# 
# df_Phi %>% ggplot(aes(x = Class, y = Phi, fill = Class)) +
#   theme_bw() +
#   geom_violin() + geom_boxplot(width = 0.1) +
#   ylim(0, 1) +
#   ggtitle("Parameter estimation for conditional outcome probabilities") +
#   ylab("Estimated P(Y=1|Class)")
# 
# ### Plot Phi with additional covariates
# df_Phi <- data.frame(pnorm(res$analysis_adj$xi_med_adj))
# colnames(df_Phi) <- c("Ref", age_categs[-1], racethnic_categs[-1], 
#                       # educ_categs[-1], 
#                       smoker_categs[-1],
#                       physactive_categs[-1])
# # colnames(df_Phi) <- c(age_categs, educ_categs, racethnic_categs, smoker_categs)
# # colnames(df_Phi) <- colnames(res$data_vars$V)
# df_Phi$Class <- factor(1:K)
# df_Phi <- df_Phi %>% gather("Covariate_Level", "Phi", -Class)
# df_Phi$Covariate <- c(rep("Ref", K),
#                       rep("Age", K*(length(age_categs) - 1)), 
#                       rep("Race/Ethnicity", K*(length(racethnic_categs) - 1)), 
#                       # rep("Education", K*(length(educ_categs) - 1)), 
#                       rep("Smoking", K*(length(smoker_categs) - 1)),
#                       rep("Physical Activity", K*(length(physactive_categs) - 1)))
# # df_Phi$Covariate <- c(rep("Age", K*6), rep("Education", K*3), 
# #                       rep("Race/Ethnicity", K*5), rep("Smoking", K*2))
# # df_Phi %>% ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
# #   theme_bw() + scale_color_brewer(palette="Set2") + 
# #   geom_point() + geom_line() + 
# #   facet_grid(~Covariate) +
# #   ggtitle("Parameter estimation for conditional outcome probabilities") +
# #   ylab("P(Y=1|-")
# p_age <- df_Phi %>% filter(Covariate == "Age") %>%
#   ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
#   theme_bw() + scale_color_brewer(palette="Set2") + 
#   geom_point() + geom_line() + ylim(0,1) + 
#   ylab("Hypertension Probability") + xlab("Age") +
#   theme(axis.text=element_text(size=10),
#         axis.title=element_text(size=10))
# # p_educ <- df_Phi %>% filter(Covariate == "Education") %>%
# #   ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
# #   theme_bw() + scale_color_brewer(palette="Set2") + 
# #   geom_point() + geom_line() +  ylim(0,1) + 
# #   ylab("Hypertension Probability") + xlab("Education") +
# #   theme(axis.text=element_text(size=10),
# #         axis.title=element_text(size=10))
# p_race <- df_Phi %>% filter(Covariate == "Race/Ethnicity") %>%
#   mutate(Covariate_Level = factor(Covariate_Level, levels = racethnic_categs)) %>%
#   ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
#   theme_bw() + scale_color_brewer(palette="Set2") + 
#   geom_point() + geom_line() +  ylim(0,1) + 
#   ylab("Hypertension Probability") + xlab("Race/Ethnicity") +
#   theme(axis.text=element_text(size=10),
#         axis.title=element_text(size=10)) +
#   ggtitle("Hypertension Risk by Dietary Pattern and Race/Ethnicity")
# p_smoker <- df_Phi %>% filter(Covariate == "Smoking") %>%
#   ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
#   theme_bw() + scale_color_brewer(palette="Set2") + 
#   geom_point() + geom_line() +  ylim(0,1) + 
#   ylab("Hypertension Probability") + xlab("Smoking") +
#   theme(axis.text=element_text(size=10),
#         axis.title=element_text(size=10))
# p_physactive <- df_Phi %>% filter(Covariate == "Physical Activity") %>%
#   ggplot(aes(x = Covariate_Level, y = Phi, col = Class, group = Class)) +
#   theme_bw() + scale_color_brewer(palette="Set2") + 
#   geom_point() + geom_line() +  ylim(0,1) + 
#   ylab("Hypertension Probability") + xlab("Physical Activity") +
#   theme(axis.text=element_text(size=10),
#         axis.title=element_text(size=10))
# ggarrange(p_age, p_race, p_smoker, p_physactive, ncol = 3, nrow=2,
#           common.legend = TRUE, legend = "top", widths = c(1,0.7,1,0.5))
# # df_Phi %>% ggplot(aes(x = Stratum, y = Phi, col = Class, group = Class)) + 
# #   theme_bw() + 
# #   geom_point() + 
# #   geom_line() + 
# #   ylim(0, 0.5) + 
# #   ggtitle("Parameter estimation for conditional outcome probabilities") +
# #   ylab("P(Y=1|Class, Stratum")
# 
# #### Plot Phi line plot OLD VERSION
# df_Phi <- data.frame(Stratum = factor(res$data_vars$s_all), 
#                      Class = factor(res$analysis_adj$c_all), 
#                      Phi = res$analysis_adj$Phi_med)
# df_Phi %>% ggplot(aes(x = Class, y = Phi, col = Stratum, group = Stratum)) + 
#   theme_bw() + 
#   geom_point() + 
#   geom_line() + 
#   ylim(0, 0.5) + 
#   ggtitle("Parameter estimation for conditional outcome probabilities") +
#   ylab("P(Y=1|Class, Stratum)")
# df_Phi %>% ggplot(aes(x = Stratum, y = Phi, col = Class, group = Class)) + 
#   theme_bw() + 
#   geom_point() + 
#   geom_line() + 
#   ylim(0, 0.5) + 
#   ggtitle("Parameter estimation for conditional outcome probabilities") +
#   ylab("P(Y=1|Class, Stratum)")
# 
# 
# #============ Plot Phi boxplots separately for each covariate ==================
# Phi_red_plot <- xi_red_plot %>%
#   mutate(value = pnorm(value))
# Phi_red_plot %>% filter(Covariate == "Age") %>%
#   ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
#   theme_bw() + scale_fill_brewer(palette="Set2") + 
#   geom_boxplot() + ylab("Hypertension Risk") + xlab("Age")
# # Phi_red_plot %>% filter(Covariate == "Education") %>%
# #   ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
# #   theme_bw() + scale_fill_brewer(palette="Set2") + 
# #   geom_boxplot() + ylab("Hypertension Risk") + xlab("Education")
# Phi_red_plot %>% filter(Covariate == "Race/Ethnicity") %>%
#   ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
#   theme_bw() + scale_fill_brewer(palette="Set2") + 
#   geom_boxplot() + ylab("Hypertension Risk") + xlab("Race/Ethnicity")
# Phi_red_plot %>% filter(Covariate == "Smoking") %>%
#   ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
#   theme_bw() + scale_fill_brewer(palette="Set2") + 
#   geom_boxplot() + ylab("Hypertension Risk") + xlab("Smoking Status")

# Testing
# results_adj <- WSOLCA_app_Rcpp(data_path = data_path, res_path = res_path,
#                                 adj_path = adj_path, stan_path = stan_path,
#                                 save_res = FALSE, n_runs = 60, burn = 30,
#                                 thin = 3)
