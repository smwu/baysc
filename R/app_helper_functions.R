#===================================================
## Helper functions for data application
## Programmer: SM Wu   
## Data: NHANES Application with Covariates  
## Date updated: 2023/07/14
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
library(kableExtra)
library(ggdendro)

#===================== Process data ============================================
# `process_data` reads in data and converts it into a processed format for the 
# 'WSOLCA_app_covs_Rcpp' function
# Input:
#   data_path: String specifying path for input data
#   covs: String vector of covariates to include in regression. Default is NULL
#   formula: String specifying formula for regression model
# Output: list 'data_vars' containing the following objects:
#   x_mat: Matrix of categorical exposure variables for all individuals; nxp
#   y_all: Vector of binary outcome for all individuals; nx1 
#   s_all: Vector of stratum indicator for all individuals; nx1
#   clus_id_all: Vector of cluster indicator for all individuals; nx1
#   sample_wt: Vector of sampling weights for all individuals; nx1
#   V: Model matrix for regression; nxq
#   V_data: Matrix of cleaned variables to use in regression; nx(num_vars)
# Example: process_data(data_path = "C:/Documents", 
#                       covs = c("age_cat", "racethnic", "smoker", "physactive"),
#                       formula = "~ age_cat + racethnic + smoker + physactive")
process_data <- function(data_path, covs = NULL, formula = NULL) {
  # Read in data
  data_vars_f_low <- read.csv(data_path)
  # Drop those with refused or unknown education data
  # Complete n = 2003
  data_vars_f_low <- data_vars_f_low %>% filter(!DMDEDUC2 %in% c(7, 9))
  # # Drop those with race/ethnic equal to Other/Mixed 
  # # Complete n = 1921
  # data_vars_f_low <- data_vars_f_low %>% filter(!RIDRETH3 == 7)
  # Drop legumes (vegs) because duplicate of legumes (proteins), just computed
  # as cup eq vs oz eq
  data_vars_f_low <- data_vars_f_low %>% select(-leg_veg)
  
  # Obtain exposure and outcome data
  x_mat <- as.matrix(data_vars_f_low %>% select(citrus:drinks))
  y_all <- data_vars_f_low$BP_flag
  
  # Get stratum IDs
  s_all <- data_vars_f_low$SDMVSTRA  # 30 strata
  # Create unique nested PSU IDs using stratum IDs
  clus_id_all <- s_all * 10 + data_vars_f_low$SDMVPSU
  # Get sampling weights
  sample_wt <- data_vars_f_low$dietwt4yr
  
  # Other covariates to be included in the probit model 
  if (is.null(covs)) {
    n <- length(y_all)
    V <- matrix(1, nrow = n) 
    q <- 1
  } else {
    V_data <- data_vars_f_low %>% 
      select(RIDAGEYR, RIDRETH3, DMDEDUC2, smoker, Phys_Active) %>%
      mutate(
        age_cat = factor(case_when(  
          20 <= RIDAGEYR & RIDAGEYR <= 39 ~ 1,
          40 <= RIDAGEYR & RIDAGEYR <= 59 ~ 2,
          RIDAGEYR >= 60 ~ 3)),
        # RIDRETH3: Race/Hispanic origin w/ NH Asian: 1=Mex_Amer, 2=Other_Hisp,
        # 3=NH_White, 4=NH_Black, 6=NH_Asian, 7=Other/Mixed
        racethnic = factor(case_when(
          RIDRETH3 == 3 ~ 1,  # NH White
          RIDRETH3 == 4 ~ 2,  # NH Black
          RIDRETH3 == 6 ~ 3,  # NH Asian
          # RIDRETH3 == 1 ~ 4,  # Mexican-American
          # RIDRETH3 == 2 ~ 5,  # Other Hispanic
          RIDRETH3 %in% c(1, 2) ~ 4,  # Mexican-American/Other Hispanic
          # RIDRETH3 == 7 ~ 6,  # Other/Mixed
          RIDRETH3 == 7 ~ 5,  # Other/Mixed
          .default = NA)),
        # DMDEDUC2: Education level for adults 20+: 1= <9th, 2=9-11th, 3=HS/GED,
        # 4=Some college/AA, 5=college grad or above, 7=refused, 9=don't know
        educ = factor(case_when(
          DMDEDUC2 %in% c(4, 5) ~ 1,  # At least some college
          DMDEDUC2 == 3 ~ 2,  # HS/GED
          DMDEDUC2 %in% c(1, 2) ~ 3,  # Less than HS
          .default = NA)),
        smoker = factor(smoker),
        physactive = factor(Phys_Active, levels=c(0,1), 
                            labels=c("Inactive", "Active")),
        .keep = "unused")
    V_data <- V_data %>% select(all_of(covs))
    
    # Regression design matrix without class assignment, nxq
    # Exclude stratifying variable as well
    if (is.null(formula)) {
      stop("Error: Need to specify a string formula for the probit regression.")
    }
    V <- model.matrix(as.formula(formula), data = V_data)
  }
  
  # Return processed data
  data_vars <- list(x_mat = x_mat, y_all = y_all, s_all = s_all, 
                    clus_id_all = clus_id_all, sample_wt = sample_wt, 
                    V = V)
  if (!is.null(covs)) {
    data_vars$V_data <- V_data
  }
  return(data_vars)
}

#============= Reorder latent classes ==========================================
reorder_classes <- function(res, model, new_order) {
  if (model == "wsOFMM") {
    res_new <- res
    res_new$analysis_adj$pi_red_adj <- res$analysis_adj$pi_red_adj[, new_order]
    res_new$analysis_adj$theta_red_adj <- res$analysis_adj$theta_red_adj[, , new_order, ]
    res_new$analysis_adj$xi_red_adj <- res$analysis_adj$xi_red_adj[, new_order, ]
    res_new$analysis_adj$pi_med_adj <- res$analysis_adj$pi_med_adj[new_order]
    res_new$analysis_adj$theta_med_adj <- res$analysis_adj$theta_med_adj[, new_order, ]
    res_new$analysis_adj$xi_med_adj <- res$analysis_adj$xi_med_adj[new_order, ]
    for (i in 1:5) {
      res_new$analysis_adj$c_all[res$analysis_adj$c_all == new_order[i]] <- i
    }
  } else {
    res_new <- res
    res_new$analysis$pi_red <- res$analysis$pi_red[, new_order]
    res_new$analysis$theta_red <- res$analysis$theta_red[, , new_order, ]
    res_new$analysis$xi_red <- res$analysis$xi_red[, new_order, ]
    res_new$analysis$pi_med <- res$analysis$pi_med[new_order]
    res_new$analysis$theta_med <- res$analysis$theta_med[, new_order, ]
    res_new$analysis$xi_med <- res$analysis$xi_med[new_order, ]
    for (i in 1:5) {
      res_new$analysis$c_all[res$analysis$c_all == new_order[i]] <- i
    }
  }
  return(res_new)
}

#============= Create demographics x latent class table ========================
# Converts from combination coding to P(Y=1|-) conditional probabilities
# Input: xi_comb: Matrix of xi parameter estimates. Kxq
# Output: List containing the following objects
#   Phi_age: Dataframe of hypertension probs for age categories. Kx3
#   Phi_age: Dataframe of hypertension probs for race/ethnicity categories. Kx5
#   Phi_smoker: Dataframe of hypertension probs for smoking categories. Kx2
#   Phi_phys: Dataframe of hypertension probs for physical activity categories. Kx2
convert_to_probs <- function(xi_comb) {
  # Assumes the following column order for xi_comb: 
  # Intercept + Age40 + Age60 + NH_Black + NH_Asian + Hispanic + Other + Smoker + Inactive
  K <- nrow(xi_comb)
  Phi_age <- data.frame(Class = factor(1:5),
                        Age20 = pnorm(xi_comb[, 1]),
                        Age40 = pnorm(xi_comb[, 1] + xi_comb[, 2]),
                        Age60 = pnorm(xi_comb[, 1] + xi_comb[, 3]))
  Phi_racethnic <- data.frame(Class = factor(1:5),
                              NH_White = pnorm(xi_comb[, 1]),
                              NH_Black = pnorm(xi_comb[, 1] + xi_comb[, 4]),
                              NH_Asian = pnorm(xi_comb[, 1] + xi_comb[, 5]),
                              Hispanic = pnorm(xi_comb[, 1] + xi_comb[, 6]),
                              Other = pnorm(xi_comb[, 1] + xi_comb[, 7]))
  Phi_smoker <- data.frame(Class = factor(1:5),
                           Non_smoker = pnorm(xi_comb[, 1]),
                           Smoker = pnorm(xi_comb[, 1] + xi_comb[, 8]))
  Phi_phys <- data.frame(Class = factor(1:5),
                         Inactive = pnorm(xi_comb[, 1]),
                         Active = pnorm(xi_comb[, 1] + xi_comb[, 9]))
  return(list(Phi_age = Phi_age, Phi_racethnic = Phi_racethnic, 
              Phi_smoker = Phi_smoker, Phi_phys = Phi_phys))
}

# Normalize a vector to sum to 1
normalize <- function(x) {
  return(x / sum(x))
}

# Obtain the row counts and percentages across latent classes for a given 
# demographic variable
# Inputs:
#   res_demog: Dataframe containing demographic vars, latent classes, and weights
#   cov: String specifying covariate of interest. `sample` gives sample counts 
# `population` gives population counts
#   level: Numeric level of interest for the factor variables. Default is NULL
get_row_props <- function(res_demog, cov, level = NULL) {
  if (cov == "sample") {
    lc_counts <- sapply(1:5, function(x) nrow(res_demog %>% filter(Class == x)))
  } else if (cov == "population") {
    lc_counts <- sapply(1:5, function(x) 
      round(sum(res_demog %>% filter(Class == x) %>% select(w)) / 1000, 0))
  } else {
    lc_counts <- sapply(1:5, function(x) 
      nrow(res_demog %>% filter(!!sym(cov) == level & Class == x)))
  }
  overall <- sum(lc_counts)
  lc_percents <- format(round(lc_counts / overall * 100, 1), nsmall = 1)
  lc_text <- sapply(1:5, function(i) paste0(lc_counts[i], " (", lc_percents[i], "%)"))
  output <- c(lc_text, overall)
  return(output)
}


# Obtain the column percentages and standard errors across latent classes for a 
# given demographic variable, calculated for the population using survey weights
# and linearization variance
# Inputs:
#   res_demog: Dataframe containing demographic vars, latent classes, and weights
#   cov: String specifying covariate of interest
#   categ_labels: String vector specifying category lavels
#   res: Results from model fit, if necessary. Default is NULL
# Output: 'output' dataframe with column percentage and standard errors for all
# levels of the input 'cov'
get_pop_props <- function(res_demog, cov, categ_labels, res = NULL) {
  # Set survey design
  nhanes_svy <- svydesign(id = ~clus_id_all,
                          weights = ~sample_wt,
                          strata=~s_all,
                          data = res_demog)
  # Number of categories
  n_levels <- length(categ_labels)
  # Initialize output table
  output <- as.data.frame(matrix(NA, nrow = n_levels, ncol = 6))
  colnames(output) <- c(paste0("Pattern", 1:5), "Total")
  # If sample, just tabulate number of individuals in each latent class. No SE
  if (cov == "sample") {
    temp_tot <- nrow(res_demog)
    temp_class <- sapply(1:5, function(x) nrow(res_demog %>% filter(Class == x)))
    output[, 6] <- 100
    output[, -6] <- format(round(temp_class / temp_tot * 100, 1), nsmall = 1)
    # If population, estimate pop mean and SE for each class using posterior pi
  } else if (cov == "population") {
    output[, 6] <- 100
    output[, -6] <- apply(res$analysis_adj$pi_red_adj, 2, function(x) 
      paste0(format(round(mean(x)*100, 1), nsmall = 1), " (", 
             format(round(sd(x)*100, 1), nsmall = 1), ")"))
    # temp_class <- as.data.frame(svymean(~Class, nhanes_svy))
    # output[, 6] <- 100
    # output[, -6] <- sapply(1:5, function(x)
    #   paste0(format(round(temp_class[x, 1]*100, 1), nsmall = 1), " (", 
    #          format(round(temp_class[x, 2]*100, 1), nsmall = 1), ")"))
    # If covariate, estimate pop mean and SE for each class and covariate level
  } else if (cov == "hei") {
    temp_tot <- as.data.frame(svymean(as.formula(paste0("~", cov)), nhanes_svy))
    temp_class <- as.data.frame(svyby(as.formula(paste0("~", cov)), ~Class, 
                                      nhanes_svy, svymean))
    output[, 6] <- paste0(format(round(temp_tot$mean, 1), nsmall = 1), " (",
                          format(round(temp_tot$hei, 1), nsmall = 1), ")")
    output[, -6] <- sapply(1:5, function(x)
      paste0(format(round(temp_class[x, 1+1:n_levels], 1), nsmall = 1), " (", 
             format(round(temp_class[x, 1+n_levels+1:n_levels], 1), nsmall = 1), 
             ")"))
  } else {
    temp_tot <- as.data.frame(svymean(as.formula(paste0("~", cov)), nhanes_svy))
    temp_class <- as.data.frame(svyby(as.formula(paste0("~", cov)), ~Class, 
                                      nhanes_svy, svymean))
    output[, 6] <- paste0(format(round(temp_tot$mean*100, 1), nsmall = 1), " (",
                          format(round(temp_tot$SE*100, 1), nsmall = 1), ")")
    output[, -6] <- sapply(1:5, function(x)
      paste0(format(round(temp_class[x, 1+1:n_levels]*100, 1), nsmall = 1), " (", 
             format(round(temp_class[x, 1+n_levels+1:n_levels]*100, 1), nsmall = 1), 
             ")"))
  }
  return(output)
}

create_demog_table_pop <- function(res_demog, age_categs, racethnic_categs,
                                   smoker_categs, physactive_categs, res) {
  demog_df <- as.data.frame(matrix(NA, nrow = 15, ncol = 8))
  colnames(demog_df) <- c("Variable", "Level", "Multicultural", 
                          "Healthy", 
                          "Western", "Restrict Veg", 
                          "Restrict Amer", "Overall")
  # demog_df[, 1] <- c("N in 1000s", "n",
  #                    paste0("Age ", age_categs), racethnic_categs,  
  #                    smoker_categs, physactive_categs)
  demog_df[, 1] <- c("Sample Size", "Population Size", "HEI-2015 Score",
                     c("Age Group", rep("", length(age_categs) - 1)),
                     c("Race/Ethnicity", rep("", length(racethnic_categs) - 1)),
                     c("Current Smoker", rep("", length(smoker_categs) - 1)),
                     c("Physical Activity", rep("", length(physactive_categs) - 1)))
  demog_df[, 2] <- c("", "", "", age_categs, racethnic_categs,  
                     smoker_categs, physactive_categs)
  demog_df[1, -c(1, 2)] <- get_pop_props(res_demog, cov = "sample", 
                                         categ_labels = "")
  demog_df[2, -c(1, 2)] <- get_pop_props(res_demog, cov = "population", 
                                         categ_labels = "", res = res)
  demog_df[3, -c(1, 2)] <- get_pop_props(res_demog, cov = "hei", categ_labels = "")
  demog_df[4:6, -c(1, 2)] <- get_pop_props(res_demog, cov = "age_cat", 
                                           categ_labels = age_categs)
  demog_df[7:11, -c(1, 2)] <- get_pop_props(res_demog, cov = "racethnic", 
                                            categ_labels = racethnic_categs)
  demog_df[12:13, -c(1, 2)] <- get_pop_props(res_demog, cov = "smoker", 
                                             categ_labels = smoker_categs)
  demog_df[14:15, -c(1, 2)] <- get_pop_props(res_demog, cov = "physactive", 
                                             categ_labels = physactive_categs)
  demog_df %>% kable(format = "latex", booktabs = TRUE)
}

# 'create_demog_table' creates a table of latent class distributions across
# sociodemographic characteristics
# Inputs:
#   res_demog: Dataframe containing demographic vars, latent classes, and weights
#   age_categs: String vector of age category names
#   racethnic_categs: String vector of race/ethnicity category names
#   educ_categs: String vector of education category names
#   smoker_categs: String vector of smoker category names
#   physactive_categs: String vector of physical activity category names
# Output: Formatted demographics table in latex
create_demog_table <- function(res_demog, age_categs, racethnic_categs,
                               smoker_categs, physactive_categs) {
  demog_df <- as.data.frame(matrix(NA, nrow = 14, ncol = 7))
  colnames(demog_df) <- c("Variable", "Overall", "Pattern 1", "Pattern 2", 
                          "Pattern 3", "Pattern 4", "Pattern 5")
  demog_df[, 1] <- c("Size in Population", "Size in Sample",
                     paste0("Age ", age_categs), racethnic_categs,  
                     smoker_categs, physactive_categs)
  demog_df[1, -1] <- get_row_props(res_demog, cov = "population", level = NULL)
  demog_df[2, -1] <- get_row_props(res_demog, cov = "sample", level = NULL)
  demog_df[3, -1] <- get_row_props(res_demog, cov = "age_cat", level = 1)
  demog_df[4, -1] <- get_row_props(res_demog, cov = "age_cat", level = 2)
  demog_df[5, -1] <- get_row_props(res_demog, cov = "age_cat", level = 3)
  demog_df[6, -1] <- get_row_props(res_demog, cov = "racethnic", level = 1)
  demog_df[7, -1] <- get_row_props(res_demog, cov = "racethnic", level = 2)
  demog_df[8, -1] <- get_row_props(res_demog, cov = "racethnic", level = 3)
  demog_df[9, -1] <- get_row_props(res_demog, cov = "racethnic", level = 4)
  demog_df[10, -1] <- get_row_props(res_demog, cov = "racethnic", level = 5)
  demog_df[11, -1] <- get_row_props(res_demog, cov = "smoker", level = 0)
  demog_df[12, -1] <- get_row_props(res_demog, cov = "smoker", level = 1)
  demog_df[13, -1] <- get_row_props(res_demog, cov = "physactive", level = "Inactive")
  demog_df[14, -1] <- get_row_props(res_demog, cov = "physactive", level = "Active")
  # demog_df[15, -1] <- get_row_props(res_demog, cov = "educ", level = 1)
  # demog_df[16, -1] <- get_row_props(res_demog, cov = "educ", level = 2)
  # demog_df[17, -1] <- get_row_props(res_demog, cov = "educ", level = 3)
  demog_df %>% kable(format = "latex", booktabs = TRUE)
}

#========================= Plotting functions ==================================
# Plot modal consumption levels
# Inputs:
#   res: Output from WSOLCA, SOLCA, or WOLCA models
#   model: String specifying model. Must be one of 'wsOFMM', 'sOFMM', or 'wOFMM'
# Output: Plots modal item consumption levels
plot_theta_modes <- function(res, model) {
  if (model == "wsOFMM") {
    est_item_probs <- res$analysis_adj$theta_med_adj
  } else {
    est_item_probs <- res$analysis$theta_med
  }
  mode_item_probs <- as.data.frame(apply(est_item_probs, c(1, 2), which.max))
  food_items <- c("Citrus/Melon/Berries", "Other Fruits", "Fruit Juice", 
                  "Dark Green Vegs", "Tomatoes", "Oth Red/Orange Vegs", 
                  "Potatoes", "Other Starchy Vegs", "Other Vegetables",  
                  "Whole Grains", "Refined Grains", "Meat", "Cured Meats", 
                  "Organ Meat", "Poultry", "Seafood (High n-3)", "Seafood (Low n-3)",
                  "Eggs", "Soybean Products", "Nuts and Seeds", "Legumes (Protein)", 
                  "Milk", "Yogurt", "Cheese", "Oils", "Solid Fats", "Added Sugar", 
                  "Alcoholic Drinks")
  # class_names <- paste0("Class ", 1:(dim(est_item_probs)[2]))
  class_names <- 1:(dim(est_item_probs)[2])
  rownames(mode_item_probs) <- food_items
  colnames(mode_item_probs) <- class_names
  mode_item_probs$Item <- rownames(mode_item_probs)
  mode_plot <- mode_item_probs %>% gather("Class", "Level", -Item) 
  mode_plot %>% ggplot(aes(x=Class, y=factor(Item, levels = rev(food_items)), 
                           fill=factor(Level))) + 
    geom_tile(color="black", linewidth = 0.3) + 
    # scale_fill_brewer(palette = "Spectral", direction = -1) + 
    scale_fill_brewer(type="seq", palette="RdYlBu", direction = -1,
                      name = "Consumption Level",
                      labels = c("None", "Low", "Med", "High")) +
    # scale_fill_manual(values = c("#89B0DF", "#607CAD", "#394C7E", "#0F2151"),
    #                   name = "Level", 
    #                   labels = c("1: None", "2: Low", "3: Med", "4: High"))+ 
    # scale_fill_gradient(trans = "reverse", high = "#89B0DF", 
    #                     low = "#0F2151") +
    # geom_text(aes(label = Level), col="white", cex=3) +
    # theme(legend.position="none") +
    xlab("Dietary Pattern") + ylab("") +
    theme_classic() +
    scale_x_discrete(labels=c("1" = "1\nMulticultural", "2" = "2\nHealthy\nAmerican",
                              "3" = "3\nWestern", "4" = "4\nRestricted\nVegetarian",
                              "5" = "5\nRestricted\nAmerican")) + 
    theme(text = element_text(size = 15),
          axis.text.x = element_text(size = 12, color = "black"), 
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 13, color = "black", face = "bold"),
          legend.title = element_text(size = 13, color = "black", face = "bold"),
          legend.text = element_text(size = 11, color = "black"),
          legend.position = "right", legend.box.spacing = unit(0.2, "pt"))
  
}

### Plot probabilities of consumption levels, grouping classes together
# Inputs:
#   res: Output from WSOLCA, SOLCA, or WOLCA models
#   model: String specifying model. Must be one of 'wsOFMM', 'sOFMM', or 'wOFMM'
# Output: Plots probability of item consumption levels
plot_theta_probs <- function(res, model) {
  if (model == "wsOFMM") {
    est_item_probs <- res$analysis_adj$theta_med_adj
  } else {
    est_item_probs <- res$analysis$theta_med
  }
  mode_item_probs <- as.data.frame(apply(est_item_probs, c(1, 2), which.max))
  food_items <- c("Citrus/Melon/Berries", "Other Fruits", "Fruit Juice", 
                  "Dark Green Vegs", "Tomatoes", "Oth Red/Orange Vegs", 
                  "Potatoes", "Other Starchy Vegs", "Other Vegetables",  
                  "Whole Grains", "Refined Grains", "Meat", "Cured Meats", 
                  "Organ Meat", "Poultry", "Seafood (High n-3)", "Seafood (Low n-3)",
                  "Eggs", "Soybean Products", "Nuts and Seeds", "Legumes (Protein)", 
                  "Milk", "Yogurt", "Cheese", "Oils", "Solid Fats", "Added Sugar", 
                  "Alcoholic Drinks")
  # class_names <- paste0("Class ", 1:(dim(est_item_probs)[2]))
  class_names <- 1:(dim(est_item_probs)[2])
  dimnames(est_item_probs)[[1]] <- food_items
  dimnames(est_item_probs)[[2]] <- class_names
  lcmodel <- reshape2::melt(est_item_probs, level=2)
  colnames(lcmodel) <- c("Item", "Class", "Level", "Probability")
  lcmodel %>%
    ggplot(aes(x = factor(Class, levels = 1:5), y = Probability,
               fill = factor(Level))) + 
    geom_bar(stat = "identity", position = "stack") + 
    facet_wrap(factor(Item, levels = food_items) ~ ., nrow = 4) + 
    scale_fill_brewer(type="seq", palette="RdYlBu", direction = -1,
                      name = "Consumption Level",
                      labels = c("None", "Low", "Med", "High")) +
    # scale_fill_viridis(discrete = TRUE, 
    #                    name = "Level", 
    #                    labels = c("1: None", "2: Low", "3: Med", "4: High")) + 
    # scale_fill_scico_d(palette = "batlow") + 
    theme_bw() + 
    labs(x = "Dietary Pattern", y="Consumption Level Probability") + 
    theme(panel.grid.major.x=element_blank(),
          axis.text.x = element_text(size = 11, color = "black"), 
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 13, color = "black", face = "bold"),
          axis.title.y = element_text(size = 13, color = "black", face = "bold"),
          legend.title = element_text(size = 13, color = "black", face = "bold"),
          legend.text = element_text(size = 11, color = "black"),
          legend.position = "right", legend.box.spacing = unit(0, "pt"),
          strip.text = element_text(size = 9),
          strip.background = element_rect(fill = "gray90"))
}


#=========== Plot Phi line plots, separately for each covariate ================
# Line plots of Phi values for all covariates
plot_Phi_line <- function(res, model, age_categs, racethnic_categs,
                          educ_categs, smoker_categs, physactive_categs) {
  if (model == "wsOFMM") {
    Phi_dfs <- convert_to_probs(res$analysis_adj$xi_med_adj)
  } else {
    Phi_dfs <- convert_to_probs(res$analysis$xi_med)
  }
  
  # Age
  p1 <- Phi_dfs$Phi_age %>% 
    pivot_longer(cols = -Class, names_to = "Age", values_to = "Phi") %>%
    ggplot(aes(x = factor(Age, levels = c("Age20", "Age40", "Age60"), 
                          labels = age_categs), 
               y = Phi, group = Class, col = Class)) + 
    theme_bw() + 
    scale_color_brewer(palette = "Set2", 
                       labels = c("Multicultural", "Healthy American", "Western", 
                                  "Restricted Vegetarian", "Restricted American")) +
    labs(col = "Dietary Pattern") +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +
    ylab("Probability of Hypertension") + xlab("Age") +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(size = 11, color = "black"), 
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 13, color = "black", face = "bold"),
          axis.title.y = element_text(size = 13, color = "black", face = "bold"),
          legend.title = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 13, color = "black")) +
    theme(plot.margin = unit(c(0.2,0,0.2,0.2), 'lines'))
  # Race/Ethnicity
  p2 <- Phi_dfs$Phi_racethnic %>% 
    pivot_longer(cols = -Class, names_to = "Race", values_to = "Phi") %>%
    ggplot(aes(x = factor(Race, levels = c("NH_White", "NH_Black", "NH_Asian", 
                                           "Hispanic", "Other"),
                          labels = racethnic_categs), 
               y = Phi, group = Class, col = Class)) + 
    theme_bw() + 
    scale_color_brewer(palette = "Set2", 
                       labels = c("Multicultural", "Healthy American", "Western", 
                                  "Restricted Vegetarian", "Restricted American")) +
    labs(col = "Dietary Pattern") +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +    
    ylab("") + xlab("Race/Ethnicity") +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(size = 11, color = "black"), 
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 13, color = "black", face = "bold"),
          axis.title.y = element_text(size = 13, color = "black", face = "bold"),
          legend.title = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 13, color = "black")) +
    theme(plot.margin = unit(c(0.2,0,0.2,0), 'lines'))
  # Smoking Status
  p3 <- Phi_dfs$Phi_smoker %>% 
    pivot_longer(cols = -Class, names_to = "Smoker", values_to = "Phi") %>%
    ggplot(aes(x = factor(Smoker, levels = c("Non_smoker", "Smoker"),
                          labels = smoker_categs), 
               y = Phi, group = Class, col = Class)) + 
    theme_bw() + 
    scale_color_brewer(palette = "Set2", 
                       labels = c("Multicultural", "Healthy American", "Western", 
                                  "Restricted Vegetarian", "Restricted American")) +
    labs(col = "Dietary Pattern") +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +    
    ylab("") + xlab("Smoking Status") +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(size = 11, color = "black"), 
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 13, color = "black", face = "bold"),
          axis.title.y = element_text(size = 13, color = "black", face = "bold"),
          legend.title = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 13, color = "black")) +
    theme(plot.margin = unit(c(0.2,0,0.2,0), 'lines'))
  # Physical Activity
  p4 <- Phi_dfs$Phi_phys %>% 
    pivot_longer(cols = -Class, names_to = "Physactive", values_to = "Phi") %>%
    ggplot(aes(x = factor(Physactive, levels = c("Inactive", "Active"),
                          labels = physactive_categs), 
               y = Phi, group = Class, col = Class)) + 
    theme_bw() + 
    scale_color_brewer(palette = "Set2", 
                       labels = c("Multicultural", "Healthy American", "Western", 
                                  "Restricted Vegetarian", "Restricted American")) +
    labs(col = "Dietary Pattern") +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +    
    ylab("") + xlab("Physical Activity") + 
    theme(text = element_text(size = 15),
          axis.text.x = element_text(size = 11, color = "black"), 
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 13, color = "black", face = "bold"),
          axis.title.y = element_text(size = 13, color = "black", face = "bold"),
          legend.title = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 13, color = "black")) +
    theme(plot.margin = unit(c(0.2,0.2,0.2,0), 'lines'))
  ggarrange(p1, p2, p3, p4, common.legend = TRUE, legend = "top", nrow = 1, 
            ncol = 4, widths = c(0.7, 1, 0.45, 0.45))
}

# Line plots of Phi values for all covariates, with confidence bands
plot_Phi_line_cis <- function(res, model, age_categs, racethnic_categs,
                              educ_categs, smoker_categs, physactive_categs) {
  if (model == "wsOFMM") {
    Phi_dfs <- convert_to_probs(res$analysis_adj$xi_med_adj)
    lbs <- apply(res$analysis_adj$xi_red_adj, c(2, 3), 
                 function(x) quantile(x, 0.025))
    ubs <- apply(res$analysis_adj$xi_red_adj, c(2, 3), 
                 function(x) quantile(x, 0.925))
    Phi_lbs <- convert_to_probs(lbs)
    Phi_ubs <- convert_to_probs(ubs)
  } else {
    Phi_dfs <- convert_to_probs(res$analysis$xi_med)
  }
  
  # Create dataframes
  # Age
  Phi_age_est <- Phi_dfs$Phi_age %>% 
    pivot_longer(cols = -Class, names_to = "Age", values_to = "Phi")
  Phi_age_lbs <- Phi_lbs$Phi_age %>% 
    pivot_longer(cols = -Class, names_to = "Age", values_to = "Phi_lbs")
  Phi_age_ubs <- Phi_ubs$Phi_age %>% 
    pivot_longer(cols = -Class, names_to = "Age", values_to = "Phi_ubs")
  Phi_age <- Phi_age_est %>%
    left_join(Phi_age_lbs, by = c("Class", "Age")) %>%
    left_join(Phi_age_ubs, by = c("Class", "Age")) %>%
    mutate(Age = factor(Age, levels = c("Age20", "Age40", "Age60"), 
                        labels = age_categs))
  # Race/Ethnicity
  Phi_racethnic_est <- Phi_dfs$Phi_racethnic %>% 
    pivot_longer(cols = -Class, names_to = "Race", values_to = "Phi")
  Phi_racethnic_lbs <- Phi_lbs$Phi_racethnic %>% 
    pivot_longer(cols = -Class, names_to = "Race", values_to = "Phi_lbs")
  Phi_racethnic_ubs <- Phi_ubs$Phi_racethnic %>% 
    pivot_longer(cols = -Class, names_to = "Race", values_to = "Phi_ubs")
  Phi_racethnic <- Phi_racethnic_est %>%
    left_join(Phi_racethnic_lbs, by = c("Class", "Race")) %>%
    left_join(Phi_racethnic_ubs, by = c("Class", "Race")) %>%
    mutate(Race = factor(Race, levels = c("NH_White", "NH_Black", "NH_Asian", 
                                         "Hispanic", "Other"),
                        labels = racethnic_categs))
  # Smoking
  Phi_smoker_est <- Phi_dfs$Phi_smoker %>% 
    pivot_longer(cols = -Class, names_to = "Smoker", values_to = "Phi")
  Phi_smoker_lbs <- Phi_lbs$Phi_smoker %>% 
    pivot_longer(cols = -Class, names_to = "Smoker", values_to = "Phi_lbs")
  Phi_smoker_ubs <- Phi_ubs$Phi_smoker %>% 
    pivot_longer(cols = -Class, names_to = "Smoker", values_to = "Phi_ubs")
  Phi_smoker <- Phi_smoker_est %>%
    left_join(Phi_smoker_lbs, by = c("Class", "Smoker")) %>%
    left_join(Phi_smoker_ubs, by = c("Class", "Smoker")) %>%
    mutate(Smoker = factor(Smoker, levels = c("Non_smoker", "Smoker"),
                        labels = smoker_categs))
  # Physical Activity
  Phi_phys_est <- Phi_dfs$Phi_phys %>% 
    pivot_longer(cols = -Class, names_to = "Physactive", values_to = "Phi")
  Phi_phys_lbs <- Phi_lbs$Phi_phys %>% 
    pivot_longer(cols = -Class, names_to = "Physactive", values_to = "Phi_lbs")
  Phi_phys_ubs <- Phi_ubs$Phi_phys %>% 
    pivot_longer(cols = -Class, names_to = "Physactive", values_to = "Phi_ubs")
  Phi_phys <- Phi_phys_est %>%
    left_join(Phi_phys_lbs, by = c("Class", "Physactive")) %>%
    left_join(Phi_phys_ubs, by = c("Class", "Physactive")) %>%
    mutate(Physactive = factor(Physactive, levels = c("Inactive", "Active"),
                        labels = physactive_categs))
  
  # Age
  p1 <- Phi_age %>% 
    ggplot(aes(x = Age, 
               y = Phi, group = Class, col = Class)) + 
    theme_bw() + 
    scale_color_brewer(palette = "Set2", 
                       labels = c("Multicultural", "Healthy American", "Western", 
                                  "Restricted Vegetarian", "Restricted American")) +
    labs(col = "Dietary Pattern") +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +
    geom_ribbon(aes(ymin = Phi_lbs, ymax = Phi_ubs, fill = Class, colour = NA), 
                alpha = 0.15, show.legend = FALSE) + 
    scale_fill_brewer(palette = "Set2", 
                      labels = c("Multicultural", "Healthy American", "Western", 
                                 "Restricted Vegetarian", "Restricted American")) +
    ylab("Probability of Hypertension") + xlab("Age") +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(size = 11, color = "black"), 
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 13, color = "black", face = "bold"),
          axis.title.y = element_text(size = 13, color = "black", face = "bold"),
          legend.title = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 13, color = "black")) +
    theme(plot.margin = unit(c(0.2,0,0.2,0.2), 'lines'))
  
  # Race/Ethnicity
  p2 <- Phi_racethnic %>%
    ggplot(aes(x = Race, 
               y = Phi, group = Class, col = Class)) + 
    theme_bw() + 
    scale_color_brewer(palette = "Set2", 
                       labels = c("Multicultural", "Healthy American", "Western", 
                                  "Restricted Vegetarian", "Restricted American")) +
    labs(col = "Dietary Pattern") +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +    
    geom_ribbon(aes(ymin = Phi_lbs, ymax = Phi_ubs, fill = Class, colour = NA), 
                alpha = 0.15, show.legend = FALSE) + 
    ylab("") + xlab("Race/Ethnicity") +
    scale_fill_brewer(palette = "Set2", 
                      labels = c("Multicultural", "Healthy American", "Western", 
                                 "Restricted Vegetarian", "Restricted American")) +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(size = 11, color = "black"), 
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 13, color = "black", face = "bold"),
          axis.title.y = element_text(size = 13, color = "black", face = "bold"),
          legend.title = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 13, color = "black")) +
    theme(plot.margin = unit(c(0.2,0,0.2,0.2), 'lines'))
  # Smoking Status
  p3 <- Phi_smoker %>% 
    ggplot(aes(x = Smoker, 
               y = Phi, group = Class, col = Class)) + 
    theme_bw() + 
    scale_color_brewer(palette = "Set2", 
                       labels = c("Multicultural", "Healthy American", "Western", 
                                  "Restricted Vegetarian", "Restricted American")) +
    labs(col = "Dietary Pattern") +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +    
    ylab("") + xlab("Smoking Status") +
    geom_ribbon(aes(ymin = Phi_lbs, ymax = Phi_ubs, fill = Class, colour = NA), 
                alpha = 0.15, show.legend = FALSE) + 
    scale_fill_brewer(palette = "Set2", 
                      labels = c("Multicultural", "Healthy American", "Western", 
                                 "Restricted Vegetarian", "Restricted American")) +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(size = 11, color = "black"), 
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 13, color = "black", face = "bold"),
          axis.title.y = element_text(size = 13, color = "black", face = "bold"),
          legend.title = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 13, color = "black")) +
    theme(plot.margin = unit(c(0.2,0,0.2,0.2), 'lines'))
  # Physical Activity
  p4 <- Phi_phys %>% 
    ggplot(aes(x = Physactive, 
               y = Phi, group = Class, col = Class)) + 
    theme_bw() + 
    scale_color_brewer(palette = "Set2", 
                       labels = c("Multicultural", "Healthy American", "Western", 
                                  "Restricted Vegetarian", "Restricted American")) +
    labs(col = "Dietary Pattern") +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +    
    ylab("") + xlab("Physical Activity") + 
    geom_ribbon(aes(ymin = Phi_lbs, ymax = Phi_ubs, fill = Class, colour = NA), 
                alpha = 0.15, show.legend = FALSE) + 
    scale_fill_brewer(palette = "Set2", 
                      labels = c("Multicultural", "Healthy American", "Western", 
                                 "Restricted Vegetarian", "Restricted American")) +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(size = 11, color = "black"), 
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.x = element_text(size = 13, color = "black", face = "bold"),
          axis.title.y = element_text(size = 13, color = "black", face = "bold"),
          legend.title = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 13, color = "black")) +
    theme(plot.margin = unit(c(0.2,0,0.2,0.2), 'lines'))
  ggarrange(p1, p2, p3, p4, common.legend = TRUE, legend = "top", nrow = 1, 
            ncol = 4, widths = c(0.7, 1, 0.45, 0.45))
  
  
} 


#### Plot xi boxplots ~ age_cat + racethnic + educ + smoker + Phys_Active
plot_xi_boxplots <- function(res, model, age_categs, racethnic_categs,
                             educ_categs, smoker_categs, physactive_categs) {
  if (model == "wsOFMM") {
    xi_dims <- dim(res$analysis_adj$xi_red_adj)
    K <- xi_dims[2]
    xi_red <- as.data.frame(t(matrix(res$analysis_adj$xi_red_adj, 
                                     nrow = xi_dims[1], 
                                     ncol = K*xi_dims[3], byrow = FALSE)))
  } else {
    xi_dims <- dim(res$analysis$xi_red)
    K <- xi_dims[2]
    xi_red <- as.data.frame(t(matrix(res$analysis$xi_red, 
                                     nrow = xi_dims[1], 
                                     ncol = K*xi_dims[3], byrow = FALSE)))
  }
  # Class variable: 1-K, 1-K, 1-K...
  xi_red$Class <- as.character(factor(c(rep(1:K, times = xi_dims[3]))))
  # Covariate level: RefxK, [40,60)xK, >=60xK... 
  xi_red$Covariate_Level <- rep(c("Ref", age_categs[-1], racethnic_categs[-1], 
                                  # educ_categs[-1], 
                                  smoker_categs[-1],
                                  physactive_categs[-1]), each = K)
  xi_red$Covariate <- c(rep("Ref", K),
                        rep("Age", K*(length(age_categs) - 1)), 
                        rep("Race/Ethnicity", K*(length(racethnic_categs) - 1)), 
                        # rep("Education", K*(length(educ_categs) - 1)), 
                        rep("Smoking", K*(length(smoker_categs) - 1)),
                        rep("Physical Activity", K*(length(physactive_categs) - 1)))
  xi_red_plot <- xi_red %>% 
    pivot_longer(cols = -c("Class", "Covariate", "Covariate_Level"), 
                 names_to = "iter", values_to = "value")
  p3_ref <- xi_red_plot %>% filter(Covariate == "Ref") %>%
    ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
    theme_bw() + scale_fill_brewer(palette="Set2") + 
    geom_boxplot() + xlab("Reference") + ylab(expression(xi))
  p3_age <- xi_red_plot %>% filter(Covariate == "Age") %>%
    ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
    theme_bw() + scale_fill_brewer(palette="Set2") + 
    geom_boxplot() + xlab("Age") + ylab(expression(xi))
  p3_race <- xi_red_plot %>% filter(Covariate == "Race/Ethnicity") %>%
    ggplot(aes(x = factor(Covariate_Level, 
                          levels = c("NH Black", "NH Asian", "Hispanic/Latino", 
                                     "Other/Mixed")), y = value, fill = Class)) +
    theme_bw() + scale_fill_brewer(palette="Set2") + 
    geom_boxplot() + xlab("Race/Ethnicity") + ylab(expression(xi))
  p3_smoke <- xi_red_plot %>% filter(Covariate == "Smoking") %>%
    ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
    theme_bw() + scale_fill_brewer(palette="Set2") + 
    geom_boxplot() + xlab("Smoking") + ylab(expression(xi))
  p3_physactive <- xi_red_plot %>% filter(Covariate == "Physical Activity") %>%
    ggplot(aes(x = Covariate_Level, y = value, fill = Class)) +
    theme_bw() + scale_fill_brewer(palette="Set2") + 
    geom_boxplot() + xlab("Physical Activity") + ylab(expression(xi))
  ggarrange(p3_ref, p3_age, p3_race, p3_smoke, p3_physactive, 
            common.legend = TRUE, legend = "right", ncol = 3, nrow = 2)
  
}

#===================== Plot pi boxplots ========================================
plot_pi_boxplots <- function(res, model) {
  if (model == "wsOFMM") {
    pi_red <- as.data.frame(res$analysis_adj$pi_red_adj)
  } else {
    pi_red <- as.data.frame(res$analysis$pi_red)
  }
  colnames(pi_red) <- paste0("pi_", 1:(dim(pi_red)[2]))
  pi_red_plot <- pi_red %>% pivot_longer(cols = everything(), names_to = "pi_comp", 
                                         values_to = "value")
  pi_red_plot %>% ggplot(aes(x = pi_comp, y = value, fill = pi_comp)) + 
    theme_bw() + scale_fill_brewer(palette="Set2") + 
    geom_boxplot() + 
    ggtitle(as.expression(bquote("Parameter estimation for "~pi~" "))) +
    xlab("Latent Class") + ylab(as.expression(bquote(~pi~"Value"))) + ylim(0,0.5)
}

#====================== Get regression coefs ===================================
get_ci <- function(x) {
  quantiles <- format(round(quantile(x, c(0.025, 0.975)), 2), nsmall = 2)
  # quantiles <- round(quantile(x, c(0.05, 0.95)), 3)
  ci <- paste0("(", quantiles[1], ", ", quantiles[2], ")")
  return(ci)
} 
get_pep <- function(x, x_med) {
  if (x_med > 0) {
    pep <- mean(x < 0)
  } else {
    pep <- mean(x > 0)
  }
  pep <- format(round(pep, 2), nsmall = 2)
  if (pep < 0.001) {
    pep <- "<0.001"
  }
  return(pep)
}

get_prob_pos <- function(x) {
  prob_pos <- format(round(mean(x > 0), 2), nsmall = 2)
  if (prob_pos < 0.01) {
    prob_pos <- "<0.01"
  }
  return(prob_pos)
}
# Converts from combination coding to reference cell coding
# Input: xi_comb: Matrix of xi parameter estimates. Kxq
# Output: List containing the following objects
#   Phi_age: Dataframe of hypertension probs for age categories. Kx3
#   Phi_age: Dataframe of hypertension probs for race/ethnicity categories. Kx5
#   Phi_smoker: Dataframe of hypertension probs for smoking categories. Kx2
#   Phi_phys: Dataframe of hypertension probs for physical activity categories. Kx2
convert_to_ref <- function(xi_med, xi_red, age_categs, racethnic_categs,
                           smoker_categs, physactive_categs, format = "latex") {
  # Assumes the following column order for xi_comb: 
  # Intercept + Age40 + Age60 + NH_Black + NH_Asian + Hispanic + Other + Smoker + Inactive
  K <- nrow(xi_med)
  q <- ncol(xi_med)
  beta <- as.data.frame(matrix(NA, nrow = 45, ncol = 4))
  # colnames(beta) <- c("Covariate", "Median Estimate", "95% Credible Interval", 
  #                     "Posterior Error Probability")
  colnames(beta) <- c("Covariate", "Median", "95% Credible Interval", 
                      "P($\\xi$ > 0)")
  cov_categs <- c(age_categs[-1], racethnic_categs[-1], smoker_categs[-1], 
                  physactive_categs[-1])
  beta[, 1] <- c("Reference", "Class2", "Class3", "Class4", "Class5",
                 cov_categs,
                 sapply(cov_categs, function(x) paste0(x, ":Class", 2:5)))
  # beta[1, -1] <- c(xi_med[1, 1], get_ci(xi_red[, 1, 1]),
  #                  get_pep(xi_red[, 1, 1], xi_med[1, 1]))
  beta[1, -1] <- c(xi_med[1, 1], get_ci(xi_red[, 1, 1]),
                   get_prob_pos(xi_red[, 1, 1]))
  for (i in 2:K) {
    # beta[i, -1] <- c(xi_med[i, 1] - xi_med[1, 1], 
    #                  get_ci(xi_red[, i, 1] - xi_med[1, 1]),
    #                  get_pep(xi_red[, i, 1] - xi_med[1, 1], xi_med[i, 1] - xi_med[1, 1]))
    beta[i, -1] <- c(median(xi_red[, i, 1] - xi_red[, 1, 1]), 
                     get_ci(xi_red[, i, 1] - xi_red[, 1, 1]),
                     get_prob_pos(xi_red[, i, 1] - xi_red[, 1, 1]))
  }
  for (i in 2:q) {
    # beta[K + (i-1), -1] <- c(xi_med[1, i], get_ci(xi_red[, 1, i]),
    #                          get_pep(xi_red[, 1, i], xi_med[1, i]))
    beta[K + (i-1), -1] <- c(xi_med[1, i], get_ci(xi_red[, 1, i]),
                             get_prob_pos(xi_red[, 1, i]))
  }
  for (i in 2:q) {
    for (j in 2:K) {
      # beta[q + (i-1)*(K-1) + (j-1), -1] <- c(xi_med[j, i] - xi_med[1, i],
      #                                        get_ci(xi_red[, j, i] - xi_med[1, i]),
      #                                        get_pep(xi_red[, j, i] - xi_med[1, i],
      #                                                xi_med[j, i] - xi_med[1, i]))
      beta[q + (i-1)*(K-1) + (j-1), -1] <- c(median(xi_red[, j, i] - xi_red[, 1, i]),
                                             get_ci(xi_red[, j, i] - xi_red[, 1, i]),
                                             get_prob_pos(xi_red[, j, i] - xi_red[, 1, i]))
    }
  }
  beta[, 2] <- round(as.numeric(beta[, 2]), 2)
  # Get significant entries
  pep_signif <- numeric(nrow(beta))
  # pep_signif <- ifelse(beta[, 4] == "<0.001", 1, 
  #                      ifelse(as.numeric(beta[, 4]) < 0.05, 1, 0))
  beta %>% 
    kable(format = format, booktabs = TRUE) %>%
    kable_styling() %>%
    column_spec(2:4, bold = ifelse(pep_signif == 1, TRUE, FALSE))
}

#================ Plot colored dendrogram ======================================
# from https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
dendro_data_k <- function(hc, k) {
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  hcdata
}

set_labels_params <- function(nbLabels, direction = c("tb", "bt", "lr", "rl"),
                              fan = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}

plot_ggdendro <- function(hcdata, direction = c("lr", "rl", "tb", "bt"),
                          fan = FALSE, scale.color = NULL, branch.size = 1,
                          label.size  = 3, nudge.label = 0.01, expand.y = 0.1) {
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(ggdendro::segment(hcdata)$y, n = 5)
  ymax      <- max(ggdendro::segment(hcdata)$y)
  ## branches
  p <- ggplot() +
    geom_segment(data         =  ggdendro::segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "square",
                 show.legend  =  FALSE,
                 size         =  branch.size)
  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }
  # labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle
  # p <- p +
  #   geom_text(data        =  label(hcdata),
  #             aes(x       =  x,
  #                 y       =  y,
  #                 label   =  label,
  #                 colour  =  factor(clust),
  #                 angle   =  angle),
  #             vjust       =  labelParams$vjust,
  #             hjust       =  labelParams$hjust,
  #             nudge_y     =  ymax * nudge.label,
  #             size        =  label.size,
  #             show.legend =  FALSE)
  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }
  # ylim <- -round(ymax * expand.y, 1)
  # p    <- p + expand_limits(y = ylim)
  p
}
