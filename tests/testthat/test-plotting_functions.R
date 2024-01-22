# Load results
data("run_nhanes_swolca_results")
res <- run_nhanes_swolca_results

# Test errors for plot_theta_modes
expect_error(plot_theta_modes(res = res, item_labels = c("test", "label")),
             "length of item_labels must equal the number of exposure items, J = 28",
             fixed = TRUE)
expect_error(plot_theta_modes(res = res, class_labels = c("test", "label")),
             "length of class_labels must equal the number of latent classes, K = 5",
             fixed = TRUE)
expect_error(plot_theta_modes(res = res, categ_labels = c("test", "label")),
             "length of categ_labels must equal the number of exposure categories, R = 4",
             fixed = TRUE)

# Test errors for plot_theta_probs
expect_error(plot_theta_probs(res = res, item_labels = c("test", "label")),
             "length of item_labels must equal the number of exposure items, J = 28",
             fixed = TRUE)
expect_error(plot_theta_probs(res = res, class_labels = c("test", "label")),
             "length of class_labels must equal the number of latent classes, K = 5",
             fixed = TRUE)
expect_error(plot_theta_probs(res = res, categ_labels = c("test", "label")),
             "length of categ_labels must equal the number of exposure categories, R = 4",
             fixed = TRUE)

# Test errors for plot_pi_boxplots
expect_error(plot_pi_boxplots(res = res, class_labels = c("test", "label")),
             "length of class_labels must equal the number of latent classes, K = 5",
             fixed = TRUE)

# Test errors for plot_Phi_line
expect_error(plot_Phi_line(res = res, cov_name = "racethnic", ci_level = 2),
             "ci_level must be between 0 and 1",
             fixed = TRUE)
expect_error(plot_Phi_line(res = res, cov_name = "racethnic", 
                           cov_labels = c("test", "label")),
             "length of cov_labels must equal the number of covariate categories: 5",
             fixed = TRUE)
expect_error(plot_Phi_line(res = res, cov_name = "racethnic", 
                           class_labels = c("test", "label")),
             "length of class_labels must equal the number of latent classes, K = 5",
             fixed = TRUE)
expect_error(plot_Phi_line(res = res, cov_name = "test"),
             "cov_name must be one of the variables specified in glm_form",
             fixed = TRUE)

# Test errors for plot_outcome_probs
expect_error(plot_outcome_probs(res = res, cov_name = c("age_cat", "racethnic"),
                                cov_labels = c("[20,40)", "[40,60)", ">=60")),
             "cov_name has length 2 but cov_labels has length 1. The two must match.",
             fixed = TRUE)
# plot_outcome_probs(res = res_wolca_svyglm, cov_name = c("age_cat", "racethnic"),
#                    cov_labels = list(age_cat_categs, racethnic_categs), 
#                    class_labels = class_labels, x_title = "Age Group", 
#                    ci_level = NULL, add_lines = TRUE)