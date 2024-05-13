# Load results
data("run_nhanes_swolca_results")
res <- run_nhanes_swolca_results

# Test errors for plot_pattern_profiles
expect_error(plot_pattern_profiles(res = res, item_labels = c("test", "label")),
             "length of item_labels must equal the number of exposure items, J = 28",
             fixed = TRUE)
expect_error(plot_pattern_profiles(res = res, class_labels = c("test", "label")),
             "length of class_labels must equal the number of latent classes, K = 5",
             fixed = TRUE)
expect_error(plot_pattern_profiles(res = res, categ_labels = c("test", "label")),
             "length of categ_labels must equal the number of exposure categories, R = 4",
             fixed = TRUE)

# Test errors for plot_pattern_probs
expect_error(plot_pattern_probs(res = res, item_labels = c("test", "label")),
             "length of item_labels must equal the number of exposure items, J = 28",
             fixed = TRUE)
expect_error(plot_pattern_probs(res = res, class_labels = c("test", "label")),
             "length of class_labels must equal the number of latent classes, K = 5",
             fixed = TRUE)
expect_error(plot_pattern_probs(res = res, categ_labels = c("test", "label")),
             "length of categ_labels must equal the number of exposure categories, R = 4",
             fixed = TRUE)

# Test errors for plot_class_dist
expect_error(plot_class_dist(res = res, class_labels = c("test", "label")),
             "length of class_labels must equal the number of latent classes, K = 5",
             fixed = TRUE)

# Test errors for plot_outcome_probs
expect_error(plot_outcome_probs(res = res, cov_name = "racethnic", ci_level = 2),
             "ci_level must be between 0 and 1",
             fixed = TRUE)
expect_error(plot_outcome_probs(res = res, cov_name = "racethnic", 
                           cov_labels = c("test", "label")),
             "length of cov_labels for covariate racethnic must equal the number of categories: 5",
             fixed = TRUE)
expect_error(plot_outcome_probs(res = res, cov_name = "racethnic", 
                           class_labels = c("test", "label")),
             "length of class_labels must equal the number of latent classes, K = 5",
             fixed = TRUE)
expect_error(plot_outcome_probs(res = res, cov_name = "test"),
             "all variables in cov_name must be specified in glm_form",
             fixed = TRUE)

# Test errors for plot_outcome_probs
expect_error(plot_outcome_probs(res = res, cov_name = c("age_cat", "racethnic"),
                                cov_labels = c("[20,40)", "[40,60)", ">=60")),
             "cov_name is a vector of length 2, while cov_labels is a list of length 1. The two must be of the same length.",
             fixed = TRUE)
# plot_outcome_probs(res = res_wolca_svyglm, cov_name = c("age_cat", "racethnic"),
#                    cov_labels = list(age_cat_categs, racethnic_categs), 
#                    class_labels = class_labels, x_title = "Age Group", 
#                    ci_level = NULL, add_lines = TRUE)