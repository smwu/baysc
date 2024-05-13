#' Plot theta modal exposure categories for each latent class
#'
#' @description
#' `plot_pattern_profiles` plots a heatmap of the latent class patterns, where the 
#' patterns are defined by the category with the highest probability 
#' (i.e., the model category) for each exposure item.
#'
#' @param res An object of class `"swolca"` or `"wolca"`, resulting from a call 
#' to [swolca()], [swolca_var_adjust()], [wolca()], or [wolca_var_adjust()]. 
#' @param item_labels String vector of names for the exposure items. Jx1. If
#' `NULL` (default), numbers from 1 to J are used.
#' @param item_title String specifying the title for the exposure items. 
#' Default is `"Item"`.
#' @param categ_labels String vector of names for the item categories. Rx1. If
#' `NULL` (default), numbers from 1 to R are used.
#' @param categ_title String specifying the title for the item categories. 
#' Default is `"Consumption Level"`.
#' @param class_labels String vector of names for the latent classes. Kx1. 
#' If `NULL` (default), numbers from 1 to K are used, where K is the final 
#' determined number of latent classes.
#' @param class_title String specifying the title for the latent classes. 
#' Default is `"Dietary Pattern"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying a heatmap of the latent class patterns.
#' 
#' @seealso [plot_pattern_probs()] [plot_class_dist()] [plot_outcome_probs()]
#' 
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_brewer labs 
#' theme_classic scale_x_discrete theme element_text
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples 
#' data(run_nhanes_swolca_results)
#' 
#' # Default labels
#' plot_pattern_profiles(res = run_nhanes_swolca_results)
#' 
#' # Specifying labels
#' item_labels <- c("Citrus/Melon/Berries", "Other Fruits", "Fruit Juice", 
#'                   "Dark Green Vegs", "Tomatoes", "Oth Red/Orange Vegs",
#'                   "Potatoes", "Other Starchy Vegs", "Other Vegetables",
#'                   "Whole Grains", "Refined Grains", "Meat", "Cured Meats",
#'                   "Organ Meat", "Poultry", "Seafood (High n-3)", "Seafood (Low n-3)",
#'                   "Eggs", "Soybean Products", "Nuts and Seeds", "Legumes (Protein)",
#'                   "Milk", "Yogurt", "Cheese", "Oils", "Solid Fats", "Added Sugar",
#'                   "Alcoholic Drinks")
#' categ_labels <- c("None", "Low", "Med", "High")
#' class_labels <- paste0("Class ", 1:5)
#' plot_pattern_profiles(res = run_nhanes_swolca_results, item_labels = item_labels,
#'                  categ_labels = categ_labels, class_labels = class_labels)
#'
plot_pattern_profiles <- function(res, item_labels = NULL, item_title = "Item",
                               categ_labels = NULL, 
                               categ_title = "Consumption Level",
                               class_labels = NULL, 
                               class_title = "Dietary Pattern", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  # Obtain theta estimates
  if (!is.null(res$estimates_adjust)) {
    # Adjusted estimates
    est_item_probs <- res$estimates_adjust$theta_med
  } else {
    # Unadjusted estimates 
    est_item_probs <- res$estimates$theta_med
  }
  
  # Obtain theta modes as the category with highest probability for each item
  mode_item_probs <- as.data.frame(apply(est_item_probs, c(1, 2), which.max))
  
  # Define item, latent class, and category names if not specified
  if (is.null(item_labels)) {
    item_labels <- 1:res$data_vars$J
  } else if (length(item_labels) != res$data_vars$J) {
    stop(paste0("length of item_labels must equal the number of exposure items, J = ",
                res$data_vars$J))
  }
  K <- dim(mode_item_probs)[2]
  if (is.null(class_labels)) {
    class_labels <- 1:K
  } else if (length(class_labels) != K) {
    stop(paste0("length of class_labels must equal the number of latent classes, K = ", K))
  }
  if (is.null(categ_labels)) {
    categ_labels <- 1:res$data_vars$R
  } else if (length(categ_labels) != res$data_vars$R) {
    stop(paste0("length of categ_labels must equal the number of exposure categories, R = ", 
                res$data_vars$R))
  }
  rownames(mode_item_probs) <- item_labels
  colnames(mode_item_probs) <- class_labels
  mode_item_probs$Item <- rownames(mode_item_probs)
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Class <- Item <- Level <- NULL
  
  # Create plot
  mode_plot <- mode_item_probs %>% 
    tidyr::gather("Class", "Level", -Item) 
  mode_plot %>% ggplot2::ggplot(ggplot2::aes(x=Class, 
                                             y=factor(Item, 
                                                      levels = rev(item_labels)), 
                                             fill = factor(Level))) + 
    ggplot2::geom_tile(color="black", linewidth = 0.3) + 
    ggplot2::scale_fill_brewer(type="seq", palette="RdYlBu", direction = -1,
                               name = categ_title, labels = categ_labels) +
    ggplot2::labs(x = class_title, y = item_title) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete() + 
    ggplot2::theme(text = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 12, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.title = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.text = ggplot2::element_text(size = 11, color = "black"),
                   legend.position = "top")
}



#' Plot probabilities of exposure categories, theta, for each latent class
#'
#' @description
#' `plot_pattern_probs` plots a grouped barplot of the probability of the exposure
#' categories, for each exposure item and each latent class.
#' 
#' @inheritParams plot_pattern_profiles
#' @param y_title String specifying the title for the y-axis. Default is 
#' `"Consumption Level Probability"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying a grouped barplot of the probability of 
#' the exposure categories, for each exposure item and each latent class
#' 
#' @seealso [plot_pattern_profiles()] [plot_class_dist()] [plot_outcome_probs()]
#' 
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap scale_fill_brewer labs 
#' theme_bw theme element_text element_blank element_rect
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples 
#' data(run_nhanes_swolca_results)
#' 
#' # Default labels
#' plot_pattern_probs(res = run_nhanes_swolca_results)
#' 
#' # Specifying labels
#' item_labels <- c("Citrus/Melon/Berries", "Other Fruits", "Fruit Juice", 
#'                   "Dark Green Vegs", "Tomatoes", "Oth Red/Orange Vegs",
#'                   "Potatoes", "Other Starchy Vegs", "Other Vegetables",
#'                   "Whole Grains", "Refined Grains", "Meat", "Cured Meats",
#'                   "Organ Meat", "Poultry", "Seafood (High n-3)", "Seafood (Low n-3)",
#'                   "Eggs", "Soybean Products", "Nuts and Seeds", "Legumes (Protein)",
#'                   "Milk", "Yogurt", "Cheese", "Oils", "Solid Fats", "Added Sugar",
#'                   "Alcoholic Drinks")
#' categ_labels <- c("None", "Low", "Med", "High")
#' class_labels <- paste0("C", 1:5)
#' plot_pattern_probs(res = run_nhanes_swolca_results, item_labels = item_labels,
#'                  categ_labels = categ_labels, class_labels = class_labels)
#' 
plot_pattern_probs <- function(res, item_labels = NULL, categ_labels = NULL, 
                               categ_title = "Consumption Level",
                               class_labels = NULL,
                               class_title = "Dietary Pattern", 
                               y_title = "Consumption Level Probability", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  # Obtain theta estimates
  if (!is.null(res$estimates_adjust)) {
    # Adjusted estimates 
    est_item_probs <- res$estimates_adjust$theta_med
  } else {
    # Unadjusted estimates
    est_item_probs <- res$estimates$theta_med
  }
  # Get number of latent classes
  K <- dim(est_item_probs)[2]
  
  # Define item, latent class, and category names if not specified
  if (is.null(item_labels)) {
    item_labels <- 1:res$data_vars$J
  } else if (length(item_labels) != res$data_vars$J) {
    stop(paste0("length of item_labels must equal the number of exposure items, J = ",
                res$data_vars$J))
  }
  if (is.null(class_labels)) {
    class_labels <- 1:K
  } else if (length(class_labels) != K) {
    stop(paste0("length of class_labels must equal the number of latent classes, K = ", K))
  }
  if (is.null(categ_labels)) {
    categ_labels <- 1:res$data_vars$R
  } else if (length(categ_labels) != res$data_vars$R) {
    stop(paste0("length of categ_labels must equal the number of exposure categories, R = ", 
                res$data_vars$R))
  }
  
  dimnames(est_item_probs)[[1]] <- item_labels
  dimnames(est_item_probs)[[2]] <- class_labels
  
  # Convert to dataframe with each row corresponding to a value in the array
  # Use base R instead of reshape2 to reduce number of package imports
      # theta_plot <- reshape2::melt(est_item_probs, level = 2)
  theta_plot <- data.frame(expand.grid(lapply(dim(est_item_probs), seq_len)), 
                           value = as.vector(est_item_probs))
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Item <- Class <- Probability <- Level <- NULL
  
  # Create plot
  colnames(theta_plot) <- c("Item", "Class", "Level", "Probability")
  theta_plot %>%
    ggplot2::ggplot(ggplot2::aes(x = factor(Class, labels = class_labels), 
                                 y = Probability,
                                 fill = factor(Level))) + 
    ggplot2::geom_bar(stat = "identity", position = "stack") + 
    ggplot2::facet_wrap(factor(Item, labels = item_labels) ~ ., nrow = 4) + 
    ggplot2::scale_fill_brewer(type="seq", palette="RdYlBu", direction = -1,
                               name = categ_title, labels = categ_labels) +
    ggplot2::theme_bw() + 
    ggplot2::labs(x = class_title, y = y_title) + 
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 11, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.title = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   legend.text = ggplot2::element_text(size = 11, color = "black"),
                   legend.position = "top",
                   strip.text = ggplot2::element_text(size = 9),
                   strip.background = ggplot2::element_rect(fill = "gray90"))
}


#' Plot distribution of classes in the population across posterior sample iterations 
#' 
#' @description
#' `plot_class_dist` plots a boxplot of the class membership probabilities, 
#' \eqn{\pi}, in the posterior samples, for each latent class.
#' 
#' @inheritParams plot_pattern_profiles
#' @param y_title String specifying the title for the y-axis. Default is 
#' `"Class Membership Probability"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying a boxplot of the distribution of class 
#' membership probabilities in the posterior samples, for each latent class.
#' 
#' @seealso [plot_pattern_probs()] [plot_pattern_profiles()] [plot_outcome_probs()]
#' 
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_fill_brewer labs 
#' theme_bw theme element_text
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
#'
#' @examples 
#' data(run_nhanes_swolca_results)
#' 
#' # Default labels
#' plot_class_dist(res = run_nhanes_swolca_results)
#' 
#' # Specifying labels
#' class_labels <- paste0("Class ", 1:5)
#' plot_class_dist(res = run_nhanes_swolca_results, class_labels = class_labels)
plot_class_dist <- function(res, class_labels = NULL, 
                             class_title = "Dietary Pattern",
                             y_title = "Class Membership Probability", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  # Obtain pi estimates
  if (!is.null(res$estimates)) {
    # Adjusted estimates 
    pi_red <- as.data.frame(res$estimates_adjust$pi_red)
  } else {
    # Unadjusted estimates 
    pi_red <- as.data.frame(res$estimates$pi_red)
  }
  # Get number of latent classes
  K <- dim(pi_red)[2]
  
  # Set class labels to 1:K if not provided
  if (is.null(class_labels)) {
    class_labels <- 1:K
  } else if (length(class_labels) != K) {
    stop(paste0("length of class_labels must equal the number of latent classes, K = ", K))
  }
  colnames(pi_red) <- class_labels
  
  # Convert to longer format for plotting
  pi_red_plot <- pi_red %>% tidyr::pivot_longer(cols = tidyselect::everything(), 
                                                names_to = "pi_comp", 
                                                values_to = "value")
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  pi_comp <- value <- NULL
  
  # Plot pi boxplots
  pi_red_plot %>% ggplot2::ggplot(ggplot2::aes(x = pi_comp, 
                                               y = value)) + 
    ggplot2::theme_bw() + ggplot2::scale_fill_brewer(palette="Set2") + 
    ggplot2::geom_boxplot() + 
    ggplot2::labs(x = class_title, y = y_title) + 
    ggplot2::theme(text = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 12, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 13, color = "black", face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 13, color = "black", face = "bold"))
}

#' Plot regression coefficients
#' 
#' @description
#' `plot_regr_coefs` plots point estimates and error bars for all the regression 
#' coefficients in reference cell coding format, using the output from the
#' [get_regr_coefs()] function. Estimates corresponding to the same latent 
#' class are displayed with the same color.  
#' 
#' @inheritParams plot_outcome_probs
#' @param regr_coefs Dataframe output from [get_regr_coefs()] containing at 
#' least the following columns:
#' \describe{
#'   \item{\code{Covariate}}{Names of covariate terms.}
#'   \item{\code{Estimate}}{Regression estimates in reference cell coding format.}
#'   \item{\code{LB}}{Regression estimate lower bound}
#'   \item{\code{UB}}{Regression estimate upper bound}
#' }
#' @param cov_labels Optional string vector of labels whose order must 
#' correspond exactly with that of the `Covariate` column in `regr_coefs`. 
#' Default is `NULL` and original covariate names are used. 
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying point estimates and error bars for the 
#' probit regression coefficients in reference cell coding format. Estimates and 
#' intervals are obtained from the [get_regr_coefs()] function using the 
#' interval level specified there.
#' 
#' @seealso [get_regr_coefs()] [plot_pattern_probs()] [plot_class_dist()] 
#' [plot_pattern_profiles()] 
#' 
#' @importFrom ggplot2 ggplot aes theme_bw geom_point theme element_text 
#' geom_hline geom_errorbar scale_color_brewer
#' @importFrom magrittr %>%
#' @export
#'
#' @examples 
#' data(run_nhanes_swolca_results)
#' res <- run_nhanes_swolca_results
#' 
#' # Get table of regression coefficients
#' regr_coefs <- get_regr_coefs(res = res, ci_level = 0.95, digits = 2)
#' 
#' # Define new vector of covariate labels
#' class_dummies <- paste0("C", 2:5)
#' reps <- length(class_dummies)
#' age_dummies <- paste0("Age", c("40_60", "60"))
#' race_dummies <- paste0("RaceEth", c("NH_Black", "NH_Asian", "Hisp", "Other"))
#' cov_labels <- c("Intercept", class_dummies, age_dummies, race_dummies,
#'                 "SmokerYes", "PhysActive",
#'                 paste0(class_dummies, ":", rep(age_dummies, each = reps)),
#'                 paste0(class_dummies, ":", rep(race_dummies, each = reps)),
#'                 paste0(class_dummies, ":", rep("SmokerYes", each = reps)),
#'                 paste0(class_dummies, ":", rep("PhysActive", each = reps)))
#' 
#' # Plot regression coefficient estimates and error bars
#' plot_regr_coefs(regr_coefs = regr_coefs, res = res, cov_labels = cov_labels)
#' 
plot_regr_coefs <- function(regr_coefs, res, cov_labels = NULL, ...) { 
  # Check object class
  if (!(class(res) %in% c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  # Check covariate names
  if (regr_coefs$Covariate[2] != "c_all2") {
    stop("regr_coefs must have latent class covariates with names of the form 
    `c_allk`, where k ranges from 2 to the total number of classes")
  }
  # If no covariate labels are specified, original labels are used
  if (is.null(cov_labels)) {
    cov_labels <- levels(as.factor(regr_coefs$Covariate))
  } else {
    # Check covariate labels 
    if ((length(cov_labels) != length(regr_coefs$Covariate)) | 
        !is.character(cov_labels)) {
      stop("cov_labels must be a string vector that exactly corresponds to 
         `regr_coefs$Covariate` and is of the same length")
    }
  }
  
  # Get number of latent classes from adjusted or unadjusted estimates
  if (!is.null(res$estimates_adjust)) {
    K <- length(res$estimates_adjust$pi_med)
  } else {
    K <- length(res$estimates$pi_med)
  }
  plot_df <- regr_coefs
  plot_df$Class <- 1
  for (k in 1:K) {
    class_str <- paste0("c_all", k)
    plot_df$Class[grepl(class_str, plot_df$Covariate)] <- k
  }
  plot_df$Class <- as.factor(plot_df$Class)
  plot_df$Covariate <- factor(plot_df$Covariate, levels = plot_df$Covariate, 
                              labels = cov_labels)
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Class <- Covariate <- Estimate <- LB <- UB <- NULL
  
  plot_df %>% 
    ggplot2::ggplot(ggplot2::aes(x = Covariate, y = Estimate, col = Class)) + 
    ggplot2::theme_bw() + 
    ggplot2::geom_point() + 
    ggplot2::scale_color_brewer(palette = "Set2", labels = Class, 
                                aesthetics = c("color")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 70, vjust = 1, 
                                                       hjust = 1)) + 
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") + 
    ggplot2::geom_errorbar(ggplot2::aes(ymin = LB, ymax = UB, color = Class, 
                                        width = 0.2))
}


#' Plot conditional probability of outcome for each latent class for one or two 
#' categorical covariates. 
#' 
#' @description
#' `plot_outcome_probs` plots the conditional probability of the outcome, 
#' ranging from 0 to 1 and obtained by transforming the probit regression 
#' coefficients to the probability scale, for up to two categorical covariates. 
#' Point estimates and error bars are colored by latent class. 
#' 
#' @inheritParams plot_pattern_profiles
#' @param cov_name String specifying the covariate to plot. To plot 
#' interactions between two covariates, specify a string vector with both 
#' covariate names. Covariates must be included in `glm_form` and `V_data`.
#' @param ci_level Optional number from 0 to 1 specifying the confidence/credible 
#' interval level. Default is `0.95`, which gives a 95\% equal-tailed interval 
#' composed of the 2.5\% and 97.5\% quantiles. For `"wolca"` objects, this must 
#' match the `ci_level` parameter in the main [wolca()] function. Set to `NULL` 
#' to exclude credible intervals from the plot. 
#' @param add_lines Boolean specifying whether lines should be added to connect 
#' the points between categories. Default is `FALSE`.
#' @param cov_labels Optional list of 1 or 2 string vectors specifying the 
#' corresponding category labels for the covariate(s) of interest. Must be the 
#' same length as `cov_name` and have the appropriate number of categories. 
#' If `NULL` (default), the covariate categories from the data are used.
#' @param x_title Optional string specifying x-axix label. If `NULL` (default), 
#' the first string in `cov_name`is used.
#' @param y_title Optional string specifying the title for the y-axis. Default is 
#' `"Probability of Outcome"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying the point estimates and error bars of 
#' the conditional probability of the outcome across categories of one or two 
#' covariates of interest, colored by latent class.
#' 
#' @seealso [get_regr_coefs()] [plot_pattern_probs()] [plot_class_dist()] [plot_pattern_profiles()] 
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_brewer labs 
#' theme_bw theme element_text position_dodge facet_grid geom_errorbar
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join
#' @importFrom stats quantile
#' @importFrom magrittr %>%
#' @export
#'
#' @examples 
#' data(run_nhanes_swolca_results)
#' res <- run_nhanes_swolca_results
#' 
#' # Default labels
#' plot_outcome_probs(res = res, cov_name = "racethnic")
#' 
#' # Specify labels
#' cov_labels <- c("NH White", "NH Black", "NH Asian", "Hispanic/Latino", 
#'                 "Other/Mixed")
#' class_labels <- paste0("Class", 1:5)
#' x_title <- "Race and Ethnicity"
#' plot_outcome_probs(res = res, cov_name = "racethnic", cov_labels = cov_labels, 
#'                    class_labels = class_labels, x_title = x_title)
#' 
#' # Two covariates, remove error bars, add lines
#' age_cat_categs <- c("(20,40)", "(40,60)", ">=60")
#' racethnic_categs <- c("NH White", "NH Black", "NH Asian", "Hispanic/Latino", 
#'                       "Other/Mixed")
#' plot_outcome_probs(res = res, cov_name = c("age_cat", "racethnic"),
#'                    cov_labels = list(age_cat_categs, racethnic_categs), 
#'                    class_labels = class_labels, x_title = "Age Group", 
#'                    ci_level = NULL, add_lines = TRUE)         
#' \dontrun{                    
#' # Multiple plots for various covariates
#' educ_categs <- c("Some College", "HS/GED", "<HS")
#' smoker_categs <- c("Non-Smoker", "Smoker")
#' physactive_categs <- c("Inactive", "Active")
#' p1 <- plot_outcome_probs(res = res, cov_name = "age_cat",
#'                          cov_labels = age_cat_categs, 
#'                          class_labels = class_labels, 
#'                          x_title = "Age Group")
#' p2 <- plot_outcome_probs(res = res, cov_name = "racethnic",
#'                          cov_labels = racethnic_categs, 
#'                          class_labels = class_labels, 
#'                          x_title = "Race and Ethnicity")
#' p3 <- plot_outcome_probs(res = res, cov_name = "smoker",
#'                          cov_labels = smoker_categs, 
#'                          class_labels = class_labels, 
#'                          x_title = "Current Smoking Status")
#' p4 <- plot_outcome_probs(res = res, cov_name = "physactive",
#'                          cov_labels = physactive_categs, 
#'                          class_labels = class_labels, 
#'                          x_title = "Physical Activity")
#' ggpubr::ggarrange(p1, p2, p3, p4, common.legend = TRUE, legend = "top", 
#'                   nrow = 1, ncol = 4, widths = c(0.7, 1, 0.45, 0.45)) 
#' }                              
#'                    
plot_outcome_probs <- function(res, cov_name, ci_level = 0.95, add_lines = FALSE, 
                              cov_labels = NULL, class_labels = NULL, 
                              class_title = "Dietary Pattern", x_title = NULL, 
                              y_title = "Probability of Outcome", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "wolca"))) {
    stop("res must be an object of class `swolca` or `wolca`, resulting 
         from a call to one of these functions")
  }
  # Check ci_level
  if (!is.null(ci_level)) {
    if (!(ci_level > 0 & ci_level < 1)) {
      stop("ci_level must be between 0 and 1")
    }
    quant_lb <- (1 - ci_level) / 2
    quant_ub <- 1 - quant_lb
  } 
  
  # Check cov_name and cov_labels
  if (!(length(cov_name) %in% c(1, 2)) | !(is.character(cov_name[[1]]))) {
    stop("cov_name must be a string vector of length 1 or 2")
  } 
  if (!is.null(cov_labels)) {
    # Convert cov_labels to list with one string vector if necessary
    if (!is.list(cov_labels)) {
      cov_labels <- list(cov_labels)
    }
    if (!(length(cov_labels) %in% c(1, 2)) | !(is.character(cov_labels[[1]]))) {
      print(cov_labels)
      stop("cov_labels must be a list of length 1 or 2 composed of string vectors")
    } else if (length(cov_name) != length(cov_labels)) {
      stop(paste0("cov_name is a vector of length ", length(cov_name), 
                  ", while cov_labels is a list of length ", length(cov_labels), 
                  ". The two must be of the same length."))
    }
  }
  
  
  # Obtain xi median and lower bound and upper bound estimates 
  # Obtain xi for `wolca()`
  if (is(res, "wolca")) {
    if (is.null(res$estimates_svyglm)) {
      stop("wolca object does not have regression estimates. Please run wolca_svyglm().")
    }
    estimates <- res$estimates_svyglm
    est_xi <- estimates$xi_est
    if (!is.null(ci_level)) {
      if (ci_level != res$data_vars$ci_level) {
        stop("ci_level must match the specified ci_level in the wolca() function")
      }
      est_lb <- estimates$xi_est_lb
      est_ub <- estimates$xi_est_ub
    }
    
  # Obtain xi for `swolca()`
  } else {  
    if (!is.null(res$estimates_adjust)) {
      # Adjusted estimates
      estimates <- res$estimates_adjust
    } else {
      # Unadjusted estimates
      estimates <- res$estimates
    }
    est_xi <- estimates$xi_med
    if (!is.null(ci_level)) {
      est_lb <- apply(estimates$xi_red, c(2, 3), 
                      function(x) stats::quantile(x, quant_lb))
      est_ub <- apply(estimates$xi_red, c(2, 3),
                      function(x) stats::quantile(x, quant_ub))
    }
  }
  # Number of latent classes
  K <- nrow(est_xi)
  
  # Get new dataframe of desired covariate levels and corresponding Phi values
  Phi_df <- convert_to_probs(est_xi = est_xi, glm_form = res$data_vars$glm_form,
                                 V_data = res$data_vars$V_data, cov_name = cov_name)
  if (!is.null(ci_level)) {
    Phi_lb <- convert_to_probs(est_xi = est_lb, glm_form = res$data_vars$glm_form,
                               V_data = res$data_vars$V_data, cov_name = cov_name)
    Phi_ub <- convert_to_probs(est_xi = est_ub, glm_form = res$data_vars$glm_form,
                               V_data = res$data_vars$V_data, cov_name = cov_name)
  }
  # Convert to long format with a single Class column for plotting
  Phi_df_long <- Phi_df %>%
    tidyr::pivot_longer(cols = 1:K, names_to = "Class", values_to = "Phi")
  if (!is.null(ci_level)) {
    Phi_lb_long <- Phi_lb %>%
      tidyr::pivot_longer(cols = 1:K, names_to = "Class", values_to = "Phi_lb")
    Phi_ub_long <- Phi_ub %>%
      tidyr::pivot_longer(cols = 1:K, names_to = "Class", values_to = "Phi_ub")
    # Get names of all columns except "Phi", "Phi_lb", or "Phi_ub"
    col_names <- colnames(Phi_df_long)[-ncol(Phi_df_long)]
    # Combine the Phi dataframes
    Phi_df_long <- Phi_df_long %>%
      dplyr::left_join(Phi_lb_long, by = col_names) %>%
      dplyr::left_join(Phi_ub_long, by = col_names)
  }
  
  
  # Default labels
  # Set x_title to cov_name if not specified
  if (is.null(x_title)) {
    x_title <- cov_name[1]
  }
  # Set cov_label to be the existing factor levels if NULL
  if (is.null(cov_labels)) {
    cov_labels <- lapply(cov_name, function(x) levels(res$data_vars$V_data[[x]]))
  # Otherwise, check that cov_label vectors are the correct length
  } else {
    for (i in 1:length(cov_labels)) {
      num_categs <- length(levels(res$data_vars$V_data[[cov_name[i]]]))
      if (length(cov_labels[[i]]) != num_categs) {
        stop(paste0("length of cov_labels for covariate ", cov_name[i], 
                    " must equal the number of categories: ", num_categs))
      }
    }
  }
  # Set class_labels to be 1:K if not specified
  if (is.null(class_labels)) {
    class_labels <- 1:K
  } else if (length(class_labels) != K) {
    stop(paste0("length of class_labels must equal the number of latent classes, K = ", K))
  }
  
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Class <- Cov1 <- Cov2 <- Phi <- NULL
  # Relabel the x-axis variable and facet labels if necessary
  Phi_df_long$Cov1 <- factor(Phi_df_long$Cov1, labels = cov_labels[[1]])
  if (length(cov_name) == 2) {
    Phi_df_long$Cov2 <- factor(Phi_df_long$Cov2, labels = cov_labels[[2]])
  }
  # Create plot for one key covariate
  g <- Phi_df_long %>%
    ggplot2::ggplot(ggplot2::aes(x = Cov1, y = Phi, group = Class, col = Class)) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_brewer(palette = "Set2", labels = class_labels) +
    ggplot2::labs(col = class_title, x = x_title, y = y_title) +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::theme(text = ggplot2::element_text(size = 15),
                   axis.text.x = ggplot2::element_text(size = 10, color = "black"),
                   axis.text.y = ggplot2::element_text(size = 10, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 12, color = "black",
                                                        face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 12, color = "black",
                                                        face = "bold"),
                   legend.title = ggplot2::element_text(size = 12, color = "black"),
                   legend.text = ggplot2::element_text(size = 11, color = "black"),
                   legend.position = "top")
    
  # Add faceting and labels if two key covariates
  if (length(cov_name) == 2) {
    g <- g + ggplot2::facet_grid(~ Cov2)
  }
  # Add error bars
  if (!is.null(ci_level)) {
    g <- g + ggplot2::geom_errorbar(ggplot2::aes(ymin = Phi_lb, ymax = Phi_ub, 
                                                 col = Class), 
                                    width = 0.5, alpha = 0.7,
                                    position = ggplot2::position_dodge(width = 0.5))
  }
  # Add lines connecting points
  if (add_lines) {
    g <- g + ggplot2::geom_line(linewidth = 0.7, alpha = 0.3,
                                position = ggplot2::position_dodge(width = 0.5))
  }
  # Return plot
  return(g)
  
  # ggarrange(p1, p2, p3, p4, common.legend = TRUE, legend = "top", nrow = 1, 
  #           ncol = 4, widths = c(0.7, 1, 0.45, 0.45))
}


#' Create traceplots
#' 
#' @param param_mcmc List output from [get_param_mcmc()] containing:
#' \describe{
#'   \item{\code{pi_mcmc}}{MxK dataframe of the class membership probability 
#'   parameters, \eqn{\pi}, where M is the number of iterations and K is the 
#'   number of latent classes.}
#'   \item{\code{theta_mcmc}}{Mx(JxKxR) dataframe of the item level probability 
#'   parameters, \eqn{\theta}, where J is the number of items and R is the 
#'   maximum number of item levels.}
#'   \item{\code{xi_mcmc}}{If output for a `"swolca"` object, Mx(KxQ) dataframe of 
#'   the regression parameters, \eqn{\xi}, where Q is the number of covariates, 
#'   excluding latent class indicators, in the regression.}
#' }
#' @param param_names String vector of parameter names to create plots for. All 
#' names must be found as column names of the dataframes in `param_mcmc`.
#' 
#' @return
#' Creates a grid of traceplots for the specified parameters
#' 
#' @seealso [get_param_mcmc()] [create_acfplot]
#' @importFrom dplyr bind_cols
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes theme_bw geom_line facet_wrap xlab ylab ggtitle
#' @export
#'
#' @examples
#' data(run_nhanes_swolca_results)
#' param_mcmc <- get_param_mcmc(res = run_nhanes_swolca_results)
#' # Specify selection of pi_1 to pi_5
#' param_names <- colnames(param_mcmc$pi_mcmc)
#' # Create traceplots
#' create_traceplot(param_mcmc = param_mcmc, param_names = param_names)
#' 
create_traceplot <- function(param_mcmc, param_names) {
  if (!is.list(param_mcmc)) {
    stop("param_mcmc must be a list outputted from the 'get_param_mcmc()' function")
  }
  if (!is.character(param_names) | length(param_names) < 1) {
    stop("param_names must be a string vector of length at least 1")
  }
  all_colnames <- unlist(lapply(param_mcmc, function(x) colnames(x)))
  if (!(all(param_names %in% all_colnames))) {
    stop("all names in 'param_names' must be found as column names of the dataframes in 'param_mcmc'")
  }
  
  # Create dataframe with all parameters
  all_params <- dplyr::bind_cols(param_mcmc)
  # Select only the specified parameters
  mcmc_df <- all_params[colnames(all_params) %in% param_names]
  # Create iteration variable
  mcmc_df$Iteration <- 1:nrow(mcmc_df)
  # Convert to long format for plotting (facet by parameter)
  mcmc_df_long <- mcmc_df %>%
    tidyr::pivot_longer(cols = -Iteration, names_to = "Parameter", 
                        values_to = "Value")
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Iteration <- Value <- Parameter <- NULL
  
  # Create traceplots for all specified parameters
  mcmc_df_long %>%
    ggplot2::ggplot(ggplot2::aes(x = Iteration, y = Value)) + 
    ggplot2::theme_bw() + 
    ggplot2::geom_line() + 
    ggplot2::facet_wrap(~Parameter) + 
    ggplot2::xlab("Iteration") + 
    ggplot2::ylab("Value") + 
    ggplot2::ggtitle("Trace Plot")
}


#' Create ACF plots
#' 
#' @inheritParams create_traceplot
#' 
#' @return
#' Creates a grid of ACF plots for the specified parameters
#' 
#' @seealso [get_param_mcmc()] [create_traceplot]
#' @importFrom dplyr bind_cols
#' @importFrom stats acf
#' @importFrom ggplot2 ggplot aes theme_bw geom_hline geom_segment facet_wrap 
#' xlab ylab ggtitle
#' @export
#'
#' @examples
#' data(run_nhanes_swolca_results)
#' param_mcmc <- get_param_mcmc(res = run_nhanes_swolca_results)
#' # Specify selection of pi_1 to pi_5
#' param_names <- colnames(param_mcmc$pi_mcmc)
#' # Create ACF plots
#' create_acfplot(param_mcmc = param_mcmc, param_names = param_names)
#' 
create_acfplot <- function(param_mcmc, param_names) {
  if (!is.list(param_mcmc)) {
    stop("param_mcmc must be a list outputted from the 'get_param_mcmc()' function")
  }
  if (!is.character(param_names) | length(param_names) < 1) {
    stop("param_names must be a string vector of length at least 1")
  }
  all_colnames <- unlist(lapply(param_mcmc, function(x) colnames(x)))
  if (!(all(param_names %in% all_colnames))) {
    stop("all names in 'param_names' must be found as column names of the dataframes in 'param_mcmc'")
  }
  
  # Create dataframe with all parameters
  all_params <- dplyr::bind_cols(param_mcmc)
  # Select only the specified parameters
  mcmc_df <- all_params[colnames(all_params) %in% param_names]
  # Get acf and lag for first parameter specified
  param_acfs <- stats::acf(mcmc_df[[param_names[1]]], plot = FALSE)
  acf_df <- data.frame(lag = param_acfs$lag, acf = param_acfs$acf, 
                       parameter = param_names[1])
  # Get acf and lag for remaining parameters
  if (length(param_names) > 1) {
    for (i in 2:length(param_names)) {
      param <- param_names[i]
      param_acfs <- stats::acf(mcmc_df[[param]], plot = FALSE)
      acf_df_param <- data.frame(lag = param_acfs$lag, acf = param_acfs$acf, 
                                 parameter = param)
      acf_df <- rbind(acf_df, acf_df_param)
    }
  }
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  lag <- acf <- parameter <- NULL
  
  # Create ACF plots for all specified parameters
  acf_df %>%
    ggplot2::ggplot(ggplot2::aes(x = lag, y = acf)) + 
    ggplot2::theme_bw() + 
    ggplot2::geom_hline(yintercept = 0, col = "blue", linetype = "dashed") + 
    ggplot2::geom_segment(mapping = ggplot2::aes(xend = lag, yend = 0)) +
    ggplot2::facet_wrap(~parameter) + 
    ggplot2::xlab("Lag") + 
    ggplot2::ylab("ACF") + 
    ggplot2::ggtitle("ACF Plot")
}


