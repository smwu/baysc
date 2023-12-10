#' Plot theta modal exposure categories for each latent class
#'
#' @description
#' `plot_theta_modes` plots a heatmap of the latent class patterns, where the 
#' patterns are defined by the category with the highest probability 
#' (i.e., the model category) for each exposure item.
#'
#' @param res An object of class `"swolca"`, `"solca"`, or `"wolca"`, resulting 
#' from a call to [swolca()] [solca()] or [wolca()]. 
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
#' @seealso [plot_theta_probs()] [plot_pi_boxplots()] [plot_Phi_line()]
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
#' plot_theta_modes(res = run_nhanes_swolca_results)
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
#' plot_theta_modes(res = run_nhanes_swolca_results, item_labels = item_labels,
#'                  categ_labels = categ_labels, class_labels = class_labels)
#'
plot_theta_modes <- function(res, item_labels = NULL, item_title = "Item",
                             categ_labels = NULL, 
                             categ_title = "Consumption Level",
                             class_labels = NULL, 
                             class_title = "Dietary Pattern", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "solca", "wolca"))) {
    stop("res must be an object of class `swolca`, `solca`, or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  # Obtain theta estimates
  if (!is.null(res$estimates)) {
    # Estimates for `solca()` or `wolca()` or adjusted estimates for `swolca()` 
    est_item_probs <- res$estimates$theta_med
  } else {
    # Unadjusted estimates for `swolca()`
    est_item_probs <- res$estimates_unadj$theta_med
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
#' `plot_theta_probs` plots a grouped barplot of the probability of the exposure
#' categories, for each exposure item and each latent class.
#' 
#' @inheritParams plot_theta_modes
#' @param y_title String specifying the title for the y-axis. Default is 
#' `"Consumption Level Probability"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying a grouped barplot of the probability of 
#' the exposure categories, for each exposure item and each latent class
#' 
#' @seealso [plot_theta_modes()] [plot_pi_boxplots()] [plot_Phi_line()]
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
#' plot_theta_probs(res = run_nhanes_swolca_results)
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
#' plot_theta_probs(res = run_nhanes_swolca_results, item_labels = item_labels,
#'                  categ_labels = categ_labels, class_labels = class_labels)
#' 
plot_theta_probs <- function(res, item_labels = NULL, categ_labels = NULL, 
                             categ_title = "Consumption Level",
                             class_labels = NULL,
                             class_title = "Dietary Pattern", 
                             y_title = "Consumption Level Probability", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "solca", "wolca"))) {
    stop("res must be an object of class `swolca`, `solca`, or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  # Obtain theta estimates
  if (!is.null(res$estimates)) {
    # Estimates for `solca()` or `wolca()` or adjusted estimates for `swolca()` 
    est_item_probs <- res$estimates$theta_med
  } else {
    # Unadjusted estimates for `swolca()`
    est_item_probs <- res$estimates_unadj$theta_med
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


#' Plot pi boxplots across posterior samplers for each latent class
#' 
#' @description
#' `plot_pi_boxplots` plots a boxplot of the class membership probabilities, 
#' \eqn{\pi}, in the posterior samples, for each latent class.
#' 
#' @inheritParams plot_theta_modes
#' @param y_title String specifying the title for the y-axis. Default is 
#' `"Class Membership Probability"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying a boxplot of the distribution of class 
#' membership probabilities in the posterior samples, for each latent class.
#' 
#' @seealso [plot_theta_probs()] [plot_theta_modes()] [plot_Phi_line()]
#' 
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_fill_brewer labs 
#' theme_bw theme element_size
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
#' plot_pi_boxplots(res = run_nhanes_swolca_results)
#' 
#' # Specifying labels
#' class_labels <- paste0("Class ", 1:5)
#' plot_pi_boxplots(res = run_nhanes_swolca_results, class_labels = class_labels)
plot_pi_boxplots <- function(res, class_labels = NULL, 
                             class_title = "Dietary Pattern",
                             y_title = "Class Membership Probability", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "solca", "wolca"))) {
    stop("res must be an object of class `swolca`, `solca`, or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  # Obtain pi estimates
  if (!is.null(res$estimates)) {
    # Estimates for `solca()` or `wolca()` or adjusted estimates for `swolca()` 
    pi_red <- as.data.frame(res$estimates$pi_red)
  } else {
    # Unadjusted estimates for `swolca()`
    pi_red <- as.data.frame(res$estimates_unadj$pi_red)
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



#' Plot conditional probability of outcome for each latent class for a specified
#' categorical covariate. 
#' 
#' @description
#' `plot_Phi_line` plots a grouped line plot of the conditional probability of 
#' the outcome, obtained by transforming the probit regression coefficients to 
#' the probability scale, for a categorical covariate, with each line 
#' corresponding to a latent class.
#' 
#' @inheritParams plot_theta_modes
#' @param cov_name String specifying the covariate to plot
#' @param ci_level Numeric from 0 to 1 specifying the credible interval level. 
#' If `NULL` (default), no credible intervals are included in the plot. 
#' Specifying 0.95 gives a 95\% equal-tailed interval composed of the 2.5\% and 
#' 97.5\% quantiles. For `wolca()` results, this must match the `ci_level` 
#' parameter in the main function. 
#' @param cov_labels String vector specifying the category labels for the 
#' covariate of interest. Must be the same length as the number of categories
#' in the covariate specified by `cov_name`. If `NULL` (default), numbers from 1
#' to the number of covariate categories are used.
#' @param x_title String specifying x-axix label. If `NULL` (default), `cov_name`
#' is used
#' @param y_title String specifying the title for the y-axis. Default is 
#' `"Probability of Outcome"`.
#' @param \dots Additional arguments passed
#' 
#' @return
#' Returns a `ggplot2` object displaying a grouped line plot of the conditional 
#' probability of the outcome across categories of the covariate, for each 
#' latent class.
#' 
#' @seealso [get_regr_coefs()] [plot_theta_probs()] [plot_pi_boxplots()] [plot_theta_modes()] 
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_brewer labs 
#' theme_bw theme element_text
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate left_join
#' @importFrom stats quantile
#' @importFrom magrittr %>%
#' @importFrom rlang .data 
#' @export
#'
#' @examples 
#' data(run_nhanes_swolca_results)
#' 
#' # Default labels
#' plot_Phi_line(res = run_nhanes_swolca_results, cov_name = "racethnic")
#' 
#' # Specifying labels
#' cov_labels <- c("NH White", "NH Black", "NH Asian", "Hispanic/Latino", 
#'                 "Other/Mixed")
#' class_labels <- paste0("Class", 1:5)
#' x_title <- "Race and Ethnicity"
#' plot_Phi_line(res = run_nhanes_swolca_results, cov_name = "racethnic",
#'              cov_labels = cov_labels, class_labels = class_labels, 
#'              x_title = x_title)
plot_Phi_line <- function(res, cov_name, ci_level = NULL,
                          cov_labels = NULL, class_labels = NULL, 
                          class_title = "Dietary Pattern", x_title = NULL, 
                          y_title = "Probability of Outcome", ...) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "solca", "wolca"))) {
    stop("res must be an object of class `swolca`, `solca`, or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  if (!is.null(ci_level)) {
    if (!(ci_level > 0 & ci_level < 1)) {
      stop("ci_level must be between 0 and 1")
    }
    quant_lb <- (1 - ci_level) / 2
    quant_ub <- 1 - quant_lb
  } 
  
  
  # Obtain xi median and lower bound and upper bound estimates 
  if (!is.null(res$estimates)) {
    if (!is.null(res$estimates$xi_med)) {
      # Estimates for `solca()` or adjusted estimates for `swolca()` 
      est_xi <- res$estimates$xi_med
      if (!is.null(ci_level)) {
        est_lb <- apply(res$estimates$xi_red, c(2, 3), 
                        function(x) quantile(x, quant_lb))
        est_ub <- apply(res$estimates$xi_red, c(2, 3),
                        function(x) quantile(x, quant_ub))
      }
    } else {
      # Estimates for `wolca()`
      est_xi <- res$estimates$xi_est
      if (ci_level != res$data_vars$ci_level) {
        stop("ci_level must match the specified ci_level in the wolca() function")
      }
      if (!is.null(ci_level)) {
        est_lb <- res$estimates$xi_est_lb
        est_ub <- res$estimates$xi_est_ub
      }
    }
  } else {
    # Unadjusted estimates for `swolca()`
    est_xi <- res$estimates_unadj$xi_med
    if (!is.null(ci_level)) {
      est_lb <- apply(res$estimates_unadj$xi_red, c(2, 3), 
                      function(x) stats::quantile(x, quant_lb))
      est_ub <- apply(res$estimates_unadj$xi_red, c(2, 3),
                      function(x) stats::quantile(x, quant_ub))
    }
  }
  
  # Convert to Phi
  Phi_df <- convert_to_probs(est_xi = est_xi, glm_form = res$data_vars$glm_form,
                             V = res$data_vars$V, cov_name = cov_name)
  if (!is.null(ci_level)) {
    Phi_lb <- convert_to_probs(est_xi = est_lb, glm_form = res$data_vars$glm_form,
                               V = res$data_vars$V, cov_name = cov_name)
    Phi_ub <- convert_to_probs(est_xi = est_ub, glm_form = res$data_vars$glm_form,
                               V = res$data_vars$V, cov_name = cov_name)
  }
  
  # Set x_title to cov_name if not specified
  if (is.null(x_title)) {
    x_title <- cov_name
  }
  # Set cov_label to be numbered from one to the number of covariate categories 
  num_categs <- ncol(Phi_df) - 1
  if (is.null(cov_labels)) {
    cov_labels <- paste0(cov_name, 1:num_categs)
  } else if (length(cov_labels) != num_categs) {
    stop(paste0("length of cov_labels must equal the number of covariate categories: ", 
                num_categs))
  }
  # Set class_labels to be 1:K if not specified
  K <- nrow(est_xi)
  if (is.null(class_labels)) {
    class_labels <- 1:K
  } else if (length(class_labels) != K) {
    stop(paste0("length of class_labels must equal the number of latent classes, K = ", K))
  }
  
  # Initialize variables to NULL to avoid global binding notes in R CMD check
  Class <- Cov <- Phi <- NULL
  
  # Convert to longer format for plotting
  Phi_df_long <- Phi_df %>%
    tidyr::pivot_longer(cols = -Class, names_to = "Cov", values_to = "Phi") %>%
    dplyr::mutate(Cov = factor(Cov, labels = cov_labels),
                  Class = as.factor(Class))
  
  if (!is.null(ci_level)) {
    Phi_lb_long <- Phi_lb %>%
      tidyr::pivot_longer(cols = -Class, names_to = "Cov", values_to = "Phi_lb") %>%
      dplyr::mutate(Cov = factor(Cov, labels = cov_labels),
                    Class = as.factor(Class))
    Phi_ub_long <- Phi_ub %>%
      tidyr::pivot_longer(cols = -Class, names_to = "Cov", values_to = "Phi_ub") %>%
      dplyr::mutate(Cov = factor(Cov, labels = cov_labels),
                    Class = as.factor(Class))
    Phi_all_long <- Phi_df_long %>%
      dplyr::left_join(Phi_lb_long, by = c("Class", "Cov")) %>%
      dplyr::left_join(Phi_ub_long, by = c("Class", "Cov"))
  }
  
  # Plot with credible intervals
  if (!is.null(ci_level)) {
    Phi_all_long %>% ggplot2::ggplot(ggplot2::aes(x = Cov, 
                                                  y = Phi, 
                                                  group = Class, 
                                                  col = Class)) + 
      ggplot2::theme_bw() + 
      ggplot2::scale_color_brewer(palette = "Set2", labels = class_labels, 
                                  aesthetics = c("color", "fill")) +
      ggplot2::labs(col = class_title, x = x_title, y = y_title) + 
      ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2) + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = Phi_lb, ymax = Phi_ub, 
                                        fill = Class, colour = NA),
                           alpha = 0.2, show.legend = FALSE) + 
      ggplot2::theme(text = ggplot2::element_text(size = 15),
                     axis.text.x = ggplot2::element_text(size = 11, color = "black"), 
                     axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                     axis.title.x = ggplot2::element_text(size = 13, color = "black", 
                                                          face = "bold"),
                     axis.title.y = ggplot2::element_text(size = 13, color = "black", 
                                                          face = "bold"),
                     legend.title = ggplot2::element_text(size = 14, color = "black"),
                     legend.text = ggplot2::element_text(size = 13, color = "black"),
                     legend.position = "top")
    # Plot without credible intervals
  } else { 
    Phi_df_long %>% ggplot2::ggplot(ggplot2::aes(x = Cov, 
                                                 y = Phi, 
                                                 group = Class, 
                                                 col = Class)) + 
      ggplot2::theme_bw() + 
      ggplot2::scale_color_brewer(palette = "Set2", labels = class_labels) +
      ggplot2::labs(col = class_title, x = x_title, y = y_title) + 
      ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2) + 
      ggplot2::theme(text = ggplot2::element_text(size = 15),
                     axis.text.x = ggplot2::element_text(size = 11, color = "black"), 
                     axis.text.y = ggplot2::element_text(size = 11, color = "black"),
                     axis.title.x = ggplot2::element_text(size = 13, color = "black", 
                                                          face = "bold"),
                     axis.title.y = ggplot2::element_text(size = 13, color = "black", 
                                                          face = "bold"),
                     legend.title = ggplot2::element_text(size = 14, color = "black"),
                     legend.text = ggplot2::element_text(size = 13, color = "black"),
                     legend.position = "top")
  }
  
  # ggarrange(p1, p2, p3, p4, common.legend = TRUE, legend = "top", nrow = 1, 
  #           ncol = 4, widths = c(0.7, 1, 0.45, 0.45))
}





#' Obtains table of regression coefficients
#' 
#' @description
#' `get_regr_coefs` produces a summary table of the regression coefficients,
#' converted to standard reference cell coding. 
#' 
#' @inheritParams plot_theta_modes
#' @inheritParams plot_Phi_line
#' @param ci_level Numeric from 0 to 1 specifying the credible interval level. 
#' Default is 0.95, which gives a 95\% equal-tailed interval composed of the 
#' 2.5\% and 97.5\% quantiles. For `wolca()` results, this must match the 
#' `ci_level` parameter in the main function. 
#' @param digits Integer indicating the number of decimal places to be used. 
#' Default is 2, which rounds to the nearest hundredth. 
#' 
#' @return
#' Returns a character vector of the source code for regression coefficients 
#' table using the format specified in `format`.
#' 
#' @seealso [plot_Phi_line()] 
#' 
#' @importFrom dplyr mutate_if
#' @importFrom stats terms as.formula median quantile
#' @export
#'
#' @examples
#' data(run_nhanes_swolca_results)
#' get_regr_coefs(res = run_nhanes_swolca_results, ci_level = 0.95, digits = 2)
#' 
get_regr_coefs <- function(res, ci_level = 0.95, digits = 2) {
  
  # Check object class
  if (!(class(res) %in% c("swolca", "solca", "wolca"))) {
    stop("res must be an object of class `swolca`, `solca`, or `wolca`, resulting 
         from a call to one of these functions")
  }
  
  if (!(ci_level > 0 & ci_level < 1)) {
    stop("ci_level must be between 0 and 1")
  }
  quant_lb <- (1 - ci_level) / 2
  quant_ub <- 1 - quant_lb
  
  # Estimates for `wolca()` can be obtained directly from the svyglm output
  if (is.null(res$estimates$xi_med)) {
    if (ci_level != res$data_vars$ci_level) {
      stop("ci_level must match the specified ci_level in the wolca() function")
    }
    beta <- as.data.frame(summary(res$estimates$fit)$coefficients)
    beta[, 4] <- ifelse(beta[, 4] < 10^(-digits), paste0("<", 10^(-digits)), 
                        format(round(beta[, 4], digits), digits))
    
    # Otherwise, for `swolca()` and `solca()`, create regression output after
    # converting from factor reference coding
  } else {
    # Obtain xi median and lower bound and upper bound estimates
    if (!is.null(res$estimates)) {
      # Estimates for `solca()` or adjusted estimates for `swolca()`
      est_xi <- res$estimates$xi_med
      est_lb <- apply(res$estimates$xi_red, c(2, 3),
                      function(x) stats::quantile(x, quant_lb))
      est_ub <- apply(res$estimates$xi_red, c(2, 3),
                      function(x) stats::quantile(x, quant_ub))
      est_red <- res$estimates$xi_red
    } else {
      # Unadjusted estimates for `swolca()`
      est_xi <- res$estimates_unadj$xi_med
      est_lb <- apply(res$estimates_unadj$xi_red, c(2, 3),
                      function(x) stats::quantile(x, quant_lb))
      est_ub <- apply(res$estimates_unadj$xi_red, c(2, 3),
                      function(x) stats::quantile(x, quant_ub))
      est_red <- res$estimates$xi_red
    }
    
    # Add outcome and latent class main and interaction terms to formula
    terms <- labels(stats::terms(stats::as.formula(res$data_vars$glm_form)))
    if (length(terms) > 0) {
      full_glm_form <- paste0("y_all ~ ", 
                              paste0("c_all * ", terms, collapse = " + ")) 
    } else {
      full_glm_form <- paste0("y_all ~ c_all") 
    }
    full_data <- data.frame(c_all = as.factor(res$estimates$c_all), 
                            y_all = res$data_vars$Y_data,
                            res$data_vars$V_data)
    model_matrix <- model.matrix(as.formula(full_glm_form), data = full_data)
    
    K <- nrow(est_xi)
    q <- ncol(est_xi)
    beta <- as.data.frame(matrix(NA, nrow = ncol(model_matrix), ncol = 4))
    beta[, 1] <- colnames(model_matrix)
    colnames(beta) <- c("Covariate", "Estimate", "95% Cred Int",
                        "P(xi > 0)")
    
    # Intercept estimates
    beta[1, -1] <- c(est_xi[1, 1], get_ci(post_samp = est_red[, 1, 1]),
                     get_prob_pos(est_red[, 1, 1]))
    
    # Latent class main effect estimates
    for (i in 2:K) {
      beta[i, -1] <- c(stats::median(est_red[, i, 1] - est_red[, 1, 1]),
                       get_ci(est_red[, i, 1] - est_red[, 1, 1], digits = digits),
                       get_prob_pos(est_red[, i, 1] - est_red[, 1, 1], digits = digits))
    }
    
    # Additional covariates main effect estimates
    for (i in 2:q) {
      beta[K + (i-1), -1] <- c(est_xi[1, i], get_ci(est_red[, 1, i]),
                               get_prob_pos(est_red[, 1, i]))
    }
    
    # Additional covariates latent class interaction terms
    for (i in 2:q) {
      for (j in 2:K) {
        beta[q + (i-1)*(K-1) + (j-1), -1] <- 
          c(stats::median(est_red[, j, i] - est_red[, 1, i]),
            get_ci(est_red[, j, i] - est_red[, 1, i], digits = digits),
            get_prob_pos(est_red[, j, i] - est_red[, 1, i], digits = digits))
      }
    }
    beta$Estimate <- as.numeric(beta$Estimate)
  }
  
  # Print output
  beta <- dplyr::mutate_if(beta, is.numeric, round, digits = digits)
  beta
}
