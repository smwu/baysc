#' Plot theta modal exposure categories for each latent class
#'
#' @description
#' `plot_theta_modes` plots a heatmap of the latent class patterns, where the 
#' patterns are defined by the category with the highest probability 
#' (i.e., the model category) for each exposure item.
#'
#' @param res Output from `swolca()`, `solca()`, or `wolca()` containing at 
#' least `estimates` (or `estimates_unadj`) and `data_vars`
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
#' @seealso [swolca()] [solca()] [wolca()]
#' 
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_brewer labs 
#' theme_classic scale_x_discrete theme element_text
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' #examples 
#' 
plot_theta_modes <- function(res, item_labels = NULL, item_title = "Item",
                             categ_labels = NULL, 
                             categ_title = "Consumption Level",
                             class_labels = NULL, 
                             class_title = "Dietary Pattern", ...) {
  
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
  }
  if (is.null(class_labels)) {
    K <- dim(mode_item_probs)[2]
    class_labels <- 1:K
  }
  if (is.null(categ_labels)) {
    categ_labels <- 1:res$data_vars$R
  }
  rownames(mode_item_probs) <- item_labels
  colnames(mode_item_probs) <- class_labels
  mode_item_probs$Item <- rownames(mode_item_probs)
  
  # Create plot
  mode_plot <- mode_item_probs %>% tidyr::gather("Class", "Level", -(rlang::.data$Item)) 
  mode_plot %>% ggplot2::ggplot(ggplot2::aes(x=rlang::.data$Class, 
                                             y=factor(rlang::.data$Item, 
                                                      levels = rev(item_labels)), 
                                             fill = factor(rlang::.data$Level))) + 
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
#' @seealso [swolca()] [solca()] [wolca()]
#' 
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap scale_fill_brewer labs 
#' theme_bw theme element_text element_blank element_rect
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' #examples 
#' 
plot_theta_probs <- function(res, item_labels = NULL, categ_labels = NULL, 
                             categ_title = "Consumption Level",
                             class_labels = NULL,
                             class_title = "Dietary Pattern", 
                             y_title = "Consumption Level Probability", ...) {
  
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
  }
  if (is.null(class_labels)) {
    class_labels <- 1:K
  }
  if (is.null(categ_labels)) {
    categ_labels <- 1:res$data_vars$R
  }
  
  dimnames(est_item_probs)[[1]] <- item_labels
  dimnames(est_item_probs)[[2]] <- class_labels
  
  # Create plot
  
  # Convert to dataframe with each row corresponding to a value in the array
  # Use base R instead of reshape2 to reduce number of package imports
      # theta_plot <- reshape2::melt(est_item_probs, level = 2)
  theta_plot <- data.frame(expand.grid(lapply(dim(est_item_probs), seq_len)), 
                           value = as.vector(est_item_probs))
  theta_plot <- est_item_probs %>% tidyr::gather("Class", "Level", -(rlang::.data$Item)) 
  
  colnames(theta_plot) <- c("Item", "Class", "Level", "Probability")
  theta_plot %>%
    ggplot2::ggplot(ggplot2::aes(x = factor(rlang::.data$Class, levels = 1:K), 
                                 y = rlang::.data$Probability,
                                 fill = factor(rlang::.data$Level))) + 
    ggplot2::geom_bar(stat = "identity", position = "stack") + 
    ggplot2::facet_wrap(factor(rlang::.data$Item, levels = item_labels) ~ ., nrow = 4) + 
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
#' @seealso [swolca()] [solca()] [wolca()]
#' 
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_fill_brewer labs 
#' theme_bw theme element_size
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
#'
#' #examples 
plot_pi_boxplots <- function(res, class_labels = NULL, 
                             class_title = "Dietary Pattern",
                             y_title = "Class Membership Probability", ...) {
  
  # Obtain pi estimates
  if (!is.null(res$estimates)) {
    # Estimates for `solca()` or `wolca()` or adjusted estimates for `swolca()` 
    pi_red <- as.data.frame(res$estimates$pi_red)
  } else {
    # Unadjusted estimates for `swolca()`
    pi_red <- as.data.frame(res$estimates_unadj$pi_red)
  }
  
  # Set class labels to 1:K if not provided
  if (is.null(class_labels)) {
    class_labels <- 1:(dim(pi_red)[2])
  } 
  colnames(pi_red) <- class_labels
  # Convert to longer format for plotting
  pi_red_plot <- pi_red %>% tidyr::pivot_longer(cols = tidyselect::everything(), 
                                                names_to = "pi_comp", 
                                                values_to = "value")
  # Plot pi boxplots
  pi_red_plot %>% ggplot2::ggplot(ggplot2::aes(x = rlang::.data$pi_comp, 
                                               y = rlang::.data$value)) + 
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
#' in the covariate specified by `cov_name`
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
#' @seealso [swolca()] [solca()] [wolca()]
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
#' #examples 
plot_Phi_line <- function(res, cov_name, ci_level = NULL,
                          cov_labels = NULL, class_labels = NULL, 
                          class_title = "Dietary Pattern", x_title = NULL, 
                          y_title = "Probability of Outcome", ...) {
  
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
  if (is.null(cov_labels)) {
    cov_labels <- paste0(cov_name, 1:(ncol(Phi_df) - 2))
  }
  # Set class_labels to be 1:K if not specified
  if (is.null(class_labels)) {
    class_labels <- 1:nrow(est_xi)
  }
  
  # Convert to longer format for plotting
  Phi_df_long <- Phi_df %>%
    tidyr::pivot_longer(cols = -rlang::.data$Class, names_to = cov_name, values_to = "Phi") %>%
    dplyr::mutate(Cov = factor(rlang::.data[[cov_name]], labels = cov_labels),
                  Class = as.factor(rlang::.data$Class))
  
  if (!is.null(ci_level)) {
    Phi_lb_long <- Phi_lb %>%
      tidyr::pivot_longer(cols = -rlang::.data$Class, names_to = cov_name, values_to = "Phi_lb") %>%
      dplyr::mutate(Cov = factor(rlang::.data[[cov_name]], labels = cov_labels),
                    Class = as.factor(rlang::.data$Class))
    Phi_ub_long <- Phi_ub %>%
      tidyr::pivot_longer(cols = -rlang::.data$Class, names_to = cov_name, values_to = "Phi_ub") %>%
      dplyr::mutate(Cov = factor(rlang::.data[[cov_name]], labels = cov_labels),
                    Class = as.factor(rlang::.data$Class))
    Phi_all_long <- Phi_df_long %>%
      dplyr::left_join(Phi_lb_long, by = c("Class", "Cov")) %>%
      dplyr::left_join(Phi_ub_long, by = c("Class", "Cov"))
  }
  
  # Plot with credible intervals
  if (!is.null(ci_level)) {
    Phi_all_long %>% ggplot2::ggplot(ggplot2::aes(x = rlang::.data$Cov, 
                                                  y = rlang::.data$Phi, 
                                                  group = rlang::.data$Class, 
                                                  col = rlang::.data$Class)) + 
      ggplot2::theme_bw() + 
      ggplot2::scale_color_brewer(palette = "Set2", labels = class_labels, 
                                  aesthetics = c("color", "fill")) +
      ggplot2::labs(col = class_title, x = x_title, y = y_title) + 
      ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2) + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = rlang::.data$Phi_lb, ymax = rlang::.data$Phi_ub, 
                                        fill = rlang::.data$Class, colour = NA),
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
    Phi_df_long %>% ggplot2::ggplot(ggplot2::aes(x = rlang::.data$Cov, 
                                                 y = rlang::.data$Phi, 
                                                 group = rlang::.data$Class, 
                                                 col = rlang::.data$Class)) + 
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


