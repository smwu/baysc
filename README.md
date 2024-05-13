# **baysc**: BAYesian Survey Clustering

An R package for running Bayesian supervised and unsupervised clustering methods on survey data.

**Maintainer**: Stephanie M. Wu ([swu\@g.harvard.edu](mailto:swu@g.harvard.edu){.email})

**Contributors**: Matthew R. Williams ([mrwilliams\@rti.org](mailto:mrwilliams@rti.org){.email}); Terrance D. Savitsky ([savitsky.terrance\@bls.gov](mailto:savitsky.terrance@bls.gov){.email}); Briana J.K. Stephenson ([bstephenson\@hsph.harvard.edu](mailto:bstephenson@hsph.harvard.edu){.email})

| Citation                                                                                                                                                                                                                                                   | Link                                              |
|------------------------------------------------|-----------------------|
| Wu, S. M., Williams, M. R., Savitsky, T. D., & Stephenson, B. J. (2023). Derivation of outcome-dependent dietary patterns for low-income women obtained from survey data using a Supervised Weighted Overfitted Latent Class Analysis. *ArXiv:2310.01575*. | [Link to paper](https://arxiv.org/abs/2310.01575) |

## Table of content

-   [1. Installation](#id-section1)
-   [2. Overview](#id-section2)
-   [3. Functions](#id-section3)
-   [2. Data](#id-section4)
-   [2. Example](#id-section5)
-   [2. Help](#id-section6)

<div id='id-section1'/>

## Installation

``` r
# Install devtools for package loading 
install.packages(devtools)
library(devtools)
# Install baysc from GitHub
devtools::install_github("smwu/baysc")
library(baysc)
```

<div id='id-section2'/>

## Overview

`baysc` is an R package for running Bayesian clustering methods on survey data. A Bayesian latent class analysis (LCA), termed the Weighted Overfitted Latent Class Analysis (WOLCA), is available for eliciting underlying cluster patterns from multivariate categorical data, incorporating survey sampling weights and other survey design elements. Options also exist for relating the patterns to a binary outcome, either by using a two-step approach that applies WOLCA and then runs a survey-weighted regression, or by utilizing a one-step supervised approach where creation of the cluster patterns is directly informed by the outcome, referred to as the Supervised Weighted Overfitted Latent Class Analysis (SWOLCA). More information about the models can be found in the paper linked above. Summary and plotting functions for visualizing output are also available, as are diagnostic functions for examining convergence of the sampler.

<div id='id-section3'/>

## Functions

Use the `wolca()` function to run an unsupervised WOLCA and obtain pattern profiles. `wolca_var_adjust()` provides a post-hoc variance adjustment that enables correct uncertainty estimation. To examine the association of pattern profiles with a binary outcome through a two-step appraoch, run `wolca_svyglm()`. Use the `swolca()` function to run a SWOLCA model that allows information about the binary outcome to directly inform the creation of the pattern profiles. `swolca_var_adjust()` provides a post-hoc variance adjustment that enbales correct uncertainty estimation. Detailed information about the functions and related statistical details can be found in the vignette, "[An introduction to the basyc package](vignettes/baysc.pdf)," and in the paper linked above.

<div id='id-section4'/>

## Data

`baysc` applies Bayesian latent class analysis using the following input data:

-   Multivariate categorical exposure: $nxJ$ matrix, where $n$ is the sample size and $J$ is the number of categorical item variables. Each item must be a categorical variable.
-   (Optional) binary outcome: $nx1$ vector
-   (Optional) additional confounders to adjust for when evaluating the exposure-outcome association: $nxQ$ dataframe, where $Q$ is the number of additional confounders.
-   (Optional) survey design elements such as stratum indicators, cluster indicators, and sampling weights: each formatted as a $nx1$ vector.

We provide an example dataset from the National Health and Nutrition Examination Survey (NHANES) that includes multivariate categorical dietary intake data as well as binary hypertension data for low-income women in the United States. Survey sampling weights and information on stratification and clustering are included to allow for adjustment for survey design when conducting estimation and inference.

A simulated dataset named `sim_data` is also provided to help with familiarization of the package, and simulation functions `simulate_pop()` and `simulate_samp()` are available for creating example data population and survey sample data for simulation studies, with more information provided in the "[Simulating Data](vignettes/simulating_data.pdf)" vignette.

<div id='id-section5'/>

## Example

``` r
library(baysc)

#==== Create data ====
# Load NHANES dataset
data("data_nhanes")

# Exposure matrix composed of food groups
x_mat <- as.matrix(data_nhanes[, 11:38])

# Outcome data on hypertension
y_all <- data_nhanes$BP_flag

# Survey stratum indicators
stratum_id <- data_nhanes$stratum_id
# Survey cluster indicators
cluster_id <- data_nhanes$cluster_id
# Survey sampling weights
sampling_wt <- data_nhanes$sample_wt

# Create dataframe of additional confounders
V_data <- data_nhanes[, c("age_cat", "racethnic", "smoker", "physactive")]
# Regression formula for additional confounders
glm_form <- "~ age_cat + racethnic + smoker + physactive"

#==== Run model ====
# Run SWOLCA
res_swolca <- swolca(x_mat = x_mat, y_all = y_all, V_data = V_data,
                     glm_form = glm_form, sampling_wt = sampling_wt,
                     cluster_id = cluster_id, stratum_id = stratum_id,
                     adapt_seed = 888, n_runs = 300, burn = 150, thin = 3,
                     update = 50, save_res = FALSE)
# Apply variance adjustment
res_swolca_adjust <- swolca_var_adjust(res = res_swolca, adjust_seed = 888,
                                       num_reps = 100, save_res = FALSE)

#==== Display results ====
# Plot derived patterns
plot_pattern_profiles(res = res_swolca_adjust)
# Plot outcome regression coefficients
regr_coefs <- get_regr_coefs(res = res_swolca_adjust, ci_level = 0.95, digits = 2)
plot_regr_coefs(regr_coefs = regr_coefs, res = res_swolca_adjust)
```

When running the variance adjustment, the following warning messages may appear: *"the number of chains is less than 1; sampling not done"* and *"In mrbweights(design\$cluster, design\$strata, design\$fpc, ...) : Design is sampled with replacement: only first stage used."* These messages do not pose an issue to the statistical validity of the methods and can be ignored.

<div id='id-section6'/>

## Contributing and Getting Help

Please report bugs by opening an [issue](https://github.com/smwu/baysc/issues/new/choose). If you wish to contribute, please make a pull request. If you have questions, you can open a [discussion thread](https://github.com/smwu/baysc/discussions).
