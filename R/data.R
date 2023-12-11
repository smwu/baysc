#' NHANES diet and hypertension for low-income women 
#' 
#' Cleaned dataset containing NHANES 2015-2018 data on dietary intake and 
#' hypertension among adult women aged 20 or over who are classified as 
#' low-income (reported household income at or below 185\% of the federal 
#' poverty level). 
#' 
#' @docType data
#' @usage data(data_nhanes)
#' @format A dataframe with 2004 rows and 71 variables:
#' \describe{
#'   \item{\code{SEQN}}{5- or 6-digit unique individual identifier}
#'   \item{\code{stratum_id}}{Survey design stratum indicator}
#'   \item{\code{cluster_id}}{Survey design cluster indicator, modified to be 
#'   unique across strata using the formula: SDMVSTRA*10 + SDMVPSU}
#'   \item{\code{sample_wt}}{Individual survey sampling weights}
#'   \item{\code{age_cat}}{Categorical age variable with three categories: 
#'   1: 20-39, 2: 40-59, 3: >=60, created from the RIAGEYR variable}
#'   \item{\code{racethnic}}{Race and Hispanic origin ethnicity with 5 categories:
#'   1: NH White, 2: NH Black, 3: NH Asian, 4: Hispanic, 5: Other/Mixed, created
#'   from the RIDRETH3 variable}
#'   \item{\code{educ}}{Education level with 3 categories: 1: at least some 
#'   college, 2: high school/GED, 3: less than high school, created from the
#'   DMDEDUC2 variable}
#'   \item{\code{smoker}}{Binary current smoking status, defined as 1 (yes)
#'   if having smoked 100 cigarettes in lifetime and currently smokes cigarettes,
#'   and 0 (no) otherwise, created from SMQ040 and SMQ020 variables}
#'   \item{\code{physactive}}{Physically active as a factor with levels "Inactive"
#'   and "Active", with active defined as having at least 150 minutes of moderate
#'   or vigorous exercise per week}
#'   \item{\code{BP_flag}}{Binary 0/1 hypertension variable, with 1 (yes) 
#'   defined as having average systolic blood pressure (BP) above 130, having 
#'   average diastolic BP above 80, currently taking BP medication, or having 
#'   been told have high BP at least two times}
#'   \item{\code{citrus}}{Consumption of citrus, melons, or berries with 4
#'   categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{oth_fruit}}{Consumption of other fruits with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{fruit_juice}}{Consumption of fruit juices with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{dark_green}}{Consumption of dark green vegetables with 4
#'   categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{tomatoes}}{Consumption of tomatoes and tomato products with 4
#'   categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{oth_red}}{Consumption of other red and orange vegetables with 4
#'   categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{potatoes}}{Consumption of white potatoes with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{oth_starchy}}{Consumption of other starchy vegetables with 4
#'   categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{oth_veg}}{Consumption of other vegetables with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{whole_grain}}{Consumption of whole grains with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{ref_grain}}{Consumption of refined grains with 4
#'   categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{meat}}{Consumption of beef, veal, pork, lamb, and game meat, 
#'   excluding organ and cured meat, with 4 categories: 1) None, 2) Low, 3) Med, 
#'   4) High}
#'   \item{\code{cured_meats}}{Consumption of cured meats with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{organ}}{Consumption of organ meat from beef, veal, pork, lamb,
#'   game, and poultry with 4 categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{poultry}}{Consumption of poultry, excluding organ and cured 
#'   meat, with 4 categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{seafood_high}}{Consumption of seafood high in n-3 fatty acids, 
#'   with 4 categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{seafood_low}}{Consumption of seafood low in n-3 fatty acids, 
#'   with 4 categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{eggs}}{Consumption of eggs and egg substitutdes, with 4 
#'   categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{soybean}}{Consumption of soy products, excluding soymilk and 
#'   products made with raw soybean, with 4 categories: 1) None, 2) Low, 3) Med, 
#'   4) High}
#'   \item{\code{nuts}}{Consumption of nuts and seeds with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{leg_protein}}{Consumption of legumes (beans, peas, lentils) 
#'   computed as protein foods, with 4 categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{milk}}{Consumption of milk and soymilk, with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{yogurt}}{Consumption of yogurt, with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{cheese}}{Consumption of cheese, with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{oils}}{Consumption of oils, with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{solid_fats}}{Consumption of solid fats, with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{add_sugars}}{Consumption of added sugars, with 4 categories: 
#'   1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{drinks}}{Consumption of alcoholic beverages and alcohol added 
#'   to foods after cooking, with 4 categories: 1: None, 2: Low, 3: Med, 4: High}
#'   \item{\code{cycle}}{NHANES two-year cycle, either "2015-2016" or "2017-2018"}
#'   \item{\code{nrecall}}{Number of 24-hr dietary recalls completed by the 
#'   individual}
#'   \item{\code{SDMVSTRA}}{Original NHANES stratum indicator}
#'   \item{\code{SDMVPSU}}{Original NHANES cluster indicator}
#'   \item{\code{RIAGENDR}}{Original NHANES gender variable. Equal to 2 for all 
#'   individuals because restricted to women.}
#'   \item{\code{RIDAGEYR}}{Original NHANES age in years. Restricted to those 
#'   age 20 or over.}
#'   \item{\code{RIDRETH3}}{Original NHANES race/ethnicty variable with 7 
#'   categories: 1=Mex_Amer, 2=Other_Hisp, 3=NH_White, 4=NH_Black, 6=NH_Asian, 
#'   7=Other/Mixed}
#'   \item{\code{DMDEDUC2}}{Original NHANES education level with categories:
#'   1= <9th, 2=9-11th, 3=HS/GED, 4=Some college/AA, 5=college grad or above, 
#'   7=refused, 9=don't know}
#'   \item{\code{SMQ020}}{Original NHANES binary variable of having smoked at 
#'   least 100 cigarettes in life}
#'   \item{\code{SMQ040}}{Original NHANES binary variable of currently smoking
#'   cigarettes}
#'   \item{\code{Mins_Active}}{Number of minutes of moderate or vigorous 
#'   exercise per week, created from NHANES variables PAQ610, PAQ615,
#'   PAQ625, PAQ630, PAQ640, PAQ645, PAQ655, PAQ660, PAQ670, PAQ675}
#'   \item{\code{SBP_avg}}{Average systolic BP over four readings, created from
#'   NHANES variables BPXSY1, BPXSY2, BPXSY3, BPXSY4}
#'   \item{\code{DBP_avg}}{Average diastolic BP over four readings, created from
#'   NHANES variables BPXDI1, BPXDI2, BPXDI3, BPXDI4}
#'   \item{\code{BPQ030}}{Original NHANES binary variable of having been told 
#'   have high BP at least two times}
#'   \item{\code{BPQ050A}}{Original NHANES binary variable of currently taking
#'   prescribed medication for high BP}
#'   \item{\code{BPQ020}}{Original NHANES binary variable of having ever been 
#'   told have high BP}
#' }
#' @source https://github.com/smwu/SWOLCA/tree/main/Data. For more information 
#' on dataset preparation, see Wu et al. (2023)
#' @references Wu, S. M., Williams, M. R., Savitsky, T. D., & Stephenson, B. J. 
#' (2023). Derivation of outcome-dependent dietary patterns for low-income women 
#' obtained from survey data using a Supervised Weighted Overfitted Latent Class 
#' Analysis. arXiv preprint arXiv:2310.01575.
#' @keywords datasets
#' @examples 
#' data(data_nhanes)
#' x_mat <- data_nhanes %>% dplyr::select(citrus:drinks)
#' y_all <- data_nhanes$BP_flag
"data_nhanes"



#' Simulated data
#' 
#' A dataset simulated using [simulate_pop()] and [simulate_samp()].
#' @docType data
#' @usage data(sim_data)
#' @format A list with 23 elements:
#' \describe{
#'   \item{\code{samp_ind}}{Vector of population indices for sampled 
#'   individuals. nx1, where n=4000.}
#'   \item{\code{sample_wt}}{Vector of sampling weights for sampled 
#'   individuals. nx1}
#'   \item{\code{N}}{Population size}
#'   \item{\code{J}}{Number of exposure items}
#'   \item{\code{R}}{Number of exposure categories}
#'   \item{\code{H}}{Number of stratum (i.e., levels of S)}
#'   \item{\code{N_s}}{Vector of population sizes for levels of S. Hx1, where H = 2.}
#'   \item{\code{true_K}}{True number of latent classes, K}
#'   \item{\code{true_Ai}}{`NULL` or vector of additional continuous variable 
#'   for sampled individuals. nx1}
#'   \item{\code{true_Bi}}{`NULL` or vector of additional binary variable for 
#'   sampled individuals. nx1}
#'   \item{\code{true_Si}}{Vector of true stratum indicators for sampled 
#'   individuals. nx1}
#'   \item{\code{true_Ci}}{Vector of true latent class indicators for sampled 
#'   individuals. nx1}
#'   \item{\code{true_pi}}{Vector of true pi values overall in the population. 
#'   Kx1, where K = 3.}
#'   \item{\code{true_pi_s}}{`NULL` or HxK matrix of true pi values within each 
#'   level of S}
#'   \item{\code{X_data}}{Matrix of multivariate categorical exposure for 
#'   sampled individuals. nxJ, where J = 30.}
#'   \item{\code{true_global_patterns}}{Matrix of true global exposure patterns 
#'   defined by modal category. JxK}
#'   \item{\code{true_global_thetas}}{Array of true thetas. JxKxR, where R = 4.}
#'   \item{\code{Y_data}}{Vector of binary outcome for sampled individuals. nx1}
#'   \item{\code{cluster_id}}{Vector of true cluster indicators for sampled 
#'   individuals. nx1}
#'   \item{\code{cluster_size}}{Cluster size}
#'   \item{\code{true_xi}}{Matrix of probit regression coefficients in mixture 
#'   reference coding. Kxq, where q is the number of covariate terms in the 
#'   regression, excluding all terms involving C; q = 2 for this dataset.}
#'   \item{\code{true_Phi}}{Vector of true outcome probabilities for sampled 
#'   individuals. nx1}
#'   \item{\code{true_Phi_mat}}{`NULL` or matrix of true outcome probabilities
#'   for individuals aggregated by C and S. KxH}
#' }
#' @source Simulated data using [simulate_pop()] and [simulate_samp()].
#' @keywords datasets
#' @examples 
#' data(sim_data)
#' x_mat <- sim_data$X_data
"sim_data"


#' NHANES SWOLCA results 
#' 
#' Results from running `swolca()` on the `data_nhanes` dataset
#' @docType data
#' @usage data(run_nhanes_swolca_results)
#' @format List `res` containing:
#' \describe{
#'   \item{\code{estimates_unadj}}{List of unadjusted posterior model results}
#'   \item{\code{runtime}}{Total runtime for model}
#'   \item{\code{data_vars}}{List of data variables used}
#'   \item{\code{MCMC_out}}{List of full MCMC output}
#'   \item{\code{post_MCMC_out}}{List of MCMC output after relabeling}
#'   \item{\code{K_fixed}}{Number of classes used for the fixed sampler}
#'   \item{\code{estimates}}{List of adjusted posterior model results with 
#'   correct uncertainty estimation}
#' }
#' @source Result from running `swolca()` on `data_nhanes`
#' @keywords datasets
#' @examples 
#' data(run_nhanes_swolca_results)
"run_nhanes_swolca_results"