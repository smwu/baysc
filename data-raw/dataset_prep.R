
#================= Code to prepare `data_nhanes` dataset =======================
library(dplyr)
# Download 'nhanes1518_adult_low_f_12jul2023.csv' file from 
# https://github.com/smwu/SWOLCA. 

# Set working directory to where the downloaded file is, and read in the dataset
data_nhanes <- read.csv("./nhanes1518_adult_low_f_12jul2023.csv")

# Clean variables
data_nhanes <- data_nhanes %>% 
  # Drop legumes (vegs) because duplicate of legumes (proteins), just computed
  # as cup eq vs oz eq
  select((-leg_veg)) %>%
  # Drop those with refused or unknown education data
  filter(!DMDEDUC2 %in% c(7, 0)) %>%
  mutate(
    # Get stratum IDs
    stratum_id = SDMVSTRA,
    # Create unique nested PSU IDs using stratum IDs
    cluster_id = stratum_id * 10 + SDMVPSU,
    # Get sampling weights
    sample_wt = dietwt4yr,
    # Categorize age
    age_cat = factor(case_when(  
      20 <= RIDAGEYR & RIDAGEYR <= 39 ~ 1,
      40 <= RIDAGEYR & RIDAGEYR <= 59 ~ 2,
      RIDAGEYR >= 60 ~ 3)),
    # Re-categorize race and Hispanic origin ethnicity
    # Original RIDRETH3 categories: 1=Mex_Amer, 2=Other_Hisp, 3=NH_White, 
    # 4=NH_Black, 6=NH_Asian, 7=Other/Mixed
    # New categories: 1=NH White, 2=NH Black, 3=NH Asian, 4=Hispanic, 5=Other/Mixed
    racethnic = factor(case_when(
      RIDRETH3 == 3 ~ 1,  # NH White
      RIDRETH3 == 4 ~ 2,  # NH Black
      RIDRETH3 == 6 ~ 3,  # NH Asian
      RIDRETH3 %in% c(1, 2) ~ 4,  # Mexican-American/Other Hispanic
      RIDRETH3 == 7 ~ 5,  # Other/Mixed
          # RIDRETH3 == 1 ~ 4,  # Mexican-American
          # RIDRETH3 == 2 ~ 5,  # Other Hispanic
          # RIDRETH3 == 7 ~ 6,  # Other/Mixed
      .default = NA)),
    # Re-categorize education level
    # Original DMDEDUC2 categories: 1= <9th, 2=9-11th, 3=HS/GED,
    # 4=Some college/AA, 5=college grad or above, 7=refused, 9=don't know
    # New DMDEDUC2 categories: 1=at least some college, 2=HS/GED, 3=less than HS
    educ = factor(case_when(
      DMDEDUC2 %in% c(4, 5) ~ 1,  # At least some college
      DMDEDUC2 == 3 ~ 2,  # HS/GED
      DMDEDUC2 %in% c(1, 2) ~ 3,  # Less than HS
      .default = NA)),
    # According to CDC Tobacco Glossary, define current smoker as having smoked 
    # 100 cigarettes in lifetime and currently smokes cigarettes
    # SMQ020: Smoked at least 100 cigarettes in life
    # SMQ040 (nested): Now smoke cigarettes
    smoker = factor(case_when(
      SMQ040 %in% c(1, 2) ~ 1,  # current smoker
      SMQ020 == 1 ~ 0,  # did not smoke 100 cigs in life
      SMQ020 == 1 & SMQ040 == 3 ~ 0, # former smoker
      .default = 0), levels = c(0, 1)),
    # Physically active = 1 if >=150 mins of moderate or vigorous exercise per week
    physactive = factor(
      case_when(
        Mins_Active >= 150 ~ 1,
        Mins_Active < 150 ~ 0,
        .default = NA), levels=c(0, 1), labels=c("Inactive", "Active")),
    # BPQ020: Ever told you had high blood pressure
    # BPQ030: Told had high blood pressure 2+ times
    # BPQ050A: Now taking prescribed medicine for HBP
    # Assign hypertension to those currently taking blood pressure medication
    # or those told they had hypertension 2+ times
    BP_flag = case_when(
      SBP_avg > 130 ~ 1,  # Systolic BP above 130
      DBP_avg > 80 ~ 1,  # Diagnostic BP above 80
      BPQ030 == 1 ~ 1,  # Told high BP 2+ times
      BPQ050A == 1 ~ 1,  # Taking BP medication
      BPQ020 == 2 ~ 0,  # Never told had high blood pressure
      .default = 0
    ))

# Keep only necessary variables
data_nhanes <- data_nhanes %>%
  select(SEQN, stratum_id, cluster_id, sample_wt, age_cat, racethnic, educ, 
         smoker, physactive, BP_flag, citrus:drinks, cycle, nrecall, SDMVSTRA, 
         SDMVPSU, RIDAGEYR, RIDRETH3, DMDEDUC2, SMQ020, SMQ040, 
         Mins_Active, SBP_avg, DBP_avg, BPQ030, BPQ050A, BPQ020)

usethis::use_data(data_nhanes, overwrite = TRUE)


#========= Code to prepare `run_nhanes_swolca_results` dataset =================
data("data_nhanes")
x_mat <- as.matrix(dplyr::select(data_nhanes, citrus:drinks))
y_all <- data_nhanes$BP_flag
stratum_id <- data_nhanes$stratum_id
cluster_id <- data_nhanes$cluster_id
sampling_wt <- data_nhanes$sample_wt
V_data <- dplyr::select(data_nhanes, age_cat, racethnic, smoker, physactive)
glm_form <- "~ age_cat + racethnic + smoker + physactive"

class_cutoff <- 0.05
K_max <- 30
alpha_adapt <- eta_adapt <- mu0_adapt <- Sig0_adapt <- NULL
alpha_fixed <- eta_fixed <- mu0_fixed <- Sig0_fixed <- K_fixed <- fixed_seed <- NULL
run_sampler <- "both"
adapt_seed = 20230225
n_runs = 20000
burn = 10000
thin = 5
save_res = TRUE
save_path = "~/data/"

library(devtools)
load_all()

# runtime 2.878507 hours
res_nhanes <- swolca(x_mat = x_mat, y_all = y_all, sampling_wt = sampling_wt,
                     cluster_id = cluster_id, stratum_id = stratum_id, 
                     V_data = V_data, run_sampler = "both", 
                     glm_form = glm_form, adapt_seed = 20230225, 
                     n_runs = 20000, burn = 19800, thin = 5, save_res = TRUE,
                     save_path = "~/Downloads/run_nhanes")

# load("~/Downloads/run_nhanes_swolca_adapt.RData")
# run_nhanes_swolca_adapt <- res
# usethis::use_data(run_nhanes_swolca_adapt, overwrite = TRUE, compress = "xz")

load("~/Downloads/run_nhanes_swolca_results.RData")
run_nhanes_swolca_results <- res
usethis::use_data(run_nhanes_swolca_results, overwrite = TRUE, compress = "xz")



#================= Code to prepare `sim_data` dataset ==========================
# Create population with default parameters
sim_pop <- simulate_pop(save_res = FALSE)
# Create sample with default parameters
sim_samp <- simulate_samp(sim_pop = sim_pop, save_res = FALSE)
sim_data <- sim_samp
usethis::use_data(sim_data, overwrite = TRUE)
