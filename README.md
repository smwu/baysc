# **baysc**: BAYesian Survey Clustering

An R package for running Bayesian supervised and unsupervised clustering methods on survey data.

**Maintainer**: Stephanie M. Wu ([swu\@g.harvard.edu](mailto:swu@g.harvard.edu){.email})

**Contributors**: Matthew R. Williams ([mrwilliams\@rti.org](mailto:mrwilliams@rti.org){.email}); Terrance D. Savitsky ([savitsky.terrance\@bls.gov](mailto:savitsky.terrance@bls.gov){.email}); Briana J.K. Stephenson ([bstephenson\@hsph.harvard.edu](mailto:bstephenson@hsph.harvard.edu){.email})

**Citation**: Wu S, Williams M, Savitsky T, Stephenson B (2025). *baysc: BAYesian Survey Clustering*. R package version 0.1.0, <https://github.com/smwu/baysc>.

## Table of contents

-   [1. Installation](#id-section1)
-   [2. Overview](#id-section2)
-   [3. Functions](#id-section3)
-   [4. Data](#id-section4)
-   [5. Example](#id-section5)
-   [6. Supplementary](#id-section6)
-   [7. Help](#id-section7)

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

During installation, the following errors may arise:

-   *No package called 'rstantools'*: Please install the `rstantools` package using `install.packages("rstantools")`.
-   *Library 'gfortran' not found*: This is a compiler configuration issue that can arise when using Rcpp on Mac computers with Apple silicon (e.g., M1-M4 chips). Users may need to install Xcode, GNU Fortran, and OpenMP, and edit the `~/.R/Makevars` file. For more details, see the "Supplementary" section below.
-   *Library 'emutls_w' not found*: This is a toolchain mismatch issue that can arise when using `Rcpp`-dependent packages on Mac computers with Apple silicon (e.g., M1-M4 chips). Users may need to install gfortran and edit the `~/.R/Makevars` file. For more details, see the "Supplementary" section below.
</div>

<div id='id-section2'/>
## Overview

`baysc` is an R package for running Bayesian clustering methods on survey data. A Bayesian latent class analysis (LCA), termed the Weighted Overfitted Latent Class Analysis (WOLCA), is available for eliciting underlying cluster patterns from multivariate categorical data, incorporating survey sampling weights and other survey design elements. Options also exist for relating the patterns to a binary outcome, either by using a two-step approach that applies WOLCA and then runs a survey-weighted regression, or by utilizing a one-step supervised approach where creation of the cluster patterns is directly informed by the outcome, referred to as the Supervised Weighted Overfitted Latent Class Analysis (SWOLCA). More information about the models can be found in the paper linked above. Summary and plotting functions for visualizing output are also available, as are diagnostic functions for examining convergence of the sampler.
</div>

<div id='id-section3'/>
## Functions

Use the `wolca()` function to run an unsupervised WOLCA and obtain pattern profiles. `wolca_var_adjust()` provides a post-hoc variance adjustment that enables correct uncertainty estimation. To examine the association of pattern profiles with a binary outcome through a two-step appraoch, run `wolca_svyglm()`. Use the `swolca()` function to run a SWOLCA model that allows information about the binary outcome to directly inform the creation of the pattern profiles. `swolca_var_adjust()` provides a post-hoc variance adjustment that enbales correct uncertainty estimation. Detailed information about the functions and related statistical details can be found in the vignette, "[An introduction to the baysc package](https://raw.githubusercontent.com/smwu/baysc/refs/heads/main/vignettes/baysc.pdf)," and in the paper linked above.
</div>

<div id='id-section4'/>
## Data

`baysc` applies Bayesian latent class analysis using the following input data:

-   Multivariate categorical exposure: $nxJ$ matrix, where $n$ is the sample size and $J$ is the number of categorical item variables. Each item must be a categorical variable.
-   (Optional) binary outcome: $nx1$ vector
-   (Optional) additional confounders to adjust for when evaluating the exposure-outcome association: $nxQ$ dataframe, where $Q$ is the number of additional confounders.
-   (Optional) survey design elements such as stratum indicators, cluster indicators, and sampling weights: each formatted as a $nx1$ vector.

We provide an example dataset from the National Health and Nutrition Examination Survey (NHANES) that includes multivariate categorical dietary intake data as well as binary hypertension data for low-income women in the United States. Survey sampling weights and information on stratification and clustering are included to allow for adjustment for survey design when conducting estimation and inference.
</div>

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
</div>

<div id='id-section6'/>
## Supplementary

### gfortran error

For users experiencing a "**library 'gfortran' not found**" error message during installation, additional steps are needed to install the `baysc` package. This is a compiler configuration issue that can arise when using Rcpp-dependent packages on Mac computers with Apple silicon (e.g., M1, M2, M3). Please follow the instructions listed below, adapted from instructions posted at [https://stackoverflow.com/questions/70638118/configuring-compilers-on-apple-silicon-m1-m2-m3-for-rcpp-and-other-tool)](https://stackoverflow.com/questions/70638118/configuring-compilers-on-apple-silicon-m1-m2-m3-for-rcpp-and-other-tool).

1.  Download an R binary from CRAN at this link: <https://cran.r-project.org/bin/macosx/>. Select the binary built for Apple silicon (M1-M3), which will typically the top link on the left under “Latest release”, and install.
2.  Go through the instructions to install R, entering in passwords if necessary.
3.  Install Xcode by opening Terminal and running the following in the command line: `sudo xcode-select --install` . This will install the latest release version of Apple’s Command Line Tools for Xcode, which includes Apple Clang.
4.  Download the GNU Fortran binary .tar.xz file at this link: <https://github.com/R-macos/gcc-12-branch/releases/tag/12.2-darwin-r0>. Install GNU fortran by running the following code in the command line:

```         
curl -LO https://github.com/R-macos/gcc-12-branch/releases/download/12.2-darwin-r0/gfortran-12.2-darwin20-r0-universal.tar.xz
sudo tar xvf gfortran-12.2-darwin20-r0-universal.tar.xz -C /
sudo ln -sfn $(xcrun --show-sdk-path) /opt/gfortran/SDK
```

5.  Check your Apple Clang version by running the following in the command line: `clang --version`.
6.  Download OpenMP at <https://mac.r-project.org/openmp/> by clicking on the Release.tar.gz link corresponding to your Apple Clang version.
7.  Install OpenMP by running the following in the command line, making sure to include the correct version. For example, for Apple clang version 15.0.0 (i.e., clang 1500.x), the corresponding OpenMP version is *16.0.4*, so the commands would be:

```         
curl -LO https://mac.r-project.org/openmp/openmp-16.0.4-darwin20-Release.tar.gz
sudo mkdir -p /opt/R/$(uname -m)
sudo tar -xvf openmp-16.0.4-darwin20-Release.tar.gz --strip-components=2 -C /opt/R/$(uname -m)
```

8.  Navigate to the R Makevars file by running in command line: `cd ~/.R/Makevars`. If the file does not exist, create it by running `mkdir ~/.R/Makevars`.
9.  In command line, run `vi`, add the below lines to the file, then save and exit out of the file by typing `ESC` followed by `:wq`.

```         
CPPFLAGS += -Xclang -fopenmp
LDFLAGS += -lomp
```

10. Retry the `baysc` package installation by running `devtools::install_github("smwu/baysc")` in R.

### emutls_w error

For users experiencing a "**library 'emutls_w' not found**" error message during installation, additional steps are needed to install the `baysc` package. This is likely due to a toolchain mismatch issue that can arise when using `Rcpp`-dependent packages on Mac computers with Apple silicon (e.g., M1, M2, M3, M4 chips). The issue occurs on these computers because R uses a Clang-based compiler toolchain that can sometimes conflict with gfortran, especially when the compilers are built for different architectures. In particular, the system may fail to locate the `emutls_w` library, which is needed by the Clang compiler to support thread-local storage in multithreaded C++ code. Please follow the instructions listed below to resolve the issue.

1.  Install gfortran by going to this link: <https://mac.r-project.org/tools/>. Under “Mandatory tools” and bullet point “GNU Fortran compiler”, download the latest .dmg (e.g., `gfortran-14.2-arm64.dmg`). Make sure the version you install matches the version that is missing.

2.  Confirm installation of gfortran: running `ls /opt/gfortran` in Terminal should yield folders including “bin”, “lib”, etc.

3.  Find out where `libemutls_w.a` is stored. In Terminal, run:

    `find /opt/gfortran -name "libemutls_w.a"`

    You should see something like `/opt/gfortran/lib/gcc/aarch64-apple-darwin23/14.2.0/libemutls_w.a`.

4.  Explicitly point R to the library where `libemutls_w.a` lives by updating your R Makevars file. In Terminal, navigate to your R Makevars file by running `cd \~/.R/` followed by `vi Makevars`. If the `.R` folder doesn’t exist, run `mkdir \~/.R/`.

5.  Add the below lines to the Makevars file, then save and exit out of the file by typing `ESC` followed by `:wq`.

```         
CC = clang 
CXX = clang++ 
FC = /opt/gfortran/bin/gfortran

CFLAGS = -O2 -arch arm64 
CXXFLAGS = -O2 -arch arm64 
FFLAGS = -O2 -arch arm64

# Linker flags to find emutls_w
LDFLAGS = -L/opt/gfortran/lib/gcc/aarch64-apple-darwin23/14.2.0 -arch arm64
```

In the last line, the text after -L should match the output from Step 3 specifying where the `libemutls_w.a` library is. For example, here, it is `/opt/gfortran/lib/gcc/aarch64-apple-darwin23/14.2.0`.

6.  Restart R and retry the `baysc` package installation.
</div>

<div id='id-section7'/>
## Contributing and Getting Help

Please report bugs by opening an [issue](https://github.com/smwu/baysc/issues/new/choose). If you wish to contribute, please make a pull request. If you have questions, you can open a [discussion thread](https://github.com/smwu/baysc/discussions).
</div>
