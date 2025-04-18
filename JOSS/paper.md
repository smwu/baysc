---
title: 'baysc: An R package for Bayesian survey clustering'
tags:
  - R
  - Bayesian
  - survey
  - model-based clustering
  - dietary patterns
authors:
  - name: Stephanie M. Wu
    corresponding: true
    orcid: 0000-0001-5110-8407
    affiliation: 1 
  - name: Matthew R. Williams
    orcid: 0000-0001-8894-1240
    affiliation: 2
  - name: Terrance D. Savitsky
    orcid: 0000-0003-1843-3106
    affiliation: 3
  - name: Briana J.K. Stephenson
    orcid: 0000-0002-6147-1039
    affiliation: 4
affiliations:
 - name: Division of Psychiatry, UCL, London, U.K.
   index: 1
 - name: RTI International, Research Triangle Park, North Carolina, U.S.A
   index: 2
 - name: Office of Survey Methods Research, U.S. Bureau of Labor Statistics, Washington, DC, U.S.A
   index: 3
 - name: Department of Biostatistics, Harvard T.H. Chan School of Public Health, Boston, Massachusetts, U.S.A
   index: 4
citation_author: Wu et. al.
date: 8 April 2024
year: 2024
bibliography: baysc.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---

# Summary

Model-based clustering methods allow a large number of correlated variables to be summarized into underlying patterns, where each pattern describes a cluster and each individual is assigned to a cluster. Example applications include identifying dietary patterns from dietary intake data [@stephenson2020empirically] and creating profiles of health and development among children [@lanza2016latent]. Bayesian formulations of such analyses allow the number of clusters to be determined by the data rather than through researcher post-hoc analyses. Clustering methods can also be extended to the hybrid supervised setting where interest lies in the association between the identified clusters and an outcome. When such clustering methods are applied to survey data, failure to account for the complex survey design and incorporate survey weights into the estimation leads to biased estimation and inference when the results are generalized to the population outside of the survey data.

The `baysc` R package provides functionality to allow for Bayesian clustering analyses, both unsupervised and supervised, to be performed while incorporating survey weights and design features that account for complex survey sampling designs. Asymptotically correct point estimates and credible intervals are produced with respect to the underlying population from which the observed sample was generated. This novel feature allows for application of latent class analysis (LCA) to datasets realized from surveys administered by government statistical agencies. The package uses methods derived from the LCA literature and focuses on clustering in the setting where the correlated variables are categorical and the outcome, where applicable, is binary. The package includes additional functions for plotting and summarizing output, and an example dataset from the National Health and Nutrition Examination Survey (NHANES) containing dietary intake and hypertension data among low-income women in the United States [@nchs2023homepage].

# Statement of Need

A number of R packages provide functionality for model-based clustering in R. Frequentist approaches include `poLCA` [@linzer2011polca] for classical LCA, `randomLCA` [@beath2017randomlca] for LCA with individual-specific random effects, and `mclust` [@scrucca2023model] and `tidyLPA` [@rosenberg2019tidylpa] for clustering of continuous variables. `BayesLCA` [@white2014bayeslca] and `BayesBinMix` [@papastamoulis2017bayesbinmix] use Bayesian approaches for categorical and binary data, respectively. `PReMiuM` [@liverani2015premium] fits a wide variety of supervised models that handle various types of discrete and continuous exposure and outcome data. However, these packages do not allow for survey weights and complex survey design to be incorporated to ensure valid estimation and inference when using survey data.

The `baysc` package implements a weighted pseudo-likelihood approach proposed in @wu2024derivation that can integrate sampling weights when creating patterns using categorical data. The models adjust for stratification, clustering, and informative sampling to provide accurate point and variance estimation. When interest lies in how clusters are related to an outcome, supervised methods are available that jointly model the exposure and outcome, capturing the exposure-outcome association with more precision than two-step approaches that perform the clustering and the regression analyses sequentially. In addition, interaction effects between the clusters and the outcome are able to be captured through a mixture reference coding scheme.

# Acknowledgement

This work was supported in part by the National Institute of Allergy and Infectious Diseases (NIAID: T32 AI007358), the National Heart, Lung, and Blood Institute (NHLBI: R25 HL105400), and the Harvard Data Science Initiative $\text{Bias}^2$ Program. The authors declare no conflicts of interest.

# References
