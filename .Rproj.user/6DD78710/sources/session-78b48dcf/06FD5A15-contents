# ============================================================
# Title:      Coverage Simulations for Confidence Intervals
# Author:     Jan-Lukas Wermuth
# Date:       2025-04-29
# Purpose:    This script simulates the empirical coverage rates
#             for all the DGPs considered in the paper.
# ============================================================
rm(list = ls())

library(devtools)
library(arrangements)
library(rstatix)
library(foreach)
library(readxl)
library(doParallel)
library(DescTools)
library(compositions)
library(mvtnorm)

# install_github("jan-lukas-wermuth/NCor")
library(NCor)
setwd("~Dropbox/Pohle Wermuth/NominalCorrelation/replication_NCor")
invisible(lapply(list.files("results/Simulations/True_gammas", pattern = "\\.RData$", full.names = TRUE), function(x) load(x, envir = globalenv())))
results_folder <- "results/Simulations/Coverage"

# Parameter Specifications ------------------------------------------------
MC <- 1000
SampleSizes <- c(50, 200, 800)
categories <- c("A", "B", "C")

# Regression Normal -------------------------------------------------------
rhos <- (seq(0, 1, length.out = 100))^2 * 4
decision_wermuth_array <- array(data = NA, dim = c(length(SampleSizes), MC, length(rhos)), dimnames = list(SampleSizes, 1:MC, rhos))

for (rho in rhos){
  for (n in SampleSizes){
    for (i in 1:MC){
      print(c(rho, n, i))
      XY <- Gen_RegressionDGP(n, rho, i, Inf) # Generates data according to specified DGP
      res <- NCor(XY[,1], XY[,2], nominal = "r", CIs = TRUE, Test = FALSE)[[1]]
      decision_wermuth_array[as.character(n),as.character(i),as.character(rho)] <- as.numeric(data.table::between(gammas_RegNorm[as.character(rho)], res[2], res[3]))
    }
  }
}

save(decision_wermuth_array, file = paste(results_folder, "RegNormal_CIs.RData", sep = "/"))

# Regression Cauchy -------------------------------------------------------
rhos <- (seq(0, 1, length.out = 100))^3 * 40
decision_wermuth_array <- array(data = NA, dim = c(length(SampleSizes), MC, length(rhos)), dimnames = list(SampleSizes, 1:MC, rhos))

for (rho in rhos){
  for (n in SampleSizes){
    for (i in 1:MC){
      print(c(rho, n, i))
      XY <- Gen_RegressionDGP(n, rho, i, 1) # Generates data according to specified DGP
      res <- NCor(XY[,1], XY[,2], nominal = "r", CIs = TRUE, Test = FALSE)[[1]]
      decision_wermuth_array[as.character(n),as.character(i),as.character(rho)] <- as.numeric(data.table::between(gammas_RegCauchy[as.character(rho)], res[2], res[3]))
    }
  }
} 

save(decision_wermuth_array, file = paste(results_folder, "RegCauchy_CIs.RData", sep = "/"))


# Multinomial Logit Normal ------------------------------------------------------
rhos <- (seq(0, 1, length.out = 100))^2 * 9
decision_wermuth_array <- array(data = NA, dim = c(length(SampleSizes), MC, length(rhos)), dimnames = list(SampleSizes, 1:MC, rhos))

for (rho in rhos){
  for (n in SampleSizes){
    for (i in 1:MC){
      print(c(rho, n, i))
      XY <- Gen_MultinomialDGP(n, rho, i, Inf) # Generates data according to specified DGP
      res <- NCor(XY[,2], XY[,1], nominal = "r", CIs = TRUE, Test = FALSE)[[1]]
      decision_wermuth_array[as.character(n),as.character(i),as.character(rho)] <- as.numeric(data.table::between(gammas_MultNorm[as.character(rho)], res[2], res[3]))
    }
  }
} 

save(decision_wermuth_array, file = paste(results_folder, "MultNormal_CIs.RData", sep = "/"))

# Multinomial Logit Cauchy ------------------------------------------------------
rhos <- (seq(0, 1, length.out = 100))^2 * 9
decision_wermuth_array <- array(data = NA, dim = c(length(SampleSizes), MC, length(rhos)), dimnames = list(SampleSizes, 1:MC, rhos))

for (rho in rhos){
  for (n in SampleSizes){
    for (i in 1:MC){
      print(c(rho, n, i))
      XY <- Gen_MultinomialDGP(n, rho, i, 1) # Generates data according to specified DGP
      res <- NCor(XY[,2], XY[,1], nominal = "r", CIs = TRUE, Test = FALSE)[[1]]
      decision_wermuth_array[as.character(n),as.character(i),as.character(rho)] <- as.numeric(data.table::between(gammas_MultCauchy[as.character(rho)], res[2], res[3]))
    }
  }
}

save(decision_wermuth_array, file = paste(results_folder, "MultCauchy_CIs.RData", sep = "/"))


# 3 x 3 Skewed Uniform ------------------------------------------------------
rhos <- (seq(0, 1, length.out = 100))^2 * 0.04
decision_wermuth_array <- array(data = NA, dim = c(length(SampleSizes), MC, length(rhos)), dimnames = list(SampleSizes, 1:MC, rhos))

for (rho in rhos){
  probabilities <- matrix(c(76/300 + 2*rho, 76/300 + 2*rho, 76/300 - 4*rho, 4/100 - rho, 4/100 - rho, 4/100 + 2*rho, 4/100 - rho, 4/100 - rho, 4/100 + 2*rho), nrow = 3)
  prob_vector <- as.vector(probabilities)
  for (n in SampleSizes){
    for (i in 1:MC){
      print(c(rho, n, i))
      contingency_matrix <- Gen_3x3DGP(n, i) # Generates data according to specified DGP
      res <- NCor(contingency_matrix, nominal = "rc", CIs = TRUE, Test = FALSE)[[1]]
      decision_wermuth_array[as.character(n),as.character(i),as.character(rho)] <- as.numeric(data.table::between(gammas_3x3SU[as.character(rho)], res[2], res[3]))
    }
  }
} 

save(decision_wermuth_array, file = paste(results_folder, "3x3SU_CIs.RData", sep = "/"))

# 3 x 3 Uniform Uniform ------------------------------------------------------
rhos <- (seq(0, 1, length.out = 100))^2 / 36
decision_wermuth_array <- array(data = NA, dim = c(length(SampleSizes), MC, length(rhos)), dimnames = list(SampleSizes, 1:MC, rhos)) # Initialize results array

for (rho in rhos){
  probabilities <- matrix(c(1/9 + 2*rho, 1/9 + 2*rho, 1/9 - 4*rho, 1/9 - rho, 1/9 - rho, 1/9 + 2*rho, 1/9 - rho, 1/9 - rho, 1/9 + 2*rho), nrow = 3)
  prob_vector <- as.vector(probabilities)
  for (n in SampleSizes){
    for (i in 1:MC){
      print(c(rho, n, i))
      contingency_matrix <- Gen_3x3DGP(n, i) # Generates data according to specified DGP
      res <- NCor(contingency_matrix, nominal = "rc", CIs = TRUE, Test = FALSE)[[1]]
      decision_wermuth_array[as.character(n),as.character(i),as.character(rho)] <- as.numeric(data.table::between(gammas_3x3UU[as.character(rho)], res[2], res[3]))
    }
  }
} 

save(decision_wermuth_array, file = paste(results_folder, "3x3UU_CIs.RData", sep = "/"))

