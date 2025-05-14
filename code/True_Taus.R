# Title:      True tau (and generalizations for the discrete case) for a Grid of Dependence Parameters
# Author:     Jan-Lukas Wermuth
# Date:       2025-05-09
# Purpose:    This script simulates the true value of tau as well 
#             as gamma in the discrete cases
#             for every value of alpha, the dependence parameter
#             in the DGPs.

library(mnorm)
library(DescTools)
library(data.table)
library(tsDyn)
library(doParallel)
library(doRNG)
library(foreach)
library(pcaPP)
library(VGAM)
library(here)
# library(devtools)
# install_github("jan-lukas-wermuth/RCor")
library(RCor)

invisible(lapply(list.files(here("code/functions"), pattern = "\\.R$", full.names = TRUE), source)) # Load all the functions

# Parameter Speicification ------------------------------------------------
MC <- 1000
Ti <- 1000
alphas_PowerGraph <- round(head(seq(-1, 1, 0.01), -1)[-1], digits = 2) # get rid of the first (-1) and last (1) element of the alpha-vector
alphas_CIs <- c(sort(-log10(seq(1,9.9999,0.0999))[-1]), log10(seq(1,9.9999,0.0999)))
alphas_fileendings <- c("", "_CIs", "_CIs_short")

# Continuous IID Processes -------------------------------------------------
## "Multivariate t4" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(-0.95, -0.59, 0, 0.59, 0.95)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_Mvt_iid(Ti, alpha, i, 4)
      RCor::RCor(XY[,1], XY[,2], method = "tau", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("taus_t4_iid", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("taus_t4_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("taus_t4_iid", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

## "Multivariate t1" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(-0.985, -0.52, 0, 0.52, 0.985)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_Mvt_iid(Ti, alpha, i, 1)
      RCor::RCor(XY[,1], XY[,2], method = "tau", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("taus_t1_iid", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("taus_t1_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("taus_t1_iid", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

## Normal Exponential DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(-0.95, -0.59, 0, 0.59, 0.95)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_NExp_iid(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "tau", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("taus_NExp_iid", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("taus_NExp_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("taus_NExp_iid", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# Discrete IID Processes -------------------------------------------------
## "Multivariate Poisson" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(0, 0.34, 0.69)
alphas_list <- list(alphas_PowerGraph[100:199], alphas_CIs[91:181], alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvPois_iid(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "gamma", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("gammas_MvPois_iid", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("gammas_MvPois_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("gammas_MvPois_iid", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

## "Multivariate Zipf" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(0, 0.174, 0.677)
alphas_list <- list(alphas_PowerGraph[100:199], alphas_CIs[91:181], alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvZipf_iid(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "gamma", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("gammas_MvZipf_iid", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("gammas_MvZipf_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("gammas_MvZipf_iid", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

## "Multivariate Skellam" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(-0.665, -0.345, 0, 0.345, 0.665)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvSkellam_iid(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "gamma", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("gammas_MvSkellam_iid", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("gammas_MvSkellam_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("gammas_MvSkellam_iid", alphas_fileendings[j]), ".RData"))))
}

# Continuous Time Series Processes -------------------------------------------------
## "Multivariate Normal" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.95, -0.59, 0, 0.59, 0.95)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_Mvt_TS(Ti, alpha, i, Inf)
      RCor::RCor(XY[,1], XY[,2], method = "tau", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("taus_norm_TS", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("taus_norm_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("taus_norm_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# "Multivariate t(4)" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.96, -0.58, 0, 0.58, 0.96)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_Mvt_TS(Ti, alpha, i, 4)
      RCor::RCor(XY[,1], XY[,2], method = "tau", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("taus_t4_TS", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("taus_t4_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("taus_t4_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# "Multivariate t(1)" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.985, -0.51, 0, 0.51, 0.985)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_Mvt_TS(Ti, alpha, i, 1)
      RCor::RCor(XY[,1], XY[,2], method = "tau", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("taus_t1_TS", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("taus_t1_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("taus_t1_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# "Multivariate TEAR(1)" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(0, 0.48, 0.845)
alphas_list <- list(alphas_PowerGraph[100:199], alphas_CIs[91:181], alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_TEAR_TS(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "tau", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("taus_TEAR_TS", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("taus_TEAR_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("taus_TEAR_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# Discrete Time Series Processes -------------------------------------------------
# "Multivariate Poisson" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(0, 0.34, 0.69)
alphas_list <- list(alphas_PowerGraph[100:199], alphas_CIs[91:181], alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvPois_TS(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "gamma", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("gammas_MvPois_TS", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("gammas_MvPois_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("gammas_MvPois_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# "Multivariate Zipf" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(0, 0.46, 0.865)
alphas_list <- list(alphas_PowerGraph[100:199], alphas_CIs[91:181], alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvZipf_TS(Ti, alpha, i, 1000)
      RCor::RCor(XY[,1], XY[,2], method = "gamma", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("gammas_MvZipf_TS", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("gammas_MvZipf_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("gammas_MvZipf_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# "Multivariate Skellam" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.65, -0.33, 0, 0.33, 0.65)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  kendall_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    kendall <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvSkellam_TS(Ti, alpha, i, 1000)
      RCor::RCor(XY[,1], XY[,2], method = "gamma", Inference = FALSE)[[1]]
    }
    kendall_array[, as.character(alpha)] <- kendall
  }
  assign(paste0("gammas_MvSkellam_TS", alphas_fileendings[j]), colMeans(kendall_array))
  save(list = paste0("gammas_MvSkellam_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_taus/", paste0(paste0("gammas_MvSkellam_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

