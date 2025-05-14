# Title:      True rho (and generalizations for the discrete case) for a Grid of Dependence Parameters
# Author:     Jan-Lukas Wermuth
# Date:       2025-05-09
# Purpose:    This script simulates the true value of rho as well 
#             as rhob in the discrete cases
#             for every value of alpha, the dependence parameter
#             in the DGPs (alpha in the paper).

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
alphas_CIs_short <- c(-0.835, -0.4, 0, 0.4, 0.835)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_Mvt_iid(Ti, alpha, i, 4)
      RCor::RCor(XY[,1], XY[,2], method = "rho", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhos_t4_iid", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhos_t4_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhos_t4_iid", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

## "Multivariate t1" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(-0.92, -0.31, 0, 0.31, 0.92)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_Mvt_iid(Ti, alpha, i, 1)
      RCor::RCor(XY[,1], XY[,2], method = "rho", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhos_t1_iid", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhos_t1_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhos_t1_iid", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

## Normal Exponential DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(-0.81, -0.42, 0, 0.42, 0.81)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_NExp_iid(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "rho", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhos_NExp_iid", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhos_NExp_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhos_NExp_iid", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# Discrete IID Processes -------------------------------------------------
## "Multivariate Poisson" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(0, 0.42, 0.81)
alphas_list <- list(alphas_PowerGraph[100:199], alphas_CIs[91:181], alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvPois_iid(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "rho_b", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhobs_MvPois_iid", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhobs_MvPois_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhobs_MvPois_iid", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

## "Multivariate Zipf" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(0, 0.255, 0.815)
alphas_list <- list(alphas_PowerGraph[100:199], alphas_CIs[91:181], alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvZipf_iid(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "rho_b", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhobs_MvZipf_iid", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhobs_MvZipf_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhobs_MvZipf_iid", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

## "Multivariate Skellam" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(-0.728, -0.365, 0, 0.365, 0.728)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvSkellam_iid(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "rho_b", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhobs_MvSkellam_iid", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhobs_MvSkellam_iid", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhobs_MvSkellam_iid", alphas_fileendings[j]), ".RData"))))
}

# Continuous Time Series Processes -------------------------------------------------
## "Multivariate Normal" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.815, -0.42, 0, 0.42, 0.815)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_Mvt_TS(Ti, alpha, i, Inf)
      RCor::RCor(XY[,1], XY[,2], method = "rho", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhos_norm_TS", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhos_norm_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhos_norm_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# "Multivariate t(4)" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.82, -0.41, 0, 0.41, 0.82)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_Mvt_TS(Ti, alpha, i, 4)
      RCor::RCor(XY[,1], XY[,2], method = "rho", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhos_t4_TS", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhos_t4_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhos_t4_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# "Multivariate t(1)" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.915, -0.305, 0, 0.305, 0.915)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, as.character(alphas))) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_Mvt_TS(Ti, alpha, i, 1)
      RCor::RCor(XY[,1], XY[,2], method = "rho", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhos_t1_TS", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhos_t1_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhos_t1_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# "Multivariate TEAR(1)" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(0, 0.4, 0.8)
alphas_list <- list(alphas_PowerGraph[100:199], alphas_CIs[91:181], alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_TEAR_TS(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "rho", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhos_TEAR_TS", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhos_TEAR_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhos_TEAR_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# Discrete Time Series Processes -------------------------------------------------
# "Multivariate Poisson" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(0, 0.42, 0.813)
alphas_list <- list(alphas_PowerGraph[100:199], alphas_CIs[91:181], alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvPois_TS(Ti, alpha, i)
      RCor::RCor(XY[,1], XY[,2], method = "rho_b", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhobs_MvPois_TS", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhobs_MvPois_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhobs_MvPois_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# "Multivariate Zipf" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(0, 0.37, 0.81)
alphas_list <- list(alphas_PowerGraph[100:199], alphas_CIs[91:181], alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvZipf_TS(Ti, alpha, i, 1000)
      RCor::RCor(XY[,1], XY[,2], method = "rho_b", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhobs_MvZipf_TS", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhobs_MvZipf_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhobs_MvZipf_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

# "Multivariate Skellam" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.745, -0.37, 0, 0.37, 0.745)
alphas_list <- list(alphas_PowerGraph, alphas_CIs, alphas_CIs_short)

cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_list)){
  alphas <- alphas_list[[j]]
  spearman_array <- array(data = NA, dim = c(MC, length(alphas)), dimnames = list(1:MC, alphas)) # Initialize results array
  for (alpha in alphas) {
    spearman <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      XY <- Gen_MvSkellam_TS(Ti, alpha, i, 1000)
      RCor::RCor(XY[,1], XY[,2], method = "rho_b", Inference = FALSE)[[1]]
    }
    spearman_array[, as.character(alpha)] <- spearman
  }
  assign(paste0("rhobs_MvSkellam_TS", alphas_fileendings[j]), colMeans(spearman_array))
  save(list = paste0("rhobs_MvSkellam_TS", alphas_fileendings[j]), file = here(paste0("results/simulations/true_rhos/", paste0(paste0("rhobs_MvSkellam_TS", alphas_fileendings[j]), ".RData"))))
}
stopCluster(cl)

