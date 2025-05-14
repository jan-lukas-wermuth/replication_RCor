# Title:      Coverage Simulations for Confidence Intervals
# Author:     Jan-Lukas Wermuth
# Date:       2025-05-09
# Purpose:    This script simulates the empirical coverage rates
#             for rho and rho_b for all the DGPs considered in the paper.
rm(list = ls())

library(mnorm)
library(DescTools)
library(data.table)
library(tsDyn)
library(doParallel)
library(doRNG)
library(foreach)
library(here)
# library(devtools)
# install_github("jan-lukas-wermuth/RCor")
library(RCor)

invisible(lapply(list.files(here("code/functions"), pattern = "\\.R$", full.names = TRUE), source)) # Load all the functions
invisible(lapply(list.files(here("results/simulations/true_rhos"), pattern = "\\.RData$", full.names = TRUE), function(x) load(x, envir = globalenv())))

# Parameter Specification ------------------------------------------------
MC <- 1000
SampleSizes <- c(50, 200, 800)
# rhos_CIs <- c(sort(-log10(seq(1,9.9999,0.0999))[-1]), log10(seq(1,9.9999,0.0999)))

# Continuous IID Processes ------------------------------------------------
## Multivariate normal DGP ------------------------------------------------
alphas_CIs_short <- c(-0.81, -0.42, 0, 0.42, 0.81)
decision_spearman_array_norm_iid <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (alpha in alphas_CIs_short){
  for (Ti in SampleSizes){
    decision_spearman_norm_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvN_iid(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = TRUE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(6/pi*asin(alpha/2), spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_norm_iid[as.character(Ti),,as.character(alpha)] <- decision_spearman_norm_iid
  }
}
stopCluster(cl)
save(decision_spearman_array_norm_iid, file = here("results/simulations/coverages/spearman/rho_cov_norm_iid.RData"))

decision_spearman_array_norm_iid_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (alpha in alphas_CIs_short){
  for (Ti in SampleSizes){
    decision_spearman_norm_iid_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvN_iid(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = TRUE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(6/pi*asin(alpha/2), spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_norm_iid_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_norm_iid_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_norm_iid_fis, file = here("results/simulations/coverages/spearman/rho_cov_norm_iid_fis.RData"))

## "Multivariate t(4)" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(-0.835, -0.4, 0, 0.4, 0.835)
decision_spearman_array_t4_iid <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_t4_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_Mvt_iid(Ti, alpha, i, 4)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = TRUE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rhos_t4_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_t4_iid[as.character(Ti),,as.character(alpha)] <- decision_spearman_t4_iid
  }
}
stopCluster(cl)
save(decision_spearman_array_t4_iid, file = here("results/simulations/coverages/spearman/rho_cov_t4_iid.RData"))

decision_spearman_array_t4_iid_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_t4_iid_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_Mvt_iid(Ti, alpha, i, 4)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = TRUE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rhos_t4_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_t4_iid_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_t4_iid_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_t4_iid_fis, file = here("results/simulations/coverages/spearman/rho_cov_t4_iid_fis.RData"))

## "Multivariate t(1)" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(-0.985, -0.52, 0, 0.52, 0.985)
decision_spearman_array_t1_iid <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_t1_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_Mvt_iid(Ti, alpha, i, 1)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = TRUE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rhos_t1_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_t1_iid[as.character(Ti),,as.character(alpha)] <- decision_spearman_t1_iid
  }
}
stopCluster(cl)
save(decision_spearman_array_t1_iid, file = here("results/simulations/coverages/spearman/rho_cov_t1_iid.RData"))

decision_spearman_array_t1_iid_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_t1_iid_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_Mvt_iid(Ti, alpha, i, 1)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = TRUE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rhos_t1_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_t1_iid_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_t1_iid_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_t1_iid_fis, file = here("results/simulations/coverages/spearman/rho_cov_t1_iid_fis.RData"))

## Normal Exponential DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(-0.95, -0.59, 0, 0.59, 0.95)
decision_spearman_array_NExp_iid <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_NExp_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_NExp_iid(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = TRUE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rhos_NExp_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_NExp_iid[as.character(Ti),,as.character(alpha)] <- decision_spearman_NExp_iid
  }
}
stopCluster(cl)
save(decision_spearman_array_NExp_iid, file = here("results/simulations/coverages/spearman/rho_cov_NExp_iid.RData"))

decision_spearman_array_NExp_iid_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_NExp_iid_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_NExp_iid(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = TRUE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rhos_NExp_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_NExp_iid_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_NExp_iid_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_NExp_iid_fis, file = here("results/simulations/coverages/spearman/rho_cov_NExp_iid_fis.RData"))

# Discrete IID Processes ------------------------------------------------
## "Multivariate Poisson" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(0, 0.34, 0.69)
decision_spearman_array_MvPois_iid <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvPois_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvPois_iid(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = TRUE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvPois_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvPois_iid[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvPois_iid
  }
}
stopCluster(cl)
save(decision_spearman_array_MvPois_iid, file = here("results/simulations/coverages/spearman/rho_b_cov_MvPois_iid.RData"))

decision_spearman_array_MvPois_iid_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvPois_iid_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvPois_iid(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = TRUE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvPois_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvPois_iid_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvPois_iid_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_MvPois_iid_fis, file = here("results/simulations/coverages/spearman/rho_b_cov_MvPois_iid_fis.RData"))

## "Multivariate Zipf" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(0, 0.174, 0.677)
decision_spearman_array_MvZipf_iid <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvZipf_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvZipf_iid(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = TRUE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvZipf_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvZipf_iid[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvZipf_iid
  }
}
stopCluster(cl)
save(decision_spearman_array_MvZipf_iid, file = here("results/simulations/coverages/spearman/rho_b_cov_MvZipf_iid.RData"))

decision_spearman_array_MvZipf_iid_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvZipf_iid_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvZipf_iid(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = TRUE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvZipf_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvZipf_iid_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvZipf_iid_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_MvZipf_iid_fis, file = here("results/simulations/coverages/spearman/rho_b_cov_MvZipf_iid_fis.RData"))

## "Multivariate Skellam" DGP (iid) -------------------------------------------------
alphas_CIs_short <- c(-0.665, -0.345, 0, 0.345, 0.665)
decision_spearman_array_MvSkellam_iid <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvSkellam_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvSkellam_iid(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = TRUE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvSkellam_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvSkellam_iid[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvSkellam_iid
  }
}
stopCluster(cl)
save(decision_spearman_array_MvSkellam_iid, file = here("results/simulations/coverages/spearman/rho_b_cov_MvSkellam_iid.RData"))

decision_spearman_array_MvSkellam_iid_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvSkellam_iid_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvSkellam_iid(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = TRUE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvSkellam_iid_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvSkellam_iid_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvSkellam_iid_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_MvSkellam_iid_fis, file = here("results/simulations/coverages/spearman/rho_b_cov_MvSkellam_iid_fis.RData"))

# Continuous Time Series Processes ------------------------------------------------
## "Multivariate Normal" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.95, -0.59, 0, 0.59, 0.95)
decision_spearman_array_norm_TS <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_norm_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_Mvt_TS(Ti, alpha, i, Inf)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = FALSE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rhos_norm_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_norm_TS[as.character(Ti),,as.character(alpha)] <- decision_spearman_norm_TS
  }
}
stopCluster(cl)
save(decision_spearman_array_norm_TS, file = here("results/simulations/coverages/spearman/rho_cov_norm_TS.RData"))

decision_spearman_array_norm_TS_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_norm_TS_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_Mvt_TS(Ti, alpha, i, Inf)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = FALSE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rhos_norm_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_norm_TS_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_norm_TS_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_norm_TS_fis, file = here("results/simulations/coverages/spearman/rho_cov_norm_TS_fis.RData"))

## "Multivariate t(4)" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.96, -0.58, 0, 0.58, 0.96)
decision_spearman_array_t4_TS <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_t4_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_Mvt_TS(Ti, alpha, i, 4)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = FALSE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rhos_t4_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_t4_TS[as.character(Ti),,as.character(alpha)] <- decision_spearman_t4_TS
  }
}
stopCluster(cl)
save(decision_spearman_array_t4_TS, file = here("results/simulations/coverages/spearman/rho_cov_t4_TS.RData"))

decision_spearman_array_t4_TS_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_t4_TS_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_Mvt_TS(Ti, alpha, i, 4)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = FALSE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rhos_t4_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_t4_TS_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_t4_TS_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_t4_TS_fis, file = here("results/simulations/coverages/spearman/rho_cov_t4_TS_fis.RData"))

## "Multivariate t(1)" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.985, -0.51, 0, 0.51, 0.985)
decision_spearman_array_t1_TS <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_t1_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_Mvt_TS(Ti, alpha, i, 1)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = FALSE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rhos_t1_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_t1_TS[as.character(Ti),,as.character(alpha)] <- decision_spearman_t1_TS
  }
}
stopCluster(cl)
save(decision_spearman_array_t1_TS, file = here("results/simulations/coverages/spearman/rho_cov_t1_TS.RData"))

decision_spearman_array_t1_TS_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_t1_TS_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_Mvt_TS(Ti, alpha, i, 1)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = FALSE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rhos_t1_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_t1_TS_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_t1_TS_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_t1_TS_fis, file = here("results/simulations/coverages/spearman/rho_cov_t1_TS_fis.RData"))

## "Multivariate TEAR(1)" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(0, 0.48, 0.845)
decision_spearman_array_TEAR_TS <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_TEAR_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_TEAR_TS(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = FALSE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rhos_TEAR_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_TEAR_TS[as.character(Ti),,as.character(alpha)] <- decision_spearman_TEAR_TS
  }
}
stopCluster(cl)
save(decision_spearman_array_TEAR_TS, file = here("results/simulations/coverages/spearman/rho_cov_TEAR_TS.RData"))

decision_spearman_array_TEAR_TS_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_TEAR_TS_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_TEAR_TS(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho", IID = FALSE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rhos_TEAR_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_TEAR_TS_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_TEAR_TS_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_TEAR_TS_fis, file = here("results/simulations/coverages/spearman/rho_cov_TEAR_TS_fis.RData"))

# Discrete Time Series Processes ------------------------------------------------
## "Multivariate Poisson" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(0, 0.34, 0.69)
decision_spearman_array_MvPois_TS <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvPois_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvPois_TS(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = FALSE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvPois_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvPois_TS[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvPois_TS
  }
}
stopCluster(cl)
save(decision_spearman_array_MvPois_TS, file = here("results/simulations/coverages/spearman/rho_b_cov_MvPois_TS.RData"))

decision_spearman_array_MvPois_TS_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvPois_TS_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvPois_TS(Ti, alpha, i)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = FALSE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvPois_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvPois_TS_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvPois_TS_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_MvPois_TS_fis, file = here("results/simulations/coverages/spearman/rho_b_cov_MvPois_TS_fis.RData"))

## "Multivariate Zipf" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(0, 0.46, 0.865)
decision_spearman_array_MvZipf_TS <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvZipf_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvZipf_TS(Ti, alpha, i, 1000)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = FALSE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvZipf_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvZipf_TS[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvZipf_TS
  }
}
stopCluster(cl)
save(decision_spearman_array_MvZipf_TS, file = here("results/simulations/coverages/spearman/rho_b_cov_MvZipf_TS.RData"))

decision_spearman_array_MvZipf_TS_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvZipf_TS_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvZipf_TS(Ti, alpha, i, 1000)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = FALSE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvZipf_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvZipf_TS_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvZipf_TS_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_MvZipf_TS_fis, file = here("results/simulations/coverages/spearman/rho_b_cov_MvZipf_TS_fis.RData"))

## "Multivariate Skellam" DGP (Time Series) -------------------------------------------------
alphas_CIs_short <- c(-0.65, -0.33, 0, 0.33, 0.65)
decision_spearman_array_MvSkellam_TS <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvSkellam_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvSkellam_TS(Ti, alpha, i, 1000)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = FALSE, Fisher = FALSE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvSkellam_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvSkellam_TS[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvSkellam_TS
  }
}
stopCluster(cl)
save(decision_spearman_array_MvSkellam_TS, file = here("results/simulations/coverages/spearman/rho_b_cov_MvSkellam_TS.RData"))

decision_spearman_array_MvSkellam_TS_fis <- array(data = NA, dim = c(length(SampleSizes), MC, length(alphas_CIs_short)), dimnames = list(SampleSizes, 1:MC, alphas_CIs_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (j in seq_along(alphas_CIs_short)){
  alpha <- alphas_CIs_short[j]
  for (Ti in SampleSizes){
    decision_spearman_MvSkellam_TS_fis <- foreach(i = 1:MC, .combine = 'c') %dopar% {
      set.seed(i)
      XY <- Gen_MvSkellam_TS(Ti, alpha, i, 1000)
      spearman <- RCor::RCor(XY[,1], XY[,2], method = "rho_b", IID = FALSE, Fisher = TRUE, Inference = TRUE)
      as.numeric(data.table::between(rho_bs_MvSkellam_TS_CIs_short[j], spearman[[2]], spearman[[3]]))
    }
    decision_spearman_array_MvSkellam_TS_fis[as.character(Ti),,as.character(alpha)] <- decision_spearman_MvSkellam_TS_fis
  }
}
stopCluster(cl)
save(decision_spearman_array_MvSkellam_TS_fis, file = here("results/simulations/coverages/spearman/rho_b_cov_MvSkellam_TS_fis.RData"))
