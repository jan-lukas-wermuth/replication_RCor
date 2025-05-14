# Title:      Size and Power Simulations
# Author:     Jan-Lukas Wermuth
# Date:       2025-05-09
# Purpose:    This script simulates the size and power values
#             for tau and gamma for all the DGPs considered in the paper.
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

# Parameter Specification ------------------------------------------------
MC <- 1000
SampleSizes <- c(50, 200, 800)
coefficients_cont <- c("r", "tau", "rho")
coefficients_disc <- c("r", "tau", "tau_b", "gamma", "rho_b")

# Continuous IID Processes ------------------------------------------------
## Multivariate normal DGP ------------------------------------------------
alphas_Test_short <- c(-0.16, -0.08, 0, 0.08, 0.16)
pval_array_norm_iid <- array(data = NA, dim = c(length(coefficients_cont), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_cont, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_cont){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_norm_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_MvN_iid(Ti, alpha, i)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = TRUE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_norm_iid[coef,as.character(Ti),,as.character(alpha)] <- pval_norm_iid
    }
  }
}
stopCluster(cl)
save(pval_array_norm_iid, file = here("results/simulations/size_power/pval_norm_iid.RData"))

## "Multivariate t(4)" DGP (iid) -------------------------------------------------
alphas_Test_short <- c(-0.14, -0.07, 0, 0.07, 0.14)
pval_array_t4_iid <- array(data = NA, dim = c(length(coefficients_cont), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_cont, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_cont){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_t4_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_Mvt_iid(Ti, alpha, i, 4)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = TRUE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_t4_iid[coef,as.character(Ti),,as.character(alpha)] <- pval_t4_iid
    }
  }
}
stopCluster(cl)
save(pval_array_t4_iid, file = here("results/simulations/size_power/pval_t4_iid.RData"))

## "Multivariate t(1)" DGP (iid) -------------------------------------------------
alphas_Test_short <- c(-0.07, -0.03, 0, 0.03, 0.07)
pval_array_t1_iid <- array(data = NA, dim = c(length(coefficients_cont), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_cont, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_cont){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_t1_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_Mvt_iid(Ti, alpha, i, 1)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = TRUE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_t1_iid[coef,as.character(Ti),,as.character(alpha)] <- pval_t1_iid
    }
  }
}
stopCluster(cl)
save(pval_array_t1_iid, file = here("results/simulations/size_power/pval_t1_iid.RData"))

## Normal Exponential DGP (iid) -------------------------------------------------
alphas_Test_short <- c(-0.16, -0.08, 0, 0.08, 0.16)
pval_array_Nexp_iid <- array(data = NA, dim = c(length(coefficients_cont), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_cont, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_cont){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_Nexp_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_NExp_iid(Ti, alpha, i)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = TRUE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_Nexp_iid[coef,as.character(Ti),,as.character(alpha)] <- pval_Nexp_iid
    }
  }
}
stopCluster(cl)
save(pval_array_Nexp_iid, file = here("results/simulations/size_power/pval_Nexp_iid.RData"))

# Discrete IID Processes ------------------------------------------------
## "Multivariate Poisson" DGP (iid) -------------------------------------------------
alphas_Test_short <- c(0, 0.04, 0.09)
pval_array_MvPois_iid <- array(data = NA, dim = c(length(coefficients_disc), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_disc, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_disc){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_MvPois_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_MvPois_iid(Ti, alpha, i)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = TRUE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_MvPois_iid[coef,as.character(Ti),,as.character(alpha)] <- pval_MvPois_iid
    }
  }
}
stopCluster(cl)
save(pval_array_MvPois_iid, file = here("results/simulations/size_power/pval_MvPois_iid.RData"))

## "Multivariate Zipf" DGP (iid) -------------------------------------------------
alphas_Test_short <- c(0, 0.008, 0.019)
pval_array_MvZipf_iid <- array(data = NA, dim = c(length(coefficients_disc), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_disc, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_disc){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_MvZipf_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_MvZipf_iid(Ti, alpha, i)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = TRUE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_MvZipf_iid[coef,as.character(Ti),,as.character(alpha)] <- pval_MvZipf_iid
    }
  }
}
stopCluster(cl)
save(pval_array_MvZipf_iid, file = here("results/simulations/size_power/pval_MvZipf_iid.RData"))

## "Multivariate Skellam" DGP (iid) -------------------------------------------------
alphas_Test_short <- c(-0.11, -0.055, 0, 0.055, 0.11)
pval_array_MvSkellam_iid <- array(data = NA, dim = c(length(coefficients_disc), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_disc, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_disc){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_MvSkellam_iid <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_MvSkellam_iid(Ti, alpha, i)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = TRUE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_MvSkellam_iid[coef,as.character(Ti),,as.character(alpha)] <- pval_MvSkellam_iid
    }
  }
}
stopCluster(cl)
save(pval_array_MvSkellam_iid, file = here("results/simulations/size_power/pval_MvSkellam_iid.RData"))

# Continuous Time Series Processes ------------------------------------------------
## "Multivariate Normal" DGP (Time Series) -------------------------------------------------
alphas_Test_short <- c(-0.16, -0.08, 0, 0.08, 0.16)
pval_array_MvNorm_TS <- array(data = NA, dim = c(length(coefficients_cont), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_cont, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_cont){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_MvNorm_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_Mvt_TS(Ti, alpha, i, Inf)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = FALSE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_MvNorm_TS[coef,as.character(Ti),,as.character(alpha)] <- pval_MvNorm_TS
    }
  }
}
stopCluster(cl)
save(pval_array_MvNorm_TS, file = here("results/simulations/size_power/pval_MvNorm_TS.RData"))

## "Multivariate t(4)" DGP (Time Series) -------------------------------------------------
alphas_Test_short <- c(-0.15, -0.08, 0, 0.08, 0.15)
pval_array_t4_TS <- array(data = NA, dim = c(length(coefficients_cont), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_cont, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_cont){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_t4_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_Mvt_TS(Ti, alpha, i, 4)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = FALSE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_t4_TS[coef,as.character(Ti),,as.character(alpha)] <- pval_t4_TS
    }
  }
}
stopCluster(cl)
save(pval_array_t4_TS, file = here("results/simulations/size_power/pval_Mvt4_TS.RData"))

## "Multivariate t(1)" DGP (Time Series) -------------------------------------------------
alphas_Test_short <- c(-0.06, -0.03, 0, 0.03, 0.06)
pval_array_t1_TS <- array(data = NA, dim = c(length(coefficients_cont), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_cont, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_cont){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_t1_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_Mvt_TS(Ti, alpha, i, 1)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = FALSE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_t1_TS[coef,as.character(Ti),,as.character(alpha)] <- pval_t1_TS
    }
  }
}
stopCluster(cl)
save(pval_array_t1_TS, file = here("results/simulations/size_power/pval_Mvt1_TS.RData"))

## "Multivariate TEAR(1)" DGP (Time Series) -------------------------------------------------
alphas_Test_short <- c(0, 0.07, 0.14)
pval_array_TEAR_TS <- array(data = NA, dim = c(length(coefficients_cont), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_cont, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_cont){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_TEAR_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_TEAR_TS(Ti, alpha, i)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = FALSE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_TEAR_TS[coef,as.character(Ti),,as.character(alpha)] <- pval_TEAR_TS
    }
  }
}
stopCluster(cl)
save(pval_array_TEAR_TS, file = here("results/simulations/size_power/pval_TEAR_TS.RData"))

# Discrete Time Series Processes ------------------------------------------------
## "Multivariate Poisson" DGP (Time Series) -------------------------------------------------
alphas_Test_short <- c(0, 0.04, 0.085)
pval_array_MvPois_TS <- array(data = NA, dim = c(length(coefficients_disc), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_disc, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_disc){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_MvPois_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_MvPois_TS(Ti, alpha, i)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = FALSE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_MvPois_TS[coef,as.character(Ti),,as.character(alpha)] <- pval_MvPois_TS
    }
  }
}
stopCluster(cl)
save(pval_array_MvPois_TS, file = here("results/simulations/size_power/pval_MvPois_TS.RData"))

## "Multivariate Zipf" DGP (Time Series) -------------------------------------------------
alphas_Test_short <- c(0, 0.047, 0.098)
pval_array_MvZipf_TS <- array(data = NA, dim = c(length(coefficients_disc), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_disc, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_disc){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_MvZipf_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_MvZipf_TS(Ti, alpha, i)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = FALSE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_MvZipf_TS[coef,as.character(Ti),,as.character(alpha)] <- pval_MvZipf_TS
    }
  }
}
stopCluster(cl)
save(pval_array_MvZipf_TS, file = here("results/simulations/size_power/pval_MvZipf_TS.RData"))

## "Multivariate Skellam" DGP (Time Series) -------------------------------------------------
alphas_Test_short <- c(-0.09, -0.045, 0, 0.045, 0.09)
pval_array_MvSkellam_TS <- array(data = NA, dim = c(length(coefficients_disc), length(SampleSizes), MC, length(alphas_Test_short)), dimnames = list(coefficients_disc, SampleSizes, 1:MC, alphas_Test_short)) # Initialize results array
cl <- makeCluster(detectCores() - 1, type = "PSOCK")
registerDoParallel(cl)
for (coef in coefficients_disc){
  for (alpha in alphas_Test_short){
    for (Ti in SampleSizes){
      pval_MvSkellam_TS <- foreach(i = 1:MC, .combine = 'c') %dopar% {
        set.seed(i)
        XY <- Gen_MvSkellam_TS(Ti, alpha, i)
        RCor::RCor(XY[,1], XY[,2], method = coef, IID = FALSE, Fisher = FALSE, Inference = TRUE)[[4]]
      }
      pval_array_MvSkellam_TS[coef,as.character(Ti),,as.character(alpha)] <- pval_MvSkellam_TS
    }
  }
}
stopCluster(cl)
save(pval_array_MvSkellam_TS, file = here("results/simulations/size_power/pval_MvSkellam_TS.RData"))
