# Title:      DGP Simulation Collection
# Author:     Jan-Lukas Wermuth
# Date:       2025-05-08
# Purpose:    This script collects functions that simulate
#             from a given DGP. It contains all DGPs considered
#             in the paper.

# Mutivariate Normal DGP (iid) -------------------------------------------------
Gen_MvN_iid <- function(Ti, rho, i){
  set.seed(i) 
  XY <- mnorm::rmnorm(Ti, mean = rep(0, 2), sigma = matrix(c(1, rho, rho, 1), ncol=2))
  return(XY)
}

# "Multivariate t" DGP (iid) -------------------------------------------------
Gen_Mvt_iid <- function(Ti, rho, i, df){
  set.seed(i) 
  X <- rt(Ti, df = df)
  Y <- rho * X + sqrt(1 - rho^2) * rt(Ti, df = df)
  return(cbind(X, Y))
}

# Normal Exponential DGP (iid) -------------------------------------------------
Gen_NExp_iid <- function(Ti, rho, i){
  set.seed(i) 
  X <- rnorm(Ti)
  Y <- qexp(pnorm(rho * X + sqrt(1 - rho^2) * rnorm(Ti)))
  return(cbind(X, Y))
}

# "Multivariate Poisson" DGP (iid) -------------------------------------------------
Gen_MvPois_iid <- function(Ti, rho, i){
  set.seed(i) 
  X <- rpois(Ti, lambda = 1)
  Y <- rep(NA, Ti)
  for (t in 1:Ti) {
    Y[t] <- rbinom(1, X[t], rho) + rbinom(1, rpois(1, lambda = 1), 1 - rho)
  }
  return(cbind(X, Y))
}

# "Multivariate Zipf" DGP (iid) -------------------------------------------------
Gen_MvZipf_iid <- function(Ti, rho, i){
  set.seed(i) 
  X <- VGAM::rzeta(Ti, shape = 1) - 1
  Y <- rep(NA, Ti)
  for (t in 1:Ti) {
    Y[t] <- rbinom(1, X[t], rho) + rbinom(1, VGAM::rzeta(1, shape = 1) - 1, 1 - rho)
  }
  return(cbind(X, Y))
}

# "Multivariate Skellam" DGP (iid) -------------------------------------------------
Gen_MvSkellam_iid <- function(Ti, rho, i){
  set.seed(i) 
  X <- rpois(Ti, lambda = 1) - rpois(Ti, lambda = 1)
  Y <- rep(NA, Ti)
  for (t in 1:Ti) {
    u <- rpois(1, lambda = 1) - rpois(1, lambda = 1)
    Y[t] <- sign(rho) * sign(X[t]) * rbinom(1, abs(X[t]), abs(rho)) + ifelse(rho >= 0, 1, -1) * sign(1 - abs(rho)) * sign(u) * rbinom(1, abs(u), abs(1 - abs(rho)))
  }
  return(cbind(X, Y))
}

# "Multivariate Normal and t" DGP (Time Series) -------------------------------------------------
Gen_Mvt_TS <- function(Ti, rho, i, df){
  set.seed(i) 
  X <- arima.sim(model = list(ar = 0.8), n = Ti, rand.gen = function(n) rt(n, df = df))
  Y <- rho * X + sqrt(1 - rho^2) * arima.sim(model = list(ar = 0.8), n = Ti, rand.gen = function(n) rt(n, df = df))
  return(cbind(X, Y))
}

# "Multivariate TEAR(1)" DGP (Time Series) -------------------------------------------------
Gen_TEAR_TS <- function(Ti, rho, i){
  set.seed(i) 
  Xt <- rexp(1)
  eps <- rexp(Ti)
  X <- rep(NA, Ti)
  for(t in 1:Ti){
    Xt <- rbinom(1, 1, 0.8) * Xt + (1 - 0.8) * eps[t]
    X[t] <- Xt
  }
  Ut <- rexp(1)
  nu <- rexp(Ti)
  U <- rep(NA, Ti)
  for(t in 1:Ti){
    Ut <- rbinom(1, 1, 0.8) * Ut + (1 - 0.8) * nu[t]
    U[t] <- Ut
  }
  B_rho <- rbinom(Ti, 1, rho)
  Y <- B_rho * X + (1 - B_rho) * U
  return(cbind(X, Y))
}

# "Multivariate Poisson" DGP (Time Series) -------------------------------------------------
Gen_MvPois_TS <- function(Ti, rho, i){
  set.seed(i) 
  Xt <- rpois(1, 1)
  eps <- rpois(Ti, 1*(1-0.8))
  X <- rep(NA, Ti)
  for(t in 1:Ti){
    Xt <- rbinom(1, Xt, 0.8) + eps[t]
    X[t] <- Xt
  }
  Ut <- rpois(1, 1)
  nu <- rpois(Ti, 1*(1-0.8))
  U <- rep(NA, Ti)
  for(t in 1:Ti){
    Ut <- rbinom(1, Ut, 0.8) + nu[t]
    U[t] <- Ut
  }
  Y <- rep(NA, Ti)
  for (t in 1:Ti) {
    Y[t] <- rbinom(1, X[t], rho) + rbinom(1, U[t], 1 - rho)
  }
  return(cbind(X, Y))
}

# "Multivariate Zipf" DGP (Time Series) -------------------------------------------------
Gen_MvZipf_TS <- function(Ti, rho, i, prerun){
  set.seed(i) 
  Xt <- 1
  eps <- VGAM::rzeta(Ti + prerun, shape = 1.5)
  X <- rep(NA, Ti + prerun)
  for(t in 1:(Ti + prerun)){
    Xt <- rbinom(1, Xt, 0.8) + eps[t]
    X[t] <- Xt
  }
  Ut <- 1
  nu <- VGAM::rzeta(Ti + prerun, shape = 1.5)
  U <- rep(NA, Ti + prerun)
  for(t in 1:(Ti + prerun)){
    Ut <- rbinom(1, Ut, 0.8) + nu[t]
    U[t] <- Ut
  }
  X <- X[1001:(Ti + prerun)]
  U <- U[1001:(Ti + prerun)]
  Y <- rep(NA, Ti)
  for (t in 1:(Ti)) {
    Y[t] <- rbinom(1, X[t], rho) + rbinom(1, U[t], 1 - rho)
  }
  return(cbind(X, Y))
}

# "Multivariate Skellam" DGP (Time Series) -------------------------------------------------
Gen_MvSkellam_TS <- function(Ti, rho, i, prerun){
  set.seed(i) 
  Xt <- 1
  eps <- rpois(Ti + prerun, 0.2) - rpois(Ti + prerun, 0.2)
  X <- rep(NA, Ti + prerun)
  for(t in 1:(Ti + prerun)){
    Xt <- sign(Xt) * rbinom(1, abs(Xt), 0.8) + eps[t]
    X[t] <- Xt
  }
  Ut <- 1
  nu <- rpois(Ti + prerun, 0.2) - rpois(Ti + prerun, 0.2)
  U <- rep(NA, Ti + prerun)
  for(t in 1:(Ti + prerun)){
    Ut <- sign(Ut) * rbinom(1, abs(Ut), 0.8) + nu[t]
    U[t] <- Ut
  }
  X <- X[1001:(Ti + prerun)]
  U <- U[1001:(Ti + prerun)]
  Y <- rep(NA, Ti)
  for (t in 1:Ti) {
    Y[t] <- sign(rho) * sign(X[t]) * rbinom(1, abs(X[t]), abs(rho)) + ifelse(rho >= 0, 1, -1) * sign(1 - abs(rho)) * sign(U[t]) * rbinom(1, abs(U[t]), abs(1 - abs(rho)))
  }
  return(cbind(X, Y))
}



