########################## Tau Independent Long-Run Variance Estimation #############################
Tau_ind_LRV <- function(X, Y, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)
  
  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")
  
  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  
  # Calculate autocovariances in a vector with row = lag
  x_autoc <- (n - 1) / n * acf((rank(X) - 0.5) / n - 0.5, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # This is the acf of the demeaned grade. Therefore, demean = FALSE
  y_autoc <- (n - 1) / n * acf((rank(Y) - 0.5) / n - 0.5, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # This is the acf of the demeaned grade. Therefore, demean = FALSE
  
  # Calculate estimator of LRV for tau under independence
  Tau_ind_LRV <- 64 * sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))
  
  return(Tau_ind_LRV)
}

########################## Rho Independent Long-Run Variance Estimation #############################
Rho_ind_LRV <- function(X, Y, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)
  
  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")
  
  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  
  # Calculate autocovariances in a vector with row = lag
  x_autoc <- acf(X, plot = FALSE, type = "correlation", demean = TRUE, lag.max = n - 1)$acf 
  y_autoc <- acf(Y, plot = FALSE, type = "correlation", demean = TRUE, lag.max = n - 1)$acf 
  
  # Calculate estimator of LRV for rho under independence
  Rho_ind_LRV <- sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))
  
  return(Rho_ind_LRV)
}

########################## Rhob Independent Long-Run Variance Estimation #############################
Rhob_ind_LRV <- function(X, Y, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)
  
  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")
  
  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  
  # Calculate autocovariances in a vector with row = lag
  x_autoc <- acf((rank(X) - 0.5) / n - 0.5, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf / cov((rank(X) - 0.5) / n - 0.5, (rank(X) - 0.5) / n - 0.5) # This is the acf of the demeaned grade. Therefore, demean = FALSE
  y_autoc <- acf((rank(Y) - 0.5) / n - 0.5, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf / cov((rank(Y) - 0.5) / n - 0.5, (rank(Y) - 0.5) / n - 0.5) # This is the acf of the demeaned grade. Therefore, demean = FALSE
  
  # Calculate estimator of LRV for rho_b under independence
  Rhob_ind_LRV <- sum(x_autoc[1] * y_autoc[1], 2 * (w * x_autoc[-1] * y_autoc[-1]))
  
  return(Rhob_ind_LRV)
}


########################## Tau Long-Run Variance Estimation #############################
Tau_LRV <- function(X, Y, kendall, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)
  
  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")
  
  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  
  # Need to define helper functions
  G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
  G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  
  # Define kernel realizations
  k_XY <- 4 * G_XY(X, Y) - 2 * (G_X(X) + G_Y(Y)) + 1 - kendall
  
  # Calculate autocovariances in a vector with row = lag
  k_XY_autoc <- (n - 1) / n * acf(k_XY, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY has mean 0. Therefore, demean = FALSE

  # Calculate estimator of LRV for tau
  Tau_LRV <- 4 * sum(k_XY_autoc[1], 2 * (w * k_XY_autoc[-1]))
  
  return(Tau_LRV)
}

########################## (Pearson) Rho Independent Long-Run Variance Estimation #############################
Rho_LRV <- function(X, Y, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)
  
  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")
  
  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  
  # Estimate values
  mean_x <- mean(X)
  mean_y <- mean(Y)
  sigma_xy <- (n - 1) / n * cov(X, Y)
  var_x <- (n - 1) / n * var(X)
  var_y <- (n - 1) / n * var(Y)
  x_autoc <- (n - 1) / n * acf(X, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  y_autoc <- (n - 1) / n * acf(Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  x2_autoc <- (n - 1) / n * acf(X^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  y2_autoc <- (n - 1) / n * acf(Y^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  xy_autoc <- (n - 1) / n * acf(X*Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  xy_crossc <- (n - 1) / n * ccf(X, Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  xx2_crossc <- (n - 1) / n * ccf(X, X^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  yx2_crossc <- (n - 1) / n * ccf(Y, X^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  xy2_crossc <- (n - 1) / n * ccf(X, Y^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  yy2_crossc <- (n - 1) / n * ccf(Y, Y^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  x2y2_crossc <- (n - 1) / n * ccf(X^2, Y^2, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  xxy_crossc <- (n - 1) / n * ccf(X, X*Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  yxy_crossc <- (n - 1) / n * ccf(Y, X*Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  x2xy_crossc <- (n - 1) / n * ccf(X^2, X*Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  y2xy_crossc <- (n - 1) / n * ccf(Y^2, X*Y, plot = FALSE, type = "covariance", demean = TRUE, lag.max = n - 1)$acf
  
  # Estimate Long-Run Variances
  x_LRVa <- sum(x_autoc[1], 2 * (w * x_autoc[-1]))
  y_LRVa <- sum(y_autoc[1], 2 * (w * y_autoc[-1]))
  x2_LRVa <- sum(x2_autoc[1], 2 * (w * x2_autoc[-1]))
  y2_LRVa <- sum(y2_autoc[1], 2 * (w * y2_autoc[-1]))
  xy_LRVa <- sum(xy_autoc[1], 2 * (w * xy_autoc[-1]))
  xy_LRVc <- sum(c(sort(w), 1, w) * xy_crossc)
  xx2_LRVc <- sum(c(sort(w), 1, w) * xx2_crossc)
  yx2_LRVc <- sum(c(sort(w), 1, w) * yx2_crossc)
  xy2_LRVc <- sum(c(sort(w), 1, w) * xy2_crossc)
  yy2_LRVc <- sum(c(sort(w), 1, w) * yy2_crossc)
  x2y2_LRVc <- sum(c(sort(w), 1, w) * x2y2_crossc)
  xxy_LRVc <- sum(c(sort(w), 1, w) * xxy_crossc)
  yxy_LRVc <- sum(c(sort(w), 1, w) * yxy_crossc)
  x2xy_LRVc <- sum(c(sort(w), 1, w) * x2xy_crossc)
  y2xy_LRVc <- sum(c(sort(w), 1, w) * y2xy_crossc)
  
  # Fill the matrix with long-Run variances
  Sigma <- diag(c(x_LRVa, y_LRVa, x2_LRVa, y2_LRVa, xy_LRVa))
  Sigma[upper.tri(Sigma)] <- c(xy_LRVc, xx2_LRVc, yx2_LRVc, xy2_LRVc, yy2_LRVc, x2y2_LRVc, xxy_LRVc, yxy_LRVc, x2xy_LRVc, y2xy_LRVc)
  Sigma <- as.matrix(forceSymmetric(Sigma, uplo = "U"))
  
  # Create Jacobian matrices for Delta method
  A <- matrix(c(-2*mean_x, -mean_y, 0, 0, -mean_x, -2*mean_y, 1, 0, 0, 0, 0, 1, 0, 1, 0), nrow = 3)
  B <- c(-0.5 * sigma_xy / sqrt(var_y) * sqrt(var_x) ^ (-3), 1 / sqrt(var_y * var_x), -0.5 * sigma_xy / sqrt(var_x) * sqrt(var_y) ^ (-3))
  
  Rho_LRV <- B %*% A %*% Sigma %*% t(A) %*% B
  
  return(Rho_LRV)
}


########################## (Spearman) SRho Long-Run Variance Estimation #############################
SRho_LRV <- function(X, Y, spearman, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)
  
  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")
  
  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  
  # Need to define helper functions
  G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
  G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  g_x <- Vectorize(function(x_val) mean(G_XY(x_val, Y)))
  g_y <- Vectorize(function(y_val) mean(G_XY(X, y_val)))
  
  # Define kernel realizations
  G_XX <- G_X(X)
  G_YY <- G_Y(Y)
  k_XY <- 4 * (g_x(X) + g_y(Y) + G_XX * G_YY - G_XX - G_YY) + 1 - spearman
  
  # Calculate autocovariances in a vector with row = lag
  k_XY_autoc <- (n - 1) / n * acf(k_XY, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY has mean 0. Therefore, demean = FALSE
  
  # Calculate estimator of LRV for srho
  SRho_LRV <- 9 * sum(k_XY_autoc[1], 2 * (w * k_XY_autoc[-1]))
  
  return(SRho_LRV)
}

########################## Gamma Long-Run Variance Estimation #############################
Gamma_LRV <- function(X, Y, kendall, tie_prob, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)
  
  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")
  
  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  
  # Need to define helper functions
  G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
  G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  x_eq_y_eq <- Vectorize(function(x_val, y_val) mean(X == x_val & Y == y_val))
  x_eq <- Vectorize(function(x_val) mean(X == x_val))
  y_eq <- Vectorize(function(y_val) mean(Y == y_val))
  
  # Define kernel realizations
  k_XY_tau <- 4 * G_XY(X, Y) - 2 * (G_X(X) + G_Y(Y)) + 1 - kendall
  k_XY_tie <- x_eq(X) + y_eq(Y) - x_eq_y_eq(X, Y) - tie_prob
  
  # Calculate autocovariances in a vector with row = lag
  k_XY_tau_autoc <- (n - 1) / n * acf(k_XY_tau, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY_tau has mean 0. Therefore, demean = FALSE
  k_XY_tie_autoc <- (n - 1) / n * acf(k_XY_tie, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY_tie has mean 0. Therefore, demean = FALSE
  k_XY_tautie_crossc <- (n - 1) / n * ccf(k_XY_tau, k_XY_tie, plot = FALSE, type = "covariance", lag.max = n - 1)$acf
  
  # Calculate estimator of LRV for gamma
  sigma_tau_sq <- 4 * sum(k_XY_tau_autoc[1], 2 * (w * k_XY_tau_autoc[-1]))
  sigma_nu_sq <- 4 * sum(k_XY_tie_autoc[1], 2 * (w * k_XY_tie_autoc[-1]))
  sigma_taunu <- 4 * sum(c(sort(w), 1, w) * k_XY_tautie_crossc)
  gamma <- kendall/(1-tie_prob)
  
  Gamma_LRV <- (sigma_tau_sq + gamma^2 * sigma_nu_sq + 2 * gamma * sigma_taunu) / (1 - tie_prob)^2
  
  return(Gamma_LRV)
}

########################## TauB Long-Run Variance Estimation #############################
TauB_LRV <- function(X, Y, kendall, kendall_X, kendall_Y, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)
  
  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")
  
  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  
  # Need to define helper functions
  G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
  G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  x_neq <- Vectorize(function(x_val) mean(X != x_val))
  y_neq <- Vectorize(function(y_val) mean(Y != y_val))
  
  # Define kernel realizations
  k_XY_tau <- 4 * G_XY(X, Y) - 2 * (G_X(X) + G_Y(Y)) + 1 - kendall
  k_X_tau <- x_neq(X) - kendall_X
  k_Y_tau <- y_neq(Y) - kendall_Y
  
  # Calculate autocovariances in a vector with row = lag
  k_XY_tau_autoc <- (n - 1) / n * acf(k_XY_tau, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY_tau has mean 0. Therefore, demean = FALSE
  k_X_tau_autoc <- (n - 1) / n * acf(k_X_tau, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_X_tau has mean 0. Therefore, demean = FALSE
  k_Y_tau_autoc <- (n - 1) / n * acf(k_Y_tau, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_Y_tau has mean 0. Therefore, demean = FALSE
  k_X_tau_crossc <- (n - 1) / n * ccf(k_XY_tau, k_X_tau, plot = FALSE, type = "covariance", lag.max = n - 1)$acf
  k_Y_tau_crossc <- (n - 1) / n * ccf(k_XY_tau, k_Y_tau, plot = FALSE, type = "covariance", lag.max = n - 1)$acf
  k_XY_tautau_crossc <- (n - 1) / n * ccf(k_X_tau, k_Y_tau, plot = FALSE, type = "covariance", lag.max = n - 1)$acf
  
  # Calculate estimator of LRV for taub
  sigma_tau_sq <- 4 * sum(k_XY_tau_autoc[1], 2 * (w * k_XY_tau_autoc[-1]))
  sigma_tauX_sq <- 4 * sum(k_X_tau_autoc[1], 2 * (w * k_X_tau_autoc[-1]))
  sigma_tauY_sq <- 4 * sum(k_Y_tau_autoc[1], 2 * (w * k_Y_tau_autoc[-1]))
  sigma_tautauX <- 4 * sum(c(sort(w), 1, w) * k_X_tau_crossc)
  sigma_tautauY <- 4 * sum(c(sort(w), 1, w) * k_Y_tau_crossc)
  sigma_tauXtauY <- 4 * sum(c(sort(w), 1, w) * k_XY_tautau_crossc)

  TauB_LRV <- (sigma_tau_sq - kendall * (sigma_tautauX / kendall_X - sigma_tautauY / kendall_Y) + kendall^2 / 4 * (sigma_tauX_sq / kendall_X^2 + sigma_tauY_sq / kendall_Y^2 + (2 * sigma_tauXtauY) / kendall_X / kendall_Y)) / (kendall_X * kendall_Y)
  return(TauB_LRV)
}

########################## Rhob Long-Run Variance Estimation #############################
Rhob_LRV <- function(X, Y, spearman, spearman_X, spearman_Y, bandwidth = "Dehling"){
  if (length(X) != length(Y)){stop("X and Y must have equal length")}
  n <- length(X)
  
  # Determine bandwidth
  if (bandwidth == "StockWatson"){b <- floor(0.75 * n^(1/3))}
  else if (bandwidth == "Dehling"){b <- floor(2 * n^(1/3))}
  else stop("Please insert a valid bandwith calculation method")
  
  # Calculate weights
  h <- 1:(n-1)
  w <- pmax(1 - abs(h) / (b + 1), 0)
  
  # Need to define helper functions
  G_XY <- Vectorize(function(x_val, y_val) (mean(X <= x_val & Y <= y_val) + mean(X <= x_val & Y < y_val) + mean(X < x_val & Y <= y_val) + mean(X < x_val & Y < y_val)) / 4)
  G_X <- Vectorize(function(x_val) (mean(X < x_val) + mean(X <= x_val)) / 2)
  F_X <- Vectorize(function(x_val) mean(X <= x_val))
  F_X_ <- Vectorize(function(x_val) mean(X < x_val))
  G_Y <- Vectorize(function(y_val) (mean(Y < y_val) + mean(Y <= y_val)) / 2)
  F_Y <- Vectorize(function(y_val) mean(Y <= y_val))
  F_Y_ <- Vectorize(function(y_val) mean(Y < y_val))
  g_x <- Vectorize(function(x_val) mean(G_XY(x_val, Y)))
  g_y <- Vectorize(function(y_val) mean(G_XY(X, y_val)))
  f_x1 <- Vectorize(function(x_val) mean(pmin(F_X(x_val), F_X(X))))
  f_x2 <- Vectorize(function(x_val) mean(pmin(F_X_(x_val), F_X(X))))
  f_x3 <- Vectorize(function(x_val) mean(pmin(F_X(x_val), F_X_(X))))
  f_x4 <- Vectorize(function(x_val) mean(pmin(F_X_(x_val), F_X_(X))))
  f_y1 <- Vectorize(function(y_val) mean(pmin(F_Y(y_val), F_Y(Y))))
  f_y2 <- Vectorize(function(y_val) mean(pmin(F_Y_(y_val), F_Y(Y))))
  f_y3 <- Vectorize(function(y_val) mean(pmin(F_Y(y_val), F_Y_(Y))))
  f_y4 <- Vectorize(function(y_val) mean(pmin(F_Y_(y_val), F_Y_(Y))))
  
  # Define kernel realizations
  G_XX <- G_X(X)
  G_YY <- G_Y(Y)
  k_XY_rho <- 4 * (g_x(X) + g_y(Y) + G_XX * G_YY - G_XX - G_YY) + 1 - spearman
  k_X_rho <- 2 * (f_x1(X) + f_x2(X) + f_x3(X) + f_x4(X) + 2 * G_XX^2 - 4 * G_XX) + 1 - spearman_X
  k_Y_rho <- 2 * (f_y1(Y) + f_y2(Y) + f_y3(Y) + f_y4(Y) + 2 * G_YY^2 - 4 * G_YY) + 1 - spearman_Y
  
  # Calculate autocovariances in a vector with row = lag
  k_XY_rho_autoc <- (n - 1) / n * acf(k_XY_rho, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_XY has mean 0. Therefore, demean = FALSE
  k_X_rho_autoc <- (n - 1) / n * acf(k_X_rho, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_X_tie has mean 0. Therefore, demean = FALSE
  k_Y_rho_autoc <- (n - 1) / n * acf(k_Y_rho, plot = FALSE, type = "covariance", demean = FALSE, lag.max = n - 1)$acf # k_Y_tie has mean 0. Therefore, demean = FALSE
  k_X_rho_crossc <- (n - 1) / n * ccf(k_XY_rho, k_X_rho, plot = FALSE, type = "covariance", lag.max = n - 1)$acf
  k_Y_rho_crossc <- (n - 1) / n * ccf(k_XY_rho, k_Y_rho, plot = FALSE, type = "covariance", lag.max = n - 1)$acf
  k_XY_rho_crossc <- (n - 1) / n * ccf(k_X_rho, k_Y_rho, plot = FALSE, type = "covariance", lag.max = n - 1)$acf
  
  # Calculate estimator of LRV for srho
  sigma_rho_sq <- 9 * sum(k_XY_rho_autoc[1], 2 * (w * k_XY_rho_autoc[-1]))
  sigma_rhoX_sq <- 9 * sum(k_X_rho_autoc[1], 2 * (w * k_X_rho_autoc[-1]))
  sigma_rhoY_sq <- 9 * sum(k_Y_rho_autoc[1], 2 * (w * k_Y_rho_autoc[-1]))
  sigma_rhorhoX <- 9 * sum(c(sort(w), 1, w) * k_X_rho_crossc)
  sigma_rhorhoY <- 9 * sum(c(sort(w), 1, w) * k_Y_rho_crossc)
  sigma_rhoXrhoY <- 9 * sum(c(sort(w), 1, w) * k_XY_rho_crossc)
  
  Rhob_LRV <- (sigma_rho_sq - spearman * (sigma_rhorhoX / spearman_X + sigma_rhorhoY / spearman_Y) + spearman^2 / 4 * (sigma_rhoX_sq / spearman_X^2 + sigma_rhoY_sq / spearman_Y^2 + (2 * sigma_rhoXrhoY) / spearman_Y / spearman_X)) / (spearman_X * spearman_Y)
  
  return(Rhob_LRV)
}

