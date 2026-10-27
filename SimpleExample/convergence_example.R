################################# FRONT MATTER #################################

rm(list=ls())
gc()

library(dplyr)
library(ggplot2)
library(tidyr)
library(mvtnorm)
library(clusterGeneration)
library(xtable)

set.seed(1234)

# X takes K + 1 values 0:K with K even
# Defiers for all X < K / 2, compliers for all X > K / 2, no switchers for
# X = K / 2.
# Assume a fixed covariate type design: n observations for each covariate value
# with a non-zero number of switchers. 

################################### FUNCTIONS ##################################

gen_data <- function(n, sigma, pZ, K, a, d, B){
  
  # Generates (D,X,Y,Z)
  # inputs
  # n: number of observations for each X value
  # sigma: covariance matrix for the unobservables
  # pZ: probability Z = 1
  # K: x values range from 0 to K
  # a,d: first stage D = 1(a + d * Z >= V)
  # B: ATE
  
  # Unobservables
  U <- rmvnorm(n * K, mean = rep(0, 3), sigma = sigma)
  U <- data.frame(U)
  names(U) <- c("U1", "U0", "V")
  
  # X and Z
  X <- rep(c(0:(K / 2 -1), ((K / 2) + 1):K), each = n)
  Z <- rbinom(n * K, 1, pZ)
  
  # Construct D and Y
  D <- (a + d[X + 1] * Z >= U$V)
  Y <- (1 - D) * U$U0 + D *(B + U$U1)
  
  return(data.frame(D = D, X = X, Y = Y, Z = Z))
  
}

# gen_data <- function(n, sigma, p_x, pZ, K, a, d, B){
# 
#   # Generates (D,X,Y,Z)
#   # inputs
#   # n: number of observations
#   # sigma: covariance matrix for the unobservables
#   # p_x: distribution for x
#   # pZ: probability Z = 1
#   # K: x values range from 0 to K
#   # a,d: first stage D = 1(a + d * Z >= V)
#   # B: ATE
# 
#   # Unobservables
#   U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
#   U <- data.frame(U)
#   names(U) <- c("U1", "U0", "V")
# 
#   # X and Z
#   X <- sample(0:K, size = n, replace = T, prob = p_x)
#   Z <- rbinom(n, 1, pZ)
# 
#   # Construct D and Y
#   D <- (a + d[X + 1] * Z >= U$V)
#   Y <- (1 - D) * U$U0 + D *(B + U$U1)
# 
#   return(data.frame(D = D, X = X, Y = Y, Z = Z))
# 
# }

gen_vars <- function(data, K, true_s){
  
  # adds in the demeaned Z, true sign function, and estimated sign function
  # inputs:
  # data: D, X, Y, Z, ordered by X, n draws for each X
  # K x values range from 0 to K
  # true_s value of the sign function for each X
  
  X_vals <- sort(unique(data$X))
  
  # Z_tilde = Z - E[Z|X]
  bar_Z <- rep(0, K + 1)
  for (k in X_vals){
    bar_Z[k + 1] <- mean(data$Z[data$X == k])
  }
  data$Z_tilde <- data$Z - bar_Z[data$X + 1]
  
  # Oracle sign function
  data$s_star <- true_s[data$X + 1] 
  
  # Estimated sign function
  est_s <- rep(0, K + 1)
  for (k in X_vals){
    
    temp <- dplyr::filter(data, X == k)
    n_1 <- sum(temp$Z)
    n_0 <- sum(1 - temp$Z)
    
    if (n_1 > 0 & n_0 > 0){
      est_s[k + 1] <- 
        2 * (sum(temp$Z * temp$D) / n_1 >= sum((1 - temp$Z) * temp$D) / n_0) - 1 
    }
  }
  data$s_hat <- est_s[data$X + 1]
  
  return(data)
  
}

# gen_vars <- function(data, K, true_s){
#   
#   # adds in the demeaned Z, true sign function, and estimated sign function
#   # inputs:
#   # data data frame containing D, X, Y, Z
#   # K x values range from 0 to K
#   # true_s value of the sign function for each X
#   
#   X_vals <- sort(unique(data$X))
#   
#   # Z_tilde = Z - E[Z|X]
#   bar_Z <- rep(0, K + 1)
#   for (k in X_vals){
#     bar_Z[k + 1] <- mean(data$Z[data$X == k])
#   }
#   data$Z_tilde <- data$Z - bar_Z[data$X + 1]
#   
#   # Oracle sign function
#   data$s_star <- true_s[data$X + 1] 
#   
#   # Estimated sign function
#   est_s <- rep(0, K + 1)
#   for (k in X_vals){
#     
#     temp <- dplyr::filter(data, X == k)
#     n_1 <- sum(temp$Z)
#     n_0 <- sum(1 - temp$Z)
#     
#     if (n_1 > 0 & n_0 > 0){
#       est_s[k + 1] <- 
#         2 * (sum(temp$Z * temp$D) / n_1 > sum((1 - temp$Z) * temp$D) / n_0) - 1 
#     }
#   }
#   data$s_hat <- est_s[data$X + 1]
#   
#   return(data)
#   
# }

calc_sd <- function(data, tau_feasible, n, K){
  
  # returns the estimate of the asymptotic SD
  
  # data: a data frame containing the data (including generated variables)
  # tau_feasible: the estimate (from the feasible estimator)
  # n: number of observations for each value of X
  # K: number of unique X values
  
  # "errors"
  e_Y <- data$Y - data$D * tau_feasible
  e_Z <- data$Z_tilde
  
  # generate the indicators (which are the implicit covariates)
  X_ind <- matrix(0, nrow(data), K)
  index <- 0
  for (k in 1:K){
    
    X_ind[index + 1:n, k] <- 1
    index <- index + n
    
  }
   
  # sigma_X = E_n[XX^T]
  X_mean <- colMeans(X_ind)
  X_temp <- X_ind - matrix(X_mean, n * K, K, byrow = T)
  sigma_X <- cov(X_temp) * ((n * K) - 1) / (n * K) + 
    matrix(X_mean, nrow = K, ncol = 1) %*% matrix(X_mean, nrow = 1, ncol = K)
  
  # C = E_n[e_Y * s_hat(X) * X]
  C <- apply(X_ind, 2, function(X_col, e_Y, s_hat){
    return(X_col * e_Y * s_hat)
  }, e_Y = e_Y, s_hat = data$s_hat)
  C <- colSums(C) / (n * K)
  
  # denominator
  denom <- (sum(data$D * data$s_hat * e_Z) / (n * K))
  
  # numerator term 1: E_n[e_Y ^ 2 * e_Z ^ 2 * s_hat]
  num_1 <- sum((e_Y ^ 2) * (e_Z ^ 2)) / (n * K)
  
  # numerator term 2: C^T %*% sigma_X^{-1} %*% E_n[e_Y * s_hat * (e_Z ^ 2) * X]
  # 
  num_2 <- apply(X_ind, 2, function(X_col, e_Y, s_hat, e_Z){
    return(e_Y * s_hat * (e_Z ^ 2) * X_col)
  }, e_Y = e_Y, s_hat = data$s_hat, e_Z = e_Z)
  num_2 <- matrix(colSums(num_2) / (n * K), nrow = K, ncol = 1)
  num_2 <- matrix(C, nrow = 1, ncol = K) %*% solve(sigma_X) %*% num_2
  
  # numerator term 3: C^T %*% sigma_X^{-1} %*% E_n[e_Z^2 XX^T] %*% sigma_X^{-1}
  # %*% C
  X_temp <- apply(X_ind, 2, function(X_col, e_Z){
    return(X_col * e_Z)
  }, e_Z = e_Z)
  X_mean <- colMeans(X_temp)
  X_temp <- X_temp - matrix(X_mean, n * K, K, byrow = T)
  num_3 <-  cov(X_temp) * ((n * K) - 1) / (n * K) + 
    matrix(X_mean, nrow = K, ncol = 1) %*% matrix(X_mean, nrow = 1, ncol = K)
  num_3 <- matrix(C, nrow = 1, ncol = K) %*% solve(sigma_X) %*%
    num_3 %*% solve(sigma_X) %*% matrix(C, nrow = K, ncol = 1)
  
  # numerator = num_1 - 2 * num_2 + num_3
  num <- num_1 - 2 * num_2 + num_3
  
  return(sqrt(num) / denom)
 
}

# Given a vector containing estimates of tau, this returns summary statistics
get_stats <- function(results, B){
  
  # results is a vector of estimates of tau.
  # B is the ATE, which equals the asymptotic estimand. 
  
  # statistics returned are: mean, sd, MSE, median, 0.25 and 0.75 quantiles
  results <- results[!is.infinite(results)]
  stats <- rep(0, 6)
  stats[1] <- mean(results, na.rm = TRUE)
  stats[2] <- sd(results, na.rm = TRUE)
  stats[3] <- (stats[1] - B) ^ 2 + stats[2] ^ 2
  stats[4:6] <- quantile(results, c(0.5, 0.25, 0.75), na.rm = TRUE)
  return(stats)

}

get_coverage <- function(tau, sigma_o, B, n, K){
  
  # returns the 90%, 95%, and 99% coverage rates of the estimators
  # all inputs should be results for a given n
  # tau: a vector of estimates (either feasible or oracle)
  # sigma_o: a vector of estimated standard deviations
  # B: the true AT (asymptotic limit of the estimator)
  # n: number of observations for each X value
  # K: number of X values
  
  coverage <- matrix(0, 3, length(tau))
  
  sigma_o <- sigma_o / sqrt(n * K)
  
  CI <- matrix(0, 3, length(sigma_o))
  CI[1, ] <- sigma_o * qnorm(0.95)
  CI[2, ] <- sigma_o * qnorm(0.975)
  CI[3, ] <- sigma_o * qnorm(0.99)
  rownames(CI) <- c("90", "95", "99")
  
  for (i in 1:3){
    
    coverage[i, ] <- as.numeric((tau + CI[i, ]) > B)
    coverage[i, ] <- coverage[i, ] * as.numeric((tau - CI[i, ]) < B)
    
  }
  
  coverage <- rowMeans(coverage)
  names(coverage) <- c("90", "95", "99")
  
  return(coverage)
  
}

# get_coverage <- function(tau_oracle, tau_feasible, sigma_o, B, n, K){
#   
#   # returns the 90%, 95%, and 99% coverage rates of the estimators
#   # all inputs should be results for a given n
#   # tau_oracle: a vector of oracle estimates
#   # tau_feasible: a vector of feasible estimates
#   # sigma_o: a vector of estimated standard deviations
#   # B: the true AT (asymptotic limit of the estimator)
#   # n: number of observations for each X value
#   # K: number of X values
#   
#   coverage_feasible <- matrix(0, 3, length(tau_oracle))
#   coverage_oracle <- matrix(0, 3, length(tau_oracle))
#   
#   sigma_o <- sigma_o / sqrt(n * K)
#   
#   CI <- matrix(0, 3, length(sigma_o))
#   CI[1, ] <- sigma_o * qnorm(0.95)
#   CI[2, ] <- sigma_o * qnorm(0.975)
#   CI[3, ] <- sigma_o * qnorm(0.99)
#   rownames(CI) <- c("90", "95", "99")
#   
#   for (i in 1:3){
#     
#     coverage_feasible[i, ] <- as.numeric(tau_feasible + CI[i, ] > 0)
#     coverage_feasible[i, ] <- coverage_feasible[i, ] * 
#       as.numeric(tau_feasible - CI[i, ] < 0)
#     
#     coverage_oracle[i, ] <- as.numeric(tau_oracle + CI[i, ] > 0)
#     coverage_oracle[i, ] <- coverage_oracle[i, ] * 
#       as.numeric(tau_oracle - CI[i, ] < 0)
#     
#   }
#   
#   coverage_feasible <- rowMeans(coverage_feasible)
#   names(coverage_feasible) <- c("90", "95", "99")
#   coverage_oracle <- rowMeans(coverage_oracle)
#   names(coverage_oracle) <- c("90", "95", "99")
#   
#   return(list(feasible = feasible, oracle = oracle))
#   
# }

make_plot <- function(tau_oracle, tau_feasible, x_lims){
  
  oracle_dens <- density(tau_oracle)
  oracle_data <- data.frame(x = oracle_dens$x, y = oracle_dens$y)
  
  feasible_dens <- density(tau_feasible)
  feasible_data <- data.frame(x = feasible_dens$x, y = feasible_dens$y)
  
  ggplot() +
    geom_line(aes(x = x, y = y, col = "Oracle"), 
              data = oracle_data, 
              linetype = "dashed") +
    geom_line(aes(x = x, y = y, col = "Feasible"), 
              data = feasible_data) +
    geom_vline(aes(xintercept = 0, col = "black"), 
               linetype = "dotted", show.legend = FALSE) + 
    theme_bw() +
    theme(legend.title = element_blank()) +
    scale_color_manual(values = c("Oracle" = "#F8766D",
                                  "Feasible" = "#619CFF"),
                       labels = c("Oracle" = expression(sqrt(n) * 
                                                          "(" * 
                                                          tau[o] - 
                                                          tau[r]^"*" 
                                                        * ")"),
                                  "Feasible" = expression(sqrt(n) * 
                                                            "(" * 
                                                            tau[r] - 
                                                            tau[r]^"*" 
                                                          * ")"))) + 
    xlim(x_lims[1], x_lims[2]) +  
    labs(x = "", y = "") + 
    theme(legend.position = c(0.2, 0.9))
  
}

# scale_color_manual(values = c("Oracle" = "#F8766D",
#                               "Feasible" = "#619CFF",
#                               "Asymptotic Limit" = "black"),
#                    labels = c("Oracle" = expression(sqrt(n) * 
#                                                       "(" * 
#                                                       tau[o] - 
#                                                       tau[r]^"*" 
#                                                     * ")"),
#                               "Feasible" = expression(sqrt(n) * 
#                                                         "(" * 
#                                                         tau[r] - 
#                                                         tau[r]^"*" 
#                                                       * ")"),
#                               "Asymptotic Limit" = "")) 

################################## PARAMETERS ##################################

# Number of X values (is K)
K <- 8

# Switcher shares at each value of X are evenly spaced between -delta and delta
delta <- 0.2
w <- delta * (2 * (0:K) / K - 1)

# First stage 1(a + dZ >= V)
a <- 0
d <- qnorm(pnorm(a) + w) - a 

# ATE (also equals the asymptotic limit of the estimand)
B <- 0 

# True value of the sign function for each value of X
true_s <- 2 * (w >= 0) - 1

################################# DISTRIBUTIONS ################################

# Randomly generate a covariance matrix for the unobservables
# A <- matrix(runif(9, -1, 1), 3, 3)
# sigma <- A %*% t(A)
# # Scale to be a correlation matrix
# for (k in 1:3){
#   
#   sd_k <- sqrt(sigma[k, k])
#   sigma[, k] <- sigma[, k] / sd_k
#   sigma[k, ] <- sigma[k, ] / sd_k
#   
# }

# Set a correlation matrix for the unobservables (order is U1, U0, V)
sigma <- diag(x = 1, 3, 3)
sigma[2, 1] <- -0.25
sigma[3, 1] <- -0.5
sigma[3, 2] <- 0.5
sigma[1, 2] <- sigma[2, 1]
sigma[1, 3] <- sigma[3, 1]
sigma[2, 3] <- sigma[3, 2]
# As V is resistance to treatment, this specification implies a positive 
# selection effect. 

# # Distribution for X
# p_x <- vector(mode = "list", length = 2)
# p_x[[1]] <- c(1 / 2, 0, 1 / 2)
# p_x[[2]] <- rep(1 / 3, 3)

# Probability that Z = 1
pZ <- 1 / 2

################################# GENERATE DATA ################################

n <- c(50, 100, 200, 400) # number of observations for each X
runs <- 1000 # number of draws for each sample size

tau_oracle <- matrix(0, length(n), runs)
tau_feasible <- matrix(0, length(n), runs)
sigma_o <- matrix(0, length(n), runs)

for (i in 1:length(n)){
  for (j in 1:runs){
    
    data <- gen_data(n = n[i],
                     sigma = sigma,
                     pZ = pZ,
                     K = K,
                     a = a,
                     d = d,
                     B = B)
    data <- gen_vars(data = data,
                     K = K,
                     true_s = true_s)
    
    # estimates
    tau_oracle[i, j] <- sum(data$Y * data$s_star * data$Z_tilde) / 
      sum(data$D * data$s_star * data$Z_tilde)
    tau_feasible[i, j] <- sum(data$Y * data$s_hat * data$Z_tilde) / 
      sum(data$D * data$s_hat * data$Z_tilde)
    
    # standard deviation (calculated using the asymptotic formula)
    sigma_o[i, j] <- calc_sd(data, tau_feasible[i, j], n[i], K)
    
  }
}

############################## SUMMARY STATISTICS ##############################

# Returns the mean, standard error, MSE, median, 0.25 and 0.75 quartiles of the 
# sampling distribution of the oracle and feasible estimators. 

feasible_sum_stats <- matrix(0, length(n), 6)
oracle_sum_stats <- matrix(0, length(n), 6)

for (i in 1:length(n)){
  
  feasible_sum_stats[i, ] <- get_stats(results = tau_feasible[i, ], B = B)
  oracle_sum_stats[i, ] <- get_stats(results = tau_oracle[i, ], B = B)
  
}

rownames(feasible_sum_stats) <- n
rownames(oracle_sum_stats) <- n

xtable(feasible_sum_stats, digits = 3)
xtable(oracle_sum_stats, digits = 3)

# Add in coverage rates 
coverage_feasible <- matrix(0, length(n), 3)
coverage_oracle <- matrix(0, length(n), 3)

# n varies over rows, columns are coverage rates for 90%, 95% and 99% CIs
for (i in 1:length(n)){
  
  coverage_feasible[i, ] <- get_coverage(tau = tau_feasible[i, ], 
                                         sigma_o = sigma_o[i, ],
                                         B = B, n = n[i], K = K)
  coverage_oracle[i, ] <- get_coverage(tau = tau_oracle[i, ], 
                                       sigma_o = sigma_o[i, ],
                                       B = B, n = n[i], K = K)
  
}

colnames(coverage_feasible) <- c("90", "95", "99")
rownames(coverage_feasible) <- n

colnames(coverage_oracle) <- c("90", "95", "99")
rownames(coverage_oracle) <- n

xtable(cbind(coverage_feasible, coverage_oracle), digits = 3)

######################## PLOT OF SAMPLING DISTRIBUTIONS ########################

# Prepare data for plotting

# convert to sqrt(n) * (tau - tau^*)
tau_oracle <- tau_oracle - B
tau_feasible <- tau_feasible - B
for (i in 1:nrow(tau_oracle)){
  
  tau_oracle[i, ] <- sqrt(K * n[i]) * tau_oracle[i, ]
  tau_feasible[i, ] <- sqrt(K * n[i]) * tau_feasible[i, ]
  
}

x_lims <- rep(0, 2)
x_lims[1] <- quantile(c(tau_feasible, tau_oracle), 0.005, na.rm = TRUE)
x_lims[2] <- quantile(c(tau_feasible, tau_oracle), 0.995, na.rm = TRUE)

make_plot(tau_oracle = tau_oracle[1, ],
          tau_feasible = tau_feasible[1, ],
          x_lims = x_lims)
make_plot(tau_oracle = tau_oracle[2, ],
          tau_feasible = tau_feasible[2, ],
          x_lims = x_lims)
make_plot(tau_oracle = tau_oracle[3, ],
          tau_feasible = tau_feasible[3, ],
          x_lims = x_lims)
make_plot(tau_oracle = tau_oracle[4, ],
          tau_feasible = tau_feasible[4, ],
          x_lims = x_lims)



