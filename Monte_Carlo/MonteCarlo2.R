#################################### SUMMARY ################################### 

# This scipt performs Monte Carlo for the feasible reordered IV estimator

# The first stage is assumed to be
# D = 1(Z *(d0 + x'd) > V)
# where the marginal distribution of V is N(0, 1), so the propensity score is
# p(x,z) = Phi(Z *(d0 + x'd))
# d0 + x'd controls the sign of the relationship between Z and D. 

# We use the maximum rank correlation estimator to obtain estimates of the first
# stage parameters. This estimator does not estimate a_0, and the remaining 
# parameters are only estimated up to scale (adding a constant to all terms 
# does not perturb ranks, and multiplying all non-constant terms by a positive
# constant would also not affect ranks)

# The second stage is 
# Y(0) = U(0)
# Y(1) = B_0 + x'B + U(1)
# Y = Y(0) + D(Y(1) - Y(0))
# where (V, U(1), U(0)) are jointly normal. 

################################# FRONT MATTER #################################

rm(list = ls())
gc()

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(boot)
library(xtable)

path <- getwd()

set.seed(0451)

################################### FUNCTIONS ##################################

# Generate parameter vectors b and B which match sd_B, sd_d, and rho_Bd
gen_Bd <- function(K, sd_B, sd_d, rho_Bd){
  
  # create orthonormal vectors e_1 and e_2
  e_1 <- runif(K, -1, 1)
  e_1 <- e_1 / sqrt(sum(e_1 ^ 2))
  u <- runif(K, -1, 1)
  e_2 <- u - (sum(e_1 * u)) * e_1 / (sum(e_1 ^ 2))
  e_2 <- e_2 / sqrt(sum(e_2 ^ 2))
  
  # a is a scaled version of e_1
  B <- sd_B * e_1
  # d is a combination of e_1 and e_2
  d <- sd_d * rho_Bd * e_1 + sd_d * sqrt(1 - rho_Bd ^ 2) * e_2
  
  return(list(B = B,
              d = d))
  
}

# Efficient calculation of the MRC objective function (as suggested by 
# Windmeijer/Abrevata) 
# * theta is a vector of estimated parameters (a, d_0, d)
# * data is a data frame containing X, Z, Z * X
# * varnames is a list of all variables used in the regression (X_1, X_2, etc)
S_theta <- function(theta, data, varnames){
  
  theta <- theta / sqrt(sum(theta ^ 2))
  index <- as.matrix(data[varnames]) %*% theta
  D_sort <- data$D[order(index)] 
  tally <- cumsum(1 - D_sort)
  val <- sum(tally[as.logical(D_sort)])
  
  return(val)
  
}

calc_sd <- function(data, tau_feasible, n, K){
  
  # returns the estimate of the asymptotic SD
  
  # data: a data frame containing the data (including generated variables)
  # tau_feasible: the estimate (from the feasible estimator)
  # n: number of observations
  # K: number of unique X values
  
  # "errors"
  e_Y <- data$Y - data$D * tau_feasible
  e_Z <- data$Z_tilde
  
  X <- data[, grep("^X", names(data))]
  
  # sigma_X = E_n[XX^T]
  sigma_X <- (t(as.matrix(X)) %*% as.matrix(X)) / n
  
  # C = E_n[e_Y * s_hat(X) * X]
  C <- apply(X, 2, function(X_col, e_Y, s_hat){
    return(X_col * e_Y * s_hat)
  }, e_Y = e_Y, s_hat = data$s_hat)
  C <- colMeans(C) 
  
  # denominator
  denom <- sum(data$D * data$s_hat * e_Z) / n
  
  # numerator term 1: E_n[e_Y ^ 2 * e_Z ^ 2 * s_hat]
  num_1 <- sum((e_Y ^ 2) * (e_Z ^ 2)) / n 
  
  # numerator term 2: C^T %*% sigma_X^{-1} %*% E_n[e_Y * s_hat * (e_Z ^ 2) * X]
  # 
  num_2 <- apply(X, 2, function(X_col, e_Y, s_hat, e_Z){
    return(e_Y * s_hat * (e_Z ^ 2) * X_col)
  }, e_Y = e_Y, s_hat = data$s_hat, e_Z = e_Z)
  num_2 <- matrix(colMeans(num_2), nrow = K, ncol = 1)
  num_2 <- matrix(C, nrow = 1, ncol = K) %*% solve(sigma_X) %*% num_2
  
  # numerator term 3: C^T %*% sigma_X^{-1} %*% E_n[e_Z^2 XX^T] %*% sigma_X^{-1}
  # %*% C
  X_temp <- apply(X, 2, function(X_col, e_Z){
    return(X_col * e_Z)
  }, e_Z = e_Z)
  X_temp <- (t(X_temp) %*% X_temp) / n
  num_3 <- matrix(C, nrow = 1, ncol = K) %*% solve(sigma_X) %*%
    X_temp %*% solve(sigma_X) %*% matrix(C, nrow = K, ncol = 1)
  
  # numerator = num_1 - 2 * num_2 + num_3
  num <- num_1 - 2 * num_2 + num_3
  
  return(sqrt(num) / denom)
  
}

# Given a vector containing estimates of tau, this returns summary statistics
get_stats <- function(results, tau){
  
  # results is a vector of estimates of tau.
  # tau is the asymptotic estimand. 
  
  # statistics returned are: mean, sd, MSE, median, 0.25 and 0.75 quantiles
  results <- results[!is.infinite(results)]
  stats <- rep(0, 6)
  stats[1] <- mean(results, na.rm = TRUE)
  stats[2] <- sd(results, na.rm = TRUE)
  stats[3] <- (stats[1] - tau) ^ 2 + stats[2] ^ 2
  stats[4:6] <- quantile(results, c(0.5, 0.25, 0.75), na.rm = TRUE)
  return(stats)
  
}

get_coverage <- function(results, se, tau, n){
  
  # returns the 90%, 95%, and 99% coverage rates of the estimators
  # all inputs should be results for a given n
  # results: a vector of estimates (either feasible or oracle)
  # se: a vector of estimated standard deviations
  # tau: the asymptotic limit of the estimator
  # n: number of observations
  
  coverage <- matrix(0, 3, length(results))
  
  se <- se / sqrt(n)
  
  CI <- matrix(0, 3, length(se))
  CI[1, ] <- se * qnorm(0.95)
  CI[2, ] <- se * qnorm(0.975)
  CI[3, ] <- se * qnorm(0.99)
  rownames(CI) <- c("90", "95", "99")
  
  for (i in 1:3){
    
    coverage[i, ] <- as.numeric((results + CI[i, ]) > tau)
    coverage[i, ] <- coverage[i, ] * as.numeric((results - CI[i, ]) < tau)
    
  }
  
  coverage <- rowMeans(coverage)
  names(coverage) <- c("90", "95", "99")
  
  return(coverage)
  
}

make_plot <- function(tau_oracle, tau_feasible, x_lims, trimming_flag = FALSE){
  
  tau_oracle <- tau_oracle[!is.infinite(tau_oracle)]
  tau_feasible <- tau_feasible[!is.infinite(tau_feasible)]
  
  if (trimming_flag == TRUE){
    tau_oracle <- tau_oracle[as.logical((tau_oracle > 2 * x_lims[1]) * 
                                          (tau_oracle < 2 * x_lims[2]))]
  }
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

##################################### SETUP ####################################

# Load in parameters
switchers_calibration <- 
  readRDS(paste0(path, "/Monte_Carlo/switchers_calibration.rds"))
parameters <- switchers_calibration$parameters[[1]]

for (i in 1:length(parameters)){
  assign(names(parameters)[i], parameters[[i]])
}

# Settings 
K <- 3 # number of X variables
n <- c(100, 500, 1000, 2000) # number of observations
pZ <- 1 / 2 # probabilty Z = 1
B0 <- 0 # ATE

# Monte Carlo settings
runs <- 1000

############################### MONTE CARLO LOOP ###############################

results <- vector(mode = "list", length = length(n))

for (j in 1:length(n)){
  
  tsls <- rep(0, runs)
  feasible <- rep(0, runs)
  oracle <- rep(0, runs)
  probit <- rep(0, runs)
  logit <- rep(0, runs)
  lpm <- rep(0, runs)
  diagnostics <- list(s_cor = rep(0, runs),
                      theta_hat = matrix(0, nrow = 1 + (2 * K), ncol = runs),
                      theta = matrix(0, nrow = 1 + (2 * K), ncol = runs))
  se <- rep(0, runs)
  
  for (r in 1:runs){
    
    #### GENERATE DATA ####
    
    # Generate parameters according to the constraints
    Bd <- gen_Bd(K, sd_B, sd_d, rho_Bd)
    B <- Bd$B
    d <- Bd$d
    
    # Generate sample
    Z <- rbinom(n[j], 1, pZ)
    
    X <- matrix(rnorm(n[j] * K), n[j], K)
    Xd <- X %*% matrix(d, nrow = K, ncol = 1)
    XB <- X %*% matrix(B, nrow = K, ncol = 1)
    
    U <- rmvnorm(n[j], mean = rep(0, K), sigma = sigma)
    colnames(U) <- c("U1", "U0", "V")
    U <- data.frame(U)
    
    D <- as.numeric(Z * (d0 + Xd) >= U$V) 
    
    Y <- U$U0 + D * (B0 + XB + U$U1)
    
    # Build a data frame
    data <- data.frame(D = D, Z = Z, Y = Y)
    data <- cbind(data, X)
    names(data) <- c("D", "Z", "Y", paste0("X_", 1:K))
    
    # Generate interactions
    interactions <- apply(X, 2, function(X, Z){return(Z * X)}, Z = Z)
    colnames(interactions) <- paste("Z_X", 1:K, sep = "_")
    data <- cbind(data, interactions)
    
    #### MRC ESTIMATOR ####
    
    # Here theta = (a, d_0, d) (scaled to have length 1)
    varnames <- c(paste("X", 1:K, sep = "_"), "Z", paste("Z_X", 1:K, sep = "_"))
    
    # Get pilot estimates by (1) probit (2) logit (3) lpm
    probit_est <- glm(as.formula(paste0("D ~ ", paste(varnames, collapse = " + "))),
                      family = binomial(link = "probit"), 
                      data = data)$coefficients
    theta_probit <- probit_est[names(probit_est) %in% varnames]
    
    logit_est <- glm(as.formula(paste0("D ~ ", paste(varnames, collapse = " + "))),
                     family = binomial(link = "logit"), 
                     data = data)$coefficients
    theta_logit <- logit_est[names(logit_est) %in% varnames]
    
    lpm_est <- lm(as.formula(paste0("D ~ ", paste(varnames, collapse = " + "))),
                  data = data)$coefficients
    theta_lpm <- lpm_est[names(lpm_est) %in% varnames]
    
    # Run the Nelder-Mead-Head simplex from each pilot (probit, logit, lpm), 
    # random perturbations of each pilot, and random points, take the maximum
    
    # 3 pilot points
    n_perturb <- 5 # n_perturb random perturbations of each pilot point
    n_random <- 20 # n_random randomly generated starting points
    candidates <- vector(mode = "list", length = 3 + 3 * n_perturb + n_random)
    
    candidates[[1]] <- optim(par = theta_probit, 
                             fn = S_theta, 
                             method = "Nelder-Mead", 
                             data = data, 
                             varnames = varnames, 
                             control = list("fnscale" = -1))
    
    candidates[[2]] <- optim(par = theta_logit, 
                             fn = S_theta, 
                             method = "Nelder-Mead", 
                             data = data, 
                             varnames = varnames, 
                             control = list("fnscale" = -1))
    
    candidates[[3]] <- optim(par = theta_lpm, 
                             fn = S_theta, 
                             method = "Nelder-Mead", 
                             data = data, 
                             varnames = varnames, 
                             control = list("fnscale" = -1))
    
    # perturbations of theta probit
    theta_probit <- theta_probit / sqrt(sum(theta_probit ^ 2))
    for (p in 1:n_perturb){
      probit_perturbed <- theta_probit * 
        sample(c(0.9, 1.1), length(theta_probit), replace = TRUE)
      candidates[[3 + p]] <- optim(par = probit_perturbed, 
                                   fn = S_theta, 
                                   method = "Nelder-Mead", 
                                   data = data, 
                                   varnames = varnames, 
                                   control = list("fnscale" = -1))
    }
    
    # perturbations of theta logit
    theta_logit <- theta_logit / sqrt(sum(theta_logit ^ 2))
    for (p in 1:n_perturb){
      
      logit_perturbed <- theta_logit * 
        sample(c(0.9, 1.1), length(theta_logit), replace = TRUE)
      
      candidates[[3 + n_perturb + p]] <- optim(par = logit_perturbed, 
                                               fn = S_theta, 
                                               method = "Nelder-Mead", 
                                               data = data, 
                                               varnames = varnames, 
                                               control = list("fnscale" = -1))
    }
    
    # perturbations of theta lpm
    theta_lpm <- theta_lpm / sqrt(sum(theta_lpm ^ 2))
    for (p in 1:n_perturb){
      
      lpm_perturbed <- theta_lpm * 
        sample(c(0.9, 1.1), length(theta_lpm), replace = TRUE)
      
      candidates[[3 + (2 * n_perturb) + p]] <- 
        optim(par = lpm_perturbed, 
              fn = S_theta, 
              method = "Nelder-Mead", 
              data = data, 
              varnames = varnames, 
              control = list("fnscale" = -1))
    }
    
    # randomly generated starting points
    for (p in 1:n_random){
      
      theta_random <- runif(length(theta_probit), -1, 1) 
      theta_random <- theta_random / sqrt(sum(theta_random ^ 2))
      
      candidates[[3 + (3 * n_perturb) + p]] <- 
        optim(par = theta_random, 
              fn = S_theta, 
              method = "Nelder-Mead", 
              data = data, 
              varnames = varnames, 
              control = list("fnscale" = -1))
      
    }
    
    # Choose the best performing value
    i <- which.max(sapply(candidates, getElement, "value"))
    theta_hat <- candidates[[i]]$par
    names(theta_hat) <- names(candidates[[1]]$par)
    theta_hat <- theta_hat / sqrt(sum(theta_hat ^ 2)) # scale to be length 1
    
    s_hat <- sign(theta_hat["Z"] + 
                    X %*% matrix(theta_hat[colnames(interactions)], K, 1))
    s_star <- sign(d0 + Xd) # oracle
    s_probit <- sign(probit_est["Z"] + 
                       X %*% matrix(probit_est[colnames(interactions)], K ,1))
    s_logit <- sign(logit_est["Z"] + X %*% 
                      matrix(logit_est[colnames(interactions)], K ,1))
    s_lpm <- sign(lpm_est["Z"] + 
                    X %*% matrix(lpm_est[colnames(interactions)], K ,1))
    
    # Estimates
    Z_tilde <- Z - mean(data$Z)
    
    oracle[r] <- sum(Y * s_star * Z_tilde) / sum(D * s_star * Z_tilde)
    feasible[r] <- sum(Y * s_hat * Z_tilde) / sum(D * s_hat * Z_tilde)
    probit[r] <- sum(Y * s_probit * Z_tilde) / sum(D * s_probit * Z_tilde)
    logit[r] <- sum(Y * s_logit * Z_tilde) / sum(D * s_logit * Z_tilde)
    lpm[r] <- sum(Y * s_lpm * Z_tilde) / sum(D * s_lpm * Z_tilde)
    tsls[r] <- cov(Y, Z_tilde) / cov(D, Z_tilde)
    
    # Diagnostics
    diagnostics$s_cor[r] <- cor(s_hat, s_star)
    diagnostics$theta[ , r] <- c(rep(0, K), d0, d) / 
      sqrt(sum(c(rep(0, K), d0, d) ^ 2))
    diagnostics$theta_hat[, r] <- theta_hat 
    
    # Standard errors 
    data$Z_tilde <- Z_tilde
    data$s_hat <- s_hat
    se[r] <- calc_sd(data = data, tau_feasible =  feasible[r], n = n[j], K = K)
      
  }
  
  run_results <- list(oracle = oracle, 
                      feasible = feasible,
                      tsls = tsls, 
                      probit = probit,
                      logit = logit,
                      lpm = lpm,
                      diagnostics = diagnostics,
                      se = se)
  
  results[[j]] <- run_results
  
}

names(results) <- n

saveRDS(results, file = paste0(path, "/Monte_Carlo/results_switchers.rds"))

############################## SUMMARY STATISTICS ##############################

# Get the mean, standard error, MSE, median, 0.25 and 0.75 quartiles of the 
# sampling distribution of the oracle and feasible estimators. 

tau <- B0 + switchers_calibration$targets[[1]]$tau_U + 
  switchers_calibration$targets[[1]]$tau_B

feasible_sum_stats <- matrix(0, length(n), 6)
oracle_sum_stats <- matrix(0, length(n), 6)
probit_sum_stats <- matrix(0, length(n), 6)

for (i in 1:length(n)){
  
  feasible_sum_stats[i, ] <- get_stats(results = results[[i]]$feasible, tau = tau)
  oracle_sum_stats[i, ] <- get_stats(results = results[[i]]$oracle, tau = tau)
  probit_sum_stats[i, ] <- get_stats(results = results[[i]]$probit, tau = tau)
  
}

rownames(feasible_sum_stats) <- n
colnames(feasible_sum_stats) <- 
  c("Mean", "Standard Error", "MSE", "Median", "First Quartile", "Third Quartile")
rownames(oracle_sum_stats) <- n
colnames(oracle_sum_stats) <- 
  c("Mean", "Standard Error", "MSE", "Median", "First Quartile", "Third Quartile")
rownames(probit_sum_stats) <- n
colnames(probit_sum_stats) <- 
  c("Mean", "Standard Error", "MSE", "Median", "First Quartile", "Third Quartile")

xtable(feasible_sum_stats, digits = 3)
xtable(oracle_sum_stats, digits = 3)
xtable(probit_sum_stats, digits = 3)

# Add in coverage rates 
coverage_feasible <- matrix(0, length(n), 3)
coverage_oracle <- matrix(0, length(n), 3)

# n varies over rows, columns are coverage rates for 90%, 95% and 99% CIs
for (i in 1:length(n)){
  
  coverage_feasible[i, ] <- get_coverage(results = results[[i]]$feasible,
                                         se = results[[i]]$se,
                                         tau = tau, n = n[i])
  coverage_oracle[i, ] <- get_coverage(results = results[[i]]$oracle,
                                       se = results[[i]]$se,
                                       tau = tau, n = n[i])
  
}

colnames(coverage_feasible) <- c("90", "95", "99")
rownames(coverage_feasible) <- n

colnames(coverage_oracle) <- c("90", "95", "99")
rownames(coverage_oracle) <- n

xtable(cbind(coverage_feasible, coverage_oracle), digits = 3)

######################## PLOT OF SAMPLING DISTRIBUTIONS ########################

# Prepare data for plotting

# convert to sqrt(n) * (tau - tau^*)
tau_feasible <- matrix(0, length(n), runs)
tau_oracle <- matrix(0, length(n), runs)

for (i in 1:length(n)){
  
  tau_feasible[i, ] <- sqrt(n[i]) * (results[[i]]$feasible - tau)
  tau_oracle[i, ] <- sqrt(n[i]) * (results[[i]]$oracle - tau)
  
}

x_lims <- rep(0, 2)
x_lims[1] <- quantile(c(tau_feasible, tau_oracle), 0.01, na.rm = TRUE)
x_lims[2] <- quantile(c(tau_feasible, tau_oracle), 0.99, na.rm = TRUE)

make_plot(tau_oracle = tau_oracle[1, ],
          tau_feasible = tau_feasible[1, ],
          x_lims = x_lims, 
          trimming_flag = TRUE)
make_plot(tau_oracle = tau_oracle[2, ],
          tau_feasible = tau_feasible[2, ],
          x_lims = x_lims)
make_plot(tau_oracle = tau_oracle[3, ],
          tau_feasible = tau_feasible[3, ],
          x_lims = x_lims)
make_plot(tau_oracle = tau_oracle[4, ],
          tau_feasible = tau_feasible[4, ],
          x_lims = x_lims)





