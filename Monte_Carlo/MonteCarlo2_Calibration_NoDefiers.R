#################################### SUMMARY ################################### 

# This script chooses parameter values for the Monte Carlo exercise

# The first stage is assumed to be
# D = 1(z * d0  >=  V)
# where the marginal distribution of V is N(0, 1), so the propensity score is
# p(x,z) = Phi(z * d0)
# The sign of d0 determines whether compliers or defiers are present

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

################################### SETTINGS ###################################

parameters <- vector(mode = "list", length = 2)
parameters[[1]] <- list()
parameters[[2]] <- list()

targets <- vector(mode = "list", length = 2)
targets[[1]] <- list(compliers = 0.1,
                     tau_U = -0.1)
targets[[2]] <- list(compliers = 0.05,
                     tau_U = -0.1)

############################## CALIBRATION BLOCK 1 #############################

# Choose d0 to attain the desired share of compliers 

for (i in 1:length(targets)){
  
  parameters[[i]]$d0 <- qnorm(targets[[i]]$compliers + pnorm(0))
    
}

############################## CALIBRATION BLOCK 2 #############################

# Given d0, target a value of the asymptotic estimand. 
# We have the formula
# tau = B0 + (rho_1v - rho_0v) * term_three
# where term_three captures the effect of U. Evaluate this and choose the 
# multiplying constant.

for (i in 1:length(targets)){
  
  term_three <- (dnorm(0) - dnorm(parameters[[i]]$d0)) / 
    (pnorm(parameters[[i]]$d0) - pnorm(0))
  
  # chi_three = sd_1 rho_1V - sd_0 rho_0V
  chi_three <- targets[[i]]$tau_U / term_three
  
  # fix correlations and adjst standard deviations to match the target
  sigma <- diag(1, 3)
  sigma[2, 1] <- -0.25
  sigma[3, 1] <- -0.5
  sigma[3, 2] <- 0.5
  sigma[1, 2] <- -0.25
  sigma[1, 3] <- -0.5
  sigma[2, 3] <- 0.5
  
  # assume sd_U is common for U(1) and U(0). Back it out.
  sd_U <- chi_three / (sigma[3, 1] - sigma[3, 2])
  # multiply by sd_U so that sigma is a covariance
  sigma[, 1:2] <- sigma[, 1:2] * sd_U
  sigma[1:2, ] <- sigma[1:2, ] * sd_U
  
  parameters[[i]]$sigma <- sigma
  
}

############################## VALIDATE PARAMETERS #############################

# Check if parameters give the correct asymptotic estimand
K <- 3
n <- 100000
runs <- 100

# other parameters
pZ <- 1 / 2
Var_Z <- pZ * (1 - pZ)
B0 <- 0
sd_B <- 1

results <- vector(mode = "list", length = 2)

for (i in 1:length(parameters)){
  
  analytical <- rep(0, 5)
  names(analytical) <- c("denom", "num_1", "num_2", "num", "tau")
  analytical[1] <- Var_Z * (pnorm(parameters[[i]]$d0) - pnorm(0))
  analytical[2] <- Var_Z * B0 * (pnorm(parameters[[i]]$d0) - pnorm(0))
  analytical[3] <- Var_Z * 
    (parameters[[i]]$sigma[3, 1] - parameters[[i]]$sigma[3, 2]) * 
    (dnorm(0) - dnorm(parameters[[i]]$d0))
  analytical[4] <- sum(analytical[2:3])
  analytical[5] <- analytical[4] / analytical[1]
  
  data <- matrix(0, runs, 5)
  colnames(data) <- c("denom", "num_1", "num_2", "num", "tau")
  
  for (j in 1:runs){
    
    B <- runif(K, min = -1, max = 1)
    B <- B * sd_B / sqrt(sum(B ^ 2))
      
    X <- matrix(rnorm(n * K), nrow = n, ncol = K)
    XB <- X %*% as.matrix(B, nrow = K, ncol = 1)
    
    Z <- rbinom(n, 1, pZ)
    Z_tilde <- Z - mean(Z)
    
    U <- rmvnorm(n, mean = rep(0, K), sigma = parameters[[i]]$sigma)
    colnames(U) <- c("U1", "U0", "V")
    U <- data.frame(U)
    
    D <- as.numeric(Z * parameters[[i]]$d0 >= U$V)
    
    Y <- (1 - D) * U$U0 + D * (B0 + XB + U$U1)
    
    # Denominator
    data[j, 1] <- mean(D * Z_tilde)
    
    # Numerator
    data[j, 2] <- B0 * mean(D * Z_tilde)
    data[j, 3] <- mean(((1 - D) * U$U0 + D * U$U1) * Z_tilde)
    data[j, 4] <- sum(data[j, 2:3])
    
    # Estimand
    data[j, 5] <- data[j, 4] / data[j, 1]
    
  }
  
  results[[i]] <- list(analytical = analytical,
                       data = data)
  
}

for (i in 1:length(parameters)){
  
  # denominator
  plot(results[[i]]$data[, 1], rep(results[[i]]$analytical[1], runs))
  mean(results[[i]]$data[, 1])
  results[[i]]$analytical[1]
  
  # numerator
  plot(results[[i]]$data[, 4], rep(results[[i]]$analytical[4], runs))
  mean(results[[i]]$data[, 4])
  results[[i]]$analytical[4]
  
  # numerator part 1
  plot(results[[i]]$data[, 2], rep(results[[i]]$analytical[2], runs))
  mean(results[[i]]$data[, 2])
  results[[i]]$analytical[2]
  
  # numerator part 2
  plot(results[[i]]$data[, 3], rep(results[[i]]$analytical[3], runs))
  mean(results[[i]]$data[, 3])
  results[[i]]$analytical[3]
  
  # tau
  plot(results[[i]]$data[, 5], rep(results[[i]]$analytical[5], runs))
  mean(results[[i]]$data[, 5])
  results[[i]]$analytical[5]
  
}

################################## SAVE OUTPUT #################################

calibration <- list(targets = targets, 
                    parameters = parameters)
saveRDS(calibration, 
        file = paste0(path, "/Monte_Carlo/no_defiers_calibration.rds"))
