#################################### SUMMARY ################################### 

# This script chooses parameter values for the Monte Carlo exercise

# The first stage is assumed to be
# D = 1(z * (d0 + x'd)  >=  V)
# where the marginal distribution of V is N(0, 1), so the propensity score is
# p(x,z) = Phi(z * (d0 + x'd))
# d0 + x'd controls the sign of the relationship between Z and D. 

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
targets[[1]] <- list(compliers = 0.11,
                     defiers = 0.09,
                     tau_B = -0.05,
                     tau_U = -0.05)
targets[[2]] <- list(compliers = 0.06,
                     defiers = 0.04,
                     tau_B = -0.05,
                     tau_U = -0.05)

e_tol <- 0.00000001 # error tolerance
# search bounds for use in calibration block 1 
bounds <- list(top =  c(0.1, 0.08), bottom = c(0.05, 0.01))

############################## CALIBRATION BLOCK 1 #############################

# Choose d0 and sqrt(d'd) (sd_d) to attain the desired shares of compliers and 
# defiers

# integrate this kernel from -d0 to Inf to get the share of compliers
# integrate -1 times this kernel from -Inf to -d0 to get the share of defiers
kernel <- function(z, d0, sd_d)(
  
  return((pnorm(d0 + z) - pnorm(0)) * dnorm(z, mean = 0, sd = sd_d))
  
)

# Loop to solve for sd_d which meets the compliers target
get_sd_d <- function(d0, e_tol, targets){
  
  # initial guess
  n_points <- 101
  sd_d <- seq(0.01, 1, length.out = n_points)
  compliers <- rep(0, n_points)
  
  for (i in 1:n_points){
    
    compliers[i] <- integrate(kernel, lower = -d0, upper = Inf, 
                              d0 = d0, sd_d = sd_d[i])$value
    
  }
  
  index <- which.min(abs(compliers - targets$compliers))
  upper <- sd_d[min(n_points, index + 1)]
  lower <- sd_d[max(1, index - 1)]
  error <- abs(compliers[index] - targets$compliers)
  
  # quick loop to improve
  while(error > e_tol){
    
    sd_d <- seq(lower, upper, length.out = n_points)
    for (i in 1:n_points){
      
      compliers[i] <- integrate(kernel, lower = -d0, upper = Inf, 
                                d0 = d0, sd_d = sd_d[i])$value
      
    }
    
    index <- which.min(abs(compliers - targets$compliers))
    upper <- sd_d[min(n_points, index + 1)]
    lower <- sd_d[max(1, index - 1)]
    error <- abs(compliers[index] - targets$compliers)
    
  }
  
  # final answer
  return(list(sd_d = sd_d[index], compliers = compliers[index]))
 
}

# Split the difference type algorithm to find d0 and sd_d 
# Assume that desired compliers share > desired defiers share
# Loop is as follows:
# * Take d0 as fixed
# * Choose sd_d to attain the desired compliers shares
# * For this (d0, sd_d) pair, evaluate the defier shares
# * If the defiers share is too small, d0 is too large -> revise upper bound
# * If the defiers share is too large, d0 is too small -> revise lower bound

for (i in 1:length(targets)){
  
  # upper and lower bounds of the search region
  top <- bounds$top[i]
  bottom <- bounds$bottom[i]
  
  count <- 0
  limit <- 1000 # maximum number of iterations
  
  while (count < limit){
    
    mid <- (top + bottom) / 2
    
    # solve for sd_d that matches compliers
    sd_d <- get_sd_d(d0 = mid, e_tol = e_tol, targets = targets[[i]])$sd_d
    
    # check defiers
    defiers <- - integrate(kernel, lower = - Inf, upper = -mid, 
                           d0 = mid, sd_d = sd_d)$value
    
    if (abs(defiers - targets[[i]]$defiers) < e_tol){
      print("converged")
      break
    }else{
      if (defiers > targets[[i]]$defiers){
        bottom <- mid
      }else{
        top <- mid
      }
    }
    
    count <- count + 1
    
  }
  
  parameters[[i]]$d0 <- mid
  parameters[[i]]$sd_d <- sd_d
  
}

############################## VALIDATE PARAMETERS #############################

# Check if these parameters give the correct shares for complier and defiers

K <- 3 # number of X variables
n <- 100000

results <- vector(mode = "list", length = length(parameters))

for (i in 1:length(parameters)){
  
  d <- runif(K, min = -1, max = 1)
  d <- d * parameters[[i]]$sd_d / sqrt(sum(d * d))
  
  X <- matrix(rnorm(n * K), nrow = n, ncol = K)
  Xd <- X %*% matrix(d, nrow = K, ncol = 1)
  
  comparisons <- matrix(0, 3, 2)
  
  comparisons[, 1] <- c(targets[[i]]$compliers + targets[[i]]$defiers,
                        targets[[i]]$compliers,
                        targets[[i]]$defiers)
  comparisons[1, 2] <- mean(abs(pnorm(parameters[[i]]$d0 + Xd) - pnorm(0)))
  comparisons[2, 2] <- mean(((parameters[[i]]$d0 + Xd) > 0) * 
                              abs(pnorm(parameters[[i]]$d0 + Xd) - pnorm(0)))
  comparisons[3, 2] <- mean(((parameters[[i]]$d0 + Xd) < 0) * 
                              abs(pnorm(parameters[[i]]$d0 + Xd) - pnorm(0)))
  
  rownames(comparisons) <- c("switchers", "compliers", "defiers")
  results[[i]] <- comparisons
  
}

results[[1]]
results[[2]]


############################## CALIBRATION BLOCK 2 #############################

# Given d0 and sd_d, choose B'd to target a value of the asymptotic estimand. 
# We have the formula
# tau = B0 + (B'd / d'd) * term_two + (rho_1v - rho_0v) * term_three
# where term_two captures the effect of XB and term three the effect of U. We 
# can simply evaluate these integrals and then choose the multiplying constants.

# kernel_one is s(X) * D * (Z - E[Z|X])
kernel_one <- function(z, d0, sd_d){
  
  return((2 * (d0 + z >= 0) - 1) * (pnorm(d0 + z) - pnorm(0)) * 
           dnorm(z, mean = 0, sd = sd_d))
  
}

# kernel_two is s(X) * D * (Z - E[Z|X]) * Xd 
kernel_two <- function(z, d0, sd_d){
  
  return((2 * (d0 + z >= 0) - 1) * (pnorm(d0 + z) - pnorm(0)) * 
           z * dnorm(z, mean = 0, sd = sd_d))
  
}

# kernel_three is s(X) * ((1 - D) * U0 + D * U1) * (Z - E[Z|X])
kernel_three <- function(z, d0, sd_d){
  
  return((2 * (d0 + z >= 0) - 1) * (dnorm(0) - dnorm(d0 + z)) * 
           dnorm(z, mean = 0, sd = sd_d))
  
}

for (i in 1:length(targets)){
  
  denom <- integrate(kernel_one, lower = -Inf, upper = Inf, 
                     d0 = parameters[[i]]$d0,
                     sd_d = parameters[[i]]$sd_d)$value
  
  term_two <- integrate(kernel_two, lower = -Inf, upper = Inf, 
                        d0 = parameters[[i]]$d0,
                        sd_d = parameters[[i]]$sd_d)$value / denom
  
  term_three <- integrate(kernel_three, lower = -Inf, upper = Inf, 
                          d0 = parameters[[i]]$d0,
                          sd_d = parameters[[i]]$sd_d)$value / denom
  
  # chi_two = ||B|| * B'd / ||d||
  chi_two <- targets[[i]]$tau_B / term_two
  rho_Bd <- - 0.5
  sd_B <- chi_two * parameters[[i]]$sd_d / rho_Bd
  
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
  
  # # Randomly draw Sigma
  # A <- matrix(rnorm(9), nrow = 3, ncol = 3)
  # sigma <- t(A) %*% A
  # sd_V <- sqrt(sigma[3, 3])
  # sigma[3, ] <- sigma[3, ] / sd_V
  # sigma[, 3] <- sigma[, 3] / sd_V
  # chi_2 <- sigma[3, 1] - sigma[3, 2]
  parameters[[i]]$sd_B <- sd_B
  parameters[[i]]$rho_Bd <- rho_Bd
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

results <- vector(mode = "list", length = 2)

for (i in 1:length(parameters)){
  
  analytical <- rep(0, 6)
  names(analytical) <- c("denom", "num_1", "num_2", "num_3", "num", "tau")
  analytical[1] <- Var_Z * 
    integrate(kernel_one, lower = -Inf, upper = Inf, 
              d0 = parameters[[i]]$d0, 
              sd_d = parameters[[i]]$sd_d)$value
  analytical[2] <- Var_Z * B0 * 
    integrate(kernel_one, lower = -Inf, upper = Inf, 
              d0 = parameters[[i]]$d0, 
              sd_d = parameters[[i]]$sd_d)$value
  analytical[3] <- Var_Z * 
    (parameters[[i]]$sd_B * parameters[[i]]$rho_Bd / parameters[[i]]$sd_d) *
    integrate(kernel_two, lower = -Inf, upper = Inf, 
              d0 = parameters[[i]]$d0, 
              sd_d = parameters[[i]]$sd_d)$value 
  analytical[4] <- Var_Z * 
    (parameters[[i]]$sigma[3, 1] - parameters[[i]]$sigma[3, 2]) * 
    integrate(kernel_three, lower = -Inf, upper = Inf, 
              d0 = parameters[[i]]$d0, 
              sd_d = parameters[[i]]$sd_d)$value    
  analytical[5] <- sum(analytical[2:4])
  analytical[6] <- analytical[5] / analytical[1]
  
  data <- matrix(0, runs, 6)
  colnames(data) <- c("denom", "num_1", "num_2", "num_3", "num", "tau")
  
  for (j in 1:runs){
    
    # Generate B and d subject to constraints
    # Construct orthogonal vectors 
    e_1 <- runif(K, min = - 1, max = 1)
    u <- runif(K, min = - 1, max = 1)
    e_2 <- u - sum(e_1 * u) * e_1 / (sum(e_1 ^ 2))
    e_1 <- e_1 / sqrt(sum(e_1 ^ 2))
    e_2 <- e_2 / sqrt(sum(e_2 ^ 2))
    
    # Use e_1 and e_2 to build d and B
    d <- e_1 * parameters[[i]]$sd_d
    B <- parameters[[i]]$sd_B * 
      (parameters[[i]]$rho_Bd * e_1 + sqrt(1 - parameters[[i]]$rho_Bd ^ 2) * e_2)
    
    X <- matrix(rnorm(n * K), nrow = n, ncol = K)
    Xd <- X %*% as.matrix(d, nrow = K, ncol = 1)
    XB <- X %*% as.matrix(B, nrow = K, ncol = 1)
    
    s <- 2 * as.numeric(parameters[[i]]$d0 + Xd >= 0) - 1
    
    Z <- rbinom(n, 1, pZ)
    Z_tilde <- Z - mean(Z)
    
    U <- rmvnorm(n, mean = rep(0, K), sigma = parameters[[i]]$sigma)
    colnames(U) <- c("U1", "U0", "V")
    U <- data.frame(U)
    
    D <- as.numeric(Z * (parameters[[i]]$d0 + Xd) >= U$V)
    
    Y <- (1 - D) * U$U0 + D * (B0 + XB + U$U1)
    
    # Denominator
    data[j, 1] <- mean(s * D * Z_tilde)
    
    # Numerator
    data[j, 2] <- B0 * mean(s * D * Z_tilde)
    data[j, 3] <- mean(s * D * XB * Z_tilde)
    data[j, 4] <- mean(s * ((1 - D) * U$U0 + D * U$U1) * Z_tilde)
    data[j, 5] <- sum(data[j, 2:4])
    
    # Estimand
    data[j, 6] <- data[j, 5] / data[j, 1]
    
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
  plot(results[[i]]$data[, 5], rep(results[[i]]$analytical[5], runs))
  mean(results[[i]]$data[, 5])
  results[[i]]$analytical[5]
  
  # numerator part 1
  plot(results[[i]]$data[, 2], rep(results[[i]]$analytical[2], runs))
  mean(results[[i]]$data[, 2])
  results[[i]]$analytical[2]
  
  # numerator part 2
  plot(results[[i]]$data[, 3], rep(results[[i]]$analytical[3], runs))
  mean(results[[i]]$data[, 3])
  results[[i]]$analytical[3]
  
  # numerator part 3
  plot(results[[i]]$data[, 4], rep(results[[i]]$analytical[4], runs))
  mean(results[[i]]$data[, 4])
  results[[i]]$analytical[4]
  
  # tau
  plot(results[[i]]$data[, 6], rep(results[[i]]$analytical[6], runs))
  mean(results[[i]]$data[, 6])
  results[[i]]$analytical[6]
  
}

############################### SWITCHER DENSITY ###############################

# Given these parameters, plot the density of switchers (compliers and defiers)
n_points <- 1001
C <- seq(-1/2, 1/2, length.out = n_points)

filenames <- paste0("switcher_dens_parameters_", 1:length(parameters))

for (i in 1:length(parameters)){
  
  dens <- rep(0, n_points)
  dist <- rep(0, n_points)
  
  for (j in 1:n_points){
    
    dist[j] <- pnorm(qnorm(pnorm(0) + C[j]) - parameters[[i]]$d0, 
                     mean = 0, sd = parameters[[i]]$sd_d) 
    dens[j] <-  dnorm(qnorm(pnorm(0) + C[j]) - parameters[[i]]$d0, 
                      mean = 0, sd = parameters[[i]]$sd_d) / 
      dnorm(qnorm(pnorm(0) + C[j]) - parameters[[i]]$d0)
    
  }
  
  plot(C, dist)
  plot(C, dens)
  
  # verify by simulation
  K <- 3
  n <- 100000
  d <- runif(K, min = - 1, max  = 1)
  d <- d * parameters[[i]]$sd_d / sqrt(sum(d ^ 2))
  X <- matrix(rnorm(n * K), nrow = n, ncol = K)
  Xd <- X %*% as.matrix(d, nrow = K, ncol = 1)
  Delta <- pnorm(Xd + parameters[[i]]$d0) - pnorm(0)
  
  emp_dist <- ecdf(Delta)
  data_dist <- rep(0, n_points)
  
  for (j in 1:n_points){
    
    data_dist[j] <- emp_dist(C[j])
    
  }
  
  plot_data <- data.frame(C = C, data = data_dist, analytical = dist)
  plot_data %>% ggplot() + 
    geom_line(aes(x = C, y = data, col = "data")) + 
    geom_line(aes(x = C, y = analytical, col = "analytical")) + 
    theme_bw() + 
    theme(legend.title = element_blank())
  
  data_dens <- density(Delta)
  plot(data_dens)
  lines(C, dens, lty = 2)
  
  dens[1] <- 0
  dens[n_points] <- 0
  
  # Nice plot of the density of switchers
  pdf(file = paste0(filenames[i], ".pdf"), width = 8.98, height = 5.60)
  plot_data <- data.frame(C = C, dens = dens)
  print(plot_data %>% ggplot() + 
    geom_line(aes(x = C, y = dens, col = "density")) + 
    geom_vline(aes(xintercept = pnorm(parameters[[i]]$d0) - pnorm(0), col = "d0"),
               linetype = "dotted") + 
    theme_bw() + 
    scale_color_manual(values = c("density" = "#F8766D",
                                  "d0" = "black")) + 
    labs(x = "p(x,1) - p(x,0)", y = "") +
    theme(legend.position = "null"))
  dev.off()
  
}

################################## SAVE OUTPUT #################################

calibration <- list(targets = targets, 
                    parameters = parameters)
saveRDS(calibration, 
        file = paste0(path, "/Monte_Carlo/switchers_calibration.rds"))
