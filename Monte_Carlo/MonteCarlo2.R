#################################### SUMMARY ################################### 

# This scipt performs Monte Carlo for the feasible reordered IV estimator

# The first stage is assumed to be
# D = 1(z * x'd - V > 0)
# where the marginal distribution of V is N(0, 1), so the propensity score is
# p(x,z) = Phi(z * x'd)
# x'd controls the sign of the relationship between Z and D. 

# We use the maximum rank correlation estimator to obtain estimates of the first
# stage parameters. This estimator does not estimate a_0, and the remaining 
# parameters are only estimated up to scale (adding a constant to all terms 
# does not perturb ranks, and multiplying all non-constant terms by a positive
# constant would also not affect ranks)

# The second stage is 
# Y(0) = U(0)
# Y(1) = B_0 + f(x'B) + U(1)
# Y = Y(0) + D(Y(1) - Y(0))
# where f() is a function which varies between specifications and
# (V, U(1), U(0)) are jointly normal. 

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

################################## CALIBRATION #################################

# Find constraints on d'd that give the desired number of compliers and defiers

kernel <- function(z, sd_d)(
  
  return((pnorm(z) - pnorm(0)) * dnorm(z, mean = 0, sd = sd_d))
  
)

# initial guess
n_points <- 101
sd_d <- seq(0.01, 5, length.out = n_points)
compliers <- rep(0, n_points)

for (i in 1:n_points){
  
  compliers[i] <- integrate(kernel, lower = 0, upper = Inf, sd_d = sd_d[i])$value
  
}

target <- 0.1
e_tol <- 0.00000001 # error tolerance

index <- which.min(abs(compliers - target))
upper <- sd_d[min(n_points, index + 1)]
lower <- sd_d[max(1, index - 1)]
error <- abs(compliers[index] - target)

# quick loop to improve
while(error > e_tol){
  
  sd_d <- seq(lower, upper, length.out = n_points)
  for (i in 1:n_points){
    
    compliers[i] <- integrate(kernel, lower = 0, upper = Inf, sd_d = sd_d[i])$value
    
  }
  
  index <- which.min(abs(compliers - target))
  upper <- sd_d[min(n_points, index + 1)]
  lower <- sd_d[max(1, index - 1)]
  error <- abs(compliers[index] - target)
  
}

# final answer
sd_d <- sd_d[index]

############################## VALIDATE PARAMETERS #############################

# Switcher share
K <- 3 # number of X variables
n <- 100000
X <- matrix(rnorm(n * K), nrow = n, ncol = K)

d <- runif(K, min = -1, max = 1)
d <- d * sd_d / sqrt(sum(d * d))

mean(abs(pnorm(X %*% matrix(d, nrow = K, ncol = 1)) - pnorm(0)))
2 * target

################################## CHECK LATES #################################

K <- 3
d <- runif(K, min = - 1, max = 1)
d <- d / sqrt(sum(d ^ 2))

# # Randomly generate a covariance matrix for the unobservables
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

# prespecified correlation matrix, variable order is (U1, U0, V)
sigma <- diag(1, K)
sigma[2, 1] <- -0.25
sigma[3, 1] <- -0.5
sigma[3, 2] <- 0.5
sigma[1, 2] <- -0.25
sigma[1, 3] <- -0.5
sigma[2, 3] <- 0.5

results <- list()

# randomly draw X
X <- rnorm(3, mean = 0, sd = 1)
Xd <- sum(X * d)
results$Xd = Xd

# LATE by simulation
n <- 100000

U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
colnames(U) <- c("U1", "U0", "V")
U <- data.frame(U)

Z <- rbinom(n, 1, 1 / 2)
D <- (Z * Xd >= U$V)
Y <- U$U0 + D * (U$U1 - U$U0)

data <- data.frame(Y = Y, D = D, Z = Z)

# data LATE
cov(Y, Z) / cov(D, Z)

D_hat <- lm(D ~ Z, data = data)$fitted.values
results$data <- lm(Y ~ D_hat)$coefficients[2]

# analytical LATE
results$analytical <- 
  (sigma[3, 1] - sigma[3, 2]) * (dnorm(0) - dnorm(Xd)) / (pnorm(Xd) - pnorm(0))

results

############################### SWITCHER DENSITY ###############################

n_points <- 101
C <- seq(0, 1/2, length.out = n_points)
dens <- rep(0, n_points)
dist <- rep(0, n_points)

sd_d <- 0.25

for (i in 1:n_points){
  
  dist[i] <- pnorm(qnorm(pnorm(0) + C[i]), mean = 0, sd = sd_d) -
    pnorm(qnorm(pnorm(0) - C[i]), mean = 0, sd = sd_d)
  dens[i] <- (1 / dnorm(qnorm(pnorm(0) + C[i]))) * 
    dnorm(qnorm(pnorm(0) + C[i]), mean = 0, sd = sd_d) + 
    (1 / dnorm(qnorm(pnorm(0) - C[i]))) *
    dnorm(qnorm(pnorm(0) - C[i]), mean = 0, sd = sd_d)
  
}

plot(C, dist)
plot(C, dens)

# verify by simulation
K <- 3
n <- 100000
d <- runif(K, min = - 1, max  = 1)
d <- d * sd_d / sqrt(sum(d ^ 2))
X <- matrix(rnorm(n * K), nrow = n, ncol = K)
Xd <- X %*% as.matrix(d, nrow = K, ncol = 1)
Delta <- abs(pnorm(Xd) - pnorm(0))

emp_dist <- ecdf(Delta)
data_dist <- rep(0, n_points)

for (i in 1:n_points){
  
  data_dist[i] <- emp_dist(C[i])
  
}

plot_data <- data.frame(C = C, data = data_dist, analytical = dist)
plot_data %>% ggplot() + 
  geom_line(aes(x = C, y = data, col = "data")) + 
  geom_line(aes(x = C, y = analytical, col = "analytical")) + 
  theme_bw() + 
  theme(legend.title = element_blank())

data_dens <- density(Delta)
plot(data_dens)
plot(C, dens)

########################## EXPECTED TREATMENT EFFECTS ##########################

# Get the expected treatment effect for each mass of switchers

# prespecified correlation matrix, variable order is (U1, U0, V)
sigma <- diag(1, 3)
sigma[2, 1] <- -0.25
sigma[3, 1] <- -0.5
sigma[3, 2] <- 0.5
sigma[1, 2] <- -0.25
sigma[1, 3] <- -0.5
sigma[2, 3] <- 0.5

rho_1V <- sigma[3, 1]
rho_0V <- sigma[3, 2]

K <- 3 # number of controls
n_points <- 101 # number of switcher values where we evaluate E[tau]

# Randomly generate parameter vectors 
d <- runif(K, -1, 1)
d <- d / sqrt(sum(d ^ 2))
B <- runif(K, -1, 1)
B <- B / sqrt(sum(B ^ 2))
B_0 <- 1 # ATE

# Values of switcher mass
C <- seq(-1 / 2, 1 / 2, length.out = n_points)
E_tau <- rep(0, n_points) # expected treatment effect

sd_d <- sqrt(sum(d ^ 2))
sd_B <- sqrt(sum(B ^ 2))
rho_Bd <- sum(B * d) / (sd_d * sd_B)

for (i in 1:n_points){
  
  # implied value of X'd from C
  Xd <- qnorm(pnorm(0) + C[i])
  
  # evaluate the expected treatment effect
  if (Xd == 0){
    E_tau[i] <- B_0
  }else{
    E_tau[i] <- B_0 + sd_d * rho_Bd * Xd  / sd_B + 
      (rho_1V - rho_0V) * (dnorm(0) - dnorm(Xd)) / (pnorm(Xd) - pnorm(0))
  }

}

# validate using data
n <- 100000

X <- matrix(rnorm(n *K), nrow = n, ncol = K)
Xd <- X %*% as.matrix(d, nrow = K, ncol = 1)
XB <- X %*% as.matrix(B, nrow = K, ncol = 1)

Z <- rbinom(n, 1, 1 / 2)

U <- rmvnorm(n, mean = rep(0, K), sigma = sigma)
colnames(U) <- c("U1", "U0", "V")
U <- data.frame(U)
Y_0 <- U$U0
Y_1 <- B_0 + XB + U$U1
D <- (Z * Xd >= U$V)
Y <- Y_0 + D * (Y_1 - Y_0)

data <- data.frame(Y = Y, Z = Z, D = D, Xd = Xd, XB = XB)
data$C <- pnorm(Xd) - pnorm(0)

# split into bins, get an average for each bin
C_bins <- (C[1:(n_points - 1)] + C[2:n_points]) / 2
C_bins <- c(-1 / 2, C_bins, 1 / 2)

tau_data <- rep(0, n_points)

for (i in 1:n_points){
  
  temp_data <- 
    data[as.logical((data$C >= C_bins[i]) * (data$C < C_bins[i + 1])), ]
  
  check <- (sum(temp_data$Z) > 0) * (sum(1 - temp_data$Z) > 0)
  
  if (check == 1){
    
    temp_data$D_hat <- lm(D ~ Z, data = temp_data)$fitted.values
    tau_data[i] <- lm(Y ~ D_hat, data = temp_data)$coefficients[2]
    
  }
  
}

plot(C, pmax(-3, pmin(tau_data, 3)))
plot(C, E_tau)

# Pick one point and evaluate
n_draws <- 100
comparison <- matrix(0, nrow = n_draws, ncol = 3)
colnames(comparison) <- c("Xd", "data", "analytical")
comparison <- data.frame(comparison)

for (i in 1:n_draws){
  
  X <- rnorm(K)
  Xd <- sum(X * d)
  comparison$Xd[i] <- Xd
  
  lower <- min(0, Xd)
  upper <- max(0, Xd)
  include <- as.logical((U$V > lower) * (U$V < upper))
  
  XB <- rnorm(n, mean = sd_B * rho_Bd * Xd / sd_d, sd = sd_B * sqrt(1 - rho_Bd ^ 2))
  
  comparison$data[i] <- B_0 + mean(XB) + mean((U$U1 - U$U0)[include])
  comparison$analytical[i] <- B_0 + sd_d * rho_Bd * Xd  / sd_B + 
    (rho_1V - rho_0V) * (dnorm(0) - dnorm(Xd)) / (pnorm(Xd) - pnorm(0)) 
  
}

plot(comparison$analytical, comparison$data)


################################### FUNCTIONS ##################################

# Generate parameter vectors a and d which match sd_a, sd_d, and rho_ad
gen_ad <- function(k, sd_a, sd_d, rho_ad){
  
  # create orthonormal vectors e_1 and e_2
  e_1 <- runif(k, -1, 1)
  e_1 <- e_1 / sqrt(sum(e_1 ^ 2))
  u <- runif(k, -1, 1)
  e_2 <- u - (sum(e_1 * u)) * e_1 / (sum(e_1 ^ 2))
  e_2 <- e_2 / sqrt(sum(e_2 ^ 2))
  
  # a is a scaled version of e_1
  a <- sd_a * e_1
  # d is a combination of e_1 and e_2
  d <- sd_d * rho_ad * e_1 + sd_d * sqrt(1 - rho_ad ^ 2) * e_2
  
  return(list(a = a,
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

# old_S_theta <- function(candidates, data, varnames){
#   
#   index <- as.matrix(data[varnames]) %*% candidates
#   D_sort <- data$D[order(index)] 
#   tally <- cumsum(1 - D_sort)
#   val <- sum(tally[as.logical(D_sort)])
#   
#   return(val)
#   
# }

# Estimation function written so that it can be passed to the boot function
# * varnames is a list of variables (X_1,..,X_k,Z,Z * X_1,...,Z * X_k)
# * data is the sample
# * indices is included to allow the boot function to operate
# The fun
est_theta <- function(data, indices, varnames){
  
  boot_data <- data[indices, ] 
  model <- as.formula(paste0("D ~ ", paste(varnames, collapse = " + ")))
  
  probit_est <- glm(model,
                    family = binomial(link = "probit"), 
                    data = boot_data)$coefficients
  theta_probit <- probit_est[names(probit_est) %in% varnames]
  
  logit_est <- glm(model,
                   family = binomial(link = "logit"), 
                   data = boot_data)$coefficients
  theta_logit <- logit_est[names(logit_est) %in% varnames]
  
  lpm_est <- lm(model,
                data = boot_data)$coefficients
  theta_lpm <- lpm_est[names(lpm_est) %in% varnames]
  
  # Run the Nelder-Mead-Head simplex from each pilot, take the maximum
  candidates <- vector(mode = "list", length = 3)
  
  candidates[[1]] <- optim(par = theta_probit, 
                           fn = S_theta, 
                           method = "Nelder-Mead", 
                           data = boot_data, 
                           varnames = varnames, 
                           control = list("fnscale" = -1))
  
  candidates[[2]] <- optim(par = theta_logit, 
                           fn = S_theta, 
                           method = "Nelder-Mead", 
                           data = boot_data, 
                           varnames = varnames, 
                           control = list("fnscale" = -1))
  
  candidates[[3]] <- optim(par = theta_lpm, 
                           fn = S_theta, 
                           method = "Nelder-Mead", 
                           data = boot_data, 
                           varnames = varnames, 
                           control = list("fnscale" = -1))
  
  # Choose the best performing value
  i <- which.max(sapply(candidates, getElement, "value"))
  theta_hat <- candidates[[i]]$par
  theta_hat <- theta_hat / sqrt(sum(theta_hat ^ 2)) # scale to be length 1
  
  return(theta_hat)
  
}

##################################### SETUP ####################################

# Load in parameters
parameters <- readRDS(paste0(path, "/Monte_Carlo/Parameters1.rds"))

for (i in 1:length(parameters)){
  assign(names(parameters)[i], parameters[[i]])
}

# Settings 
k <- 3 # number of parameters
N <- c(100, 500, 1000, 2000) # number of observations
pZ <- 1 / 2 # probabilty Z = 1

# Monte Carlo settings
runs <- 1000

############################### MONTE CARLO LOOP ###############################

results <- vector(mode = "list", length = length(N))

for (j in 1:length(N)){
  
  tsls <- rep(0, runs)
  corrected <- rep(0, runs)
  oracle <- rep(0, runs)
  probit <- rep(0, runs)
  logit <- rep(0, runs)
  lpm <- rep(0, runs)
  diagnostics <- list(s_cor = rep(0, runs),
                      theta_hat = matrix(0, nrow = 1 + (2 * k), ncol = runs),
                      theta = matrix(0, nrow = 1 + (2 * k), ncol = runs))
  
  n <- N[j]
  
  for (r in 1:runs){
    
    #### GENERATE DATA ####
    
    # Generate parameters according to the constraints
    ad <- gen_ad(k, sd_a, sd_d, rho_ad)
    a <- ad$a
    d <- ad$d
    
    g <- a + d
    g_0 <- a_0 + d_0
    
    B <- runif(k, min = -1, max = 1)
    B <- B * sd_B / sqrt(sum(B ^ 2))
    
    # Generate sample
    Z <- rbinom(n, 1, pZ)
    X <- matrix(rnorm(n * k), n, k)
    V <- rnorm(n)
    
    index <- a_0 + X %*% as.matrix(a, k, 1) + Z * (d_0 + X %*% as.matrix(d, k, 1))
    P <- pnorm(index)
    D <- as.numeric(index + V > 0) 
    
    U <- rmvnorm(n, mean = rep(0, k), sigma = sigma)
    colnames(U) <- c("V", "U1", "U0")
    U <- data.frame(U)
    
    Y_0 <- B_0 + X %*% as.matrix(B, k, 1) + U$U0
    Y_1 <- B_0 + zeta_0 + X %*% as.matrix(B, k, 1) + U$U1
    Y <- Y_0 + D * (Y_1 - Y_0)
    
    # Build a data frame
    data <- data.frame(D = D, Z = Z, Y = Y)
    data <- cbind(data, X)
    names(data) <- c("D", "Z", "Y", paste0("X_", 1:k))
    
    # Generate interactions
    interactions <- apply(X, 2, function(X, Z){return(Z * X)}, Z = Z)
    colnames(interactions) <- paste("Z_X", 1:k, sep = "_")
    data <- cbind(data, interactions)
    
    #### MRC ESTIMATOR ####
    
    # Here theta = (a, d_0, d) (scaled to have length 1)
    varnames <- c(paste("X", 1:k, sep = "_"), "Z", paste("Z_X", 1:k, sep = "_"))
    
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
    
    s_hat <- sign(theta_hat["Z"] + X %*% 
                    matrix(theta_hat[colnames(interactions)], k, 1))
    s_star <- sign(d_0 + as.matrix(X) %*% as.matrix(d, k, 1)) # oracle
    s_probit <- sign(probit_est["Z"] + X %*% 
                       matrix(probit_est[colnames(interactions)], k ,1))
    s_logit <- sign(logit_est["Z"] + X %*% 
                      matrix(logit_est[colnames(interactions)], k ,1))
    s_lpm <- sign(lpm_est["Z"] + X %*% 
                    matrix(lpm_est[colnames(interactions)], k ,1))
    
    # Estimates
    Z_hat <- lm(as.formula(paste0("Z ~ ", paste(paste0("X_", 1:k), collapse = " + "))),
                data = data)$fitted.values
    Z_tilde <- Z - Z_hat
    # corrected TSLS
    oracle[r] <- sum(Y * s_star * Z_tilde) / sum(D * s_star * Z_tilde)
    corrected[r] <- sum(Y * s_hat * Z_tilde) / sum(D * s_hat * Z_tilde)
    probit[r] <- sum(Y * s_probit * Z_tilde) / sum(D * s_probit * Z_tilde)
    logit[r] <- sum(Y * s_logit * Z_tilde) / sum(D * s_logit * Z_tilde)
    lpm[r] <- sum(Y * s_lpm * Z_tilde) / sum(D * s_lpm * Z_tilde)
    tsls[r] <- cov(Y, Z_tilde) / cov(D, Z_tilde)
    
    # Diagnostics
    diagnostics$s_cor[r] <- cor(s_hat, s_star)
    diagnostics$theta[ , r] <- c(a, d_0, d) / sqrt(sum(c(a, d_0, d) ^ 2))
    diagnostics$theta_hat[, r] <- theta_hat 
    
  }
  
  run_results <- list(oracle = oracle, 
                      corrected = corrected,
                      tsls = tsls, 
                      probit = probit,
                      logit = logit,
                      lpm = lpm,
                      diagnostics = diagnostics)
  
  results[[j]] <- run_results
  
}

names(results) <- N

saveRDS(results, file = paste0(path, "/Monte_Carlo/results.rds"))

############################# VISUALISING RESULTS ############################## 

sims <- length(N)

for (i in 1:sims){
  
  x_lims <- rep(0, 2)
  x_lims[1] <- quantile(results[[i]]$tsls[!is.na(results[[1]]$corrected)], 0.05)
  x_lims[2] <- quantile(results[[i]]$tsls[!is.na(results[[1]]$corrected)], 0.95)
  plot_points <- list()
  
  plot_data <- data.frame(tsls = results[[i]]$tsls, 
                          oracle = results[[i]]$oracle, 
                          corrected = results[[i]]$corrected,
                          probit = results[[i]]$probit,
                          logit = results[[i]]$logit,
                          lpm = results[[i]]$lpm)
  
  # plot_data %>% ggplot() + 
  #   geom_density(aes(x = tsls, col = "TSLS")) + 
  #   geom_density(aes(x = corrected, col = "corrected")) + 
  #   geom_density(aes(x = oracle, col = "oracle"))
  
  plot_data %>% ggplot() + 
    geom_density(aes(x = corrected, col = "corrected")) + 
    geom_density(aes(x = oracle, col = "oracle")) +
    geom_density(aes(x = probit, col = "probit")) + 
    geom_density(aes(x = logit, col = "logit")) +
    geom_density(aes(x = lpm, col = "lpm")) + 
    theme_bw() +
    theme(legend.title = element_blank()) +
    labs(x = "Estimate", y = "")
  
  plot_data %>% ggplot() + 
    geom_density(aes(x = corrected, col = "corrected")) + 
    geom_density(aes(x = oracle, col = "oracle")) +
    geom_density(aes(x = probit, col = "probit")) + 
    theme_bw() +
    theme(legend.title = element_blank()) +
    labs(x = "Estimate", y = "")
  
}


# Summary statistics
table <- matrix(0, 6 * sims, 7)

get_stats <- function(results){
  
  stats <- rep(0, 7)
  results <- results[!is.infinite(results)]
  stats[1] <- mean(results, na.rm = TRUE)
  stats[2] <- sd(results, na.rm = TRUE)
  stats[3:7] <- quantile(results, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  return(stats)
  
}

rownames(table) <- rep(N, 6)
table[1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "tsls"))},
           results = results))
table[sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "corrected"))},
           results = results))
table[2 * sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "oracle"))},
           results = results))
table[3 * sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "probit"))},
           results = results))
table[4 * sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "logit"))},
           results = results))
table[5 * sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "lpm"))},
           results = results))


xtable(table[1:sims, ], digits = 3) # tsls
xtable(table[sims + 1:sims, ], digits = 3) # corrected
xtable(table[2 * sims + 1:sims, ], digits = 3) # oracle
xtable(table[3 * sims + 1:sims, ], digits = 3) # probit
xtable(table[4 * sims + 1:sims, ], digits = 3) # logit
xtable(table[5 * sims + 1:sims, ], digits = 3) # lpm

# 
# ########################### BOOTSTRAP STANDARD ERRORS ##########################
# 
# # Obtain standard errors by non-parametric bootstrap
# temp <- boot(data = data, statistic = est_theta, varnames = varnames, R = 10)
