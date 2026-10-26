#################################### SUMMARY ################################### 

# This scipt performs Monte Carlo for the feasible reordered IV estimator

# The first stage is assumed to be
# D = 1(z * x'd - V > 0)
# where the marginal distribution of V is N(0, 1), so the propensity score is
# p(x,z) = Phi(z * x'd)
# x'd controls the sign of the relationship between Z and D. 

# The second stage is 
# Y(0) = U(0)
# Y(1) = B_0 + psi(x'B) + U(1)
# Y = Y(0) + D(Y(1) - Y(0))
# where psi() is a function which varies between specifications and
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

target <- 0.05 # share of compliers (half the share of switchers)
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

############################### SWITCHER DENSITY ###############################

n_points <- 101
C <- seq(0, 1/2, length.out = n_points)
dens <- rep(0, n_points)
dist <- rep(0, n_points)

# sd_d <- 1

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

###################### EFFECT OF ESTIMATION ON DENOMINATOR #####################

kernel_one <- function(z, sd_d){
  
  return((pnorm(z) - pnorm(0)) * dnorm(z, mean = 0, sd = sd_d))
  
}

kernel_two <- function(z, sd_d, rho_dd){

  return((1 - pnorm(-rho_dd * z / (sd_d * sqrt(1 - rho_dd ^ 2)))) * 
    (pnorm(z) - pnorm(0)) * dnorm(z, mean = 0, sd = sd_d))
}

K <- 3

pZ <- 1 / 2
Var_Z <- pZ * (1 - pZ)

n <- 100000

runs <- 100
data_diff <- rep(0, runs)
analytical_diff <- rep(0, runs)
data_term_1 <- rep(0, runs)
data_term_2 <- rep(0, runs)
analytical_term_1 <- rep(0, runs)
analytical_term_2 <- rep(0, runs)
rho_record <- rep(0, runs)
 
for (i in 1:runs){
  
  d <- runif(K, min = - 1, max  = 1)
  d <- d / sqrt(sum(d ^ 2))
  d_hat <- d + rnorm(K, mean = 0, sd = 1 / 4)
  
  sd_d <- sqrt(sum(d ^ 2))
  sd_d_hat <- sqrt(sum(d_hat ^ 2))
  rho_dd <- sum(d * d_hat) / (sd_d * sd_d_hat)
  rho_record[i] <- rho_dd

  X <- matrix(rnorm(n * K), nrow = n, ncol = K)
  Xd <- X %*% as.matrix(d, nrow = K, ncol = 1)
  Xd_hat <- X %*% as.matrix(d_hat, nrow = K, ncol = 1)
  
  s <- 2 * as.numeric(Xd >= 0) - 1
  s_hat <- 2 * as.numeric(Xd_hat >= 0) - 1
  
  Z <- rbinom(n, 1, pZ)
  Z_tilde <- Z - mean(Z)
  
  V <- rnorm(n)
  D <- as.numeric(Z * Xd >= V)
  
  data_term_1[i] <- mean(s_hat * D * Z_tilde)
  data_term_2[i] <- mean(s * D * Z_tilde)
  
  data_diff[i] <- data_term_1[i] - data_term_2[i]
  
  analytical_term_1[i] <- 2 * Var_Z * 
    integrate(kernel_two, 
              lower = -Inf, 
              upper = Inf, 
              sd_d = sd_d, 
              rho_dd = rho_dd)$value
  
  analytical_term_2[i] <- 2 * Var_Z * 
    integrate(kernel_one, lower = 0, upper = Inf, sd_d = sd_d)$value
  
  analytical_diff[i] <- analytical_term_1[i] - analytical_term_2[i]

}

plot(data_diff, analytical_diff)
plot(data_term_1, analytical_term_1)
plot(data_term_2, analytical_term_2)

plot(rho_record, data_diff)
plot(rho_record, analytical_diff)

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

#################### EXPECTED TREATMENT EFFECTS (PSI = X'B) ####################

# Assume psi(X'B) = X'B. Get the expected treatment effect. Note that, in this
# case, treatment effects are largest (in magnitude) for strata with the most
# switchers. 
 
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

# rnadomly draws points and evaluate
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
plot(comparison$Xd, comparison$analytical)
plot(comparison$Xd, comparison$data)

######################### VARIANCE OF TREATMENT EFFECTS ########################

kernel_one <- function(z, sd_d){
  
  return((dnorm(0) - dnorm(z)) * z * dnorm(z, mean = 0, sd = sd_d) / 
           (pnorm(z) - pnorm(0)))
  
}

kernel_two <- function(z, sd_d){
  
  return((dnorm(0) - dnorm(z)) ^ 2 * dnorm(z, mean = 0, sd = sd_d) / 
           ((pnorm(z) - pnorm(0)) ^ 2))
  
}

sigma <- diag(1, 3)
sigma[2, 1] <- -0.25
sigma[3, 1] <- -0.5
sigma[3, 2] <- 0.5
sigma[1, 2] <- -0.25
sigma[1, 3] <- -0.5
sigma[2, 3] <- 0.5

rho_1V <- sigma[3, 1]
rho_0V <- sigma[3, 2]

K <- 3

runs <- 100
analytical_var <- rep(0, 100)
data_var <- rep(0, 100)

B_0 <- 1
n <- 100000

for (i in 1:runs){
  
  d <- runif(K, -1, 1)
  # d <- d / sqrt(sum(d ^ 2))
  B <- runif(K, -1, 1)
  # B <- B / sqrt(sum(B ^ 2))
  sd_d <- sqrt(sum(d ^ 2))
  sd_B <- sqrt(sum(B ^ 2))
  rho_Bd <- sum(B * d) / (sd_d * sd_B)
  
  analytical_var[i] <- sd_B ^ 2 + 
    2 * (rho_1V - rho_0V) * (sd_B * rho_Bd / sd_d) * 
    integrate(kernel_one, lower = -Inf, upper = Inf, sd_d = sd_d)$value +
    2 * (rho_1V - rho_0V) ^ 2 * 
    integrate(kernel_two, lower = 0, upper = Inf, sd_d = sd_d)$value
  
  # Validate by simulation
  X <- matrix(rnorm(n * K), nrow = n, ncol = K)
  Xd <- X %*% matrix(d, nrow = K, ncol = 1)
  XB <- X %*% matrix(B, nrow = K, ncol = 1)
  
  tau_X <- B_0 + XB + (rho_1V - rho_0V) * (dnorm(0) - dnorm(Xd)) / 
    (pnorm(Xd) - pnorm(0))
  data_var[i] <- var(tau_X)
  
}

plot(data_var, analytical_var)

####################### EFFECT OF ESTIMATION ON NUMERATOR ######################

kernel_one <- function(z, d0, d0_hat, rho_dd, sd_d){
  
  return(
    (1 - (d0 + z >= 0) - 
       pnorm(-(d0_hat + rho_dd * z) / (sd_d * sqrt(1 - rho_dd ^ 2)))) * 
      (pnorm(d0 + z) - pnorm(0)) * dnorm(z, mean = 0, sd = sd_d))

}

kernel_two <- function(z, d0, d0_hat, rho_dd, sd_d){
  
  return(dnorm(-(d0_hat + rho_dd * z) / (sd_d * sqrt(1 - rho_dd ^ 2))) * 
      (pnorm(d0 + z) - pnorm(0)) * dnorm(z, mean = 0, sd = sd_d))
  
}

kernel_three <- function(z, d0, d0_hat, rho_dd, sd_d){
  
  return(
    (1 - (d0 + z >= 0) - 
       pnorm(-(d0_hat + rho_dd * z) / (sd_d * sqrt(1 - rho_dd ^ 2)))) * 
      (pnorm(d0 + z) - pnorm(0)) * z * dnorm(z, mean = 0, sd = sd_d))
  
}

kernel_four <- function(z, d0, d0_hat, rho_dd, sd_d){
  
  return(
    (1 - (d0 + z >= 0) - 
       pnorm(-(d0_hat + rho_dd * z) / (sd_d * sqrt(1 - rho_dd ^ 2)))) * 
      (dnorm(0) - dnorm(d0 + z)) * dnorm(z, mean = 0, sd = sd_d))
  
}

K <- 3
n <- 100000

sigma <- diag(1, 3)
sigma[2, 1] <- -0.25
sigma[3, 1] <- -0.5
sigma[3, 2] <- 0.5
sigma[1, 2] <- -0.25
sigma[1, 3] <- -0.5
sigma[2, 3] <- 0.5

rho_1V <- sigma[3, 1]
rho_0V <- sigma[3, 2]

pZ <- 1 / 2
Var_Z <- pZ * (1 - pZ)

# constants
d0 <- 0.2
B0 <- 0.75

# Specify lengths and inner product for d and B
sd_d <- 0.75
sd_B <- 1.25
rho_Bd <- 0.5

runs <- 100
data_term_1 <- rep(0, runs)
data_term_2 <- rep(0, runs)
data_term_3 <- rep(0, runs)
analytical_term_1 <- rep(0, runs)
analytical_term_2 <- rep(0, runs)
analytical_term_3 <- rep(0, runs)

for (i in 1:runs){
  
  # Generate B and d subject to constraints
  # Construct orthogonal vectors 
  e_1 <- runif(K, min = - 1, max = 1)
  u <- runif(K, min = - 1, max = 1)
  e_2 <- u - sum(e_1 * u) * e_1 / (sum(e_1 ^ 2))
  e_1 <- e_1 / sqrt(sum(e_1 ^ 2))
  e_2 <- e_2 / sqrt(sum(e_2 ^ 2))
  
  # Use e_1 and e_2 to build d and B
  d <- e_1 * sd_d
  B <- sd_B * (rho_Bd * e_1 + sqrt(1 - rho_Bd ^ 2) * e_2)
  
  # randomly perturbed "estimates"
  d_hat <- d + rnorm(K, mean = 0, sd = 1 / 4)
  d0_hat <- d0 + rnorm(1, mean = 0, sd = 1 / 4)
  sd_d_hat <- sqrt(sum(d_hat ^ 2))
  rho_dd <- sum(d * d_hat) / (sd_d * sd_d_hat)
  rho_Bd_hat <- sum(B * d_hat) / (sd_B * sd_d_hat)
  
  X <- matrix(rnorm(n * K), nrow = n, ncol = K)
  Xd <- X %*% as.matrix(d, nrow = K, ncol = 1)
  Xd_hat <- X %*% as.matrix(d_hat, nrow = K, ncol = 1)
  XB <- X %*% as.matrix(B, nrow = K, ncol = 1)
  
  s <- 2 * as.numeric(d0 + Xd >= 0) - 1
  s_hat <- 2 * as.numeric(d0_hat + Xd_hat >= 0) - 1
  
  Z <- rbinom(n, 1, pZ)
  Z_tilde <- Z - mean(Z)
  
  U <- rmvnorm(n, mean = rep(0, K), sigma = sigma)
  colnames(U) <- c("U1", "U0", "V")
  U <- data.frame(U)
  
  D <- as.numeric(Z * (d0 + Xd) >= U$V)
  
  Y <- (1 - D) * U$U0 + D * (B0 + XB + U$U1)
  
  # Term 1 (B_0)
  data_term_1[i] <- B0 * mean((s_hat - s)* D * Z_tilde)
  analytical_term_1[i] <- 2 * Var_Z * B0 * 
    integrate(kernel_one, lower = -Inf, upper = Inf,
              d0 = d0, d0_hat = d0_hat, rho_dd = rho_dd, sd_d = sd_d)$value
   
  
  # Term 2 (XB) 
  data_term_2[i] <- mean((s_hat - s) * D * XB * Z_tilde)
  analytical_term_2[i] <- 2 * Var_Z * (
    (sd_B * (rho_Bd_hat - rho_dd * rho_Bd) / sqrt(1 - rho_dd ^ 2)) *
      integrate(kernel_two, lower = -Inf, upper = Inf, 
                d0 = d0, d0_hat = d0_hat, rho_dd = rho_dd, sd_d = sd_d)$value +
      (sd_B * rho_Bd / sd_d) *
      integrate(kernel_three, lower = -Inf, upper = Inf, 
                d0 = d0, d0_hat = d0_hat, rho_dd = rho_dd, sd_d = sd_d)$value
  )
  
  # Term 3 (unobservables)
  data_term_3[i] <- mean((s_hat - s) * ((1 - D) * U$U0 + D * U$U1) * Z_tilde)
  
  analytical_term_3[i] <- 2 * (rho_1V - rho_0V) * Var_Z * 
    integrate(kernel_four, lower = -Inf, upper = Inf,
              d0 = d0, d0_hat = d0_hat, rho_dd = rho_dd, sd_d = sd_d)$value
 
}

data_total <- data_term_1 + data_term_2 + data_term_3
analytical_total <- analytical_term_1 + analytical_term_2 + analytical_term_3

plot(data_term_1, analytical_term_1)
plot(data_term_2, analytical_term_2)
plot(data_term_3, analytical_term_3)
plot(data_total, analytical_total)

################################# TSLS ESTIMAND ################################

kernel_one <- function(z, d0, sd_d){
  
  return((2 * (d0 + z >= 0) - 1) * (pnorm(d0 + z) - pnorm(0)) * 
           dnorm(z, mean = 0, sd = sd_d))
  
}

kernel_two <- function(z, d0, sd_d){
  
  return((2 * (d0 + z >= 0) - 1) * (pnorm(d0 + z) - pnorm(0)) * 
           z * dnorm(z, mean = 0, sd = sd_d))
  
}

kernel_three <- function(z, d0, sd_d){
  
  return((2 * (d0 + z >= 0) - 1) * (dnorm(0) - dnorm(d0 + z)) * 
           dnorm(z, mean = 0, sd = sd_d))
  
}

K <- 3
n <- 100000

sigma <- diag(1, 3)
sigma[2, 1] <- -0.25
sigma[3, 1] <- -0.5
sigma[3, 2] <- 0.5
sigma[1, 2] <- -0.25
sigma[1, 3] <- -0.5
sigma[2, 3] <- 0.5

rho_1V <- sigma[3, 1]
rho_0V <- sigma[3, 2]

pZ <- 1 / 2
Var_Z <- pZ * (1 - pZ)

# constants
d0 <- 0.2
B0 <- 0.75

# Specify lengths and inner product for d and B
sd_d <- 0.75
sd_B <- 1.25
rho_Bd <- 0.5

runs <- 100
data_denom <- rep(0, runs)
data_num_1 <- rep(0, runs)
data_num_2 <- rep(0, runs)
data_num_3 <- rep(0, runs)

analytical_denom  <- Var_Z * 
  integrate(kernel_one, lower = -Inf, upper = Inf, d0 = d0, sd_d = sd_d)$value
analytical_num_1 <- Var_Z * B0 * 
  integrate(kernel_one, lower = -Inf, upper = Inf, 
            d0 = d0, sd_d = sd_d)$value
analytical_num_2 <- Var_Z * (sd_B * rho_Bd / sd_d) *
  integrate(kernel_two, lower = -Inf, upper = Inf, 
            d0 = d0, sd_d = sd_d)$value 
analytical_num_3 <- Var_Z * (rho_1V - rho_0V) * 
  integrate(kernel_three, lower = -Inf, upper = Inf, 
            d0 = d0, sd_d = sd_d)$value    
analytical_num <- analytical_num_1 + analytical_num_2 + analytical_num_3

for (i in 1:runs){
  
  # Generate B and d subject to constraints
  # Construct orthogonal vectors 
  e_1 <- runif(K, min = - 1, max = 1)
  u <- runif(K, min = - 1, max = 1)
  e_2 <- u - sum(e_1 * u) * e_1 / (sum(e_1 ^ 2))
  e_1 <- e_1 / sqrt(sum(e_1 ^ 2))
  e_2 <- e_2 / sqrt(sum(e_2 ^ 2))
  
  # Use e_1 and e_2 to build d and B
  d <- e_1 * sd_d
  B <- sd_B * (rho_Bd * e_1 + sqrt(1 - rho_Bd ^ 2) * e_2)
  
  X <- matrix(rnorm(n * K), nrow = n, ncol = K)
  Xd <- X %*% as.matrix(d, nrow = K, ncol = 1)
  XB <- X %*% as.matrix(B, nrow = K, ncol = 1)
  
  s <- 2 * as.numeric(d0 + Xd >= 0) - 1

  Z <- rbinom(n, 1, pZ)
  Z_tilde <- Z - mean(Z)
  
  U <- rmvnorm(n, mean = rep(0, K), sigma = sigma)
  colnames(U) <- c("U1", "U0", "V")
  U <- data.frame(U)
  
  D <- as.numeric(Z * (d0 + Xd) >= U$V)
  
  Y <- (1 - D) * U$U0 + D * (B0 + XB + U$U1)
  
  # Denominator
  data_denom[i] <- mean(s * D * Z_tilde)
 
  # Numerator
  data_num_1[i] <- B0 * mean(s * D * Z_tilde)
  data_num_2[i] <- mean(s * D * XB * Z_tilde)
  data_num_3[i] <- mean(s * ((1 - D) * U$U0 + D * U$U1) * Z_tilde)
  
}

data_num <- data_num_1 + data_num_2 + data_num_3

plot(data_denom, rep(analytical_denom, runs))
plot(data_num, rep(analytical_num, runs))
mean(data_num)
analytical_num

plot(data_num_1, rep(analytical_num_1, runs))
mean(data_num_1)
analytical_num_1
plot(data_num_2, rep(analytical_num_2, runs))
mean(data_num_2)
analytical_num_2
plot(data_num_3, rep(analytical_num_3, runs))
mean(data_num_3)
analytical_num_3

data_tau <- data_num / data_denom
analytical_tau <- analytical_num / analytical_denom
plot(data_tau, rep(analytical_tau, runs))
