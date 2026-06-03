########################### FRONT MATTER ###########################

rm(list=ls())
gc()

library(dplyr)
library(ggplot2)
library(tidyr)
library(mvtnorm)
library(clusterGeneration)

####################### CHECKING CALCULATIONS ######################

# Specify process for the unobservables.

# Generate sigma
A <- matrix(runif(9, -1, 1), 3, 3)
# Scale so that all diagonal elements will be 1
for (i in 1:3){
  A[i, ] <- A[i, ] / sqrt(sum(A[i, ] ^ 2))
}
sigma <- A %*% t(A)

# Generate unobservables
n <- 100000
U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
U <- data.frame(U)
names(U) <- c("V", "U_0", "U_1")

# pick random propensity scores
p_0 <- runif(1, 0, 1)
p_1 <- p_0 + runif(1, 0, 1 - p_0)
q_0 <- qnorm(p_0)
q_1 <- qnorm(p_1)

# check expectations
temp <- dplyr::filter(U, V >= q_0 & V <= q_1)
mean(temp$U_0)
sigma[2, 1] * (dnorm(q_0) - dnorm(q_1)) / (p_1 - p_0)
mean(temp$U_1)
sigma[3, 1] * (dnorm(q_0) - dnorm(q_1)) / (p_1 - p_0)

########################### CHECKING LATE ##########################

# Simple case: no covariates, binary instrument, binary treatment
# Check that analytically calculated and estimated LATEs match

rm(list = ls())

n <- 1000000
Z <- rbinom(n, 1, 1 / 2)
tau <- 0.05 # average treatment effect is 2 * tau
a_1 <- qnorm(0.6) # pnorm(a_1) is propensity score when Z = 1
a_0 <- qnorm(0.45) # pnorm(a_0) is propensity score when Z = 0

# covariance matrix for the unobservables
A <- matrix(runif(9, -1, 1), 3, 3)
# scale so that all diagonal elements will be 1
for (i in 1:3){
  A[i, ] <- A[i, ] / sqrt(sum(A[i, ] ^ 2))
}
sigma <- A %*% t(A)

# generate unobservables
U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
U <- data.frame(U)
names(U) <- c("V", "U_0", "U_1")

# generate observables
D <- as.numeric(a_1 * Z + a_0 * (1 - Z) >= U$V)
Y_1 <- tau + U$U_1
Y_0 <- - tau + U$U_0
Y <- D * Y_1 + (1 - D) * Y_0

# estimated LATE
data <- data.frame(Y = Y, D = D, Z = Z)
D_hat <- lm(D ~ Z, data = data)$fitted.values
data$D_hat <- D_hat
lm(Y ~ D_hat, data = data)

# alternative calculation
Z_tilde <- Z - mean(Z)
cov(Y, Z_tilde) / cov(D, Z_tilde)

# calculated LATE
p_0 <- pnorm(a_0)
p_1 <- pnorm(a_1)
q_0 <- a_0
q_1 <- a_1

# check expectations
2 * tau + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(q_0) - dnorm(q_1)) / (p_1 - p_0)

########################## SIMPLE EXAMPLE ##########################

# Single binary X -> only two LATES to worry about
# Sign of the first-stage is different for each X value

rm(list = ls())

n <- 1000000
Z <- rbinom(n, 1, 1 / 2)
X <- rbinom(n, 1, 1 / 2)

X_ord <- order(X)
Z <- Z[X_ord]
X <- X[X_ord]

n_1 <- sum(X)
n_0 <- n - n_1

tau <- 0.05 # average treatment effect is 2 * tau
a_1 <- qnorm(0.6) # pnorm(a_1) is the higher propensity score 
a_0 <- qnorm(0.45) # pnorm(a_0) is the lower propensity score

# covariance matrix for the unobservables, let this be common
A <- matrix(runif(9, -1, 1), 3, 3)
# scale so that all diagonal elements will be 1
for (i in 1:3){
  A[i, ] <- A[i, ] / sqrt(sum(A[i, ] ^ 2))
}
sigma <- A %*% t(A)

# generate unobservables
U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
U <- data.frame(U)
names(U) <- c("V", "U_0", "U_1")

# generate D
D <- as.numeric((Z * (a_1 * (1 - X) + a_0 * X) +
  (1 - Z) * (a_0 * (1 - X) + a_1 * X)) >= U$V)

# check D is correct
mean(D[Z == 1 & X == 1]) - pnorm(a_0)
mean(D[Z == 0 & X == 1]) - pnorm(a_1)
mean(D[Z == 1 & X == 0]) - pnorm(a_1)
mean(D[Z == 0 & X == 0]) - pnorm(a_0)

# generate Y
Y_1 <- tau + U$U_1
Y_0 <- - tau + U$U_0
Y <- D * Y_1 + (1 - D) * Y_0

#### Estimate LATE in subsamples #### 

data <- data.frame(Y = Y, D = D, Z = Z, X = X)
data$Z_tilde <- Z - mean(Z)

p_0 <- pnorm(a_0)
p_1 <- pnorm(a_1)
q_0 <- a_0
q_1 <- a_1

# Analytical Calculation
2 * tau + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(q_0) - dnorm(q_1)) / (p_1 - p_0)

# X = 0 subsample
data_X0 <- dplyr::filter(data, X == 0)
data_X0$Z_tilde <- data_X0$Z - mean(data_X0$Z)
D_hat <- lm(D ~ Z, data = data_X0)$fitted.values
data_X0$D_hat <- D_hat
lm(Y ~ D_hat, data = data_X0)

# alternative calculation
cov(data_X0$Y, data_X0$Z_tilde) / cov(data_X0$D, data_X0$Z_tilde)

# X = 1 subsample
data_X1 <- dplyr::filter(data, X == 1)
data_X1$Z_tilde <- data_X1$Z - mean(data_X1$Z)
D_hat <- lm(D ~ Z, data = data_X1)$fitted.values
data_X1$D_hat <- D_hat
lm(Y ~ D_hat, data = data_X1)

# alternative calculation
cov(data_X1$Y, data_X1$Z_tilde) / cov(data_X1$D, data_X1$Z_tilde)

#### LATE in the full sample ####

# analytical calculation 
# as Z and X are independent and the propensity score change is
# constant, the LATE for the full sample equals the LATE in each subsample
2 * tau + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(q_0) - dnorm(q_1)) / (p_1 - p_0)

# estimated
D_hat <- lm(D ~ Z, data = data)$fitted.values
data$D_hat <- D_hat
lm(Y ~ D_hat, data = data)

# alternative method
cov(data$Y, data$Z_tilde) / cov(data$D, data$Z_tilde)

# Calculation as a weighted average (separate Z_tilde)
denom_0 <- cov(data_X0$D, data_X0$Z_tilde)
denom_1 <- cov(data_X1$D, data_X1$Z_tilde)
denom_all <- cov(data$D, data$Z_tilde)
theta_0 <- denom_0 / denom_all
theta_1 <- denom_1 / denom_all
LATE_0 <- cov(data_X0$Y, data_X0$Z_tilde) / 
  cov(data_X0$D, data_X0$Z_tilde)
LATE_1 <- cov(data_X1$Y, data_X1$Z_tilde) / 
  cov(data_X1$D, data_X1$Z_tilde)
(theta_0 * n_0 * LATE_0 + theta_1 * n_1 * LATE_1) / n
cov(data$Y, data$Z_tilde) / cov(data$D, data$Z_tilde)

# common Z_tilde
temp_X0 <- dplyr::filter(data, X == 0)
temp_X1 <- dplyr::filter(data, X == 1)
denom_0 <- cov(temp_X0$D, temp_X0$Z_tilde)
denom_1 <- cov(temp_X1$D, temp_X1$Z_tilde)
denom_all <- cov(data$D, data$Z_tilde)
theta_0 <- denom_0 / denom_all
theta_1 <- denom_1 / denom_all
LATE_0 <- cov(temp_X0$Y, temp_X0$Z_tilde) / 
  cov(temp_X0$D, temp_X0$Z_tilde)
LATE_1 <- cov(temp_X1$Y, temp_X1$Z_tilde) / 
  cov(temp_X1$D, temp_X1$Z_tilde)
(theta_0 * n_0 * LATE_0 + theta_1 * n_1 * LATE_1) / n
cov(data$Y, data$Z_tilde) / cov(data$D, data$Z_tilde)

# Relabeled data
data$Z_mod <- data$Z_tilde * (2 * data$X - 1)
cov(data$Y, data$Z_mod) / cov(data$D, data$Z_mod)

########################## BIGGER EXAMPLE ##########################

# Two covariates: X_1 and X_2. Sign varies over their range
rm(list = ls())

n <- 1000000
Z <- rbinom(n, 1, 1 / 2)
X <- rbinom(n, 1, 1 / 2)

X_ord <- order(X)
Z <- Z[X_ord]
X <- X[X_ord]

n_1 <- sum(X)
n_0 <- n - n_1

tau <- 0.05 # average treatment effect is 2 * tau
a_1 <- qnorm(0.6) # pnorm(a_1) is the higher propensity score 
a_0 <- qnorm(0.45) # pnorm(a_0) is the lower propensity score

# covariance matrix for the unobservables, let this be common
A <- matrix(runif(9, -1, 1), 3, 3)
# scale so that all diagonal elements will be 1
for (i in 1:3){
  A[i, ] <- A[i, ] / sqrt(sum(A[i, ] ^ 2))
}
sigma <- A %*% t(A)

# generate unobservables
U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
U <- data.frame(U)
names(U) <- c("V", "U_0", "U_1")

# generate D
D <- as.numeric((Z * (a_1 * (1 - X) + a_0 * X) +
                   (1 - Z) * (a_0 * (1 - X) + a_1 * X)) >= U$V)

# check D is correct
mean(D[Z == 1 & X == 1]) - pnorm(a_0)
mean(D[Z == 0 & X == 1]) - pnorm(a_1)
mean(D[Z == 1 & X == 0]) - pnorm(a_1)
mean(D[Z == 0 & X == 0]) - pnorm(a_0)

# generate Y
Y_1 <- tau + U$U_1
Y_0 <- - tau + U$U_0
Y <- D * Y_1 + (1 - D) * Y_0

#### Estimate LATE in subsamples #### 

data <- data.frame(Y = Y, D = D, Z = Z, X = X)
data$Z_tilde <- Z - mean(Z)

p_0 <- pnorm(a_0)
p_1 <- pnorm(a_1)
q_0 <- a_0
q_1 <- a_1

# Analytical Calculation
2 * tau + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(q_0) - dnorm(q_1)) / (p_1 - p_0)

# X = 0 subsample
data_X0 <- dplyr::filter(data, X == 0)
data_X0$Z_tilde <- data_X0$Z - mean(data_X0$Z)
D_hat <- lm(D ~ Z, data = data_X0)$fitted.values
data_X0$D_hat <- D_hat
lm(Y ~ D_hat, data = data_X0)

# alternative calculation
cov(data_X0$Y, data_X0$Z_tilde) / cov(data_X0$D, data_X0$Z_tilde)

# X = 1 subsample
data_X1 <- dplyr::filter(data, X == 1)
data_X1$Z_tilde <- data_X1$Z - mean(data_X1$Z)
D_hat <- lm(D ~ Z, data = data_X1)$fitted.values
data_X1$D_hat <- D_hat
lm(Y ~ D_hat, data = data_X1)

# alternative calculation
cov(data_X1$Y, data_X1$Z_tilde) / cov(data_X1$D, data_X1$Z_tilde)

#### LATE in the full sample ####

# analytical calculation 
# as Z and X are independent and the propensity score change is
# constant, the LATE for the full sample equals the LATE in each subsample
2 * tau + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(q_0) - dnorm(q_1)) / (p_1 - p_0)

# estimated
D_hat <- lm(D ~ Z, data = data)$fitted.values
data$D_hat <- D_hat
lm(Y ~ D_hat, data = data)

# alternative method
cov(data$Y, data$Z_tilde) / cov(data$D, data$Z_tilde)

# Calculation as a weighted average (separate Z_tilde)
denom_0 <- cov(data_X0$D, data_X0$Z_tilde)
denom_1 <- cov(data_X1$D, data_X1$Z_tilde)
denom_all <- cov(data$D, data$Z_tilde)
theta_0 <- denom_0 / denom_all
theta_1 <- denom_1 / denom_all
LATE_0 <- cov(data_X0$Y, data_X0$Z_tilde) / 
  cov(data_X0$D, data_X0$Z_tilde)
LATE_1 <- cov(data_X1$Y, data_X1$Z_tilde) / 
  cov(data_X1$D, data_X1$Z_tilde)
(theta_0 * n_0 * LATE_0 + theta_1 * n_1 * LATE_1) / n
cov(data$Y, data$Z_tilde) / cov(data$D, data$Z_tilde)

# common Z_tilde
temp_X0 <- dplyr::filter(data, X == 0)
temp_X1 <- dplyr::filter(data, X == 1)
denom_0 <- cov(temp_X0$D, temp_X0$Z_tilde)
denom_1 <- cov(temp_X1$D, temp_X1$Z_tilde)
denom_all <- cov(data$D, data$Z_tilde)
theta_0 <- denom_0 / denom_all
theta_1 <- denom_1 / denom_all
LATE_0 <- cov(temp_X0$Y, temp_X0$Z_tilde) / 
  cov(temp_X0$D, temp_X0$Z_tilde)
LATE_1 <- cov(temp_X1$Y, temp_X1$Z_tilde) / 
  cov(temp_X1$D, temp_X1$Z_tilde)
(theta_0 * n_0 * LATE_0 + theta_1 * n_1 * LATE_1) / n
cov(data$Y, data$Z_tilde) / cov(data$D, data$Z_tilde)

# Relabeled data
data$Z_mod <- data$Z_tilde * (2 * data$X - 1)
cov(data$Y, data$Z_mod) / cov(data$D, data$Z_mod)









X <- sample(1:K, N, replace = T)
X <- sort(X)
threshold <- floor(K / 2)
X_minus <- X <= threshold

a_1 <- 0.1 # complier share
a_0 <- runif(K, min = 0.3, max = 0.7) # reference share
U_D <- runif(N, 0 , 1)

D <- rep(0, N)
index <- 0
for (k in 1:K){
  
  N_k <- sum(X == k)
  D[index + 1:N_k] <- as.numeric(
    (a_0[k] + a_1 * Z[index + 1:N_k] * (1 - 2 * (k <= threshold))) >= 
      U_D[index + 1:N_k])
  
  index <- index + N_k
  
}

# cell-by-cell regressions
data <- data.frame(D = D,
                   Z = Z,
                   X = X)
estimates <- matrix(0, K, 2)
 
for (k in 1:K){
  
  temp_data <- dplyr::filter(data, X == k)
  estimates[k, ] <- lm("D ~ Z", data = temp_data)$coefficients
  
}

estimates <- data.frame(estimates)
colnames(estimates) <- c("a_0", "a_1")

plot(estimates[, 2])
lines(x = 1:threshold, y = rep(-a_1, threshold))
lines(x = (threshold + 1):N, y = rep(a_1, N - threshold))

plot(a_0, estimates[, 1])

# Add in the outcome