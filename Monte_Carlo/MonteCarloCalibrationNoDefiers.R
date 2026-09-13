#################################### SUMMARY ###################################

# A script used to tune the Monte Carlo processes in the "no defiers" case

# The first stage is assumed to be
# p(x,z) = Phi(a_0 + x'a + z * d_0)
# x is assumed to be a vector of independent N(0, 1) variables. Hence x'a will
# also be a normally distributed variable. 
# We have 3 quantities to calibrate: a_0, d_0, and sd_a^2

# Target
# 1) The share of compliers (E[p(x,Z = 1)] - E[p(x,Z = 0)])
# 2) The share of always takers (E[p(x,Z = 0)])
# 3) The standard deviation of p(x,Z = 0)

# The strategy is as follows:
# * Given sd_a and a_0, the first condition can be used to obtain d_0
# * Given sd_a, the second condition can be used to obtain a_0
# * The third condition depends only on a_0 and sd_a. 
# * We can reduce this to an equation in sd_a alone by using the second. 
# * Hence, we choose sd_a to satisfy the third condition, the back out the 
# remaining parameters. 

################################# FRONT MATTER #################################

rm(list = ls())
gc()

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)

################################### FUNCTIONS ##################################

# Function to calculate the sd of P(D = 1|X, Z = 0)
# Ep_Z is the average propensity score given Z = 0
# sd_Xt is sd_a
get_sd_Z <- function(Ep_Z, sd_Xt){
  
  mu <- qnorm(Ep_Z) * sqrt(1 + sd_Xt ^ 2)
  return(sqrt(pmvnorm(upper = rep(0, 2), 
                      mean = rep(- mu, 2),
                      sigma = rbind(c(1 + sd_Xt ^ 2, sd_Xt ^ 2),
                                    c(sd_Xt ^ 2, 1 + sd_Xt ^ 2)),
                      keepAttr = FALSE) - Ep_Z ^ 2))
  
}

#################################### TUNING ####################################

targets <- list(Ep_Z0 = 0.4, # E[p(X, Z = 0)]
                sd_Z0 = 0.15, # sd(p(X, Z = 0))
                C_share = 0.2) # share of compliers

e_tol <- 0.00000001

#### a_0 and sd_a ####

Ep_Z0 <- targets$Ep_Z0
sd_Z0 <- targets$sd_Z0

# initial guess
sd_a <- seq(0, 5, length.out = 101)
sd_Z <- sapply(sd_a, get_sd_Z, Ep_Z = Ep_Z0)

index <- which.min(abs(sd_Z - sd_Z0))
upper <- sd_a[min(length(sd_a), index + 1)]
lower <- sd_a[max(1, index - 1)]
error <- abs(sd_Z[index] - sd_Z0)

# quick loop to improve
while(error > e_tol){
  
  sd_a <- seq(lower, upper, length.out = 101)
  sd_Z <- sapply(sd_a, get_sd_Z, Ep_Z = Ep_Z0)
  
  index <- which.min(abs(sd_Z - sd_Z0))
  upper <- sd_a[min(length(sd_a), index + 1)]
  lower <- sd_a[max(1, index - 1)]
  error <- abs(sd_Z[index] - sd_Z0)
  
}

sd_a <- sd_a[index]
a_0 <-  qnorm(Ep_Z0) * sqrt(1 + sd_a ^ 2)

#### d_0 ####

C_share <- targets$C_share

# use that the complier share is E[p(x, Z = 1)] - E[p(x, Z = 0)]
d_0 <- sqrt(1 + sd_a ^ 2) * qnorm(Ep_Z0 + C_share) - a_0
g_0 <- a_0 + d_0

############################ TEST PARAMETER VALUES #############################

# Check targets against their sample counterparts for a DGP with parameters set 
# to the tuned values

k <- 2
n <- 100000 # sample size

# Generate parameters
a <- runif(k, -1, 1)
a <- sd_a * a / sqrt(sum(a ^ 2))

# Generate a sample to check calculations
X <- matrix(rnorm(n * k), n, k)

P_x0 <- pnorm(a_0 + X %*% as.matrix(a, k, 1))
P_x1 <- pnorm(g_0 + X %*% as.matrix(a, k, 1))

# Mean propensity score for Z = 0 observations
mean(P_x0)
targets$Ep_Z0

# Mean propensity score for Z = 1 observations
mean(P_x1)
targets$Ep_Z0 + targets$C_share

# Standard deviation for Z = 0 propensity score
sqrt(pmvnorm(upper = c(0, 0), 
             mean = rep(-a_0, 2), 
             sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                           c(sd_a ^ 2, 1 + sd_a ^ 2)),
             keepAttr = FALSE) - 
       pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2)) ^ 2)
targets$sd_Z0
sd(P_x0)

###################### BLOCK TWO APPROXIMATE CALIBRATION #######################

zeta_0 <- 0.1
B_0 <- 0
sd_B <- 0.45
pZ <- 0.5

B <- runif(k, -1, 1)
B <- B * sd_B / sqrt(sum(B ^ 2))

sigma <- diag(1, 3, 3)
sigma[1, 2:3] <- c(0.6, -0.4)
sigma[2:3, 1] <- sigma[1, 2:3]
colnames(sigma) <- c("V", "U_1", "U_0")

# Obtain the uncorrected TSLS estimate
Z <- rbinom(n, 1, pZ)
U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
colnames(U) <- c("V", "U1", "U0")
U <- data.frame(U)

index <- a_0 + X %*% as.matrix(a, k, 1) + Z * d_0
D <- as.numeric(index  > U[, 1])

# Outcomes
Y_0 <- B_0 + U$U0
Y_1 <- B_0 + zeta_0 + U$U1
Y <- (1 - D) * Y_0 + D * Y_1

# OLS bias
data <- data.frame(Y = Y, D = D)
OLS <- lm("Y ~ D", data = data)
zeta_hat <- OLS$coefficients[2]
bias <- zeta_hat - zeta_0
bias

# TSLS estimate
Z_tilde <- Z - mean(Z)
cov(Y, Z_tilde) / cov(D, Z_tilde)

################################## SAVE OUTPUT #################################

path <- getwd()
save_path <- paste0(path, "/Monte_Carlo")

parameters <- list(a_0 = a_0, 
                   d_0 = d_0,
                   sd_a = sd_a,
                   B_0 = B_0,
                   zeta_0 = zeta_0, 
                   sd_B = sd_B,
                   sigma = sigma)

saveRDS(parameters, file = paste0(path, "/Monte_Carlo/ParametersNoDefiers.rds"))

# Generate a table of parameter values
# temp <- unlist(parameters[!names(parameters) %in% "sigma"])
# 
# table <- data.frame(par = temp)
# 
# xtable(table, digits = 3)

