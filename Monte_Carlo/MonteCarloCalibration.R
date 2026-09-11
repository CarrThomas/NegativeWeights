#################################### SUMMARY ###################################

# A script used to calibrate the Monte Carlo processes

# The first stage is assumed to be
# p(x,z) = Phi(a_0 + x'a + z(d_0 + x'd))
# x is assumed to be a vector of independent N(0, 1) variables. Hence x'a and 
# x'd will be jointly normally distributed variables.
# We have 5 quantities to calibrate: a_0, d_0, and the covariance matrix 
# (diagonal elements sd_a^2, sd_d^2 off-diagonal element sd_a * sd_d * rho_ad)
# Let g_0 = a_0 + d_0, g = a + d, sd_g be sqrt(g'g), and rho_dg be g'd / 
# sd_g * sd_d. 

# The strategy is as follows
# 1) Choose a_0 and sd_a to target the average propensity score and its variance
# for Z = 0 observations
# 2) Choose g_0 and sd_g to target the average propensity score and its variance
# for Z = 1 observations
# 3) Given these parameters, choose sd_d to match the share of compliers.

# Note that the average propensity score given Z = 0 identifies the sum of 
# always takers and defiers, and the average given Z = 1 identifies the sum of
# always takers and compliers. As we target the share of compliers explicitly,
# we have three targets which pin down always-takers, defiers, and compliers. 
# Never-takers are then given by 1 less their sum. 

################################# FRONT MATTER #################################

rm(list = ls())
gc()

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

################################### FUNCTIONS ##################################

# Function to calculate the sd of P(D = 1|X, Z = z)
# Ep_Z is the average propensity score given Z = 1 or Z = 0
# sd_Xt is either sd_a or sd_d
get_sd_Z <- function(Ep_Z, sd_Xt){
  
  mu <- qnorm(Ep_Z) * sqrt(1 + sd_Xt ^ 2)
  return(sqrt(pmvnorm(upper = rep(0, 2), 
                      mean = rep(- mu, 2),
                      sigma = rbind(c(1 + sd_Xt ^ 2, sd_Xt ^ 2),
                                    c(sd_Xt ^ 2, 1 + sd_Xt ^ 2)),
                      keepAttr = FALSE) - Ep_Z ^ 2))
  
}

# Function to calculate the share of compliers.
# Compliers are calculated as the gap between p(x, z = 1) and p(x, z = 0) for
# x such that z encourages treatment. 
get_compliers <- function(sd_d, a_0, g_0, sd_a, sd_g){
  
  # derived parameters
  rho_ad <- (sd_g ^ 2 - (sd_a ^ 2 + sd_d ^ 2)) / (2 * sd_a * sd_d)
  rho_dg <- (sd_a * rho_ad + sd_d) / sd_g
  d_0 <- g_0 - a_0
  
  # always takers (integral of p(x, z = 0))
  integral_1 <- integrate(kernel, lower = - d_0, upper = Inf,
                          t_0 = a_0,
                          sd_d = sd_d,
                          sd_t = sd_a,
                          rho_td = rho_ad)$value
  
  # always takers and compliers (integral of p(x, z = 1))
  integral_2 <- integrate(kernel, lower = - d_0, upper = Inf,
                          t_0 = g_0,
                          sd_d = sd_d,
                          sd_t = sd_g,
                          rho_td = rho_dg)$value
  
  return(integral_2 - integral_1)
}

# function to be integrated
# consists of the density for Xd multiplied by either
# 1) Pr(Xa + V < - a_0|Xd) (set t_0 = a_0, sd_t = sd_a, rho_td = rho_ad)
# 2) Pr(Xg + V < - g_0|Xd) (set t_0 = g_0, sd_t = sd_g, rho_td = rho_gd)
kernel <- function(Xd, t_0, sd_d, sd_t, rho_td){
  
  sd_temp <- sqrt(1 + (sd_t ^ 2) * (1 - rho_td ^ 2))
  bound <- t_0 + rho_td * sd_t * Xd / sd_d 
  
  return(pnorm(q = bound, mean = 0, sd = sd_temp) * dnorm(x = Xd, sd = sd_d))
  
}


# Calculates shares of always takers, never takers, defiers and compliers
# always_takers_1, never_takers_1 -> the region where Z encourages treatment.
# always_takers_2, never_takers_2 -> the region where Z discourages treatment.
get_shares <- function(a_0, d_0, g_0, sd_a, sd_d, sd_g, rho_ad, rho_dg){
  
  always_takers_1 <- integrate(kernel, lower = - d_0, upper = Inf,
                               t_0 = a_0,
                               sd_d = sd_d,
                               sd_t = sd_a,
                               rho_td = rho_ad)$value
  
  never_takers_1 <- 1 - pnorm(- d_0 / sd_d) -
    integrate(kernel, lower = - d_0, upper = Inf,
              t_0 = g_0,
              sd_d = sd_d,
              sd_t = sd_g,
              rho_td = rho_dg)$value
  
  compliers <- 1 - pnorm(q = -d_0, sd = sd_d) - always_takers_1 - never_takers_1
  
  always_takers_2 <- integrate(kernel, lower = - Inf, upper = - d_0,
                               t_0 = g_0,
                               sd_d = sd_d,
                               sd_t = sd_g,
                               rho_td = rho_dg)$value
  
  never_takers_2 <- pnorm(- d_0 / sd_d) -
    integrate(kernel, lower = - Inf, upper = - d_0,
              t_0 = a_0,
              sd_d = sd_d,
              sd_t = sd_a,
              rho_td = rho_ad)$value
  
  defiers <- pnorm(q = -d_0, sd = sd_d) - always_takers_2 - never_takers_2
  
  shares <- c(always_takers_1,
              always_takers_2,
              never_takers_1,
              never_takers_2,
              compliers,
              defiers)
  names(shares) <- c("always_takers_1",
                     "always_takers_2",
                     "never_takers_1",
                     "never_takers_2",
                     "compliers",
                     "defiers")
  
  return(shares)
  
}


################################## CALIBRATION ################################# 

#List of calibration targets
targets <- list(Ep_Z0 = 0.49, # E[p(x, z = 0)] (AT + D)
                sd_Z0 = 0.15, # sd(p(x, z = 0))
                Ep_Z1 = 0.51, # E[p(x, z = 1)] (AT + C)
                sd_Z1 = 0.1, # sd(p(x, z = 1)
                C_share = 0.11) # share of compliers

e_tol <- 0.00000001 # error tolerance

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

#### g_0 and sd_g ####

Ep_Z1 <- targets$Ep_Z1
sd_Z1 <- targets$sd_Z1

# initial guess
sd_g <- seq(0, 5, length.out = 101)
sd_Z <- sapply(sd_g, get_sd_Z, Ep_Z = Ep_Z1)

index <- which.min(abs(sd_Z - sd_Z1))
upper <- sd_g[min(length(sd_g), index + 1)]
lower <- sd_g[max(1, index - 1)]
error <- abs(sd_Z[index] - sd_Z1)

while(error > e_tol){
  
  sd_g <- seq(lower, upper, length.out = 101)
  sd_Z <- sapply(sd_g, get_sd_Z, Ep_Z = Ep_Z1)
  
  index <- which.min(abs(sd_Z - sd_Z1))
  upper <- sd_g[min(length(sd_g), index + 1)]
  lower <- sd_g[max(1, index - 1)]
  error <- abs(sd_Z[index] - sd_Z1)
  
}

sd_g <- sd_g[index]
g_0 <-  qnorm(Ep_Z1) * sqrt(1 + sd_g ^ 2)

#### Compliers ####

C_share <- targets$C_share

# initial guess
sd_d <- seq(abs(sd_g - sd_a), sd_g + sd_a, length.out = 101)
compliers <- sapply(sd_d, 
                    get_compliers, 
                    a_0 = a_0, 
                    g_0 = g_0, 
                    sd_a = sd_a, 
                    sd_g = sd_g)

# initial guess 
index <- which.min(abs(compliers - C_share))
upper <- sd_d[min(length(compliers), index + 1)]
lower <- sd_d[max(1, index - 1)]
error <- abs(compliers[index] - C_share)

if (length(unique(sign(compliers - C_share))) == 1){
  
  print("targets are impossible!")
  
}else{
  
  while(error > e_tol){
    
    sd_d <- seq(lower, upper, length.out = 101)
    
    compliers <- sapply(sd_d, 
                        get_compliers, 
                        a_0 = a_0, 
                        g_0 = g_0, 
                        sd_a = sd_a, 
                        sd_g = sd_g)
    
    index <- which.min(abs(compliers - C_share))
    upper <- sd_d[min(length(compliers), index + 1)]
    lower <- sd_d[max(1, index - 1)]
    error <- abs(compliers[index] - C_share)
    
  }
  
  sd_d <- sd_d[index]
  rho_ad <- (sd_g ^ 2 - (sd_a ^ 2 + sd_d ^ 2)) / (2 * sd_a * sd_d)
  rho_dg <- (sd_a * rho_ad + sd_d) / sd_g
  d_0 <- g_0 - a_0
  
}

############################ TEST PARAMETER VALUES #############################
 
# Calculate shares of different groups by numerical integration
shares <- get_shares(a_0, d_0, g_0, sd_a, sd_d, sd_g, rho_ad, rho_dg)
shares[5] - targets$C_share # compliers
shares[1] + shares[2] - (targets$Ep_Z1- targets$C_share) # always-takers
shares[6] - (targets$Ep_Z0 + targets$C_share - targets$Ep_Z1)

# Check targets against their sample counterparts for a DGP with parameters set 
# to the tuned values

k <- 2 # number of x variables
n <- 1000000 # sample size

# Generate a and d which match sd_a, sd_d, and rho_ad
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
g <- a + d

# verify that a and d match the targets
sd_a - sqrt(sum(a ^ 2))
sd_d - sqrt(sum(d ^ 2))
sd_a * sd_d * rho_ad - sum(a * d)

# draw X and calculate propensity scores
X <- matrix(rnorm(n * k), n, k)
P_x0 <- pnorm(a_0 + X %*% as.matrix(a, k, 1))
P_x1 <- pnorm(g_0 + X %*% as.matrix(g, k, 1))
ind <- (P_x1 > P_x0)

# Shares
shares[1]
sum(P_x0[ind]) / n

shares[2]
sum(P_x1[!ind]) / n

shares[3]
sum((1 - P_x1[ind])) / n

shares[4]
sum((1 - P_x0[!ind])) / n

shares[5]
sum((P_x1[ind] - P_x0[ind])) / n

shares[6]
sum((P_x0[!ind] - P_x1[!ind])) / n

# Treatment rates for Z = 0 and Z = 1 observations
targets$Ep_Z0
mean(P_x0)

targets$Ep_Z1
mean(P_x1)

# Variance of the Z = 0 and Z = 1 treatment rates
sqrt(pmvnorm(upper = c(0, 0), 
             mean = rep(-a_0, 2), 
             sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                           c(sd_a ^ 2, 1 + sd_a ^ 2)),
             keepAttr = FALSE) - 
       pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2)) ^ 2)
sd(P_x0)
targets$sd_Z0

sqrt(pmvnorm(upper = c(0, 0), 
             mean = rep(-g_0, 2), 
             sigma = rbind(c(1 + sd_g ^ 2, sd_g ^ 2),
                           c(sd_g ^ 2, 1 + sd_g ^ 2)),
             keepAttr = FALSE) - 
       pnorm(0, mean = -g_0, sd = sqrt(1 + sd_g ^ 2)) ^ 2)
sd(P_x1)
targets$sd_Z1

# Complier share
sum((P_x1 - P_x0)[P_x1 > P_x0]) / n
targets$C_share

# ################################### BLOCK TWO ##################################
# 
# # Calibrate outcome parameters to target the OLS bias and the corrected TSLS 
# # estimate
# 
# # Parameters to be calibrated are the elements of the covariance matrix for the
# # unobservables. 
# 
# # Calibrate sigma_1 rho_1v and sigma_0 rho_0v to match the corrected TSLS 
# # estimate and the OLS bias. These conditions only identify those products. The
# # individual parameters can be set freely (e.g. sd = 1 if possible)
# 
# # Our two targets give a linear system with unknowns sigma_1 rho_1v and
# # sigma_0 rho_0v. We simply obtain the coefficients of this system and solve
# # for the unknown. 
# 
# zeta_0 <- 0.1 # ATE
# pZ <- 1 / 2
# 
# block_two_targets <- c(0.05, zeta_0 + 0.01) 
# names(block_two_targets) <- c("bias", "estimand")
# 
# # Constants for the TSLS equation
# C_1 <- zeta_0 
# C_2 <- (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
#           (1 - pnorm(-d_0,
#                      mean = -a_0 * sd_d * sd_a * rho_ad / 
#                        (1 + sd_a ^ 2),
#                      sd = sd_d * 
#                        sqrt(1 - (rho_ad ^ 2) * (sd_a ^ 2) / 
#                               (1 + sd_a ^ 2)))) - 
#           dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
#           (1 - pnorm(-d_0,
#                      mean = -g_0 * sd_d * sd_g * rho_dg / 
#                        (1 + sd_g ^ 2),
#                      sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
#                                         (1 + sd_g ^ 2)))) -
#           dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
#           pnorm(-d_0,
#                 mean = -a_0 * sd_d * sd_a * rho_ad / 
#                   (1 + sd_a ^ 2),
#                 sd = sd_d * sqrt(1 - rho_ad ^ 2 * sd_a ^ 2 / 
#                                    (1 + sd_a ^ 2))) + 
#           dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
#           pnorm(-d_0,
#                 mean = -g_0 * sd_d * sd_g * rho_dg / 
#                   (1 + sd_g ^ 2),
#                 sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
#                                    (1 + sd_g ^ 2)))) / (shares[5] + shares[6])
# 
# # constants for the bias equation
# integrals <- rep(0, 4)
# 
# integrals[1] <- integrate(
#   f = Vectorize(function(xg, g_0, sd_g){
#     if(pnorm(g_0 + xg) > 0){
#       return(dnorm(g_0 + xg) * dnorm(xg, sd = sd_g) / pnorm(g_0 + xg))
#     }else{
#       return(0)
#     }}, SIMPLIFY = TRUE),
#   lower = -Inf,
#   upper = Inf,
#   g_0 = g_0, 
#   sd_g = sd_g)$value
# 
# integrals[2] <- integrate(
#   f = Vectorize(function(xa, a_0, sd_a){
#     if(pnorm(a_0 + xa) > 0){
#       return(dnorm(a_0 + xa) * dnorm(xa, sd = sd_a) / pnorm(a_0 + xa))
#     }else{
#       return(0)
#     }}, SIMPLIFY = TRUE),
#   lower = -Inf,
#   upper = Inf,
#   a_0 = a_0, 
#   sd_a = sd_a)$value
# 
# integrals[3] <- integrate(
#   f = Vectorize(function(xg, g_0, sd_g){
#     if((1 - pnorm(g_0 + xg)) > 0){
#       return(dnorm(g_0 + xg) * dnorm(xg, sd = sd_g) / (1 - pnorm(g_0 + xg)))
#     }else{
#       return(0)
#     }}, SIMPLIFY = TRUE),
#   lower = -Inf,
#   upper = Inf,
#   g_0 = g_0, 
#   sd_g = sd_g)$value
# 
# integrals[4] <- integrate(
#   f = Vectorize(function(xa, a_0, sd_a){
#     if((1 - pnorm(a_0 + xa) > 0)){
#       return(dnorm(a_0 + xa) * dnorm(xa, sd = sd_a) / (1 - pnorm(a_0 + xa)))
#     }else{
#       return(0)
#     }}, SIMPLIFY = TRUE),
#   lower = -Inf,
#   upper = Inf,
#   a_0 = a_0, 
#   sd_a = sd_a)$value
# 
# # analytically calculated bias
# C_3 <- - pZ * integrals[1] - (1 - pZ) * integrals[2]
# C_4 <- - pZ * integrals[3] - (1 - pZ) * integrals[4]
# 
# # set up the linear system
# A <- cbind(c(C_2, C_3), c(-C_2, C_4))
# B <- matrix(block_two_targets - c(C_1, 0), nrow = 2, ncol = 1)
# 
# intermediate <- solve(A, B) 
# 
# sigma <- diag(1, 3, 3)
# sigma[1, 2:3] <- intermediate
# sigma[2:3, 1] <- intermediate
# colnames(sigma) <- c("V", "U1", "U0")
# 
# lambda <- sum(intermediate[1:2]^2) + 0.1
# diag(sigma)[2:3] <- diag(sigma)[2:3] * lambda
# det(sigma)
# 
# ############################ TEST PARAMETER VALUES #############################
# 
# # Obtain the uncorrected TSLS estimate
# Z <- rbinom(n, 1, pZ)
# U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
# colnames(U) <- c("V", "U1", "U0")
# U <- data.frame(U)
# 
# index <- a_0 + X %*% as.matrix(a, k, 1) + Z * (d_0 + X %*% as.matrix(d, k, 1))
# P <- pnorm(index)
# D <- as.numeric(index  > U[, 1]) 
# 
# B_0 <- 0
# 
# # Outcomes
# Y_0 <- B_0 + U$U0
# Y_1 <- B_0 + zeta_0 + U$U1
# Y <- (1 - D) * Y_0 + D * Y_1
# 
# # OLS bias
# data <- data.frame(Y = Y, D = D)
# OLS <- lm("Y ~ D", data = data)
# zeta_hat <- OLS$coefficients[2]
# bias <- zeta_hat - zeta_0
# 
# integrals <- rep(0, 4)
# 
# integrals[1] <- - sigma[1, 2] * 
#   integrate(
#     f = Vectorize(function(xg, g_0, sd_g){
#       if(pnorm(g_0 + xg) > 0){
#         return(dnorm(g_0 + xg) * dnorm(xg, sd = sd_g) / pnorm(g_0 + xg))
#       }else{
#         return(0)
#       }}, SIMPLIFY = TRUE),
#     lower = -Inf,
#     upper = Inf,
#     g_0 = g_0, 
#     sd_g = sd_g)$value
# 
# integrals[2] <- - sigma[1, 2] * 
#   integrate(
#     f = Vectorize(function(xa, a_0, sd_a){
#       if(pnorm(a_0 + xa) > 0){
#         return(dnorm(a_0 + xa) * dnorm(xa, sd = sd_a) / pnorm(a_0 + xa))
#       }else{
#         return(0)
#       }}, SIMPLIFY = TRUE),
#     lower = -Inf,
#     upper = Inf,
#     a_0 = a_0, 
#     sd_a = sd_a)$value
# 
# integrals[3] <- sigma[1, 3] *  
#   integrate(
#   f = Vectorize(function(xg, g_0, sd_g){
#     if((1 - pnorm(g_0 + xg)) > 0){
#       return(dnorm(g_0 + xg) * dnorm(xg, sd = sd_g) / (1 - pnorm(g_0 + xg)))
#     }else{
#       return(0)
#     }}, SIMPLIFY = TRUE),
#   lower = -Inf,
#   upper = Inf,
#   g_0 = g_0, 
#   sd_g = sd_g)$value
# 
# integrals[4] <- sigma[1, 3] *  
#   integrate(
#   f = Vectorize(function(xa, a_0, sd_a){
#     if((1 - pnorm(a_0 + xa) > 0)){
#       return(dnorm(a_0 + xa) * dnorm(xa, sd = sd_a) / (1 - pnorm(a_0 + xa)))
#     }else{
#       return(0)
#     }}, SIMPLIFY = TRUE),
#   lower = -Inf,
#   upper = Inf,
#   a_0 = a_0, 
#   sd_a = sd_a)$value
# 
# mean(U$U1[D == 1]) - mean(U$U0[D == 0])
# bias
# 
# # analytically calculated bias
# pZ * (integrals[1] - integrals[3]) + (1 - pZ) * (integrals[2] - integrals[4])
# 
# # TSLS numerator
# Z_tilde <- Z - mean(Z)
# cov(Y, Z_tilde)
# 
# # TSLS estimate
# cov(Y, Z_tilde) / cov(D, Z_tilde)
# (zeta_0 * (pnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) - 
#              pnorm(a_0, sd = sqrt(1 + sd_a ^ 2))) + 
#     (sigma[1, 2] - sigma[1, 3]) * 
#     (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) - 
#        dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)))) / 
#   (pnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) - 
#      pnorm(a_0, sd = sqrt(1 + sd_a ^ 2)))
# 
# cov(Y - B_0 - zeta_0 * D, Z_tilde) / cov(D, Z_tilde)
# ((sigma[1, 2] - sigma[1, 3]) * 
#     (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) - 
#        dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)))) / 
#   (pnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) - 
#      pnorm(a_0, sd = sqrt(1 + sd_a ^ 2)))
# 
# # Corrected TSLS numerator
# s <- sign(d_0 + X %*% as.matrix(d, k, 1))
# cov(D, s * Z_tilde) / (pZ * (1 - pZ)) # should be compliers + defiers
# cov(D, Z_tilde) / (pZ * (1 - pZ)) # should be compliers - defiers
# 
# cov(Y, s * Z_tilde)
# 
# rho_dg <- sum(d * g) / (sd_d * sd_g)
# 
# # analytical 
# pZ * (1 - pZ) * (zeta_0 * 0.2 + 
#                    (sigma[1, 2] - sigma[1, 3]) *
#                    (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
#                       (1 - pnorm(-d_0,
#                                  mean = -a_0 * sd_d * sd_a * rho_ad / 
#                                    (1 + sd_a ^ 2),
#                                  sd = sd_d * 
#                                    sqrt(1 - (rho_ad ^ 2) * (sd_a ^ 2) / 
#                                           (1 + sd_a ^ 2)))) - 
#                       dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
#                       (1 - pnorm(-d_0,
#                                  mean = -g_0 * sd_d * sd_g * rho_dg / 
#                                    (1 + sd_g ^ 2),
#                                  sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
#                                                     (1 + sd_g ^ 2)))) -
#                       dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
#                       pnorm(-d_0,
#                             mean = -a_0 * sd_d * sd_a * rho_ad / 
#                               (1 + sd_a ^ 2),
#                             sd = sd_d * sqrt(1 - rho_ad ^ 2 * sd_a ^ 2 / 
#                                                (1 + sd_a ^ 2))) + 
#                       dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
#                       pnorm(-d_0,
#                             mean = -g_0 * sd_d * sd_g * rho_dg / 
#                               (1 + sd_g ^ 2),
#                             sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
#                                                (1 + sd_g ^ 2)))))
# 
# cov(Y - B_0 - D * zeta_0, s * Z_tilde)
# # analytical 
# pZ * (1 - pZ) * ((sigma[1, 2] - sigma[1, 3]) *
#                    (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
#                       (1 - pnorm(-d_0,
#                                  mean = -a_0 * sd_d * sd_a * rho_ad / 
#                                    (1 + sd_a ^ 2),
#                                  sd = sd_d * 
#                                    sqrt(1 - (rho_ad ^ 2) * (sd_a ^ 2) / 
#                                           (1 + sd_a ^ 2)))) - 
#                       dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
#                       (1 - pnorm(-d_0,
#                                  mean = -g_0 * sd_d * sd_g * rho_dg / 
#                                    (1 + sd_g ^ 2),
#                                  sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
#                                                     (1 + sd_g ^ 2)))) -
#                       dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
#                       pnorm(-d_0,
#                             mean = -a_0 * sd_d * sd_a * rho_ad / 
#                               (1 + sd_a ^ 2),
#                             sd = sd_d * sqrt(1 - rho_ad ^ 2 * sd_a ^ 2 / 
#                                                (1 + sd_a ^ 2))) + 
#                       dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
#                       pnorm(-d_0,
#                             mean = -g_0 * sd_d * sd_g * rho_dg / 
#                               (1 + sd_g ^ 2),
#                             sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
#                                                (1 + sd_g ^ 2)))))
# 
# # corrected estimate
# cov(Y, s * Z_tilde) / cov(D, s * Z_tilde)
# 
# zeta_0 + ((sigma[1, 2] - sigma[1, 3]) *
#             (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
#                (1 - pnorm(-d_0,
#                           mean = -a_0 * sd_d * sd_a * rho_ad / 
#                             (1 + sd_a ^ 2),
#                           sd = sd_d * 
#                             sqrt(1 - (rho_ad ^ 2) * (sd_a ^ 2) / 
#                                    (1 + sd_a ^ 2)))) - 
#                dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
#                (1 - pnorm(-d_0,
#                           mean = -g_0 * sd_d * sd_g * rho_dg / 
#                             (1 + sd_g ^ 2),
#                           sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
#                                              (1 + sd_g ^ 2)))) -
#                dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
#                pnorm(-d_0,
#                      mean = -a_0 * sd_d * sd_a * rho_ad / 
#                        (1 + sd_a ^ 2),
#                      sd = sd_d * sqrt(1 - rho_ad ^ 2 * sd_a ^ 2 / 
#                                         (1 + sd_a ^ 2))) + 
#                dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
#                pnorm(-d_0,
#                      mean = -g_0 * sd_d * sd_g * rho_dg / 
#                        (1 + sd_g ^ 2),
#                      sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
#                                         (1 + sd_g ^ 2))))) / 0.2

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

index <- a_0 + X %*% as.matrix(a, k, 1) + Z * (d_0 + X %*% as.matrix(d, k, 1))
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

# Corrected TSLS numerator
s <- sign(d_0 + X %*% as.matrix(d, k, 1))
cov(D, s * Z_tilde) / (pZ * (1 - pZ)) # should be compliers + defiers
cov(D, Z_tilde) / (pZ * (1 - pZ)) # should be compliers - defiers

cov(Y, s * Z_tilde) / cov(D, s * Z_tilde)

################################## SAVE OUTPUT #################################

path <- getwd()
save_path <- paste0(path, "/Monte_Carlo")

parameters <- list(a_0 = a_0, 
                   d_0 = d_0,
                   sd_a = sd_a,
                   sd_d = sd_d,
                   rho_ad = rho_ad,
                   B_0 = B_0,
                   zeta_0 = zeta_0, 
                   sd_B = sd_B,
                   sigma = sigma)

saveRDS(parameters, file = paste0(path, "/Monte_Carlo/Parameters1.rds"))

# Generate a table of parameter values
# temp <- unlist(parameters[!names(parameters) %in% "sigma"])
# 
# table <- data.frame(par = temp)
# 
# xtable(table, digits = 3)

#################################### ERRATA #################################### 

# Alternative method for the final `block':
# 1) consider rho_ad between 0 and 1
# 2) given rho_ad and the other parameters calculated earlier (a_0, g_0, sd_a, 
# sd_g), calculate sd_d and rho_dg
# 3) evaluate the complier share
# This achieves the same end, but it is a bit easier to set sd_d than rho_ad. If
# we set rho_ad, we have to solve a quadratic for sd_d. 

# get_compliers <- function(rho_ad, a_0, g_0, sd_a, sd_g){
#   
#   # calculate the gap between p(x,z = 1) and p(x, z = 0) in the region where Z 
#   # encourages treatment 
#   
#   # back out other parameters
#   sd_d <- - sd_a * rho_ad + sqrt((sd_a ^ 2) *(rho_ad ^ 2 - 1) + sd_g ^ 2)
#   rho_dg <- (sd_a * rho_ad + sd_d) / sd_g
#   d_0 <- g_0 - a_0
#   
#   # always takers (integral of p(x, z = 0))
#   integral_1 <- integrate(kernel, lower = - d_0, upper = Inf,
#                              t_0 = a_0,
#                              sd_d = sd_d,
#                              sd_t = sd_a,
#                              rho_td = rho_ad)$value
#   
#   # always takers and compliers (integral of p(x, z = 1))
#   integral_2 <- integrate(kernel, lower = - d_0, upper = Inf,
#                           t_0 = g_0,
#                           sd_d = sd_d,
#                           sd_t = sd_g,
#                           rho_td = rho_dg)$value
#   
#   return(integral_2 - integral_1)
# }
# 
# if (sd_g > sd_a){
#   
#   rho_ad <- seq(-1, 1, length.out = 101)
#   
# }else{
#   
#   rho_ad <- seq(-1, -sqrt(1 - sd_g ^ 2 / sd_a ^ 2), length.out = 101)
#   
# }
# 
# compliers <- sapply(rho_ad, 
#                     get_compliers, 
#                     a_0 = a_0, 
#                     g_0 = g_0, 
#                     sd_a = sd_a, 
#                     sd_g = sd_g)
# 
# # initial guess 
# index <- which.min(abs(compliers - C))
# upper <- rho_ad[min(length(compliers), index + 1)]
# lower <- rho_ad[max(1, index - 1)]
# error <- abs(compliers[index] - C)
# 
# while(error > 0.0001){
#   
#   rho_ad <- seq(lower, upper, length.out = 101)
#   compliers <- sapply(rho_ad, 
#                       get_compliers, 
#                       a_0 = a_0, 
#                       g_0 = g_0, 
#                       sd_a = sd_a, 
#                       sd_g = sd_g)
#   
#   index <- which.min(abs(compliers - C))
#   upper <- rho_ad[min(length(compliers), index + 1)]
#   lower <- rho_ad[max(1, index - 1)]
#   error <- abs(compliers[index] - C)
#   
# }
# 
# rho_ad <- rho_ad[index]