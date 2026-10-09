################################# FRONT MATTER #################################

rm(list=ls())
gc()

library(dplyr)
library(ggplot2)
library(tidyr)
library(mvtnorm)
library(clusterGeneration)
library(xtable)

############################# CHECKING CALCULATIONS ############################

# Let X take K + 1 values, where K is an even number.
# Choose values of p(X,1) - p(X,0) that are evenly spaced between -delta and
# delta. Call this vector w. Let its limits be -delta and delta. We have
# w_0 = -delta
# w_K = delta
# w_k = delta * (2 * k / K - 1)

K <- 2
delta <- 0.2
w <- delta * (2 * (0:K) / K - 1)

# The selection model is set as
# 1{a + d(X) > V}
# We need to back out the d(X) that achieves a given p(X,1) - p(X,0). As  V is
# normally distributed
# p(X,0) = Phi(a) 
# p(X,1) = Phi(a + d(X))
# Hence,
# p(X,1) - p(X,0) = Phi(a + d(X)) - Phi(a)
# solving for d(X)
# d(X) = Phi^{-1}(Phi(a) + p(X,1) - p(X,0)) - a
# To simplify, we set a = 0 (so that Phi(a) = 0.5)
a <- 0
d <- qnorm(pnorm(a) + w) - a 

# Set the ATE
B <- 1 

#### Specify distributions ####

# Randomly generate a covariance matrix for the unobservables
A <- matrix(runif(9, -1, 1), 3, 3)
sigma <- A %*% t(A)
# Scale to be a correlation matrix
for (k in 1:3){
  
  sd_k <- sqrt(sigma[k, k])
  sigma[, k] <- sigma[, k] / sd_k
  sigma[k, ] <- sigma[k, ] / sd_k
  
}
rho_V1 <- sigma[3, 1]
rho_V0 <- sigma[3, 2]

# Distribution for X
p_x <- runif(K + 1, 0, 1)
p_x <- p_x / sum(p_x)

# Probability that Z = 1
pZ <- 1 / 2

#### Generate Variables ####

n <- 1000000

# Unobservables
U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
U <- data.frame(U)
names(U) <- c("U1", "U0", "V")

# X and Z
X <- sample(0:K, size = n, replace = T, prob = p_x)
Z <- rbinom(n, 1, pZ)

# Construct D and Y
D <- (a + d[X + 1] * Z >= U$V)
Y <- (1 - D) * U$U0 + D *(B + U$U1)

data <- data.frame(Y = Y, X = X, D = D, Z = Z)

#### Check Propensity Scores ####

w_hat <- rep(0, K + 1)
for (k in 0:K){
  
  temp <- dplyr::filter(data, X == k)
  w_hat[k + 1] <- mean(temp$D[temp$Z == 1]) - mean(temp$D[temp$Z == 0])
  
}
w_hat
w

#### Construct Z_tilde ####

# Construct E[Z|X =  x] as the sample average of Z with X = x
bar_Z <- rep(0, K + 1)
for (k in 0:K){
  
  bar_Z[k + 1] <- mean(data$Z[data$X == k])
  
}

Z_tilde <- rep(0, n)
for (i in 1:n){
  
  Z_tilde[i] <- Z[i] - bar_Z[X[i] + 1]
  
}

s_star <- rep(0, n)
for (i in 1:n){
  s_star[i] <- 2 * (w[X[i] + 1] >= 0) - 1
}

#### Check Denominator ####

mean(D * s_star * Z_tilde)
pZ * (1 - pZ) * sum(p_x * abs(w))
pZ * (1 - pZ) * sum(p_x * abs(pnorm(a + d) - pnorm(a)))
pZ * (1 - pZ) * sum(p_x * ((w >= 0) - (w < 0)) * (pnorm(a + d) - pnorm(a)))

#### Check Numerator ####

mean(Y * s_star * Z_tilde) 
B * pZ * (1 - pZ) * sum(p_x * ((w < 0) * (pnorm(a) - pnorm(a + d)) + 
                                 (w >= 0) * (pnorm(a + d) - pnorm(a)))) +
  pZ * (1 - pZ) * (rho_V1 - rho_V0) *
  sum(p_x * ((w >= 0) * (dnorm(a) - dnorm(a + d)) + 
               (w < 0) * (dnorm(a + d) - dnorm(a))))
B * pZ * (1 - pZ) * 
  sum(p_x * ((w >= 0) - (w < 0)) * (pnorm(a + d) - pnorm(a))) + 
  pZ * (1 - pZ) * (rho_V1 - rho_V0) *
  sum(p_x * ((w >= 0) - (w < 0)) * (dnorm(a) - dnorm(a + d)))


#### Check Numerator for Individual X Values ####

data$Z_tilde <- Z_tilde
data$s_star <- s_star

data_num <- rep(0, K + 1)
true_num <- rep(0, K + 1)

for (k in 0:K){
  
  temp <- dplyr::filter(data, X == k)
  data_num[k + 1] <- mean(temp$Y * temp$s_star * temp$Z_tilde)
  true_num[k + 1] <- pZ * (1 - pZ) * (
    B * ((w[k + 1] < 0) * (pnorm(a) - pnorm(a + d[k + 1])) + 
           (w[k + 1] >= 0) * (pnorm(a + d[k + 1]) - dnorm(a))) +
      (rho_V1 - rho_V0) * 
      ((w[k + 1] >= 0) * (dnorm(a) - dnorm(a + d[k + 1])) + 
         (w[k + 1] < 0) * (dnorm(a + d[k + 1]) - dnorm(a))))
  
}

data_num
true_num


#### Check the whole thing ####

mean(Y * s_star * Z_tilde) / mean(D * s_star * Z_tilde)
B + (rho_V1 - rho_V0) *
  sum(p_x * ((w >= 0) * (dnorm(a) - dnorm(a + d)) + 
               (w < 0) * (dnorm(a + d) - dnorm(a)))) / sum(p_x * abs(w))
B + (rho_V1 - rho_V0) *
  sum(p_x * (((w >= 0) - (w < 0)) * (dnorm(a) - dnorm(a + d)))) / 
  sum(p_x * abs(w))

data_tau <- mean(Y * s_star * Z_tilde) / mean(D * s_star * Z_tilde)
true_tau <- B + (rho_V1 - rho_V0) *
  sum(p_x * ((w >= 0) * (dnorm(a) - dnorm(a + d)) + 
               (w < 0) * (dnorm(a + d) - dnorm(a)))) / sum(p_x * abs(w))

data_tau
true_tau

############################ SIMPLE EXAMPLE TEST RUN ###########################

# X takes three values (-1, 0, 1) 
# Defiers at X = -1, compliers at X = 1, no switchers at X = 0
# Compare oracle and feasible estimators in two cases
# 1) P(X = 0) = 0, P(X = - 1) = P(X = 1) = 1 / 2
# 2) P(X = x) = 1 / 3 for x = -1,0,1

rm(list = ls())
gc()

# First stage parameters
K <- 2
delta <- 0.2
w <- delta * (2 * (0:K) / K - 1)
a <- 0
d <- qnorm(pnorm(a) + w) - a 

# Set the ATE
B <- 1 

#### Specify distributions ####

# Randomly generate a covariance matrix for the unobservables
A <- matrix(runif(9, -1, 1), 3, 3)
sigma <- A %*% t(A)
# Scale to be a correlation matrix
for (k in 1:3){
  
  sd_k <- sqrt(sigma[k, k])
  sigma[, k] <- sigma[, k] / sd_k
  sigma[k, ] <- sigma[k, ] / sd_k
  
}
rho_V1 <- sigma[3, 1]
rho_V0 <- sigma[3, 2]

# Distribution for X
p_x <- vector(mode = "list", length = 2)
p_x[[1]] <- c(1 / 2, 0, 1 / 2)
p_x[[2]] <- rep(1 / 3, 3)

# Probability that Z = 1
pZ <- 1 / 2

#### Generate Variables ####

l <- 2 # which distribution over X do you want to use
n <- 1000000

# Unobservables
U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
U <- data.frame(U)
names(U) <- c("U1", "U0", "V")

# X and Z
X <- sample(0:K, size = n, replace = T, prob = p_x[[l]])
Z <- rbinom(n, 1, pZ)

# Construct D and Y
D <- (a + d[X + 1] * Z >= U$V)
Y <- (1 - D) * U$U0 + D *(B + U$U1)

data <- data.frame(Y = Y, X = X, D = D, Z = Z)

#### Check Propensity Scores ####

w_hat <- rep(0, K + 1)
for (k in 0:K){
  
  temp <- dplyr::filter(data, X == k)
  w_hat[k + 1] <- mean(temp$D[temp$Z == 1]) - mean(temp$D[temp$Z == 0])
  
}
w_hat
w

#### Construct Z_tilde ####

# Construct E[Z|X =  x] as the sample average of Z with X = x
bar_Z <- rep(0, K + 1)
for (k in 0:K){
  
  bar_Z[k + 1] <- mean(data$Z[data$X == k])
  
}

Z_tilde <- rep(0, n)
for (i in 1:n){
  
  Z_tilde[i] <- Z[i] - bar_Z[X[i] + 1]
  
}

s_star <- rep(0, n)
for (i in 1:n){
  s_star[i] <- 2 * (w[X[i] + 1] >= 0) - 1
}

#### Check Denominator ####

mean(D * s_star * Z_tilde)
pZ * (1 - pZ) * sum(p_x[[l]] * abs(w))
pZ * (1 - pZ) * sum(p_x[[l]] * abs(pnorm(a + d) - pnorm(a)))
pZ * (1 - pZ) * sum(p_x[[l]] * ((w >= 0) - (w < 0)) * (pnorm(a + d) - pnorm(a)))

#### Check Numerator ####

mean(Y * s_star * Z_tilde) 
B * pZ * (1 - pZ) * sum(p_x[[l]] * ((w < 0) * (pnorm(a) - pnorm(a + d)) + 
                                      (w >= 0) * (pnorm(a + d) - pnorm(a)))) +
  pZ * (1 - pZ) * (rho_V1 - rho_V0) *
  sum(p_x[[l]] * ((w >= 0) * (dnorm(a) - dnorm(a + d)) + 
                    (w < 0) * (dnorm(a + d) - dnorm(a))))
B * pZ * (1 - pZ) * 
  sum(p_x[[l]] * ((w >= 0) - (w < 0)) * (pnorm(a + d) - pnorm(a))) + 
  pZ * (1 - pZ) * (rho_V1 - rho_V0) *
  sum(p_x[[l]] * ((w >= 0) - (w < 0)) * (dnorm(a) - dnorm(a + d)))


#### Check Numerator for Individual X Values ####

data$Z_tilde <- Z_tilde
data$s_star <- s_star

data_num <- rep(0, K + 1)
true_num <- rep(0, K + 1)

for (k in 0:K){
  
  temp <- dplyr::filter(data, X == k)
  data_num[k + 1] <- mean(temp$Y * temp$s_star * temp$Z_tilde)
  true_num[k + 1] <- pZ * (1 - pZ) * (
    B * ((w[k + 1] < 0) * (pnorm(a) - pnorm(a + d[k + 1])) + 
           (w[k + 1] >= 0) * (pnorm(a + d[k + 1]) - dnorm(a))) +
      (rho_V1 - rho_V0) * 
      ((w[k + 1] >= 0) * (dnorm(a) - dnorm(a + d[k + 1])) + 
         (w[k + 1] < 0) * (dnorm(a + d[k + 1]) - dnorm(a))))
  
}

data_num
true_num


#### Check the whole thing ####

mean(Y * s_star * Z_tilde) / mean(D * s_star * Z_tilde)
B + (rho_V1 - rho_V0) *
  sum(p_x[[l]] * ((w >= 0) * (dnorm(a) - dnorm(a + d)) + 
                    (w < 0) * (dnorm(a + d) - dnorm(a)))) / sum(p_x[[l]] * abs(w))
B + (rho_V1 - rho_V0) *
  sum(p_x[[l]] * (((w >= 0) - (w < 0)) * (dnorm(a) - dnorm(a + d)))) / 
  sum(p_x[[l]] * abs(w))

data_tau <- mean(Y * s_star * Z_tilde) / mean(D * s_star * Z_tilde)
true_tau <- B + (rho_V1 - rho_V0) *
  sum(p_x[[l]] * ((w >= 0) * (dnorm(a) - dnorm(a + d)) + 
                    (w < 0) * (dnorm(a + d) - dnorm(a)))) / sum(p_x[[l]] * abs(w))

data_tau
true_tau

##################### PROBABILITY OF CORRECT CLASSIFICATION ####################

# Works out the probability that we correctly classify switchers for X = x given
# that we use the difference in sample conditional means to classify.

# First stage D = 1{a + dZ >= V} with V normally distributed
# p(x, 1) = pnorm(a + d)
# p(x, 0) = pnorm(a)
# Fix n (number of observations with X = x), Pr(Z = 1|X = x), a and d. Work out 
# the probability that we correctly classify switchers

# fix parameters
delta <- 0
a <- 0 
d <- qnorm(pnorm(a) + delta) - a
p_Z1 <- 1 / 2 # Pr(Z = 1|X)

get_prob_correct <- function(k, n_Z0, n_Z1, a, d){
  
  # returns the probability that d_Z1 / n_Z1 >= d_Z0 / n_Z0
  # (i.e. the sample average of D with Z = 1 weakly exceeds the Z = 0 average)
  # rearranging, we need that
  # d_Z1 >= d_Z0 * n_Z1 / n_Z0 -> d_Z1 >= floor(d_Z0 * n_Z1 / n_Z0)
  # as all numbers are integers.
  
  # k is the number of D = 1 observations with Z = 0
  # n_Z0, n_Z1 are the numbers of observations with Z = 0 and Z = 1
  # a and d are propensity score parameters (p(x,z) = pnorm(a + dZ))
  
  return(
    pbinom(floor((k - 0.00001) * n_Z1 / n_Z0), n_Z1, pnorm(a + d), 
           lower.tail = FALSE) * 
      dbinom(k, n_Z0, pnorm(a)))
  
}

n <- seq(2, 50) # possible sample sizes
P_correct <- rep(0, length(n))

for (i in 1:length(n)){
  
  # loop over sample sizes
  
  # enumerate possible combinations of Z = 1 and Z = 0 
  n_Z1 <- seq(0, n[i]) # observations with Z = 1
  n_Z0 <- n[i] - n_Z1 # observations with Z = 0
  P_temp <- rep(0, n[i] - 1)
  
  # convention: if we have 0 observations with Z = 0 or Z = 1, then the cell is
  # zeroed out. 
  
  for (j in 2:n[i]){
    
    # i fixes n_Z0 and n_Z1
    # loop over the number of D = 1 draws for Z = 0
    P_temp[j - 1] <- sum(sapply(0:n_Z0[j],
                                get_prob_correct,
                                n_Z0 = n_Z0[j],
                                n_Z1 = n_Z1[j],
                                a = a,
                                d = d)) * dbinom(n_Z1[j], n[i], p_Z1)
    
  }
  
  P_correct[i] <- sum(P_temp)
  
}

plot(n, P_correct)

# validate by simulation
sims <- 10000
sim_correct <- rep(0, length(n))

for (i in 1:length(n)){
  
  for (s in 1:sims){
    
    Z <- rbinom(n[i], 1, p_Z1)
    V <- rnorm(n[i])
    D <- (a + d * Z >= V) 
    
    if (sum(Z) > 0 & sum(1 - Z > 0)){
      
      sim_correct[i] <- sim_correct[i] + 
        as.numeric(sum(D * Z) / sum(Z) >= sum(D * (1 - Z)) / sum(1 - Z))
      
    }
    
  }
  
  sim_correct[i] <- sim_correct[i] / sims
  
}

sim_correct

plot(n, P_correct)
points(n, sim_correct, col = "red")


#### BUILD A PLOT ####

# Plot the probability that switchers are classified correctly as a function of
# n for different values of the switcher mass

delta <- c(0, 0.05, 0.1, 0.15, 0.2)
a <- 0 
d <- qnorm(pnorm(a) + delta) - a
p_Z1 <- 1 / 2 # Pr(Z = 1|X)

n <- seq(2, 100) # possible sample sizes
P_correct <- matrix(0, length(n), length(d))

for (l in 1:length(d)){
  
  for (i in 1:length(n)){
    
    # loop over sample sizes
    
    # enumerate possible combinations of Z = 1 and Z = 0 
    n_Z1 <- seq(0, n[i]) # observations with Z = 1
    n_Z0 <- n[i] - n_Z1 # observations with Z = 0
    P_temp <- rep(0, n[i] - 1)
    
    # convention: if we have 0 observations with Z = 0 or Z = 1, then the cell is
    # zeroed out. 
    
    for (j in 2:n[i]){
      
      # i fixes n_Z0 and n_Z1
      # loop over the number of D = 1 draws for Z = 0
      P_temp[j - 1] <- sum(sapply(0:n_Z0[j],
                                  get_prob_correct,
                                  n_Z0 = n_Z0[j],
                                  n_Z1 = n_Z1[j],
                                  a = a,
                                  d = d[l])) * dbinom(n_Z1[j], n[i], p_Z1)
      
    }
    
    P_correct[i, l] <- sum(P_temp)
    
  }
  
}

plot(n, P_correct[, 5])

plot_data <- matrix(0, length(n) * length(d), 3)
plot_data[, 1] <- c(P_correct)
plot_data[, 2] <- rep(n, 5)
plot_data[, 3] <- rep(delta, each = length(n))
plot_data <- data.frame(plot_data)
names(plot_data) <- c("P", "n", "switchers")

plot_data %>% ggplot() +
  geom_point(aes(x = n, 
                 y = P, 
                 col = as.factor(switchers), 
                 shape = as.factor(switchers))) +
  theme_bw() +
  labs(color = expression(delta),
       shape = expression(delta))


