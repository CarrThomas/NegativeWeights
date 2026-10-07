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
# X = K / 2
# Compare oracle and feasible estimators in two cases
# 1) P(X = x) = 1 / K for all x != K / 2, P(X = K / 2) = 0
# 2) P(X = x) = 1 / (K + 1) for all x

################################### FUNCTIONS ##################################

gen_data <- function(n, sigma, p_x, pZ, K, a, d, B){
  
  # Generates (D,X,Y,Z)
  # inputs
  # n: number of observations
  # sigma: covariance matrix for the unobservables
  # p_x: distribution for x
  # pZ: probability Z = 1
  # K: x values range from 0 to K
  # a,d: first stage D = 1(a + d * Z >= V)
  # B: ATE
  
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
  
  return(data.frame(D = D, X = X, Y = Y, Z = Z))
  
}

gen_vars <- function(data, K, true_s){
  
  # adds in the demeaned Z, true sign function, and estimated sign function
  # inputs:
  # data data frame containing D, X, Y, Z
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
        2 * (sum(temp$Z * temp$D) / n_1 > sum((1 - temp$Z) * temp$D) / n_0) - 1 
    }
  }
  data$s_hat <- est_s[data$X + 1]
  
  return(data)
  
}

################################## PARAMETERS ##################################

# Number of X values
K <- 2

# Switcher shares at each value of X are evenly spaced between -delta and delta
delta <- 0.2
w <- delta * (2 * (0:K) / K - 1)

# First stage 1(a + dZ >= V)
a <- 0
d <- qnorm(pnorm(a) + w) - a 

# ATE
B <- 0 

# True value of the sign function for each value of X
true_s <- 2 * (w >= 0) - 1

################################# DISTRIBUTIONS ################################

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

################################# GENERATE DATA ################################

n <- c(100, 200, 500, 1000) # sample sizes
runs <- 10 # number of draws for each sample size

tau_oracle <- matrix(0, length(n), runs)
tau_feasible <- matrix(0, length(n), runs)

for (i in 1:length(n)){
  for (j in 1:runs){
    
    data <- gen_data(n = n[i],
                     sigma = sigma,
                     p_x = p_x[[1]],
                     pZ = pZ,
                     K = K,
                     a = a,
                     d = d,
                     B = B)
    data <- gen_vars(data = data,
                     K = K,
                     true_s = true_s)
    
    tau_oracle[i, j] <- sum(data$Y * data$s_star * data$Z_tilde) / 
      sum(data$D * data$s_star * data$Z_tilde)
    tau_feasible[i, j] <- sum(data$Y * data$s_hat * data$Z_tilde) / 
      sum(data$D * data$s_hat * data$Z_tilde)
    
    
  }
}

# summary statistics




# set first stage parameters
a_0 <- c(qnorm(0.2), qnorm(0.6)) # pnorm(a_0) is the lower propensity score
a_1 <- c(qnorm(0.6) - a_0[1], qnorm(0.2) - a_0[2])
# pnorm(a_0 + a_1) is the higher propensity score

# covariance matrix for the unobservables
# order is (V, U(0), U(1))
A <- matrix(runif(9, -1, 1), 3, 3)
sigma <- A %*% t(A)
# scale the first diagonal element to be one
sigma <- sigma / sigma[1, 1]

# specify a covariance
# to simplify, set diagonal to one -> we directly specify a correlation matrix
sigma <- diag(1, 3)
# want positive selection -> those with the highest treatment effects are more 
# likely to select treatment.
# V is `resistance to treatment` -> high V means less likely to take treatment
# to generate positive selection, we want V to be positive related to U(0) and
# negatively related to U(1)
# 
sigma[1, 2] <- 0.5
sigma[1, 3] <- -0.5
sigma[, 1] <- sigma[1 ,]

# generate unobservables
U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
U <- data.frame(U)
names(U) <- c("V", "U_0", "U_1")

# check HTE are correctly specified
mean(U$U_0)
mean(U$U_1)
# always takers -> V below a_0[1]
alwaystake <- U$V <= a_0[1]
mean((U$U_1 - U$U_0)[alwaystake])
# compliers -> a_0[1] < V < a_0[1] + a_1[1]
complier <- as.logical((a_0[1] < U$V ) * (U$V <= a_0[1] + a_1[1]))
mean((U$U_1 - U$U_0)[complier])
# never takers -> pnrom(a_0[1] + a_1[1]) < V
nevertake <- a_0[1] + a_1[1]  < U$V 
mean((U$U_1 - U$U_0)[nevertake])

(sigma[3, 1] - sigma[2, 1]) * (dnorm(a_0[1]) - dnorm(a_1[1] + a_0[1])) / 
  (pnorm(a_1[1] + a_0[1]) - pnorm(a_0[1]))

# generate D
D <- as.numeric(((1 - X) * (a_0[1] + a_1[1] * Z) + 
                   X * (a_0[2] + a_1[2] * Z)) >= U[, 1])

# check D is correct
mean(D[Z == 1 & X == 0]) - pnorm(a_0[1] + a_1[1])
mean(D[Z == 0 & X == 0]) - pnorm(a_0[1])
mean(D[Z == 1 & X == 1]) - pnorm(a_0[2] + a_1[2])
mean(D[Z == 0 & X == 1]) - pnorm(a_0[2])

# Epsilon chosen to induce a sign flip
mu_0 <- rep(0, 2)
mu_1 <- rep(0, 2) 
tau_0 <- mu_1[1] - mu_0[1] + (sigma[3, 1] - sigma[2, 1]) *
  (dnorm(a_0[1]) - dnorm(a_1[1] + a_0[1])) / 
  (pnorm(a_1[1] + a_0[1]) - pnorm(a_0[1]))
A <- (1 / 2 - gamma) / (2 * gamma)
# threshold for the TSLS estimate to be zero
tau_0 / A

epsilon <- 0.5
mu_1[2] <- mu_1[2] + epsilon

# generate Y
Y_1 <- mu_1[X + 1] + U$U_1
Y_0 <- mu_0[X + 1] + U$U_0
Y <- D * Y_1 + (1 - D) * Y_0

#### Estimate LATE in subsamples #### 

data <- data.frame(Y = Y, D = D, Z = Z, X = X)
data$Z_tilde <- data$Z - mean(data$Z)

# X = 0 subsample
data_X0 <- dplyr::filter(data, X == 0)
data_X0$Z_tilde <- data_X0$Z - mean(data_X0$Z)
D_hat <- lm(D ~ Z, data = data_X0)$fitted.values
data_X0$D_hat <- D_hat
lm(Y ~ D_hat, data = data_X0)

# alternative calculation
cov(data_X0$Y, data_X0$Z_tilde) / cov(data_X0$D, data_X0$Z_tilde)

# Analytical Calculation
mu_1 - mu_0 + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(a_0[1]) - dnorm(a_0[1] + a_1[1])) /
  (pnorm(a_0[1] + a_1[1]) - pnorm(a_0[1]))

# X = 1 subsample
data_X1 <- dplyr::filter(data, X == 1)
data_X1$Z_tilde <- data_X1$Z - mean(data_X1$Z)
D_hat <- lm(D ~ Z, data = data_X1)$fitted.values
data_X1$D_hat <- D_hat
lm(Y ~ D_hat, data = data_X1)

# alternative calculation
cov(data_X1$Y, data_X1$Z_tilde) / cov(data_X1$D, data_X1$Z_tilde)

# Analytical Calculation
# X = 0
mu_1[1] - mu_0[1] + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(a_0[2] + a_1[2]) - dnorm(a_0[2])) / 
  (pnorm(a_0[2]) - pnorm(a_0[2] + a_1[2]))
# X = 1
mu_1[2] - mu_0[2] + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(a_0[2] + a_1[2]) - dnorm(a_0[2])) / 
  (pnorm(a_0[2]) - pnorm(a_0[2] + a_1[2]))

#### LATE in the full sample ####

# Analytical Calculation 
# Z and X are independent -> the weights simplify to complier shares multiplied
# by the sign function.

# complier shares and sign function
compliers <- rep(0, 2)
compliers[1] <- pnorm(a_0[1] + a_1[1]) - pnorm(a_0[1])
compliers[2] <- pnorm(a_0[2] + a_1[2]) - pnorm(a_0[2])
S <- (compliers > 0) - (compliers < 0)
compliers <- abs(compliers)

denom = sum(c(1 / 2 + gamma, 1 / 2 - gamma) * compliers * S)

# conditional LATEs
tau <- rep(0, 2)
tau[1] <- mu_1[1] - mu_0[1] + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(a_0[1]) - dnorm(a_0[1] + a_1[1])) /
  (pnorm(a_0[1] + a_1[1]) - pnorm(a_0[1]))
tau[2] <- mu_1[2] - mu_0[2] + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(a_0[2] + a_1[2]) - dnorm(a_0[2])) / 
  (pnorm(a_0[2]) - pnorm(a_0[2] + a_1[2]))

num <- sum(c(1 / 2 + gamma, 1 / 2 - gamma) * compliers * S * tau)

num / denom  

# estimated
D_hat <- lm(D ~ Z, data = data)$fitted.values
data$D_hat <- D_hat
lm(Y ~ D_hat, data = data)

# alternative method
cov(data$Y, data$Z_tilde) / cov(data$D, data$Z_tilde)

#### LATE with sign corrected data #### 

cor_denom = sum((1 / 2) * compliers)
cor_num <- sum((1 / 2) * compliers * tau)

cor_num / cor_denom

# estimated
data$S <- (data$X == 0) - (data$X == 1)
data$Z_mod <- data$Z_tilde * data$S 
D_hat <- lm(D ~ Z_mod, data = data)$fitted.values
data$D_hat <- D_hat
lm(Y ~ D_hat, data = data)

# alternative method
cov(data$Y, data$Z_mod) / cov(data$D, data$Z_mod)

# Data for plot
plot_points <- list()
plot_points$tau <- tau
plot_points$tsls_limit <- num / denom
plot_points$corrected_limit <- cor_num / cor_denom

#### Checking if we can obtain weights ####

data$ind_0 <- as.numeric(data$X == 0)
data$ind_1 <- as.numeric(data$X == 1)

# check if we can obtain weights easily
w_0 <- cov(data$ind_0 * data$D, data$Z_tilde)
w_1 <- cov(data$ind_1 * data$D, data$Z_tilde) 

lm(D ~ Z, data = data)
(w_0 + w_1) / var(data$Z_tilde)

tau_0 <- lm(Y ~ D_hat, data = data_X0)$coefficients[2]
tau_1 <- lm(Y ~ D_hat, data = data_X1)$coefficients[2]
((tau_0 * w_0) + (tau_1 * w_1)) / (w_0 + w_1)

lm(Y ~ D_hat, data = data)

# Weights -> first stage numerators!

#### Mini Monte Carlo #### 

# Sample size 200, 800, 3200

n <- c(100, 400, 1600, 6400)
sims <- length(n)
B <- 1000

results <- vector(mode = "list", length = sims)

for (i in 1:sims){
  
  sim_results <- list(tsls = rep(0, B),
                      corrected = rep(0, B))
  
  for (b in 1:B){
    
    Z <- rbinom(n[i], 1, 1 / 2)
    X <- rbinom(n[i], 1, 1 / 2 - gamma)
    
    X_ord <- order(X)
    Z <- Z[X_ord]
    X <- X[X_ord]
    
    n_1 <- sum(X)
    n_0 <- n[i] - n_1
    
    # generate unobservables
    U <- rmvnorm(n[i], mean = rep(0, 3), sigma = sigma)
    U <- data.frame(U)
    names(U) <- c("V", "U_0", "U_1")
    
    # generate D
    D <- as.numeric(((1 - X) * (a_0[1] + a_1[1] * Z) + 
                       X * (a_0[2] + a_1[2] * Z)) >= U[, 1])
    
    # generate Y
    Y_1 <- mu_1[X + 1] + U$U_1
    Y_0 <- mu_0[X + 1] + U$U_0
    Y <- D * Y_1 + (1 - D) * Y_0
    
    data <- data.frame(Y = Y, D = D, Z = Z, X = X)
    data$Z_tilde <- data$Z - mean(data$Z)
    
    #### LATE in the full sample ####
    
    # estimated
    D_hat <- lm(D ~ Z, data = data)$fitted.values
    data$D_hat <- D_hat
    sim_results$tsls[b] <- lm(Y ~ D_hat, data = data)$coefficients[2]
    
    # estimated with sign correction
    S <- rep(0, 2)
    temp <- dplyr::filter(data, X == 0)
    S[1] <- (2 * (cor(temp$Z, temp$D) > 0)) - 1
    temp <- dplyr::filter(data, X == 1)
    S[2] <- (2 * (cor(temp$Z, temp$D) > 0)) - 1
    
    data$Z_mod <- data$Z_tilde * ((1 - data$X) * S[1] + data$X * S[2])
    
    D_hat <- lm(D ~ Z_mod, data = data)$fitted.values
    data$D_hat <- D_hat
    sim_results$corrected[b] <- lm(Y ~ D_hat, data = data)$coefficients[2]
    
  }
  
  results[[i]] <- sim_results
  
}

# Prepare data for plotting
x_lims <- rep(0, 2)
x_lims[1] <- quantile(results[[1]]$tsls, 0.05, na.rm = TRUE)
x_lims[2] <- quantile(results[[1]]$tsls, 0.95, na.rm = TRUE)

make_plot <- function(results, plot_points, x_lims){
  
  tsls <- results$tsls[!is.na(results$tsls)]
  tsls_dens <- density(tsls)
  tsls_1 <- quantile(tsls, 0.05)
  tsls_99 <- quantile(tsls, 0.95)
  tsls_data <- data.frame(x = tsls_dens$x, y = tsls_dens$y)
  tsls_data <- dplyr::filter(tsls_data, x > tsls_1 & x < tsls_99)
  
  corrected <- results$corrected
  corrected_dens <- density(corrected)
  corrected_data <- data.frame(x = corrected_dens$x, y = corrected_dens$y)
  
  ggplot() +
    geom_line(aes(x = x, y = y, col = "TSLS"), data = tsls_data, linetype = "dashed") +
    geom_line(aes(x = x, y = y, col = "Corrected"), data = corrected_data) +
    geom_vline(aes(xintercept = plot_points$tsls_limit, col = "TSLS Limit"), 
               linetype = "dotted") + 
    geom_vline(aes(xintercept = plot_points$corrected_limit, col = "Corrected Limit"),
               linetype = "dotted") + 
    xlim(x_lims[1], x_lims[2]) + 
    theme_bw() +
    theme(legend.title = element_blank()) +
    scale_color_manual(values = c("TSLS" = "#F8766D",
                                  "Corrected" = "#619CFF",
                                  "TSLS Limit" = "#F8766D",
                                  "Corrected Limit" = "#619CFF")) + 
    labs(x = "", y = "") + 
    theme(legend.position = c(0.18, 0.8))
  
}

make_plot(results[[1]], plot_points, x_lims)
make_plot(results[[2]], plot_points, x_lims)
make_plot(results[[3]], plot_points, x_lims)
make_plot(results[[4]], plot_points, x_lims)

# Summary statistics
table <- matrix(0, 2 * sims, 7)

get_stats <- function(results){
  
  stats <- rep(0, 7)
  stats[1] <- mean(results, na.rm = TRUE)
  stats[2] <- sd(results, na.rm = TRUE)
  stats[3:7] <- quantile(results, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  return(stats)
  
}

rownames(table) <- rep(n, 2)
table[1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "tsls"))},
           results = results))
table[sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "corrected"))},
           results = results))

xtable(table[1:sims, ], digits = 3)
xtable(table[sims + 1:sims, ], digits = 3)

# tsls <- results$tsls[!is.na(results$tsls)]
# tsls_dens <- density(tsls)
# tsls_1 <- quantile(tsls, 0.05)
# tsls_99 <- quantile(tsls, 0.95)
# tsls_data <- data.frame(x = tsls_dens$x, y = tsls_dens$y)
# tsls_data <- dplyr::filter(tsls_data, x > tsls_1 & x < tsls_99)
# 
# corrected <- results$corrected
# corrected_dens <- density(corrected)
# corrected_data <- data.frame(x = corrected_dens$x, y = corrected_dens$y)
# 
# ggplot() +
#   geom_line(aes(x = x, y = y, col = "TSLS"), data = tsls_data) +
#   geom_line(aes(x = x, y = y, col = "corrected"), data = corrected_data) +
#   geom_vline(aes(xintercept = tau[1], col = "tau_0")) + 
#   geom_vline(aes(xintercept = tau[2], col = "tau_1")) +
#   geom_vline(aes(xintercept = plot_points$tsls_limit, col = "tsls limit")) + 
#   geom_vline(aes(xintercept = plot_points$corrected_limit, col = "corrected limit")) + 
#   theme_bw() +
#   theme(legend.title = element_blank()) +
#   labs(x = "", y = "") + 
#   theme(legend.position = "bottom")


################################ MULTTIVALUED ################################

# X is discrete and takes an even number of values, half positive and half 
# negative
# Defiers at negative X, compliers at positive X

rm(list = ls())
gc()

gamma <- 0.1
X_vals <- 1:8
X_vals <- c(-X_vals, X_vals)
N_x <- length(X_vals)
probs <- c(rep((1 + gamma) / N_x, N_x / 2),rep((1 - gamma) / N_x, N_x / 2)) 

# set first stage parameters
a_0 <- c(qnorm(0.2), qnorm(0.6)) # pnorm(a_0) is the lower propensity score
a_1 <- c(qnorm(0.6) - a_0[1], qnorm(0.2) - a_0[2])
# pnorm(a_0 + a_1) is the higher propensity score

# specify a covariance matrix
# to simplify, set diagonal to one -> we directly specify a correlation matrix
sigma <- diag(1, 3)
# want positive selection -> those with the highest treatment effects are more 
# likely to select treatment.
# V is `resistance to treatment` -> high V means less likely to take treatment
# to generate positive selection, we want V to be positive related to U(0) and
# negatively related to U(1)
# 
sigma[1, 2] <- 0.5
sigma[1, 3] <- -0.5
sigma[, 1] <- sigma[1 ,]

# Epsilon chosen to induce a sign flip
mu_0 <- rep(0, 2)
mu_1 <- rep(0, 2) 
tau_0 <- mu_1[1] - mu_0[1] + (sigma[3, 1] - sigma[2, 1]) *
  (dnorm(a_0[1]) - dnorm(a_1[1] + a_0[1])) / 
  (pnorm(a_1[1] + a_0[1]) - pnorm(a_0[1]))
A <- (1 / 2 - gamma) / (2 * gamma)
# threshold for the TSLS estimate to be zero
tau_0 / A

epsilon <- 0.5
mu_1[2] <- mu_1[2] + epsilon

#### Mini Monte Carlo #### 

n <- c(100, 400, 1600, 6400)
sims <- length(n)
B <- 1000

results <- vector(mode = "list", length = sims)

for (i in 1:sims){
  
  sim_results <- list(tsls = rep(0, B),
                      corrected = rep(0, B))
  
  for (b in 1:B){
    
    Z <- rbinom(n[i], 1, 1 / 2)
    X <- sample(X_vals, n[i], replace = T, prob = probs)
    
    X_ord <- order(X)
    Z <- Z[X_ord]
    X <- X[X_ord]
    
    boot_X_vals <- unique(X)
    boot_N_X <- length(boot_X_vals)
    counts <- sapply(boot_X_vals, 
                     function(data, val){return(sum(data == val))}, 
                     data = X)
    
    X_pos <- (boot_X_vals > 0)
    X_pos <- unlist(lapply(1:boot_N_X, 
                           function(i, X_pos, counts){rep(X_pos[i], counts[i])},
                           X_pos = X_pos,
                           counts = counts))
    
    # generate unobservables
    U <- rmvnorm(n[i], mean = rep(0, 3), sigma = sigma)
    U <- data.frame(U)
    names(U) <- c("V", "U_0", "U_1")
    
    # generate D
    D <- as.numeric(((1 - X_pos) * (a_0[1] + a_1[1] * Z) + 
                       X_pos * (a_0[2] + a_1[2] * Z)) >= U[, 1])
    
    # generate Y
    Y_1 <- mu_1[X_pos + 1] + U$U_1
    Y_0 <- mu_0[X_pos + 1] + U$U_0
    Y <- D * Y_1 + (1 - D) * Y_0
    
    data <- data.frame(Y = Y, D = D, Z = Z, X = X)
    data$Z_tilde <- data$Z - mean(data$Z)
    
    #### LATE in the full sample ####
    
    # estimated
    D_hat <- lm(D ~ Z, data = data)$fitted.values
    data$D_hat <- D_hat
    sim_results$tsls[b] <- lm(Y ~ D_hat, data = data)$coefficients[2]
    
    # estimated with sign correction
    Z_mod <- rep(0, n[i])
    counts <- cumsum(counts)
    counts <- c(0, counts)
    
    for (j in 1:boot_N_X){
      
      temp <- data[(1 + counts[j]):counts[j + 1], ]
      # S <- (2 * (cor(temp$Z, temp$D) > 0)) - 1 
      S <- sign(mean(temp$Z * temp$D) - mean(temp$D) * mean(temp$Z))
      Z_mod[(1 + counts[j]):counts[j + 1]] <- temp$Z_tilde * S
      
    }
    
    data$Z_mod <- Z_mod
    
    D_hat <- lm(D ~ Z_mod, data = data)$fitted.values
    data$D_hat <- D_hat
    sim_results$corrected[b] <- lm(Y ~ D_hat, data = data)$coefficients[2]
    
  }
  
  results[[i]] <- sim_results
  
}

# Prepare data for plotting
x_lims <- rep(0, 2)
x_lims[1] <- quantile(results[[1]]$tsls[!is.na(results[[1]]$tsls)], 0.05)
x_lims[2] <- quantile(results[[1]]$tsls[!is.na(results[[1]]$tsls)], 0.95)
plot_points <- list()

# calculate the limits
compliers <- rep(0, N_x)
X_pos <- X_vals > 0
compliers[!X_pos] <- pnorm(a_0[1] + a_1[1]) - pnorm(a_0[1])
compliers[X_pos] <- pnorm(a_0[2] + a_1[2]) - pnorm(a_0[2])
S <- (compliers > 0) - (compliers < 0)
compliers <- abs(compliers)

denom = sum(probs * compliers * S)

tau <- rep(0, N_x)
tau[!X_pos] <- mu_1[1] - mu_0[1] + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(a_0[1]) - dnorm(a_0[1] + a_1[1])) /
  (pnorm(a_0[1] + a_1[1]) - pnorm(a_0[1]))
tau[X_pos] <- mu_1[2] - mu_0[2] + (sigma[3, 1] - sigma[2, 1]) * 
  (dnorm(a_0[2] + a_1[2]) - dnorm(a_0[2])) / 
  (pnorm(a_0[2]) - pnorm(a_0[2] + a_1[2]))

num <- sum(probs * compliers * S * tau)

plot_points$tsls_limit <- num / denom  

cor_denom = sum((1 / 2) * compliers)
cor_num <- sum((1 / 2) * compliers * tau)

plot_points$corrected_limit <- cor_num / cor_denom

make_plot <- function(results, plot_points, x_lims){
  
  tsls <- results$tsls[!is.na(results$tsls)]
  tsls_dens <- density(tsls)
  tsls_1 <- quantile(tsls, 0.05)
  tsls_99 <- quantile(tsls, 0.95)
  tsls_data <- data.frame(x = tsls_dens$x, y = tsls_dens$y)
  tsls_data <- dplyr::filter(tsls_data, x > tsls_1 & x < tsls_99)
  
  corrected <- results$corrected
  corrected_dens <- density(corrected)
  corrected_data <- data.frame(x = corrected_dens$x, y = corrected_dens$y)
  
  ggplot() +
    geom_line(aes(x = x, y = y, col = "TSLS"), data = tsls_data, linetype = "dashed") +
    geom_line(aes(x = x, y = y, col = "Corrected"), data = corrected_data) +
    geom_vline(aes(xintercept = plot_points$tsls_limit, col = "TSLS Limit"), 
               linetype = "dotted") + 
    geom_vline(aes(xintercept = plot_points$corrected_limit, col = "Corrected Limit"),
               linetype = "dotted") + 
    xlim(x_lims[1], x_lims[2]) + 
    theme_bw() +
    theme(legend.title = element_blank()) +
    scale_color_manual(values = c("TSLS" = "#F8766D",
                                  "Corrected" = "#619CFF",
                                  "TSLS Limit" = "#F8766D",
                                  "Corrected Limit" = "#619CFF")) + 
    labs(x = "", y = "") + 
    theme(legend.position = c(0.18, 0.8))
  
}

make_plot(results[[1]], plot_points, x_lims)
make_plot(results[[2]], plot_points, x_lims)
make_plot(results[[3]], plot_points, x_lims)
make_plot(results[[4]], plot_points, x_lims)

# Summary statistics
table <- matrix(0, 2 * sims, 7)

get_stats <- function(results){
  
  stats <- rep(0, 7)
  stats[1] <- mean(results, na.rm = TRUE)
  stats[2] <- sd(results, na.rm = TRUE)
  stats[3:7] <- quantile(results, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  return(stats)
  
}

rownames(table) <- rep(n, 2)
table[1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "tsls"))},
           results = results))
table[sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "corrected"))},
           results = results))

xtable(table[1:sims, ], digits = 3)
xtable(table[sims + 1:sims, ], digits = 3)

################################ BIGGER EXAMPLE ################################

# Two covariates: X_1 and X_2. Sign varies over their range
rm(list = ls())

r_1 <- 0.5
r_2 <- 0.25

N_X <- 101
X_1 <- seq(0, 1, length.out = N_X)
X_2 <- seq(0, 1, length.out = N_X)

data <- data.frame(X_1 = rep(X_1, each = N_X), 
                   X_2 = rep(X_2, N_X))
data$S <- sign((data$X_1 - r_1) *  (data$X_2 - r_2))

data %>% ggplot() + 
  geom_point(aes(x = X_1, y = X_2, col = as.factor(S))) +
  theme_bw() +
  theme(legend.position = "null")

r_3 <- 0.8
r_4 <- 0.1

data$S <- sign((data$X_1 - r_1) *  (data$X_2 - r_2) * 
                 (data$X_1 - r_3) *  (data$X_2 - r_4))

data %>% ggplot() + 
  geom_point(aes(x = X_1, y = X_2, col = as.factor(S))) +
  theme_bw() +
  labs(x = expression(X[1]), y = expression(X[2])) +
  guides(color = guide_legend(title = "s(X)")) +
  theme(legend.position = "bottom")

r_5 <- 0.7
r_6 <- 0.9 

data$S <- sign((data$X_1 - r_1) *  (data$X_2 - r_2) * 
                 (data$X_1 - r_3) *  (data$X_2 - r_4) *
                 (data$X_1 - r_5) *  (data$X_2 - r_6))

data %>% ggplot() + 
  geom_point(aes(x = X_1, y = X_2, col = as.factor(S))) +
  theme_bw() +
  labs(x = expression(X[1]), y = expression(X[2])) +
  guides(color = guide_legend(title = "s(X)")) +
  scale_color_manual(values = c("1" = "black", "-1" = "grey", "0" = "grey")) + 
  theme(legend.position = "bottom")

# Monte Carlo in this case

# n <- 1000000
# Z <- rbinom(n, 1, 1 / 2)
# X_1 <- runif(n, 0, 1)
# X_2 <- runif(n, 0, 1)
# 
# a_1 <- qnorm(0.6) # pnorm(a_1) is the higher propensity score 
# a_0 <- qnorm(0.45) # pnorm(a_0) is the lower propensity score
# 
# # covariance matrix for the unobservables, let this be common
# A <- matrix(runif(9, -1, 1), 3, 3)
# # scale so that all diagonal elements will be 1
# for (i in 1:3){
#   A[i, ] <- A[i, ] / sqrt(sum(A[i, ] ^ 2))
# }
# sigma <- A %*% t(A)
# 
# # generate unobservables
# U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
# U <- data.frame(U)
# names(U) <- c("V", "U_0", "U_1")
# 
# # generate D
# # Define regions such that D(1) > D(0) and D(0) < D(1)
# # Simple case: D(1) > D(0) if X_1 < 1/2 and vice versa
# X_plus <- X_1 < 1 / 2
# X_minus <- X_plus == FALSE
# D <- as.numeric((Z * (a_1 * X_plus + a_0 * X_minus) +
#                    (1 - Z) * (a_0 * X_plus + a_1 * X_minus)) >= U$V)
# 
# # Plot of this
# plot_data <- data.frame(xmin = c(0, 1 / 2),
#                         xmax = c(1 / 2, 1),
#                         ymin = c(0, 0),
#                         ymax = c(1, 1),
#                         region = c("positive", "negative"))
# 
# plot_data %>% ggplot() + 
#   geom_rect(aes(xmin = xmin, 
#                 xmax = xmax, 
#                 ymin = ymin, 
#                 ymax = ymax,
#                 fill = region)) +
#   theme_bw() + 
#   labs(x = "X1", y = "X2") + 
#   theme(legend.title = element_blank())
# 
# 
# # check D is correct
# mean(D[Z == 1 & X_plus == FALSE]) - pnorm(a_0)
# mean(D[Z == 0 & X_plus == FALSE]) - pnorm(a_1)
# mean(D[Z == 1 & X_plus == TRUE]) - pnorm(a_1)
# mean(D[Z == 0 & X_plus == TRUE]) - pnorm(a_0)
# 
# # plot of covariance evaluated in small tiles
# X_1_disc <- ceiling(data$X_1 * 10)
# X_2_disc <- ceiling(data$X_2 * 10)
# 
# cells <- data.frame(X_1 = rep(1:10, each = 10),
#                     X_2 = rep(1:10, length = 100))
# cells$cov <- rep(0, nrow(cells))
# for (i in 1:nrow(cells)){
#   
#   ind <- 
#     as.logical((X_1_disc == cells$X_1[i]) * (X_2_disc == cells$X_2[i]))
#   cells$cov[i] <- cov(data[ind, ]$D, data[ind, ]$Z)
#   
# }
# cells$sign <- 2 * (cells$cov > 0) - 1
# 
# cells %>% ggplot() + 
#   geom_point(aes(x = X_1 / 10, 
#                  y = X_2 / 10, 
#                  col = as.character(sign))) +
#   theme_bw() + 
#   theme(legend.title = element_blank())
# 
# # generate Y
# Y_1 <- tau + U$U_1
# Y_0 <- - tau + U$U_0
# Y <- D * Y_1 + (1 - D) * Y_0
# 
# #### Estimate LATE in subsamples #### 
# 
# data <- data.frame(Y = Y, D = D, Z = Z, X_1 = X_1, X_2 = X_2)
# data$Z_tilde <- Z - lm(Z ~ X_1 + X_2, data = data)$fitted.values
# 
# p_0 <- pnorm(a_0)
# p_1 <- pnorm(a_1)
# q_0 <- a_0
# q_1 <- a_1
# 
# # Analytical Calculation
# 2 * tau + (sigma[3, 1] - sigma[2, 1]) * 
#   (dnorm(q_0) - dnorm(q_1)) / (p_1 - p_0)
# 
# # X_plus subsample
# data_X_plus <- data[X_plus, ]
# data_X_plus$Z_tilde <- data_X_plus$Z - mean(data_X_plus$Z)
# D_hat <- lm(D ~ Z, data = data_X_plus)$fitted.values
# data_X_plus$D_hat <- D_hat
# lm(Y ~ D_hat, data = data_X_plus)
# 
# # alternative calculation
# cov(data_X_plus$Y, data_X_plus$Z_tilde) / 
#   cov(data_X_plus$D, data_X_plus$Z_tilde)
# 
# # X_minus subsample
# data_X_minus <- data[X_minus, ]
# data_X_minus$Z_tilde <- data_X_minus$Z - mean(data_X_minus$Z)
# D_hat <- lm(D ~ Z, data = data_X_minus)$fitted.values
# data_X_minus$D_hat <- D_hat
# lm(Y ~ D_hat, data = data_X_minus)
# 
# # alternative calculation
# cov(data_X_minus$Y, data_X_minus$Z_tilde) / 
#   cov(data_X_minus$D, data_X_minus$Z_tilde)
# 
# #### LATE in the full sample ####
# 
# # analytical calculation 
# # as Z and X are independent and the propensity score change is
# # constant, the LATE for the full sample equals the LATE in each subsample
# 2 * tau + (sigma[3, 1] - sigma[2, 1]) * 
#   (dnorm(q_0) - dnorm(q_1)) / (p_1 - p_0)
# 
# # estimated
# D_hat <- lm(D ~ Z, data = data)$fitted.values
# data$D_hat <- D_hat
# lm(Y ~ D_hat, data = data)
# 
# # alternative method
# cov(data$Y, data$Z_tilde) / cov(data$D, data$Z_tilde)
# 
# # Calculation as a weighted average (separate Z_tilde)
# n_plus <- sum(X_plus)
# n_minus <- sum(X_minus)
# denom_plus <- cov(data_X_plus$D, data_X_plus$Z_tilde)
# denom_minus <- cov(data_X_minus$D, data_X_minus$Z_tilde)
# denom_all <- cov(data$D, data$Z_tilde)
# theta_plus <- denom_plus / denom_all
# theta_minus <- denom_minus / denom_all
# LATE_plus <- cov(data_X_plus$Y, data_X_plus$Z_tilde) / 
#   cov(data_X_plus$D, data_X_plus$Z_tilde)
# LATE_minus <- cov(data_X_minus$Y, data_X_minus$Z_tilde) / 
#   cov(data_X_minus$D, data_X_minus$Z_tilde)
# (theta_plus * n_plus * LATE_plus + theta_minus * n_minus * LATE_minus) / n
# cov(data$Y, data$Z_tilde) / cov(data$D, data$Z_tilde)
# 
# # common Z_tilde
# temp_plus <- data[X_plus, ]
# temp_minus <- data[X_minus, ]
# denom_plus <- cov(temp_plus$D, temp_plus$Z_tilde)
# denom_minus <- cov(temp_minus$D, temp_minus$Z_tilde)
# denom_all <- cov(data$D, data$Z_tilde)
# theta_plus <- denom_plus / denom_all
# theta_minus <- denom_minus / denom_all
# LATE_plus <- cov(temp_plus$Y, temp_plus$Z_tilde) / 
#   cov(temp_plus$D, temp_plus$Z_tilde)
# LATE_minus <- cov(temp_minus$Y, temp_minus$Z_tilde) / 
#   cov(temp_minus$D, temp_minus$Z_tilde)
# (theta_plus * n_plus * LATE_plus + theta_minus * n_minus * LATE_minus) / n
# cov(data$Y, data$Z_tilde) / cov(data$D, data$Z_tilde)
# 
# # Relabeled data
# data$Z_mod <- data$Z_tilde * (2 * as.numeric(X_plus) - 1)
# cov(data$Y, data$Z_mod) / cov(data$D, data$Z_mod)
# 
# 
# #####################################################################
# 
# # CHECK THE STRUCTURE OF THE MODEL
# 
# rm(list=ls())
# gc()
# 
# 
# 
# 
# 
# X <- sample(1:K, N, replace = T)
# X <- sort(X)
# threshold <- floor(K / 2)
# X_minus <- X <= threshold
# 
# a_1 <- 0.1 # complier share
# a_0 <- runif(K, min = 0.3, max = 0.7) # reference share
# U_D <- runif(N, 0 , 1)
# 
# D <- rep(0, N)
# index <- 0
# for (k in 1:K){
#   
#   N_k <- sum(X == k)
#   D[index + 1:N_k] <- as.numeric(
#     (a_0[k] + a_1 * Z[index + 1:N_k] * (1 - 2 * (k <= threshold))) >= 
#       U_D[index + 1:N_k])
#   
#   index <- index + N_k
#   
# }
# 
# # cell-by-cell regressions
# data <- data.frame(D = D,
#                    Z = Z,
#                    X = X)
# estimates <- matrix(0, K, 2)
#  
# for (k in 1:K){
#   
#   temp_data <- dplyr::filter(data, X == k)
#   estimates[k, ] <- lm("D ~ Z", data = temp_data)$coefficients
#   
# }
# 
# estimates <- data.frame(estimates)
# colnames(estimates) <- c("a_0", "a_1")
# 
# plot(estimates[, 2])
# lines(x = 1:threshold, y = rep(-a_1, threshold))
# lines(x = (threshold + 1):N, y = rep(a_1, N - threshold))
# 
# plot(a_0, estimates[, 1])
# 
# # Add in the outcome
# 

# 
# ######################## TSLS WITH MULTIPLE INSTRUMENTS ########################
# 
# #### CHECKING MEAN INDEPENDENCE ####
# 
# # Maintain that Z is independent of all unobservables (conditional on X). Say X 
# # is endogenous. What about XZ? 
# 
# rm(list=ls())
# gc()
# 
# # X independent of Z
# n <- 1000000
# X <- rnorm(n, mean = 1)
# Z <- rbinom(n, 1, 1/3)
# XZ <- X * Z
# 
# data <- data.frame(X = X, Z = Z, XZ = XZ)
# 
# data$Z_tilde <- lm(Z ~ X, data = data)$fitted.values
# data$XZ_tilde <- lm(XZ ~ X, data = data)$fitted.values
# 
# temp <- data$X * data$Z_tilde
# cor(temp, data$XZ_tilde)
# temp - data$XZ_tilde
# mean(abs(temp - data$XZ_tilde))
# 
# # E[Z|X] is a linear function of X
# X <- runif(n, 0, 1)
# a_0 <- 1 / 3
# a_1 <- 2 / 3
# Z <- rbinom(n, 1, a_0 + a_1 * X)
# XZ <- X * Z
# 
# data <- data.frame(X = X, X_square = X ^ 2, Z = Z, XZ = XZ)
# 
# data$Z_tilde <- lm(Z ~ X, data = data)$fitted.values
# data$XZ_tilde <- lm(XZ ~ X + X_square, data = data)$fitted.values
# 
# 
# temp <- data$X * data$Z_tilde
# cor(temp, data$XZ_tilde)
# temp - data$XZ_tilde
# mean(abs(temp - data$XZ_tilde))
