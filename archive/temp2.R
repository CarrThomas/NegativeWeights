rm(list = ls())
gc()

# Check bias calculations in the case with covariates

n <- 100000 # sample size

# Specify standard deviations and correlations for x'd and x'a, generate 
# d and a to satisfy these expectations
path <- getwd()
parameters <- readRDS(paste0(path, "/Monte_Carlo/Parameters1.rds"))
for (i in 1:length(parameters)){
  assign(names(parameters)[i], parameters[[i]])
}

g_0 <- a_0 + d_0
sd_g <- sqrt(sd_a ^ 2 + 2 * rho_ad * sd_a * sd_d + sd_d ^ 2)
rho_dg <- (sd_a * rho_ad + sd_d) / sd_g

k <- 4

vals <- matrix(0, 8, 5)

for (i in 1:5){
  
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
  
  X <- matrix(rnorm(n * k), n, k)
  Xa <- X %*% as.matrix(a, k, 1)
  Xd <- X %*% as.matrix(d, k, 1)
  Xg <- Xa + Xd
  
  ind <- ((d_0 + Xd) > 0)
  
  # E[1(d_0 + x'd > 0)phi(a_0 + x'a)]
  vals[1, i] <- mean(ind * dnorm(a_0 + Xa))
  vals[2, i] <- dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
    (1 - pnorm(-d_0,
               mean = -a_0 * sd_d * sd_a * rho_ad / (1 + sd_a ^ 2),
               sd = sd_d * sqrt(1 - (rho_ad ^ 2) * (sd_a ^ 2) / (1 + sd_a ^ 2))))
  
  # kernel <- function(xd, sd_a, sd_d, rho_ad, a_0){
  #   return(dnorm(a_0, 
  #         mean = -xd * rho_ad * sd_a / sd_d, 
  #         sd = sqrt(1 + sd_a ^ 2 * (1 - rho_ad ^ 2))) *
  #     dnorm(xd, sd = sd_d))
  # }
  # integrate(kernel, lower = - d_0, upper = Inf,
  #           sd_a = sd_a, sd_d = sd_d, rho_ad = rho_ad, a_0 = a_0)
  
  # kernel <- function(xa, xd, rho_ad, sd_d, sd_a, a_0){
  #   
  #   return(dnorm(xa, mean = -a_0) *
  #            dnorm(xa, 
  #                  mean = rho_ad * sd_a * xd / sd_d, 
  #                  sd = sqrt(sd_a ^ 2 * (1 - rho_ad ^ 2))))
  # }
  # 
  # integrate(kernel, lower = -Inf, upper = Inf,
  #           rho_ad = rho_ad, sd_d = sd_d, sd_a = sd_a, a_0 = a_0,
  #           xd = 2)
  # dnorm(a_0, 
  #       mean = -2 * rho_ad * sd_a / sd_d, 
  #       sd = sqrt(1 + sd_a ^ 2 * (1 - rho_ad ^ 2)))
  
  
  # E[1(d_0 + x'd > 0)phi(g_0 + x'g)]
  vals[3, i] <- mean(ind * dnorm(g_0 + Xg))
  vals[4, i] <- dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
    (1 - pnorm(-d_0,
               mean = -g_0 * sd_d * sd_g * rho_dg / (1 + sd_g ^ 2),
               sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / (1 + sd_g ^ 2))))
  
  # E[1(d_0 + x'd < 0)phi(a_0 + x'a)]
  vals[5, i] <- mean(as.numeric(!ind) * dnorm(a_0 + Xa))
  vals[6, i] <- dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
    pnorm(-d_0,
          mean = -a_0 * sd_d * sd_a * rho_ad / (1 + sd_a ^ 2),
          sd = sd_d * sqrt(1 - rho_ad ^ 2 * sd_a ^ 2 / (1 + sd_a ^ 2)))
  
  # E[1(d_0 + x'd < 0)phi(g_0 + x'g)]
  vals[7, i] <- mean(as.numeric(!ind) * dnorm(g_0 + Xg))
  vals[8, i] <- dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
    pnorm(-d_0,
          mean = -g_0 * sd_d * sd_g * rho_dg / (1 + sd_g ^ 2),
          sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / (1 + sd_g ^ 2)))
}

vals

################################################################################

# Figuring out the selection model: what can we actually estimate from a given
# data set?
# Binary Z, Z is a valid instrument
# Binary X, X satisfies neither exogeneity or monotonicity
# Resistance to treatment V, correlated with X, indepedent of Z conditional on X

rm(list=ls())
gc()

pX <- 1 / 2
pZ <- 1 / 2

n <- 10000
X <- rbinom(n, 1, pX)
Z <- rbinom(n, 1, pZ)

# "Resistance to treatment" is a vector with two elements. 
# Each element must be mean 0 variance 1 unconditionally
# The conditional distribution can differ

temp <- unlist(parameters[!names(parameters) %in% "sigma"])

table <- data.frame(par = temp)

xtable(table, digits = 3)


################################################################################

# EXPERIMENT WITH A PRINCIPAL COMPONENTS PLOT
temp <- prcomp(data[unlist(varnames[1:i])], retx = TRUE)

plot(temp$x[, 1:2], pch = 3)

plot_data <- cbind(PC, probit_s_hat, logit_s_hat, lpm_s_hat)
plot_data <- data.frame(plot_data)
names(plot_data)[length(controls) + 1:3] <- c("probit", "logit", "lpm")

# reduce to unique values of the first two principal values, take an average
combined <- apply(plot_data[c("PC1", "PC2")], 2, paste, sep = "")[1]
apply()

plot_data %>% ggplot() +
  geom_jitter(aes(x = PC1, y = PC2, col = as.factor(probit)))

plot_data %>% ggplot() +
  geom_jitter(aes(x = PC1, y = PC2, col = as.factor(logit)))

plot_data %>% ggplot() +
  geom_jitter(aes(x = PC1, y = PC2, col = as.factor(lpm)))

plot_data$ind_1 <- plot_data$probit - plot_data$logit
plot_data$ind_2 <- plot_data$probit - plot_data$lpm

boot_mrc_theta
sapply(boot_mrc_theta, length)
n_vars <- length(boot_mrc_theta[[1]])
store <- matrix(unlist(boot_mrc_theta), B, n_vars, byrow = T)
colnames(store) <- names(boot_mrc_theta[[1]])
store <- data.frame(store)
Sigma <- cov(store)


temp <- try(1 - pchisq(matrix(mrc_est[R], 1, r) %*%
                 solve(Sigma) %*%
                 matrix(mrc_est[R], r, 1),
               r))

store <- t(apply(store, 1, function(x)(return(x / x[1]))))

################################################################################

rm(list=ls())
gc()

path <- getwd()
data <- read.csv(file = paste0(path,"/usa_00001.csv"))

S <- 1000
n <- 10000
pZ <- 1 / 2
pD <- 1 / 2
vals <- rep(0, S)

for (s in 1:S){
  
  Z <- rbinom(n, 1, pZ)
  D <- rbinom(n, 1, pD)
  Z_tilde <- Z - pZ
  
  vals[s] <- ((sum(D * Z_tilde)) ^ 2) / n
  
}

mean(vals)

################################################################################

K
vals <- 0:K
vals <- vals[vals != K / 2]
for (k in vals){
  
  temp <- dplyr::filter(data, X == k)
  print(sum(temp$D * (1 - temp$Z)) / sum(1 - temp$Z))
  print(sum(temp$D * temp$Z) / sum(temp$Z))
  
}

n <- seq(2, 100, by = 1)
trimming_1 <- n ^ (- 1 / 4)
trimming_2 <- n ^ (- 1 / 4) * (log(n) ^ - 1)
plot(n, trimming_1, type = "l", ylim = c(0, 1))
lines(n, trimming_2, lty = 2)

################### SUMMARY STATS FOR THE NON-CONVERGENT CASE ##################

store <- data
sum((data$D * data$s_star * data$Z_tilde) ^ 2)
sum((data$D * data$s_hat * data$Z_tilde) ^ 2)


############################## NON-CONVERGENT CASE #############################

rm(list = ls())
gc()

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

# D = 1{X(Z * delta)  + V \geq 0} 
# V ~ U(-1,1)
# X ~ U(-1,1) with probability (1 - pi), X = 0 with probability pi
delta <- 0.4

B <- 0 # ATE

# Randomly generate a covariance matrix for the outcome unobservables
A <- matrix(runif(4, -1, 1), 2, 2)
sigma <- A %*% t(A)
# Scale to be a correlation matrix
for (k in 1:2){
  
  sd_k <- sqrt(sigma[k, k])
  sigma[, k] <- sigma[, k] / sd_k
  sigma[k, ] <- sigma[k, ] / sd_k
  
}

# Probability that Z = 1
pZ <- 1 / 2

# generate data
n <- c(500, 1000, 2000)

runs <- 1000 # number of draws for each sample size

tau_oracle <- matrix(0, length(n), runs)
tau_feasible <- matrix(0, length(n), runs)
s_check <- matrix(0, length(n), runs)
coef_check <- vector(mode = "list", length = length(n))
break_flag <- 0

for (i in 1:length(n)){
  
  coef_check[[i]] <- vector(mode = "list", length = runs)
  
  for (j in 1:runs){
    
    X <- runif(n[i], -1, 1)
    V <- runif(n[i], -1, 1)
    U <- rmvnorm(n[i], mean = rep(0, nrow(sigma)), sigma = sigma)
    U <- data.frame(U)
    names(U) <- c("U1", "U0")
    Z <- rbinom(n[i], 1, pZ)
    
    D <- (delta * X * Z >= V)
    Y <- (1 - D) * U$U0 + D *(B + U$U1)
    
    data <- data.frame(D = D, X = X, Z = Z, XZ = X * Z, Y = Y)
    
    Z_hat <- lm(Z ~ X, data = data)$fitted.values
    data$Z_tilde <- Z - Z_hat
    
    data$s_star <- 2 * as.numeric(X >= 0) - 1
    
    p_score <- lm(D ~ X + Z + XZ, data = data)
    coef_check[[i]][[j]] <- p_score$coefficients
    Z_coef <- p_score$coefficients[names(p_score$coefficients) %in% "Z"]
    XZ_coef <- p_score$coefficients[names(p_score$coefficients) %in% "XZ"]
    data$s_hat <- 2 * as.numeric(Z_coef + XZ_coef * X >= 0) - 1
    
    s_check[i, j] <- cor(data$s_star, data$s_hat)
    
    tau_oracle[i, j] <- sum(data$Y * data$s_star * data$Z_tilde) / 
      sum(data$D * data$s_star * data$Z_tilde)
    tau_feasible[i, j] <- sum(data$Y * data$s_hat * data$Z_tilde) / 
      sum(data$D * data$s_hat * data$Z_tilde)
    
    if (abs(tau_feasible[i, j]) > 10){
      break_flag <- 1
      break
    }

  }
  
  if (break_flag == 1){
    
    break
    
  }
}

############################## SUMMARY STATISTICS ##############################

# Returns the mean, variance, MSE, median, 0.25 and 0.75 quartiles of the 
# sampling distribution of the oracle and feasible estimators. 

feasible_sum_stats <- matrix(0, length(n), 6)
oracle_sum_stats <- matrix(0, length(n), 6)

for (i in 1:length(n)){
  
  feasible_sum_stats[i, ] <- get_stats(results = tau_feasible[i, ], B = B)
  oracle_sum_stats[i, ] <- get_stats(results = tau_oracle[i, ], B = B)
  
}

rownames(feasible_sum_stats) <- n
rownames(oracle_sum_stats) <- n

xtable(feasible_sum_stats, digits = 3)
xtable(oracle_sum_stats, digits = 3)

######################## PLOT OF SAMPLING DISTRIBUTIONS ########################

# Prepare data for plotting

# convert to sqrt(n) * (tau - tau^*)
tau_oracle <- tau_oracle - B
tau_feasible <- tau_feasible - B
for (i in 1:nrow(tau_oracle)){
  
  tau_oracle[i, ] <- sqrt(n[i]) * tau_oracle[i, ]
  tau_feasible[i, ] <- sqrt(n[i]) * tau_feasible[i, ]
  
}

x_lims <- rep(0, 2)
x_lims[1] <- quantile(c(tau_feasible, tau_oracle), 0.01, na.rm = TRUE)
x_lims[2] <- quantile(c(tau_feasible, tau_oracle), 0.99, na.rm = TRUE)

make_plot(tau_oracle = tau_oracle[1, ],
          tau_feasible = tau_feasible[1, ],
          x_lims = x_lims,
          B = B)
make_plot(tau_oracle = tau_oracle[2, ],
          tau_feasible = tau_feasible[2, ],
          x_lims = x_lims,
          B = B)
make_plot(tau_oracle = tau_oracle[3, ],
          tau_feasible = tau_feasible[3, ],
          x_lims = x_lims,
          B = B)
make_plot(tau_oracle = tau_oracle[4, ],
          tau_feasible = tau_feasible[4, ],
          x_lims = x_lims,
          B = B)
