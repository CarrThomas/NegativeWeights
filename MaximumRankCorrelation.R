######################## IMPLEMENTING THE MRC ESTIMATOR ########################

# The maximum rank correlation was proposed in Han (1986). The idea is to choose 
# coefficients to maximise the rank correlation between a linear index and the 
# outcome. 

# Assume the model
#
# D = 1{\Phi(X'a + Z * X'd) >  U}
#
# where U is an iid uniform (0, 1) random variable.
#
# Intuitively, if we have observations i and j with 
# X_i'a + Z_i * X_i'd > X_j'a + Z_j * X_j'd
# then we would expect that D_i > D_j on average. Hence, we can estimate a and
# d by maximising the rank correlation between D and the linear index function 
# X'a + Z * X'd. We can do this by maximising
#
# sum_{i \neq j} (1{D_i > D_j}1{X_i'a + Z_i * X_i'd > X_j'a + Z_j * X_j'd}
#                   + 1{D_j > D_i}1{X_j'a + Z_j * X_j'd > X_i'a + Z_i * X_i'd})
#
# where the sum is over pairs (i, j) where i \neq j. If the sample size is n,
# then it can be written as the sum over i = 1:(n - 1) and j = (i + 1):n

################################ FRONT MATTER ##################################

rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(mvtnorm)
library(xtable)

################################ SIMPLE EXAMPLE ################################

# Adopt the MTE assumptions for simplicity.
# Let X consist of two normally distributed variables. 
# X is independent of Z. 

# Specify the selection equation as
# D = 1(a_0 + a_1 * X_1 + a_2 * X_2 + Z(d_0 + d_1 * X_1 + d_2 * X_2) > V)
# Here V is resistance to treatment, so is negatively related to D. Z has a 
# positive effect on D (conditional on X) if d_0 + d_1 * X_1 + d_2 * X_2 > 0. 

n <- 100000 # sample size
k <- 2
pZ <- 1 / 2 # probabilty Z = 1

# Unobservables (V, U_1, U_0) are jointly normally distributed.
# Specify the covariance matrix so there is a positive selection effect -> V 
# is negatively correlated with U(1) and positively correlated with U(0)
sigma <- diag(1, 3)
sigma[1, 2] <- 0.5
sigma[1, 3] <- -0.5
sigma[, 1] <- sigma[1 ,]

Z <- rbinom(n, 1, pZ)
X <- matrix(rnorm(n * k), n, k)
U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
colnames(U) <- c("V", "U1", "U0")

# Propensity score function
d <- runif(k, -1, 1)
a <- runif(k, -1, 1)
a_0 <- runif(1, 0, 1)
d_0 <- runif(1, 0 , 1)
g_0 <- a_0 + d_0
g <- a + d

index <- a_0 + X %*% as.matrix(a, k, 1) + Z * (d_0 + X %*% as.matrix(d, k, 1))
P <- pnorm(index)
D <- as.numeric(index + U[, 1] > 0) 
ind <- (d_0 + X %*% as.matrix(d, k, 1) > 0)
P_x0 <- pnorm(a_0 + X %*% as.matrix(a, k, 1))
P_x1 <- pnorm(g_0 + X %*% as.matrix(g, k, 1))

# Covariance matrix for (X'a, X'd, X'g)
sd_a <- sqrt(sum(a ^ 2))
sd_d <- sqrt(sum(d ^ 2))
sd_g <- sqrt(sum(g ^ 2))
rho_ad <- sum(a * d) / (sd_a * sd_d)
rho_ag <- sum(a * g) / (sd_a * sd_g)
rho_dg <- sum(g * d) / (sd_g * sd_d)

############### SHARES OF SWITCHERS, ALWAYS TAKERS, NEVER TAKERS ###############

# Proportion where Z has a postive effect on the propensity score
1 - pnorm(-d_0, mean = 0, sd = sd_d)
sum((P_x1 - P_x0) > 0) / n

#### Plot of the Z effect positive and negative regions ####

x_1 <- c(-3, 3)
x_2 <- -(d_0 + d[1] * x_1) / d[2]

plot(x_1, x_2, type = "l")

ind_1 <- (d[1] / d[2]) > 0

# region "below" the line
if (ind_1 == TRUE){
  poly_1 <- data.frame(x_1 = c(-3, 3, -3),
                       x_2 = c(min(x_2), min(x_2), max(x_2)))
}else{
  poly_1 <- data.frame(x_1 = c(-3, 3, 3),
                       x_2 = c(min(x_2), max(x_2), min(x_2)))
}

# region "above" the line
if (ind_1 == TRUE){
  poly_2 <- data.frame(x_1 = c(-3, 3, 3),
                       x_2 = c(max(x_2), max(x_2), min(x_2)))
}else{
  poly_2 <- data.frame(x_1 = c(-3, -3, 3),
                       x_2 = c(min(x_2), max(x_2), max(x_2)))
}

# Plots of positive and negative regions
# If d[2] > 0 then Z has a positive effect in the region above the line (poly 2) 
# If d[2] < 0 then Z has a positive effect in the region below the line (poly 1)  

if (d[2] > 0){
  positive <- poly_2
  negative <- poly_1
}else{
  positive <- poly_1
  negative <- poly_2
}

ggplot() + 
  geom_polygon(aes(x = x_1, 
                   y = x_2, 
                   fill = "positive"), 
               data = positive) +
  geom_polygon(aes(x = x_1, 
                   y = x_2, 
                   fill = "negative"),
               data = negative) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  theme_bw() +
  theme(legend.title = element_blank())

#### Region where Z has a positive effect on the propensity score

# Always Takers
kernel <- function(Xd, a_0, sd_d, sd_a, rho_ad){
  
  sd_temp <- sqrt(1 + (sd_a ^ 2) * (1 - rho_ad ^ 2))
  bound <- a_0 + rho_ad * sd_a * Xd / sd_d 
  
  return(pnorm(q = bound, mean = 0, sd = sd_temp) * dnorm(x = Xd, sd = sd_d))
  
}

always_takers_1 <- integrate(kernel, lower = - d_0, upper = Inf,
                             a_0 = a_0,
                             sd_d = sd_d,
                             sd_a = sd_a,
                             rho_ad = rho_ad)$value

# check against the sample (divide by n as we are not conditioning on being in 
# this region)
always_takers_1
sum(P_x0[ind]) / n

# Never takers
kernel <- function(Xd, g_0, sd_d, sd_g, rho_dg){
  
  sd_temp <- sqrt(1 + (sd_g ^ 2) * (1 - rho_dg ^ 2))
  bound <- g_0 + rho_dg * sd_g * Xd / sd_d 
  
  return(pnorm(q = bound, mean = 0, sd = sd_temp) * dnorm(x = Xd, sd = sd_d))
  
}

never_takers_1 <- 1 - pnorm(- d_0 / sd_d) -
  integrate(kernel, lower = - d_0, upper = Inf,
            g_0 = g_0,
            sd_d = sd_d,
            sd_g = sd_g,
            rho_dg = rho_dg)$value

# check against the sample (divide by n as we are not conditioning on being in 
# this region)
never_takers_1
sum((1 - P_x1[ind])) / n

# switchers in this region (compliers)
compliers <- 1 - pnorm(q = -d_0, sd = sd_d) - always_takers_1 - never_takers_1
compliers

#### Region where Z has a negative effect on the propensity score

# Always Takers
kernel <- function(Xd, g_0, sd_d, sd_g, rho_dg){
  
  sd_temp <- sqrt(1 + (sd_g ^ 2) * (1 - rho_dg ^ 2))
  bound <- g_0 + rho_dg * sd_g * Xd / sd_d 
  
  return(pnorm(q = bound, mean = 0, sd = sd_temp) * dnorm(x = Xd, sd = sd_d))
  
}

always_takers_2 <- integrate(kernel, lower = - Inf, upper = - d_0,
                             g_0 = g_0,
                             sd_d = sd_d,
                             sd_g = sd_g,
                             rho_dg = rho_dg)$value

# check against the sample (divide by n as we are not conditioning on being in 
# this region)
always_takers_2
sum(P_x1[!ind]) / n

# Never takers
kernel <- function(Xd, a_0, sd_d, sd_a, rho_ad){
  
  sd_temp <- sqrt(1 + (sd_a ^ 2) * (1 - rho_ad ^ 2))
  bound <- a_0 + rho_ad * sd_a * Xd / sd_d 
  
  return(pnorm(q = bound, mean = 0, sd = sd_temp) * dnorm(x = Xd, sd = sd_d))
  
}

never_takers_2 <- pnorm(- d_0 / sd_d) -
  integrate(kernel, lower = - Inf, upper = - d_0,
            a_0 = a_0,
            sd_d = sd_d,
            sd_a = sd_a,
            rho_ad = rho_ad)$value

# check against the sample (divide by n as we are not conditioning on being in 
# this region)
never_takers_2
sum((1 - P_x0[!ind])) / n

# switchers in this region (defiers)
defiers <- pnorm(q = -d_0, sd = sd_d) - always_takers_2 - never_takers_2
defiers

#### TSLS First Stage

# Value should be equal to compliers - defiers
data <- data.frame(X)
names(data) <- paste("X", 1:k, sep = "_")
X_names <- names(data)
data$D <- D
data$Z <- Z
lm(as.formula(paste("D ~ Z", paste(X_names, collapse = " + "), sep = " + ")),
   data = data)
compliers - defiers

#### Distribution of `switchers'

# Population distribution
kernel <- function(Xa, sd_a, sd_g, rho_ag, a_0, g_0, Delta){
  
  bound <- qnorm(Delta + pnorm(a_0 + Xa)) - g_0
  temp_mean <- sd_g * rho_ag * Xa / sd_a
  temp_sd <- sd_g * sqrt(1 - rho_ag ^ 2)
  
  return(dnorm(x = Xa, sd = sd_a) * pnorm(q = bound, mean = temp_mean, sd = temp_sd))
  
}

N <- 101
Delta <- seq(-1, 1, length.out = N)
switchers <- rep(0, N)

for (i in 2:(N - 1)){
  
  if (Delta[i] < 0){
    lower <- qnorm(-Delta[i]) - a_0
  }else{
    lower <- -Inf
  }
  
  if (Delta[i] > 0){
    upper <- qnorm(1 - Delta[i]) - a_0
  }else{
    upper <- Inf
  }
  
  switchers[i] <- integrate(kernel, lower = lower, upper = upper,
                            sd_a = sd_a,
                            sd_g = sd_g, 
                            rho_ag = rho_ag,
                            a_0 = a_0,
                            g_0 = g_0,
                            Delta = Delta[i])$value
  
  switchers[i] <- switchers[i] + 1 - pnorm(q = upper, mean = 0, sd = sd_a)
  
}
switchers[N] <- 1

# Sample distribution
Delta_p <- P_x1 - P_x0
Delta_p <- sort(Delta_p)

eps <- Delta[2] - Delta[1]
Delta_p_bin <- ceiling((Delta_p + 1) / eps)

Delta_p_dens <- tabulate(Delta_p_bin, N) / n
Delta_p_dist <- cumsum(Delta_p_dens)

# Compare the two
plot_data <- data.frame(Delta = Delta, pop = switchers, data = Delta_p_dist)

plot_data %>% ggplot() + 
  geom_line(aes(x = Delta, y = pop, col = "analytical")) +
  geom_line(aes(x = Delta, y = data, col = "data")) +
  theme_bw() + 
  theme(legend.position = "bottom")

#### Distribution of the Propensity score

N <- 1001
p <- seq(0, 1, length.out = N)
q <- qnorm(p)

# Analytical

P_dist <-  data.frame("p" = p)
P_dens <-  data.frame("p" = p)

P_dist$Z1 <- sapply(q, pnorm, mean = g_0, sd = sd_g)
P_dist$Z0 <- sapply(q, pnorm, mean = a_0, sd = sd_a)
P_dist$uncon <- pZ * P_dist$Z1 + (1 - pZ) * P_dist$Z0

P_dens$Z1 <- sapply(q, dnorm, mean = g_0, sd = sd_g) / dnorm(q)
P_dens$Z0 <- sapply(q, dnorm, mean = a_0, sd = sd_a) / dnorm(q)
P_dens$uncon <- pZ * P_dens$Z1 + (1 - pZ) * P_dens$Z0

# Sample

data_dist <- data.frame("p" = p)
data_dens <- data.frame("p" = p)
eps <- p[2] - p[1]

P_x1 <- sort(P_x1)
P_x1_bin <- ceiling(P_x1 / eps)
data_dist$Z1 <- cumsum(tabulate(P_x1_bin, N)) / n

P_x0 <- sort(P_x0)
P_x0_bin <- ceiling(P_x0 / eps)
data_dist$Z0 <- cumsum(tabulate(P_x0_bin, N)) / n

data_dist$uncon <- (sum(Z) / n) * data_dist$Z1 +
  (1 - sum(Z) / n) * data_dist$Z0

# Compare distributions

# Z = 1
ggplot() + 
  geom_line(aes(x = p, y = Z1, col = "analytical"), data = P_dist) +
  geom_line(aes(x = p, y = Z1, col = "data"), data = data_dist) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Z = 0
ggplot() + 
  geom_line(aes(x = p, y = Z0, col = "analytical"), data = P_dist) +
  geom_line(aes(x = p, y = Z0, col = "data"), data = data_dist) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank())

# unconditional 
ggplot() + 
  geom_line(aes(x = p, y = uncon, col = "analytical"), data = P_dist) +
  geom_line(aes(x = p, y = uncon, col = "data"), data = data_dist) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Densities
P_scores <- data.frame(Z1 = P_x1, Z0 = P_x0)

ggplot() +
  geom_line(aes(x = p, y = Z1, col = "analytical"), data = P_dens) +
  geom_density(aes(x = Z1, col = "data"), data = P_scores) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank())

ggplot() +
  geom_line(aes(x = p, y = Z0, col = "analytical"), data = P_dens) +
  geom_density(aes(x = Z0, col = "data"), data = P_scores) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank())

###################### MAXIMUM RANK CORRELATION ESTIMATOR ######################

# The true model is 
# p(X,Z) = Phi(a_0 + X'a + Z(d_0 + X'd))
# Let theta = (a, d_0, d),
# a_0 is dropped as when looking at ranks, the constant is not identified
# The first element of a (i.e. the coefficient for X_1) is normalised to 1 

# Efficient implementation of the objective function:
S_theta <- function(theta_hat, data, varnames){
  
  # index <- as.matrix(data[varnames]) %*% c(1, theta_hat)
  theta_hat <- theta_hat / sqrt(sum(theta_hat ^ 2))
  index <- as.matrix(data[varnames]) %*% theta_hat
  D_sort <- data$D[order(index)] 
  tally <- cumsum(1 - D_sort)
  val <- sum(tally[as.logical(D_sort)])
  
  return(val)
  
}

old_S_theta <- function(theta_hat, data, varnames){
  
  index <- as.matrix(data[varnames]) %*% theta_hat
  D_sort <- data$D[order(index)] 
  tally <- cumsum(1 - D_sort)
  val <- sum(tally[as.logical(D_sort)])
  
  return(val)
  
}

# Implement with the polytype method
n_1 <- sum(data$D)
n_0 <- sum(1 - data$D)

# Generate interactions
interactions <- apply(X, 2, function(X, Z){return(Z * X)}, Z = Z)
colnames(interactions) <- paste("Z_X", 1:k, sep = "_")
data <- cbind(data, interactions)

# Test of the optimisation function
theta_guess <- runif(2 * k + 1, -1, 1)

S_theta(theta_guess, data, varnames)
old_S_theta(theta_guess, data, varnames)

# Get pilot parameter estimates
varnames <- c(paste("X", 1:k, sep = "_"), "Z", paste("Z_X", 1:k, sep = "_"))
probit_est <- glm(as.formula(paste0("D ~ ", paste(varnames, collapse = " + "))),
                  family = binomial(link = "probit"), 
                  data = data)$coefficients
theta_probit <- probit_est[names(probit_est) %in% varnames]

old_S_theta(theta_probit, data, varnames)
S_theta(theta_probit, data, varnames)

# Solution via the Nelder-Mead-Head simplex method
optim(par = theta_probit, 
      fn = S_theta, 
      method = "Nelder-Mead", 
      data = data, 
      varnames = varnames, 
      control = list("fnscale" = -1))

############################## CHOOSING PARAMETERS #############################

# Choose parameter values to attain certain targets 
# Specify the propensity score as
# p(X,Z) = \Phi(a_0 + X'a + Z(d_0 + X'd))
# With X ~ N(0, I), the parameters we need values for are a_0, d_0, and the 
# covariance matrix of X'a and X'd (sd_a, sd_d and rho_ad)

rm(list = ls())
gc()

#### Functions ####

kernel <- function(Xd, t_0, sd_d, sd_t, rho_td){
  
  # function to be integrated
  # consists of the density for Xd multiplied by either
  # 1) Pr(Xa + V < - a_0|Xd) -> set t_0 = a_0, sd_t = sd_a, rho_td = rho_ad
  # 2) Pr(Xg + V < - g_0|Xd) -> set t_0 = g_0, sd_t = sd_g, rho_td = rho_gd
  sd_temp <- sqrt(1 + (sd_t ^ 2) * (1 - rho_td ^ 2))
  bound <- t_0 + rho_td * sd_t * Xd / sd_d 
  
  return(pnorm(q = bound, mean = 0, sd = sd_temp) * dnorm(x = Xd, sd = sd_d))
  
}

# Calculates shares of always takers, never takers, defiers and compliers
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

#### Figuring Out Stuff ####

a_0 <- qnorm(0.45)
sd_a <- 3
sd_d <- 3
rho_ad <- -0.9

# # d_0 is chosen so that a_0 + d_0 = - a_0
# d_0 <- - 2 * a_0

# Calculate derived parameters 
g_0 <- a_0 + d_0
sd_g <- sqrt(sd_a ^ 2 + 2 * sd_a * sd_d * rho_ad + sd_d ^ 2)
rho_dg <- (sd_d + sd_a * rho_ad) / sd_g

shares <- get_shares(a_0, d_0, g_0, sd_a, sd_d, sd_g, rho_ad, rho_dg)
shares

# Generate a and d subject to these constraints
pZ <- 1 / 2
k <- 2

# construct orthonormal vectors
e_1 <- runif(k, -1, 1)
e_1 <- e_1 / sqrt(sum(e_1 ^ 2))
u <- runif(k, -1, 1)
# create a vector which is orthogonal to e_1
e_2 <- u - (sum(e_1 * u)) * e_1 / (sum(e_1 ^ 2))
e_2 <- e_2 / sqrt(sum(e_2 ^ 2))

# a is a scaled version of e_1
a <- sd_a * e_1

# d is a combination of e_1 and e_2
d <- sd_d * rho_ad * e_1 + sd_d * sqrt(1 - rho_ad ^ 2) * e_2

sd_a - sqrt(sum(a ^ 2))
sd_d - sqrt(sum(d ^ 2))
sd_a * sd_d * rho_ad - sum(a * d)

g <- a + d

# Generate a sample to check calculations
n <- 100000
X <- matrix(rnorm(n * k), n, k)

P_x0 <- pnorm(a_0 + X %*% as.matrix(a, k, 1))
P_x1 <- pnorm(g_0 + X %*% as.matrix(g, k, 1))
ind <- (P_x1 > P_x0)

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
pnorm(0, mean = - a_0, sd = sqrt(1 + sd_a ^ 2))
pnorm(a_0, mean = 0, sd = sqrt(1 + sd_a ^ 2))
pnorm(a_0 / sqrt(1 + sd_a ^ 2), mean = 0, sd = 1)
mean(P_x0)

pnorm(0, mean = - g_0, sd = sqrt(1 + sd_g ^ 2))
pnorm(g_0, mean = 0, sd = sqrt(1 + sd_g ^ 2))
pnorm(g_0 / sqrt(1 + sd_g ^ 2), mean = 0, sd = 1)
mean(P_x1)

# Expectation of the squared Z = 0 and Z = 1 treatment rates
pmvnorm(upper = c(0, 0), 
        mean = rep(-a_0, 2), 
        sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                      c(sd_a ^ 2, 1 + sd_a ^ 2)))
mean(P_x0 ^ 2)

pmvnorm(upper = c(0, 0), 
        mean = rep(-g_0, 2), 
        sigma = rbind(c(1 + sd_g ^ 2, sd_g ^ 2),
                      c(sd_g ^ 2, 1 + sd_g ^ 2)))
mean(P_x1 ^ 2)

# Variance of the Z = 0 and Z = 1 treatment rates
pmvnorm(upper = c(0, 0), 
        mean = rep(-a_0, 2), 
        sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                      c(sd_a ^ 2, 1 + sd_a ^ 2)),
        keepAttr = FALSE) - pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2)) ^ 2
var(P_x0)

pmvnorm(upper = c(0, 0), 
        mean = rep(-g_0, 2), 
        sigma = rbind(c(1 + sd_g ^ 2, sd_g ^ 2),
                      c(sd_g ^ 2, 1 + sd_g ^ 2)),
        keepAttr = FALSE) - pnorm(0, mean = -g_0, sd = sqrt(1 + sd_g ^ 2)) ^ 2
var(P_x1)

#### Tuning ####

# Targets are
# 1) The fraction of observations where the effect of the instrument is positive
# 2) The share of compliers
# 3) The share of defiers 
# 4) The overall treatment take-up rate (always-takers + compliers with Z = 1 + 
# defiers with Z = 0)
# 5) The standard deviation of the treatment rate with Z = 0
# The idea is to keep the total number of "switchers" (compliers + defiers) 
# fixed, but to vary the individual shares

# Calculate the targets
get_targets <- function(a_0, d_0, g_0, sd_a, sd_d, sd_g, rho_ad, rho_dg, pZ){
  
  # target 1: share of observations where X encourages treatment
  X_plus <- 1 - pnorm(0, mean = d_0, sd = sd_d)
  
  # target 2: share of compliers
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
  
  # target 3: share of defiers
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
  
  # target 4: overall treatment take-up rate
  pD <- always_takers_1 + always_takers_2 + pZ * compliers + (1 - pZ) * defiers
  
  # target 5: standard deviation of P(D = 1|X, Z = 0)
  sd_P_x0 <- sqrt(pmvnorm(upper = c(0, 0), 
                          mean = rep(-a_0, 2), 
                          sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                                        c(sd_a ^ 2, 1 + sd_a ^ 2)),
                          keepAttr = FALSE) - 
                    pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2)) ^ 2)
  
  targets <- c(X_plus,
               compliers, 
               defiers,
               pD,
               sd_P_x0)
  names(targets) <- c("X_plus", "compliers", "defiers", "pD", "sd_P_x0")
  
  return(targets)
  
}

get_targets(a_0, d_0, g_0, sd_a, sd_d, sd_g, rho_ad, rho_dg, pZ)

# Let theta = c(a_0, d_0, sd_a, sd_d, rho_ad) be our vector of parameter values
# to be calibrated. 
# Function to be optimised takes as inputs theta, pZ, and the target values
# Given theta, it calculates the derived parameters g_0, sd_g and rho_dg
# For these parameters, it calls the function `get_targets'
# Returns the distance between the output of the function and the target values

discrepancy <- function(theta, pZ, targets){
  
  a_0 <- theta[1]
  d_0 <- theta[2]
  sd_a <- theta[3]
  sd_d <- theta[4]
  rho_ad <- theta[5]
  
  g_0 <- a_0 + d_0
  sd_g <- sqrt(sd_a ^ 2 + 2 * sd_a * sd_d * rho_ad + sd_d ^ 2)
  rho_dg <- (sd_d + sd_a * rho_ad) / sd_g
  
  vals <- get_targets(a_0 = a_0,
                      d_0 = d_0,
                      g_0 = g_0,
                      sd_a = sd_a,
                      sd_d = sd_d,
                      sd_g = sd_g,
                      rho_ad = rho_ad,
                      rho_dg = rho_dg,
                      pZ = pZ)
  
  return(sqrt(sum((vals - targets) ^ 2)))
  
}

# initial targets
# 1) effect is positive 2 / 3 of the time
# 2) compliers 0.15, defiers 0.05
# 3) take up rate 0.55
# 4) standard deviation of the Z = 0 take up rate is 0.1

targets <- c(0.66, 0.15, 0.05, 0.55, 0.1)
theta <- c(a_0, d_0, sd_a, sd_d, rho_ad)

discrepancy(theta = theta, pZ = pZ, targets = targets)

# Initial guess: suppose there are no defiers. Then the only parameters to be 
# set are a_0, d_0 and sd_a. Let the targets be
# 1) compliers 0.2
# 2) take up rate 0.55
# 3) standard deviation of the Z = 0 take up rate is 0.1
# We proceed as follows. a_0 and sd_a are chosen to target
# 1) take up rate - P(Z = 1) * compliers
# 2) standard deviation of the Z = 0 take up rate
# Given a_0 and sd_a, d_0 can be backed out from the third condition

simple_discrepancy <- function(theta, targets){
  
  a_0 <- theta[1]
  sd_a <- theta[2]
  
  vals <- rep(0, 2)
  vals[1] <- pnorm(a_0 / sqrt(1 + sd_a ^ 2))
  
  vals[2] <- sqrt(pmvnorm(upper = rep(0, 2),
                          mean = rep(-a_0, 2),
                          sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                                        c(sd_a ^ 2, 1 + sd_a ^ 2)),
                          keepAttr = FALSE) -
                    vals[1] ^ 2)
  
  # return(sqrt(sum((targets - vals) ^ 2)))
  return(max(abs(targets - vals)))
  
}

theta <- c(a_0, sd_a)
targets <- rep(0, 2)
targets[1] <- 0.55 - pZ * 0.2 # share of always-takers
targets[2] <- 0.1 # standard deviation

simple_discrepancy(theta, targets)

# even simpler... given sd_a we can solve for a_0, and reduce this to a 
# univariate problem
# a_0 = sqrt(1 + sd_a ^ 2) * qnorm(targets[1])

AT <- 0.55 - pZ * 0.2
get_sd_Z0 <- function(AT, sd_a){
  
  mu <- qnorm(AT) * sqrt(1 + sd_a ^ 2)
  return(sqrt(pmvnorm(upper = rep(0, 2), 
          mean = rep(-mu, 2),
          sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                        c(sd_a ^ 2, 1 + sd_a ^ 2)),
          keepAttr = FALSE) - AT ^ 2))
  
}

temp <- data.frame(sd_a = seq(0, 2, length.out = 101))
temp$sd_Z0 <- sapply(temp$sd_a, get_sd_Z0, AT = AT)

temp %>% ggplot() + 
  geom_line(aes(x = sd_a, y = sd_Z0, col = "sd_Z0")) + 
  theme_bw() + 
  theme(legend.position = "null")

# choose the approximate minimiser
index <- which.min(abs(temp$sd_Z0 - targets[2]))

# refine
temp <- data.frame(sd_a = seq(temp$sd_a[index - 1], temp$sd_a[index + 1], length.out = 101))
temp$sd_Z0 <- sapply(temp$sd_a, get_sd_Z0, AT = AT)

index <- which.min(abs(temp$sd_Z0 - targets[2]))

sd_a <- temp$sd_a[index]
a_0 <- qnorm(AT) * sqrt(1 + sd_a ^ 2)
d_0 <- sqrt(1 + sd_a ^ 2) * qnorm(AT + 0.2) - a_0
g_0 <- a_0 + d_0

# Check
k <- 2
a <- runif(k, -1, 1)
a <- sd_a * a / sqrt(sum(a ^ 2))

# Generate a sample to check calculations
n <- 100000
X <- matrix(rnorm(n * k), n, k)

P_x0 <- pnorm(a_0 + X %*% as.matrix(a, k, 1))
P_x1 <- pnorm(g_0 + X %*% as.matrix(a, k, 1))

sd(P_x0)
targets[2]

mean(P_x0)
AT

mean(P_x0) + pZ * mean(P_x1 - P_x0)
pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2)) + 
  pZ * (pnorm(0, mean = -g_0, sd = sqrt(1 + sd_a ^ 2)) - 
          pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2))) 

theta_no_def <- c(a_0, d_0, sd_a)
names(theta_no_def) <- c("a_0", "d_0", "sd_a")

#### Main Calibration ####

# Start with the no defiers case, gradually increase the proportion of defiers
# With each new iteration, let the initial guess be the solution from the 
# previous iteration.

theta <- matrix(0, nrow = 5, ncol = 10)
row.names(theta) <- c("a_0", "d_0", "sd_a", "sd_d", "rho_ad")
colnames(theta) <- paste0("defiers = ", seq(0, 0.09, length.out = 10))
theta[1:3, 1] <- theta_no_def
theta[4:5, 1] <- 0.001

# Targets are
# 1) The fraction of observations where the effect of the instrument is positive
# (starts close to 1 and declines)
# 2) The share of compliers
# 3) The share of defiers (0.2 - compliers)
# 4) overall treatment take-up rate 0.55
# 5) The standard deviation of the treatment rate with Z = 0 is 0.1

defiers <- seq(0, 0.09, length.out = 10)
compliers <- 0.2 - defiers
X_plus <- 1 - seq(0, 0.45, length.out = 10)

runs <- 10
for (i in 2:runs){
  
  targets <- c(X_plus[i], compliers[i], defiers[i], 0.55, 0.1) 
  guess <- as.numeric(theta[, i - 1])
  discrepancy(guess, pZ, targets)
  
}

temp <- optim(par = guess, 
      fn = discrepancy, 
      method = "Nelder-Mead", 
      pZ = pZ, 
      targets = targets)$par

temp <- optim(par = temp, 
      fn = discrepancy, 
      method = "Nelder-Mead", 
      pZ = pZ, 
      targets = targets)$par

temp_guess <- guess
temp_guess[4:5] <- c(0.1, -0.8)
temp <- optim(par = temp_guess, 
      fn = discrepancy, 
      method = "Nelder-Mead", 
      pZ = pZ, 
      targets = targets)$par



a_0 <- temp[1]
d_0 <- temp[2]
sd_a <- temp[3]
sd_d <- temp[4]
rho_ad <- temp[5]
g_0 <- a_0 + d_0
sd_g <- sqrt(sd_a ^ 2 + 2 * sd_a * sd_d * rho_ad + sd_d ^ 2)
rho_dg <- (sd_d + sd_a * rho_ad) / sd_g

get_targets(a_0 = a_0,
            d_0 = d_0,
            g_0 = g_0,
            sd_a = sd_a,
            sd_d = sd_d,
            sd_g = sd_g,
            rho_ad = rho_ad,
            rho_dg = rho_dg,
            pZ = pZ)
targets

# check 
e_1 <- runif(k, -1, 1)
e_1 <- e_1 / sqrt(sum(e_1 ^ 2))
u <- runif(k, -1, 1)
# create a vector which is orthogonal to e_1
e_2 <- u - (sum(e_1 * u)) * e_1 / (sum(e_1 ^ 2))
e_2 <- e_2 / sqrt(sum(e_2 ^ 2))

# a is a scaled version of e_1
a <- sd_a * e_1

# d is a combination of e_1 and e_2
d <- sd_d * rho_ad * e_1 + sd_d * sqrt(1 - rho_ad ^ 2) * e_2

sd_a - sqrt(sum(a ^ 2))
sd_d - sqrt(sum(d ^ 2))
sd_a * sd_d * rho_ad - sum(a * d)

g <- a + d

# Generate a sample to check calculations
n <- 100000
X <- matrix(rnorm(n * k), n, k)

P_x0 <- pnorm(a_0 + X %*% as.matrix(a, k, 1))
P_x1 <- pnorm(g_0 + X %*% as.matrix(g, k, 1))
ind <- (P_x1 > P_x0)
mean(ind)

shares <- get_shares(a_0, d_0, g_0, sd_a, sd_d, sd_g, rho_ad, rho_dg)

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
pnorm(0, mean = - a_0, sd = sqrt(1 + sd_a ^ 2))
mean(P_x0)

pnorm(0, mean = - g_0, sd = sqrt(1 + sd_g ^ 2))
mean(P_x1)

# Variance of the Z = 0 and Z = 1 treatment rates
pmvnorm(upper = c(0, 0), 
        mean = rep(-a_0, 2), 
        sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                      c(sd_a ^ 2, 1 + sd_a ^ 2)),
        keepAttr = FALSE) - pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2)) ^ 2
var(P_x0)

pmvnorm(upper = c(0, 0), 
        mean = rep(-g_0, 2), 
        sigma = rbind(c(1 + sd_g ^ 2, sd_g ^ 2),
                      c(sd_g ^ 2, 1 + sd_g ^ 2)),
        keepAttr = FALSE) - pnorm(0, mean = -g_0, sd = sqrt(1 + sd_g ^ 2)) ^ 2
var(P_x1)

# Alternative calibration. 
# 1) Fix a_0 and sd_a by jointly targeting the Z = 0 take up rate and its 
# standard deviation
# 2) Fix g_0 and sd_g by jointly targeting the Z = 1 take up rate and its 
# standard deviation
# 3) Choose rho_ad to match the compliers share

# Gunction to calculate the sd of P(D = 1|X, Z = z) given E[P(D = 1|X, Z = z)]
# and the standard deviation of the linear combination of X (sd_Xt)
get_sd_Z <- function(E, sd_Xt){
  
  mu <- qnorm(E) * sqrt(1 + sd_Xt ^ 2)
  return(sqrt(pmvnorm(upper = rep(0, 2), 
                      mean = rep(- mu, 2),
                      sigma = rbind(c(1 + sd_Xt ^ 2, sd_Xt ^ 2),
                                    c(sd_Xt ^ 2, 1 + sd_Xt ^ 2)),
                      keepAttr = FALSE) - E ^ 2))
  
}

targets <- list(E_Z0 = 0.45,
                sd_Z0 = 0.15,
                E_Z1 = 0.55,
                sd_Z1 = 0.1,
                C = 0.15)

# Block One
E_Z0 <- targets$E_Z0
sd_Z0 <- targets$sd_Z0

# quick loop to back out the SD
sd_a <- seq(0, 5, length.out = 101)
sd_Z <- sapply(sd_a, get_sd_Z, E = E_Z0)

# initial guess 
index <- which.min(abs(sd_Z - sd_Z0))
upper <- sd_a[min(length(sd_a), index + 1)]
lower <- sd_a[max(1, index - 1)]
error <- abs(sd_Z[index] - sd_Z0)

while(error > 0.0001){
  
  sd_a <- seq(lower, upper, length.out = 101)
  sd_Z <- sapply(sd_a, get_sd_Z, E = E_Z0)
  
  index <- which.min(abs(sd_Z - sd_Z0))
  upper <- sd_a[min(length(sd_a), index + 1)]
  lower <- sd_a[max(1, index - 1)]
  error <- abs(sd_Z[index] - sd_Z0)
  
}

sd_a <- sd_a[index]
a_0 <-  qnorm(E_Z0) * sqrt(1 + sd_a ^ 2)

# Block Two
E_Z1 <- targets$E_Z1
sd_Z1 <- targets$sd_Z0

# quick loop to back out the SD
sd_g <- seq(0, 5, length.out = 101)
sd_Z <- sapply(sd_g, get_sd_Z, E = E_Z1)

# initial guess 
index <- which.min(abs(sd_Z - sd_Z1))
upper <- sd_g[min(length(sd_g), index + 1)]
lower <- sd_g[max(1, index - 1)]
error <- abs(sd_Z[index] - sd_Z1)

while(error > 0.0001){
  
  sd_g <- seq(lower, upper, length.out = 101)
  sd_Z <- sapply(sd_g, get_sd_Z, E = E_Z1)
  
  index <- which.min(abs(sd_Z - sd_Z1))
  upper <- sd_g[min(length(sd_g), index + 1)]
  lower <- sd_g[max(1, index - 1)]
  error <- abs(sd_Z[index] - sd_Z1)
  
}

sd_g <- sd_g[index]
g_0 <-  qnorm(E_Z1) * sqrt(1 + sd_g ^ 2)

# Block Three
C <- targets$C

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


alt_get_compliers <- function(sd_d, a_0, g_0, sd_a, sd_g){
  
  # calculate the gap between p(x,z = 1) and p(x, z = 0) in the region where Z 
  # encourages treatment 
  
  # back out other parameters
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

sd_d <- seq(abs(sd_g - sd_a), sd_g + sd_a, length.out = 101)

compliers <- sapply(sd_d, 
                    alt_get_compliers, 
                    a_0 = a_0, 
                    g_0 = g_0, 
                    sd_a = sd_a, 
                    sd_g = sd_g)

# initial guess 
index <- which.min(abs(compliers - C))
upper <- sd_d[min(length(compliers), index + 1)]
lower <- sd_d[max(1, index - 1)]
error <- abs(compliers[index] - C)

while(error > 0.0001){
  
  sd_d <- seq(lower, upper, length.out = 101)
  
  compliers <- sapply(sd_d, 
                      alt_get_compliers, 
                      a_0 = a_0, 
                      g_0 = g_0, 
                      sd_a = sd_a, 
                      sd_g = sd_g)

  index <- which.min(abs(compliers - C))
  upper <- sd_d[min(length(compliers), index + 1)]
  lower <- sd_d[max(1, index - 1)]
  error <- abs(compliers[index] - C)
  
}

sd_d <- sd_d[index]
rho_ad <- (sd_g ^ 2 - (sd_a ^ 2 + sd_d ^ 2)) / (2 * sd_a * sd_d)
rho_dg <- (sd_a * rho_ad + sd_d) / sd_g
d_0 <- g_0 - a_0

#### Test Parameter Values ####

# Generate parameters and a sample to check calibration
e_1 <- runif(k, -1, 1)
e_1 <- e_1 / sqrt(sum(e_1 ^ 2))
u <- runif(k, -1, 1)
# create a vector which is orthogonal to e_1
e_2 <- u - (sum(e_1 * u)) * e_1 / (sum(e_1 ^ 2))
e_2 <- e_2 / sqrt(sum(e_2 ^ 2))

# a is a scaled version of e_1
a <- sd_a * e_1

# d is a combination of e_1 and e_2
d <- sd_d * rho_ad * e_1 + sd_d * sqrt(1 - rho_ad ^ 2) * e_2

sd_a - sqrt(sum(a ^ 2))
sd_d - sqrt(sum(d ^ 2))
sd_a * sd_d * rho_ad - sum(a * d)

g <- a + d

n <- 100000
X <- matrix(rnorm(n * k), n, k)

P_x0 <- pnorm(a_0 + X %*% as.matrix(a, k, 1))
P_x1 <- pnorm(g_0 + X %*% as.matrix(g, k, 1))
ind <- (P_x1 > P_x0)
mean(ind)

shares <- get_shares(a_0, d_0, g_0, sd_a, sd_d, sd_g, rho_ad, rho_dg)

shares[1]
shares[2]
shares[3]
shares[4]
shares[5]
shares[6]

# Treatment rates for Z = 0 and Z = 1 observations
pnorm(0, mean = - a_0, sd = sqrt(1 + sd_a ^ 2))
mean(P_x0)

pnorm(0, mean = - g_0, sd = sqrt(1 + sd_g ^ 2))
mean(P_x1)

# Variance of the Z = 0 and Z = 1 treatment rates
sqrt(pmvnorm(upper = c(0, 0), 
        mean = rep(-a_0, 2), 
        sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                      c(sd_a ^ 2, 1 + sd_a ^ 2)),
        keepAttr = FALSE) - pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2)) ^ 2)
sd(P_x0)

sqrt(pmvnorm(upper = c(0, 0), 
        mean = rep(-g_0, 2), 
        sigma = rbind(c(1 + sd_g ^ 2, sd_g ^ 2),
                      c(sd_g ^ 2, 1 + sd_g ^ 2)),
        keepAttr = FALSE) - pnorm(0, mean = -g_0, sd = sqrt(1 + sd_g ^ 2)) ^ 2)
sd(P_x1)
