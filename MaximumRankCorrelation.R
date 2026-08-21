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
  bound <- a_0 + rho_dg * sd_g * Xd / sd_d 
  
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

#### Region where Z has a negative effect on the propensity score

# Always Takers
kernel <- function(Xd, g_0, sd_d, sd_g, rho_dg){
  
  sd_temp <- sqrt(1 + (sd_g ^ 2) * (1 - rho_dg ^ 2))
  bound <- a_0 + rho_dg * sd_g * Xd / sd_d 
  
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




       
