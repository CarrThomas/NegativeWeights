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

################################ ADDING OUTCOMES ###############################

#### NO COVARIATES ####

# Check bias calculations in the no covariates case
rm(list = ls())
gc()

n <- 100000 # sample size
pZ <- 1 / 2 # probabilty Z = 1

# Unobservables (V, U_1, U_0) are jointly normally distributed.
# Specify the covariance matrix 
sigma <- diag(1, 3)
sigma[1, 2] <- -0.5
sigma[1, 3] <- 0.5
sigma[, 1] <- sigma[1 ,]
colnames(sigma) = c("V", "U1", "U0")

Z <- rbinom(n, 1, pZ)
U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
colnames(U) <- c("V", "U1", "U0")
U <- data.frame(U)

# Propensity score function
# a_0 <- qnorm(0.4)
# d_0 <- qnorm(0.55) - qnorm(0.4)
#a_0 <- qnorm(0.55)
#d_0 <- qnorm(0.4) - qnorm(0.55)
a_0 <- runif(1, -1, 1)
d_0 <- runif(1, -1, 1)
g_0 <- a_0 + d_0

index <- a_0 + Z * d_0
P <- pnorm(index)
D <- as.numeric(index > U$V) 
# D <- rbinom(n, 1, 1 / 2)

B_0 <- 0 # average Y(0)
zeta_0 <- 1 # ATE (E[Y(1) - Y(0)])

# Outcomes
Y_0 <- B_0 + U$U0
Y_1 <- B_0 + zeta_0 + U$U1
Y <- (1 - D) * Y_0 + D * Y_1

# OLS bias
data <- data.frame(Y = Y, D = D)
OLS <- lm("Y ~ D", data = data)
zeta_hat <- OLS$coefficients[2]
bias <- zeta_hat - zeta_0

mean(U$U1[as.logical((D == 1) * (Z == 1))])
- sigma[2, 2] * sigma[1, 2] * dnorm(g_0) / pnorm(g_0)

mean(U$U1[as.logical((D == 1) * (Z == 0))])
- sigma[2, 2] * sigma[1, 2] * dnorm(a_0) / pnorm(a_0)

mean(U$U0[as.logical((D == 0) * (Z == 1))])
sigma[3, 3] * sigma[1, 3] * dnorm(g_0) / (1 - pnorm(g_0))

mean(U$U0[as.logical((D == 0) * (Z == 0))])
sigma[3, 3] * sigma[1, 3] * dnorm(a_0) / (1 - pnorm(a_0))

mean(U$U1[D == 1])
- pZ * (sigma[2, 2] * sigma[1, 2] * dnorm(g_0) / pnorm(g_0)) -
  pZ *(sigma[2, 2] * sigma[1, 2] * dnorm(a_0) / pnorm(a_0))

mean(U$U0[D == 0])
pZ * sigma[3, 3] * sigma[1, 3] * dnorm(g_0) / (1 - pnorm(g_0)) + 
  pZ *sigma[3, 3] * sigma[1, 3] * dnorm(a_0) / (1 - pnorm(a_0))

mean(U$U1[D == 1]) - mean(U$U0[D == 0])

# analytically calculated bias
- pZ * dnorm(g_0) * ((sigma[2, 2] * sigma[1, 2]) / pnorm(g_0) + 
                       (sigma[3, 3] * sigma[1, 3]) / (1 - pnorm(g_0))) -
  (1 - pZ) * dnorm(a_0) * ((sigma[2, 2] * sigma[1, 2]) / pnorm(a_0) + 
                             (sigma[3, 3] * sigma[1, 3]) / (1 - pnorm(a_0)))
bias

# TSLS numerator
Z_tilde <- Z - mean(Z)
cov(Y, Z_tilde)

# analytical
pZ * (1 -pZ) * (zeta_0 * (pnorm(g_0) - pnorm(a_0)) + 
                  (sigma[2, 2] * sigma[1, 2] - sigma[3, 3] * sigma[1, 3]) * 
  (dnorm(a_0) - dnorm(g_0)))

#### WITH COVARIATES #### 

rm(list = ls())
gc()

# Check bias calculations in the case with covariates

n <- 100000 # sample size
k <- 4
pZ <- 1 / 2 # probabilty Z = 1

# Unobservables (V, U_1, U_0) are jointly normally distributed.
sigma <- diag(1, 3)
sigma[1, 2] <- 0.5
sigma[1, 3] <- -0.5
sigma[, 1] <- sigma[1 ,]
colnames(sigma) <- c("V", "U1", "U0")

Z <- rbinom(n, 1, pZ)
X <- matrix(rnorm(n * k), n, k)
U <- rmvnorm(n, mean = rep(0, 3), sigma = sigma)
colnames(U) <- c("V", "U1", "U0")
U <- data.frame(U)

# Specify standard deviations and correlations for x'd and x'a, generate 
# d and a to satisfy these expectations
path <- getwd()
parameters <- readRDS(paste0(path, "/Monte_Carlo/Parameters1.rds"))
for (i in 1:length(parameters)){
  assign(names(parameters)[i], parameters[[i]])
}

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

g_0 <- a_0 + d_0
sd_g <- sqrt(sum(g ^ 2))

index <- a_0 + X %*% as.matrix(a, k, 1) + Z * (d_0 + X %*% as.matrix(d, k, 1))
P <- pnorm(index)
D <- as.numeric(index  > U[, 1]) 

B_0 <- 0 # average Y(0)
zeta_0 <- 1 # ATE (E[Y(1) - Y(0)])

# Outcomes
Y_0 <- B_0 + U$U0
Y_1 <- B_0 + zeta_0 + U$U1
Y <- (1 - D) * Y_0 + D * Y_1

# OLS bias
data <- data.frame(Y = Y, D = D)
OLS <- lm("Y ~ D", data = data)
zeta_hat <- OLS$coefficients[2]
bias <- zeta_hat - zeta_0

xg <- -1
xa <- -1
mean(U$U1[g_0 + xg > U$V])
-sigma[1, 2] * dnorm(g_0 + xg) / pnorm(g_0 + xg)
mean(U$U1[a_0 + xa > U$V])
-sigma[1, 2] * dnorm(a_0 + xa) / pnorm(a_0 + xa)
mean(U$U0[g_0 + xg < U$V])
sigma[1, 3] * dnorm(g_0 + xg) / (1 - pnorm(g_0 + xg))
mean(U$U0[a_0 + xa < U$V])
sigma[1, 3] * dnorm(a_0 + xa) / (1 - pnorm(a_0 + xa))

integrals <- rep(0, 4)

integrals[1] <- - sigma[2, 2] * sigma[1, 2] * 
  integrate(
    f = Vectorize(function(xg, g_0, sd_g){
      if(pnorm(g_0 + xg) > 0){
        return(dnorm(g_0 + xg) * dnorm(xg, sd = sd_g) / pnorm(g_0 + xg))
      }else{
        return(0)
      }}, SIMPLIFY = TRUE),
    lower = -Inf,
    upper = Inf,
    g_0 = g_0, 
    sd_g = sd_g)$value
integrals[1]
mean(U$U1[as.logical((D == 1) * (Z == 1))])

integrals[2] <- - sigma[2, 2] * sigma[1, 2] * 
  integrate(
    f = Vectorize(function(xa, a_0, sd_a){
      if(pnorm(a_0 + xa) > 0){
        return(dnorm(a_0 + xa) * dnorm(xa, sd = sd_a) / pnorm(a_0 + xa))
      }else{
        return(0)
      }}, SIMPLIFY = TRUE),
    lower = -Inf,
    upper = Inf,
    a_0 = a_0, 
    sd_a = sd_a)$value
integrals[2]
mean(U$U1[as.logical((D == 1) * (Z == 0))])

integrals[3] <- sigma[3, 3] * sigma[1, 3] *  integrate(
  f = Vectorize(function(xg, g_0, sd_g){
    if((1 - pnorm(g_0 + xg)) > 0){
      return(dnorm(g_0 + xg) * dnorm(xg, sd = sd_g) / (1 - pnorm(g_0 + xg)))
    }else{
      return(0)
    }}, SIMPLIFY = TRUE),
  lower = -Inf,
  upper = Inf,
  g_0 = g_0, 
  sd_g = sd_g)$value
integrals[3]
mean(U$U0[as.logical((D == 0) * (Z == 1))])

integrals[4] <- sigma[3, 3] * sigma[1, 3] *  integrate(
  f = Vectorize(function(xa, a_0, sd_a){
    if((1 - pnorm(a_0 + xa) > 0)){
      return(dnorm(a_0 + xa) * dnorm(xa, sd = sd_a) / (1 - pnorm(a_0 + xa)))
    }else{
      return(0)
    }}, SIMPLIFY = TRUE),
  lower = -Inf,
  upper = Inf,
  a_0 = a_0, 
  sd_a = sd_a)$value
integrals[4]
mean(U$U0[as.logical((D == 0) * (Z == 0))])

mean(U$U1[D == 1]) - mean(U$U0[D == 0])
bias

# analytically calculated bias
pZ * (integrals[1] - integrals[3]) + (1 - pZ) * (integrals[2] - integrals[4])

# TSLS numerator
Z_tilde <- Z - mean(Z)
cov(Y, Z_tilde)

# analytical
pZ * (1 - pZ) * 
  (zeta_0 * (pnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) - 
                     pnorm(a_0, sd = sqrt(1 + sd_a ^ 2))) + 
               (sigma[2, 2] * sigma[1, 2] - sigma[3, 3] * sigma[1, 3]) * 
                  (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) - 
                     dnorm(g_0, sd = sqrt(1 + sd_g ^ 2))))

# TSLS estimate
cov(Y, Z_tilde) / cov(D, Z_tilde)
(zeta_0 * (pnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) - 
             pnorm(a_0, sd = sqrt(1 + sd_a ^ 2))) + 
    (sigma[2, 2] * sigma[1, 2] - sigma[3, 3] * sigma[1, 3]) * 
    (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) - 
       dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)))) / 
  (pnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) - 
     pnorm(a_0, sd = sqrt(1 + sd_a ^ 2)))

cov(Y - B_0 - zeta_0 * D, Z_tilde) / cov(D, Z_tilde)
((sigma[2, 2] * sigma[1, 2] - sigma[3, 3] * sigma[1, 3]) * 
    (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) - 
       dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)))) / 
  (pnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) - 
     pnorm(a_0, sd = sqrt(1 + sd_a ^ 2)))



# Corrected TSLS numerator
s <- sign(d_0 + X %*% as.matrix(d, k, 1))
cov(D, s * Z_tilde) / (pZ * (1 - pZ)) # should be compliers + defiers
cov(D, Z_tilde) / (pZ * (1 - pZ)) # should be compliers - defiers

cov(Y, s * Z_tilde)

rho_dg <- sum(d * g) / (sd_d * sd_g)

# analytical 
pZ * (1 - pZ) * (zeta_0 * 0.2 + 
                   (sigma[2, 2] * sigma[1, 2] - sigma[3, 3] * sigma[1, 3]) *
                   (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
                      (1 - pnorm(-d_0,
                                 mean = -a_0 * sd_d * sd_a * rho_ad / 
                                   (1 + sd_a ^ 2),
                                 sd = sd_d * 
                                   sqrt(1 - (rho_ad ^ 2) * (sd_a ^ 2) / 
                                          (1 + sd_a ^ 2)))) - 
                      dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
                      (1 - pnorm(-d_0,
                                 mean = -g_0 * sd_d * sd_g * rho_dg / 
                                   (1 + sd_g ^ 2),
                                 sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
                                                    (1 + sd_g ^ 2)))) -
                      dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
                      pnorm(-d_0,
                            mean = -a_0 * sd_d * sd_a * rho_ad / 
                              (1 + sd_a ^ 2),
                            sd = sd_d * sqrt(1 - rho_ad ^ 2 * sd_a ^ 2 / 
                                               (1 + sd_a ^ 2))) + 
                      dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
                      pnorm(-d_0,
                            mean = -g_0 * sd_d * sd_g * rho_dg / 
                              (1 + sd_g ^ 2),
                            sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
                                               (1 + sd_g ^ 2)))))

cov(Y - B_0 - D * zeta_0, s * Z_tilde)
# analytical 
pZ * (1 - pZ) * ((sigma[2, 2] * sigma[1, 2] - sigma[3, 3] * sigma[1, 3]) *
                   (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
                      (1 - pnorm(-d_0,
                                 mean = -a_0 * sd_d * sd_a * rho_ad / 
                                   (1 + sd_a ^ 2),
                                 sd = sd_d * 
                                   sqrt(1 - (rho_ad ^ 2) * (sd_a ^ 2) / 
                                          (1 + sd_a ^ 2)))) - 
                      dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
                      (1 - pnorm(-d_0,
                                 mean = -g_0 * sd_d * sd_g * rho_dg / 
                                   (1 + sd_g ^ 2),
                                 sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
                                                    (1 + sd_g ^ 2)))) -
                      dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
                      pnorm(-d_0,
                            mean = -a_0 * sd_d * sd_a * rho_ad / 
                              (1 + sd_a ^ 2),
                            sd = sd_d * sqrt(1 - rho_ad ^ 2 * sd_a ^ 2 / 
                                               (1 + sd_a ^ 2))) + 
                      dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
                      pnorm(-d_0,
                            mean = -g_0 * sd_d * sd_g * rho_dg / 
                              (1 + sd_g ^ 2),
                            sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
                                               (1 + sd_g ^ 2)))))

# corrected estimate
cov(Y, s * Z_tilde) / cov(D, s * Z_tilde)

zeta_0 + ((sigma[2, 2] * sigma[1, 2] - sigma[3, 3] * sigma[1, 3]) *
            (dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
               (1 - pnorm(-d_0,
                          mean = -a_0 * sd_d * sd_a * rho_ad / 
                            (1 + sd_a ^ 2),
                          sd = sd_d * 
                            sqrt(1 - (rho_ad ^ 2) * (sd_a ^ 2) / 
                                   (1 + sd_a ^ 2)))) - 
               dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
               (1 - pnorm(-d_0,
                          mean = -g_0 * sd_d * sd_g * rho_dg / 
                            (1 + sd_g ^ 2),
                          sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
                                             (1 + sd_g ^ 2)))) -
               dnorm(a_0, sd = sqrt(1 + sd_a ^ 2)) * 
               pnorm(-d_0,
                     mean = -a_0 * sd_d * sd_a * rho_ad / 
                       (1 + sd_a ^ 2),
                     sd = sd_d * sqrt(1 - rho_ad ^ 2 * sd_a ^ 2 / 
                                        (1 + sd_a ^ 2))) + 
               dnorm(g_0, sd = sqrt(1 + sd_g ^ 2)) * 
               pnorm(-d_0,
                     mean = -g_0 * sd_d * sd_g * rho_dg / 
                       (1 + sd_g ^ 2),
                     sd = sd_d * sqrt(1 - rho_dg ^ 2 * sd_g ^ 2 / 
                                        (1 + sd_g ^ 2))))) / 0.2


############################## CHECKING INTEGRALS ############################## 

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

