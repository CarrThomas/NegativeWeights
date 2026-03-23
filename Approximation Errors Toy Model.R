# Toy model to understand more about negative weights
library(clusterGeneration)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvtnorm)

rm(list = ls())
gc()

############################## SIMPLEST CASE ##############################

# Binary D, Binary Z, binary X

# In this case if the sign flips between X values then we should always
# be outside the convex hull of treatment effects

N <- 100000 #sample size
Z <- rbinom(N, 1, 1 / 2)
p <- 2 / 3 # Pr(X = 1)
X <- rbinom(N, 1, p)

# variance-covariance matrix for the unobservables
sig <- rbind(c(1, 0.7), c(0.7, 1))
# generate unobservables
U <- rmvnorm(N, mean = rep(0, 2), sigma = sig)
V <- U[, 1]
Uy <- U[, 2]
# first row is the selection unobservable, second is the outcome unobservable

#### FIRST STAGE ####

# positive effect at X = 1, negative effect at X = 0

# specify parameters so that at X = 1
# Pr(D = 1|Z = 1) = 0.6, Pr(D = 1|Z = 0) = 0.4
# and at X = 0
# Pr(D = 1|Z = 1) = 0.4, Pr(D = 1|Z = 0) = 0.6

a1 <- qnorm(0.6)
a0 <- qnorm(0.4)

fZ <- X * (Z * a1 + (1 - Z) * a0) + 
  (1 - X) * (Z * a0 + (1 - Z) * a1)
D <- as.numeric(fZ + V > 0)

# check the specification is correct
data_X1 <- data.frame(D = D[X == 1], Z = Z[X == 1])
mean(data_X1$D)
data_X1Z1 <- dplyr::filter(data_X1, Z == 1)
data_X1Z0 <- dplyr::filter(data_X1, Z == 0)
mean(data_X1Z1$D)
mean(data_X1Z0$D)

data_X0 <- data.frame(D = D[X == 0], Z = Z[X == 0])
mean(data_X0$D)
data_X0Z1 <- dplyr::filter(data_X0, Z == 1)
data_X0Z0 <- dplyr::filter(data_X0, Z == 0)
mean(data_X0Z1$D)
mean(data_X0Z0$D)

# Estimate the first stage
first_stage <- lm(D ~ X + Z)
summary(first_stage)

# coefficient on Z should be given by E[w(X)s(X)a(X)] where
# w(X) = Var(Z|X) / Var(Z)
# s(X) = sign of the effect of Z at X
# a(X) = (absolute value of) the effect of Z on propensity score at X
# (might need to increase the sample size!)
first_stage$coefficients[3]
p * (pnorm(a1) - pnorm(a0)) + (1 - p) * (pnorm(a0) - pnorm(a1))  

#### SECOND STAGE (CONSTANT TREATMENT EFFECT) ####

# estimate should equal the treatment effect

tau <- 1
beta <- 1
# introduce an unobservable that is correlated with D
Y <- X * beta + D * tau + Uy

# fit a model
D_hat <- first_stage$fitted.values
lm(Y ~ X + D_hat)

# should get the same by FWL
D_tilde <- lm(D_hat ~ X)$residuals
lm(Y ~ D_tilde - 1)

# result should be the same in subsamples
data_X1 <- data.frame(Y = Y[X == 1], D = D[X == 1], Z = Z[X == 1])
X1_first_stage <- lm(D ~ Z, data = data_X1)
data_X1$D_tilde <- X1_first_stage$fitted.values
X1_second_stage <- lm(Y ~ D_tilde, data = data_X1)
X1_second_stage$coefficients

data_X0 <- data.frame(Y = Y[X == 0], D = D[X == 0], Z = Z[X == 0])
X0_first_stage <- lm(D ~ Z, data = data_X0)
data_X0$D_tilde <- X0_first_stage$fitted.values
X0_second_stage <- lm(Y ~ D_tilde, data = data_X0)
X0_second_stage$coefficients

#### SECOND STAGE (non-constant treatment effect)

# Estimate should always be wrong

tau <- 1
beta_1 <- 1
beta_0 <- 0.5

# treatment effect for X = 0: tau
# treatment effect for X = 1: tau + (beta_1 - beta_0)
Y <- X * (D * beta_1 + (1 - D) * beta_0) + D * tau + Uy

# fit a model
D_hat <- first_stage$fitted.values
lm(Y ~ X + D_hat)

# coefficient on D_hat should be E[w(x)t(x)] where
# t(x) is the treatment effect at x
# w(x) is the weight on x given by
# w(x) = var(z|x) s(x) a(x) / E[var(z|x) s(x) a(x)]
# with
# s(x) the sign at x
# a(x) the mass of compliers at x

# "numerators" for both weights
o_1 <- pnorm(a_1) - pnorm(a_0)
o_0 <- pnorm(a_1) - pnorm(a_0)
# denominator
w_denom <- p * o_1 - (1 - p) * o_0
w_1 <- o_1 / w_denom
w_0 <- -o_0 / w_denom
# check weights sum to 1
p * w_1 + (1 - p) * w_0
# implied value of the estimate
p * w_1 * (tau + beta_1 - beta_0) + (1 - p) * w_0 * tau

# should get the same with FWL
D_tilde <- lm(D_hat ~ X)$residuals
lm(Y ~ D_tilde - 1)

# the treatment effect should be estimated correctly in subsamples
data_X1 <- data.frame(Y = Y[X == 1], D = D[X == 1], Z = Z[X == 1])
X1_first_stage <- lm(D ~ Z, data = data_X1)
data_X1$D_tilde <- X1_first_stage$fitted.values
X1_second_stage <- lm(Y ~ D_tilde, data = data_X1)
X1_second_stage$coefficients

data_X0 <- data.frame(Y = Y[X == 0], D = D[X == 0], Z = Z[X == 0])
X0_first_stage <- lm(D ~ Z, data = data_X0)
data_X0$D_tilde <- X0_first_stage$fitted.values
X0_second_stage <- lm(Y ~ D_tilde, data = data_X0)
X0_second_stage$coefficients

########################### THREE VALUED X CASE ###########################

# Case where X is a categorical variable that takes three values

# If the sign flips between X values we can be inside or outside the 
# convex hull depending on the distribution of treatment effects for
# Z = 1 and Z = 0

N <- 100000 # sample size
Z <- rbinom(N, 1, 1 / 2)

# specify the number of values of X and its distribution
X_vals <- 3
X_dist <- rep(1 / 3, 3)
p <- X_dist
# generate X
X <- sample(c(0, 1, 2), N, replace = TRUE, prob = X_dist)

# variance-covariance matrix for the unobservables
sig <- rbind(c(1, 0.7), c(0.7, 1))
# generate unobservables
U <- rmvnorm(N, mean = rep(0, 2), sigma = sig)
V <- U[, 1]
Uy <- U[, 2]
# first column is the selection unobservable
# second column is the outcome unobservable

#### FIRST STAGE ####

# positive effect at most values of X, negative effect at X 

# specify parameters so that at X = 1, 2
# Pr(D = 1|Z = 1) = 0.6, Pr(D = 1|Z = 0) = 0.4
# and at X = 0
# Pr(D = 1|Z = 1) = 0.4, Pr(D = 1|Z = 0) = 0.6

a1 <- qnorm(0.6)
a0 <- qnorm(0.4)

# define an indicator for "postive relationship"
ind <- as.numeric(X != 0)
fZ <- ind * (Z * a1 + (1 - Z) * a0) + 
  (1 - ind) * (Z * a0 + (1 - Z) * a1)
D <- as.numeric(fZ + V > 0)

# check the specification is correct
data_X2 <- data.frame(D = D[X == 2], Z = Z[X == 2])
mean(data_X2$D)
data_X2Z1 <- dplyr::filter(data_X2, Z == 1)
data_X2Z0 <- dplyr::filter(data_X2, Z == 0)
mean(data_X2Z1$D)
mean(data_X2Z0$D)

data_X1 <- data.frame(D = D[X == 1], Z = Z[X == 1])
mean(data_X1$D)
data_X1Z1 <- dplyr::filter(data_X1, Z == 1)
data_X1Z0 <- dplyr::filter(data_X1, Z == 0)
mean(data_X1Z1$D)
mean(data_X1Z0$D)

data_X0 <- data.frame(D = D[X == 0], Z = Z[X == 0])
mean(data_X0$D)
data_X0Z1 <- dplyr::filter(data_X0, Z == 1)
data_X0Z0 <- dplyr::filter(data_X0, Z == 0)
mean(data_X0Z1$D)
mean(data_X0Z0$D)

# Estimate the first stage
first_stage <- lm(D ~ X + Z)
summary(first_stage)

# coefficient on Z should be given by E[w(X)s(X)a(X)] where
# w(X) = Var(Z|X) / Var(Z)
# s(X) = sign of the effect of Z at X
# a(X) = (absolute value of) the effect of Z on propensity score at X
# (might need to increase the sample size!)
first_stage$coefficients[3]
(p[2] + p[3]) * (pnorm(a1) - pnorm(a0)) + 
  p[1] * (pnorm(a0) - pnorm(a1))  

#### SECOND STAGE (CONSTANT TREATMENT EFFECT) ####

# estimate should equal the treatment effect

tau <- 1
beta <- 1
# introduce an unobservable that is correlated with D
Y <- X * beta + D * tau + Uy

# fit a model
D_hat <- first_stage$fitted.values
lm(Y ~ X + D_hat)

# should get the same by FWL
D_tilde <- lm(D_hat ~ X)$residuals
lm(Y ~ D_tilde - 1)

# result should be the same in subsamples
data_X2 <- data.frame(Y = Y[X == 2], D = D[X == 2], Z = Z[X == 2])
X2_first_stage <- lm(D ~ Z, data = data_X2)
data_X2$D_tilde <- X2_first_stage$fitted.values
X2_second_stage <- lm(Y ~ D_tilde, data = data_X2)
X2_second_stage$coefficients

data_X1 <- data.frame(Y = Y[X == 1], D = D[X == 1], Z = Z[X == 1])
X1_first_stage <- lm(D ~ Z, data = data_X1)
data_X1$D_tilde <- X1_first_stage$fitted.values
X1_second_stage <- lm(Y ~ D_tilde, data = data_X1)
X1_second_stage$coefficients

data_X0 <- data.frame(Y = Y[X == 0], D = D[X == 0], Z = Z[X == 0])
X0_first_stage <- lm(D ~ Z, data = data_X0)
data_X0$D_tilde <- X0_first_stage$fitted.values
X0_second_stage <- lm(Y ~ D_tilde, data = data_X0)
X0_second_stage$coefficients

#### SECOND STAGE (non-constant treatment effect)

# write a generic model 
# Y = \sum_x 1{X = x}(B0(x) + D(B1(x) - B0(x))

# Estimate can be "right" or "wrong"
# If the treatment effect for the "defiers" is between the treatment 
# effects for the "compliers", we should be OK. Alternatively, if the
# treatment effects for the defiers are at an extreme, we will always
# be in trouble.

# put all coefficients in a matrix
# each row corresponds to one value of X (X = 0, X = 1, X = 2)
# column 1: no treatment
# column 2: treatment
B <- rbind(c(-1, 4),
           c(1, 1),
           c(2, 3))

Y <- B[X + 1, 1] + D * B[X + 1, 2] + Uy 

# coefficient on D_hat should be E[w(x)t(x)] where
# t(x) is the treatment effect at x
# w(x) is the weight on x given by
# w(x) = var(z|x) s(x) a(x) / E[var(z|x) s(x) a(x)]
# with
# s(x) the sign at x
# a(x) the mass of compliers at x

# "numerators"
w_num <- rep(pnorm(a_1) - pnorm(a_0), 3)
w_num[1] <- - w_num[1]

# denominator
w_denom <- sum(p * w_num)

# weights
w <- w_num / w_denom

# check weights sum to 1
sum(p * w)

# implied value of the estimate
sum(p * w * B[, 2])

# fit a model
D_hat <- first_stage$fitted.values
lm(Y ~ X + D_hat)

# should get the same with FWL
D_tilde <- lm(D_hat ~ X)$residuals
lm(Y ~ D_tilde - 1)

# the treatment effect should be estimated correctly in subsamples
data_X2 <- data.frame(Y = Y[X == 2], D = D[X == 2], Z = Z[X == 2])
X2_first_stage <- lm(D ~ Z, data = data_X2)
data_X2$D_tilde <- X2_first_stage$fitted.values
X2_second_stage <- lm(Y ~ D_tilde, data = data_X2)
X2_second_stage$coefficients

data_X1 <- data.frame(Y = Y[X == 1], D = D[X == 1], Z = Z[X == 1])
X1_first_stage <- lm(D ~ Z, data = data_X1)
data_X1$D_tilde <- X1_first_stage$fitted.values
X1_second_stage <- lm(Y ~ D_tilde, data = data_X1)
X1_second_stage$coefficients

data_X0 <- data.frame(Y = Y[X == 0], D = D[X == 0], Z = Z[X == 0])
X0_first_stage <- lm(D ~ Z, data = data_X0)
data_X0$D_tilde <- X0_first_stage$fitted.values
X0_second_stage <- lm(Y ~ D_tilde, data = data_X0)
X0_second_stage$coefficients

########################### MANY VALUED X CASE ###########################

# Case where X is a categorical variable that takes many values

# If the sign flips between X values we can be inside or outside the 
# convex hull depending on the distribution of treatment effects for
# Z = 1 and Z = 0

N <- 100000 # sample size
Z <- rbinom(N, 1, 1 / 2)

# specify the number of values of X and its distribution
X_vals <- 20
X_dist <- rep(1 / 3, 3)
p <- X_dist
# generate X
X <- sample(c(0, 1, 2), N, replace = TRUE, prob = X_dist)

# variance-covariance matrix for the unobservables
sig <- rbind(c(1, 0.7), c(0.7, 1))
# generate unobservables
U <- rmvnorm(N, mean = rep(0, 2), sigma = sig)
V <- U[, 1]
Uy <- U[, 2]
# first column is the selection unobservable
# second column is the outcome unobservable

#### FIRST STAGE ####

# positive effect at most values of X, negative effect at X 

# specify parameters so that at X = 1, 2
# Pr(D = 1|Z = 1) = 0.6, Pr(D = 1|Z = 0) = 0.4
# and at X = 0
# Pr(D = 1|Z = 1) = 0.4, Pr(D = 1|Z = 0) = 0.6

a1 <- qnorm(0.6)
a0 <- qnorm(0.4)

# define an indicator for "postive relationship"
ind <- as.numeric(X != 0)
fZ <- ind * (Z * a1 + (1 - Z) * a0) + 
  (1 - ind) * (Z * a0 + (1 - Z) * a1)
D <- as.numeric(fZ + V > 0)

# check the specification is correct
data_X2 <- data.frame(D = D[X == 2], Z = Z[X == 2])
mean(data_X2$D)
data_X2Z1 <- dplyr::filter(data_X2, Z == 1)
data_X2Z0 <- dplyr::filter(data_X2, Z == 0)
mean(data_X2Z1$D)
mean(data_X2Z0$D)

data_X1 <- data.frame(D = D[X == 1], Z = Z[X == 1])
mean(data_X1$D)
data_X1Z1 <- dplyr::filter(data_X1, Z == 1)
data_X1Z0 <- dplyr::filter(data_X1, Z == 0)
mean(data_X1Z1$D)
mean(data_X1Z0$D)

data_X0 <- data.frame(D = D[X == 0], Z = Z[X == 0])
mean(data_X0$D)
data_X0Z1 <- dplyr::filter(data_X0, Z == 1)
data_X0Z0 <- dplyr::filter(data_X0, Z == 0)
mean(data_X0Z1$D)
mean(data_X0Z0$D)

# Estimate the first stage
first_stage <- lm(D ~ X + Z)
summary(first_stage)

# coefficient on Z should be given by E[w(X)s(X)a(X)] where
# w(X) = Var(Z|X) / Var(Z)
# s(X) = sign of the effect of Z at X
# a(X) = (absolute value of) the effect of Z on propensity score at X
# (might need to increase the sample size!)
first_stage$coefficients[3]
(p[2] + p[3]) * (pnorm(a1) - pnorm(a0)) + 
  p[1] * (pnorm(a0) - pnorm(a1))  
gamma_z <- (p[2] + p[3]) * (pnorm(a1) - pnorm(a0)) + 
  p[1] * (pnorm(a0) - pnorm(a1))  

#### SECOND STAGE (CONSTANT TREATMENT EFFECT) ####

# estimate should equal the treatment effect

tau <- 1
beta <- 1
# introduce an unobservable that is correlated with D
Y <- X * beta + D * tau + Uy

# fit a model
D_hat <- first_stage$fitted.values
lm(Y ~ X + D_hat)

# should get the same by FWL
D_tilde <- lm(D_hat ~ X)$residuals
lm(Y ~ D_tilde - 1)

# result should be the same in subsamples
data_X2 <- data.frame(Y = Y[X == 2], D = D[X == 2], Z = Z[X == 2])
X2_first_stage <- lm(D ~ Z, data = data_X2)
data_X2$D_tilde <- X2_first_stage$fitted.values
X2_second_stage <- lm(Y ~ D_tilde, data = data_X2)
X2_second_stage$coefficients

data_X1 <- data.frame(Y = Y[X == 1], D = D[X == 1], Z = Z[X == 1])
X1_first_stage <- lm(D ~ Z, data = data_X1)
data_X1$D_tilde <- X1_first_stage$fitted.values
X1_second_stage <- lm(Y ~ D_tilde, data = data_X1)
X1_second_stage$coefficients

data_X0 <- data.frame(Y = Y[X == 0], D = D[X == 0], Z = Z[X == 0])
X0_first_stage <- lm(D ~ Z, data = data_X0)
data_X0$D_tilde <- X0_first_stage$fitted.values
X0_second_stage <- lm(Y ~ D_tilde, data = data_X0)
X0_second_stage$coefficients

#### SECOND STAGE (non-constant treatment effect)

# write a generic model 
# Y = \sum_x 1{X = x}(B0(x) + D(B1(x) - B0(x))

# Estimate can be "right" or "wrong"
# If the treatment effect for the "defiers" is between the treatment 
# effects for the "compliers", we should be OK. Alternatively, if the
# treatment effects for the defiers are at an extreme, we will always
# be in trouble.

# put all coefficients in a matrix
# each row corresponds to one value of X (X = 0, X = 1, X = 2)
# column 1: no treatment
# column 2: treatment
B <- rbind(c(-1, 0.5),
           c(1, 1),
           c(2, -2))

Y <- B[X + 1, 1] + D * B[X + 1, 2] + Uy 

# coefficient on D_hat should be E[w(x)t(x)] where
# t(x) is the treatment effect at x
# w(x) is the weight on x given by
# w(x) = var(z|x) s(x) a(x) / E[var(z|x) s(x) a(x)]
# with
# s(x) the sign at x
# a(x) the mass of compliers at x

# "numerators"
w_num <- rep(pnorm(a_1) - pnorm(a_0), 3)
w_num[1] <- - w_num[1]

# denominator
w_denom <- sum(p * w_num)

# weights
w <- w_num / w_denom

# check weights sum to 1
sum(p * w)

# implied value of the estimate
sum(p * w * B[, 2])

# fit a model
D_hat <- first_stage$fitted.values
lm(Y ~ X + D_hat)

# should get the same with FWL
D_tilde <- lm(D_hat ~ X)$residuals
lm(Y ~ D_tilde - 1)

# the treatment effect should be estimated correctly in subsamples
data_X2 <- data.frame(Y = Y[X == 2], D = D[X == 2], Z = Z[X == 2])
X2_first_stage <- lm(D ~ Z, data = data_X2)
data_X2$D_tilde <- X2_first_stage$fitted.values
X2_second_stage <- lm(Y ~ D_tilde, data = data_X2)
X2_second_stage$coefficients

data_X1 <- data.frame(Y = Y[X == 1], D = D[X == 1], Z = Z[X == 1])
X1_first_stage <- lm(D ~ Z, data = data_X1)
data_X1$D_tilde <- X1_first_stage$fitted.values
X1_second_stage <- lm(Y ~ D_tilde, data = data_X1)
X1_second_stage$coefficients

data_X0 <- data.frame(Y = Y[X == 0], D = D[X == 0], Z = Z[X == 0])
X0_first_stage <- lm(D ~ Z, data = data_X0)
data_X0$D_tilde <- X0_first_stage$fitted.values
X0_second_stage <- lm(Y ~ D_tilde, data = data_X0)
X0_second_stage$coefficients

# implied value of the estimate
sum(p * w * B[, 2])
# true estimate
sum(p * B[, 2])
# difference
(1 - gamma_z) * B[1, 2]
