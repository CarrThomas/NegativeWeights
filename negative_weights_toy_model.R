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
N <- 10000 #sample size
Z <- rbinom(N, 1, 1 / 2)
p <- 2 / 3 # Pr(X = 1)
X <- rbinom(N, 1, p)

a <- 0.1
PD1 <- X * (1 / 2 + a * (2 * Z - 1)) + 
  (1 - X) * (1 / 2 - a * (2 * Z - 1))
D <- rbinom(N, 1, PD1)

temp <- lm(D ~ X + Z)


# linear regression of D on Z, X should return a 
# coefficient of approximately a * 2 / 3 on Z

# constant treatment effect: estimates should be unbiased
tau <- 1
beta <- 1
# introduce an unobservable that is correlated with D
U <- 0.5 * D + rnorm(N)
Y <- X * beta + D * tau + U

# Binary D, Binary Z, X takes 3 values

# Pr(Z = 1) = Pr(Z = 0) = 1 / 2
# Pr(X = -1) = Pr(X = 0) = Pr(X = 1) =  1 / 3
# Z is positively related with D for X = 1, -1 
# Z is negatively related with D for X = 0
# Pr(D = 1|X = 1) = Pr(D = 1|X = -1) =  1 / 2 + a
# Pr(D = 1|X = 0) = 1 / 2 - a
# with a > 0

N <- 10000 #sample size
Z <- rbinom(N, 1, 1 / 2)
X <- rbinom(N, 1, 1 / 2)
# X <- sample(x = c(-1, 0, 1),
#             size = N,
#             replace = T,
#             prob = rep(1 / 3, 3))
a <- 0.1
ind <- as.numeric(X == 0)
p <- 1 / 2 + (1 - ind) * a * (2 * Z - 1) + ind * a * (1 - 2 * Z)
D <- rbinom(N, 1, p)

# linear regression of D on Z, X should return a 
# coefficient of approximately a * 2 / 3 on Z

# constant treatment effect: estimates should be unbiased
tau <- 1
beta <- 1
# introduce an unobservable that is correlated with D
U <- 0.5 * D + rnorm(N)
Y <- X * beta + D * tau + U

summary(lm(D ~ X + Z))
D_hat <- lm(D ~ X + Z)$fitted.values
iv <- lm(Y ~ X + D_hat)
summary(iv)

# heterogeneous treatment effect: estimates should be unbiased 
# as long as we are `right' enough
# tau is now a vector containing treatment effects for the 
# cases X = 0, X = 1
# we should be OK as long as the treatment effect for X = 1 is larger
tau <- c(1, 0, 1)

ind_0 <- as.numeric(X == 0)
ind_1 <- as.numeric(X == 1)
ind_2 <- 1 - ind_0 - ind_1

Y <-  X * beta + D * (tau[1] * ind_2 + tau[2] * ind_0 + tau[3] * ind_1) + U
iv <- lm(Y ~ X + D_hat)
summary(iv)


# generic specifications
N <- 10000
K_x <- 3 # number of controls
# correlation matrix
Sig <- genPositiveDefMat(K_x)

# Outcome equations
a_0 <- runif(1, 0, 1)
a_1 <- runif(1, 0, 1)
B_0 <- runif(K_x, 0, 1)
B_1 <- runif(K_x, 0, 1)

# Probit
g_0 <- runif(1, 0, 1)
g_1 <- runif(1, 0, 1)
g_x <- runif(K_x, 0, 1)

# Probit plot
temp <- data.frame(x = seq(-5, 5, length.out = 1001))
temp$y <- pnorm(temp$x)
temp %>% ggplot() + 
  geom_line(aes(x = x, y = y, col = "probit"))

# Generate data
# for now assume normally distributed X and independent Z
X <- matrix(rnorm(N * K_x), N, K_x)
Z <- rbinom(N, 1, 0.5)
index <- X %*% as.matrix(g_x, K_x, 1) + Z * g_1 + (1 - Z) * g_0
p_score <- pnorm(index)
# generate (U(0), U(1), V)
U <- rmvnorm(N, sigma = Sig$Sigma)
D <- as.numeric(index + U[, 3] > 0)
Y_0 <- a_0 + X %*% as.matrix(B_0, K_x, 1) + U[, 1]
Y_1 <- a_1 + X %*% as.matrix(B_1, K_x, 1) + U[, 2]
Y <- D * Y_1 + (1 - D) * Y_0

# True Expectations
E_Y_1 <- a_1 + X %*% as.matrix(B_1, K_x, 1) + 
  Sig$Sigma[2, 3] * dnorm(index)/ Sig$Sigma[3, 3]
E_Y_0 <- a_0 + X %*% as.matrix(B_0, K_x, 1) - 
  Sig$Sigma[1, 3] * dnorm(index)/ Sig$Sigma[3, 3]

# Plot of E[U_1|U_D < p(X,Z)]
p <- seq(0, 1, length.out = 101)
plot(p, Sig$Sigma[2, 3] * dnorm(qnorm(p))/ Sig$Sigma[3, 3])
plot(p, -Sig$Sigma[1, 3] * dnorm(qnorm(p))/ Sig$Sigma[3, 3])


# Linearized propensity score
p_lin <- lm(D ~ X + Z)
approx_error <- p_score - p_lin$fitted.values
temp <- ecdf(p_lin$residuals)

true_resid <- D - p_score
temp_2 <- ecdf(true_resid)

# Linearized E[Y(0)|X,Z,D = 0] and E[Y(1)|X,Z,D = 1]
data <- data.frame(Y = Y, D = D, Z = Z)
data <- cbind(data, X)
names(data)[4:(3 + K_x)] <- paste0("X", 1:K_x)

lin_E_Y0 <- 
  lm(as.formula(paste0("Y ~ Z + ",
                       paste(paste0("X", 1:K_x), collapse = " + "))),
     data = data[data$D == 0, ])

lin_E_Y1 <- 
  lm(as.formula(paste0("Y ~ Z + ",
                       paste(paste0("X", 1:K_x), collapse = " + "))),
     data = data[data$D == 1, ])

# comparison
temp_Y1 <- data.frame(true = E_Y_1[D == 1], approx = lin_E_Y1$fitted.values)
temp_Y1 %>% ggplot() + 
  geom_line(aes(x = true, y = approx, col = "Y1"))

true_resid <- (Y_1 - E_Y_1)[D == 1]
plot(ecdf(true_resid))

temp_Y0 <- data.frame(true = E_Y_0[D == 0], approx = lin_E_Y0$fitted.values)
temp_Y0 %>% ggplot() + 
  geom_line(aes(x = true, y = approx, col = "Y0"))

