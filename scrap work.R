# script for scrap work

# Does TSLS work even when X is endogeneous? (in the LATE setting)

K <- 3
B_1 <- runif(K, -1, 1)
B_0 <- runif(K, -1, 1)


rm(list = ls())
gc()

library(clusterGeneration)
library(ggplot2)
library(tidyr)
library(dplyr)

N <- 100
X <- runif(N, 0 , 2)
X_2 <- X ^ 2
X_3 <- X ^ 3

results <- vector(mode = "list", length = 5)
for (i in 1:5){
  
  Y <- X + X_2 + X_3 + rnorm(N)
  results[[i]] <- lm(Y ~ X + X_2 + X_3)$coefficients
  
}

# estimates functions
X <- seq(0, 2, length.out = 101)
Y <- matrix(0, length(X), 5)

for (i in 1:5){
  
  B <- results[[i]]
  Y[, i] <- B[1] + B[2] * X + B[3] * X ^ 2 + B[4] * X ^ 3
  
}

data <- data.frame("X" = rep(X, 5), 
                   "Y" = c(Y), 
                   "run" = c(rep("1", length(X)),
                             rep("2", length(X)),
                             rep("3", length(X)),
                             rep("4", length(X)),
                             rep("5", length(X))))

data %>% ggplot(aes(x = X, y = Y, col = run)) +
  geom_line() +
  theme_bw() 



N <- 1000
B <- 10000
X_means <- rep(0, B)
Y_means <- rep(0, B)

for (i in 1:B){
  
  X <- rbinom(N, 1, 1 / 2)
  Y <- rbinom(N, 1, 1 / 3)
  
  X_means[i] <- mean(X)
  Y_means[i] <- mean(Y)
  
}

plot(density(X_means - 1 / 2))
plot(density(Y_means - 1 / 3))

plot(density(sqrt(N) * (X_means - 1 / 2) * (Y_means - 1 / 3)))
mean(sqrt(N) * (X_means - 1 / 2) * (Y_means - 1 / 3))
sd(sqrt(N) * (X_means - 1 / 2) * (Y_means - 1 / 3))
sqrt(1 / 2 * (1 - 1 / 2))
sqrt(1 / 3 * (1 - 1 / 3))

sqrt(1 / 2 * (1 - 1 / 2)) * sqrt(1 / 3 * (1 - 1 / 3))

N <- 10000
X1 <- rbinom(N, 1, 1/2)
X2 <- rbinom(N, 1, 1/2)

d0 <- qnorm(0.3)
d1 <- qnorm(0.4) - d0
d2 <- qnorm(0.4) - d0
d3 <- qnorm(0.6) - d0 - d1 - d2

Z <- d0 + d1 * X1 + d2 * X2 + d3 * X1 * X2 + rnorm(N) >= 0

Z00 <- Z[as.logical((1 - X1) * (1 - X2))]
Z10 <- Z[as.logical(X1 * (1 - X2))]
Z01 <- Z[as.logical((1 - X1) * X2)]
Z11 <- Z[as.logical(X1 * X2)]

data <- data.frame(Z = Z, X1 = X1, X2 = X2, X1X2 = X1 * X2)

true <- lm(Z ~ X1 + X2 + X1X2, data = data)
false <- lm(Z ~ X1 + X2, data = data)

# approximation error should be postive at (0, 0) and (1, 1), and negative
# at (1, 0) and (0, 1)
mean(Z00) - false$coefficients[1]
mean(Z10) - sum(false$coefficients[1:2])
mean(Z01) - sum(false$coefficients[c(1,3)])
mean(Z11) - sum(false$coefficients)

# propensity score regression

library(mvtnorm)
# variance-covariance matrix for the unobservables
sig <- rbind(c(1, 0.7), c(0.7, 1))
# generate unobservables
U <- rmvnorm(N, mean = rep(0, 2), sigma = sig)
V <- U[, 1]
Uy <- U[, 2]

g0 <- qnorm(0.3)
g1 <- qnorm(0.4) - g0
g2 <- qnorm(0.4) - g0
g3 <- qnorm(0.6) - g0 - g1 - g2
a <- qnorm(0.4) - g0

D <- g0 + g1 * X1 + g2 * X2 + g3 * X1 * X2 + a * Z + V >= 0

data <- data.frame(D = D, X1 = X1, X2 = X2, X1X2 = X1 * X2, Z = Z)

false_pscore <- lm(D ~ X1 + X2, data = data)
true_pscore <- lm(D ~ X1 + X2 + X1X2, data = data)


runs <- 10000
N <- 10000

var_1 <- rep(0, runs)
var_2 <- rep(0, runs)
var_3 <- rep(0, runs)
var_4 <- rep(0, runs)

for (i in 1:runs){
  
  X <- rbinom(N, 1, 1/2)
  Y <- rbinom(N, 1, 1/3)

  mu_X <- sum(X) / N
  mu_Y <- sum(Y) / N
  mu_XY <- sum(X * Y) / N

  var_1[i] <- sqrt(N) * (mu_X - 1 / 2) * mu_Y
  var_2[i] <- sqrt(N) * mu_X * (mu_Y - 1 / 3)
  var_3[i] <- sqrt(N) * (mu_XY - 1 / 6)
  var_4[i] <- sqrt(N) * mu_X
  
}


