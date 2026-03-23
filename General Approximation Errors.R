# Approximation errors take two

# Consider the least restrictive possible model given the assumptions of
# Imbens and Angrist

rm(list=ls())
gc()

library(clusterGeneration)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvtnorm)
library(AER)

############################## SINGLE X CASE ##############################

# check your calculations
N <- 100000

Z <- rbinom(N, 1, 1/2)

a0 <- qnorm(0.4)
a1 <- qnorm(0.8)

sigma <- genPositiveDefMat(3)$Sigma
sigma <- sigma / sigma[3, 3]
# pull out all the elements of sigma
o0 <- sqrt(sigma[1, 1])
o1 <- sqrt(sigma[2, 2])
oV <- sqrt(sigma[3, 3])
rho_01 <- sigma[2, 1] / (o0 * o1)
rho_0V <- sigma[3, 1] / (o0 * oV)
rho_1V <- sigma[3, 2] / (o1 * oV)
 
U <- rmvnorm(N, mean = rep(0, 3), sigma = sigma)
D <- as.numeric(a0 * (1 - Z) + a1 * Z - U[, 3] > 0)

mu_0 <- 0.3
mu_1 <- 0.7
Y <- (1 - D) * (mu_0 + U[, 1]) + D * (mu_1 + U[, 2])

data <- data.frame(Z = Z, D = D, Y = Y)
mean(dplyr::filter(data, Z == 1)$D) - pnorm(a1)
mean(dplyr::filter(data, Z == 0)$D) - pnorm(a0)

# estimate LATE
ivreg(Y ~ D|Z, data = data)

# LATE
mu_1 - mu_0 + (o1 * rho_1V - o0 * rho_0V) * 
  (dnorm(a0) - dnorm(a1)) / 
  (pnorm(a1) - pnorm(a0))

####################### MANY X CASE (INDEPENDENT Z) #######################

# X is a discrete variable, and there is a separate model for each value
# of X.

N <- 500000
N_X <- 20 # number of values X can take
X <- sample(1:N_X, N, replace = TRUE)

Z <- rbinom(N, 1, 1/2)

# choose a0 and a1 for each X
a <- matrix(0, N_X, 2) #a0 in column 1, a1 in column 2
for (i in 1:N_X){
  
  draw <- runif(2, 0, 1)
  a[i, 1] <- qnorm(min(draw))
  a[i, 2] <- qnorm(max(draw))
  
}

# generate a sigma for each X
sigma <- array(0, dim = c(3, 3, N_X))

for (i in 1:N_X){
  
  sigma[, , i] <- genPositiveDefMat(3)$Sigma
  sigma[, , i] <- sigma[, , i] / sigma[3, 3, i]
  
}

# to speed up calculations, sort everything 
X <- sort(X)
indices <- (2:N)[X[2:N] != X[1:(N - 1)]]
first <- c(1, indices)
last <- c(indices - 1, N)

U <- matrix(0, N, 3)

for (i in 1:N_X){
  
  U[first[i]:last[i], ] <- rmvnorm(last[i] - first[i] + 1, 
                                   mean = rep(0, 3),
                                   sigma = sigma[,, i])
  
}

D <- rep(0, N)

for (i in 1:N_X){
  
  D[first[i]:last[i]] <- 
    as.numeric(
      a[i, 1] * (1 - Z[first[i]:last[i]]) + 
        a[i, 2] * Z[first[i]:last[i]] -
        U[first[i]:last[i], 3] > 0)
  
}

mu <- matrix(runif(2 * N_X, 0, 1), N_X, 2) 
# matrix of average treatment effects, column 1 is D = 0, 
# column 2 is D = 1

Y <- rep(0, N)

for (i in 1:N_X){
  
  Y[first[i]:last[i]] <- 
    (1 - D[first[i]:last[i]]) * (mu[i, 1] + U[first[i]:last[i], 1]) +
    D[first[i]:last[i]] * (mu[i, 2] + U[first[i]:last[i], 2])
  
}

data <- data.frame(Z = Z, D = D, Y = Y, X = X)

# check the first stage for each X
checks <- matrix(0, N_X, 2)
for (i in 1:N_X){
  
  temp <- data[first[i]:last[i], ]
  checks[i, 1] <- mean(dplyr::filter(temp, Z == 1)$D) - pnorm(a[i, 2])
  checks[i, 2] <- mean(dplyr::filter(temp, Z == 0)$D) - pnorm(a[i, 1])

}

checks

# estimate LATE
ivreg(Y ~ D|Z, data = data)

# individual LATES
# estimated
LATEs <- rep(0, length = N_X)
# true
actual <- rep(0, length = N_X)

for (i in 1:N_X){
  
  temp <- data[first[i]:last[i], ]
  LATEs[i] <- ivreg(Y ~ D|Z, data = temp)$coefficients[2]
  
  o0 <- sqrt(sigma[1, 1, i])
  o1 <- sqrt(sigma[2, 2, i])
  oV <- sqrt(sigma[3, 3, i])
  rho_01 <- sigma[2, 1, i] / (o0 * o1)
  rho_0V <- sigma[3, 1, i] / (o0 * oV)
  rho_1V <- sigma[3, 2, i] / (o1 * oV)
  
  actual[i] <- mu[i, 2] - mu[i, 1] + (o1 * rho_1V - o0 * rho_0V) * 
    (dnorm(a[i, 1]) - dnorm(a[i, 2])) / 
    (pnorm(a[i, 2]) - pnorm(a[i, 1]))
  
}

LATEs
actual

# weighted average
weights <- pnorm(a[, 2]) - pnorm(a[, 1])
sum(weights * LATEs) / sum(weights)
sum(weights * actual) / sum(weights)
ivreg(Y ~ D|Z, data = data)$coefficients[2]

# incorrect specification (should still be OK as Z is indepedent)
ivreg(Y ~ D + X|X + Z, data = data)

######################## MANY X CASE (DEPENDENT Z) ########################

# X is a discrete variable, and there is a separate model for each value
# of X.

N <- 500000
N_X <- 100 # number of values X can take
X <- sample(1:N_X, N, replace = TRUE)

# to speed up calculations, sort everything 
X <- sort(X)
indices <- (2:N)[X[2:N] != X[1:(N - 1)]]
first <- c(1, indices)
last <- c(indices - 1, N)




# Probability that Z = 1 depends on X
g <- qnorm(runif(N_X, 0.4, 0.6))
# pnorm(g) gives the probability that Z = 1 at each value of X 

Z <- rep(0, N)
U_Z <- rnorm(N)

for (i in 1:N_X){
  
  Z[first[i]:last[i]] <- as.numeric(g[i] - U_Z[first[i]:last[i]] > 0)
  
}

# check that this is correct
checks <- rep(0, N_X)

for (i in 1:N_X){
  
  checks[i] <- mean(Z[first[i]:last[i]]) - pnorm(g[i])
  
}
checks

# choose a0 and a1 for each X
a <- matrix(0, N_X, 2) #a0 in column 1, a1 in column 2
for (i in 1:N_X){
  
  draw <- runif(2, 0, 1)
  a[i, 1] <- qnorm(min(draw))
  a[i, 2] <- qnorm(max(draw))
  
}

# generate a sigma for each X
sigma <- array(0, dim = c(3, 3, N_X))

for (i in 1:N_X){
  
  sigma[, , i] <- genPositiveDefMat(3)$Sigma
  sigma[, , i] <- sigma[, , i] / sigma[3, 3, i]
  
}

U <- matrix(0, N, 3)

for (i in 1:N_X){
  
  U[first[i]:last[i], ] <- rmvnorm(last[i] - first[i] + 1, 
                                   mean = rep(0, 3),
                                   sigma = sigma[,, i])
  
}

D <- rep(0, N)

for (i in 1:N_X){
  
  D[first[i]:last[i]] <- 
    as.numeric(
      a[i, 1] * (1 - Z[first[i]:last[i]]) + 
        a[i, 2] * Z[first[i]:last[i]] -
        U[first[i]:last[i], 3] > 0)
  
}

mu <- matrix(runif(2 * N_X, 0, 1), N_X, 2) 
# matrix of average treatment effects, column 1 is D = 0, 
# column 2 is D = 1

Y <- rep(0, N)

for (i in 1:N_X){
  
  Y[first[i]:last[i]] <- 
    (1 - D[first[i]:last[i]]) * (mu[i, 1] + U[first[i]:last[i], 1]) +
    D[first[i]:last[i]] * (mu[i, 2] + U[first[i]:last[i], 2])
  
}

data <- data.frame(Z = Z, D = D, Y = Y, X = X)

# check the first stage for each X
checks <- matrix(0, N_X, 2)
for (i in 1:N_X){
  
  temp <- data[first[i]:last[i], ]
  checks[i, 1] <- mean(dplyr::filter(temp, Z == 1)$D) - pnorm(a[i, 2])
  checks[i, 2] <- mean(dplyr::filter(temp, Z == 0)$D) - pnorm(a[i, 1])
  
}

checks

# individual LATES
# estimated
LATEs <- rep(0, length = N_X)
# true
actual <- rep(0, length = N_X)

for (i in 1:N_X){
  
  temp <- data[first[i]:last[i], ]
  LATEs[i] <- ivreg(Y ~ D|Z, data = temp)$coefficients[2]
  
  o0 <- sqrt(sigma[1, 1, i])
  o1 <- sqrt(sigma[2, 2, i])
  oV <- sqrt(sigma[3, 3, i])
  rho_01 <- sigma[2, 1, i] / (o0 * o1)
  rho_0V <- sigma[3, 1, i] / (o0 * oV)
  rho_1V <- sigma[3, 2, i] / (o1 * oV)
  
  actual[i] <- mu[i, 2] - mu[i, 1] + (o1 * rho_1V - o0 * rho_0V) * 
    (dnorm(a[i, 1]) - dnorm(a[i, 2])) / 
    (pnorm(a[i, 2]) - pnorm(a[i, 1]))
  
}

# LATEs
# actual

# weighted average
weights <- (pnorm(a[, 2]) - pnorm(a[, 1])) * (pnorm(d) * (1 - pnorm(d)))
sum(weights * LATEs) / sum(weights)
sum(weights * actual) / sum(weights)

# linear specification
iv_lin <- ivreg(Y ~ D + X|Z + X, data = data)

# saturate and control
dummies <- matrix(0, N, N_X)
for (i in 1:N_X){
  
  dummies[, i] <- as.numeric(X == i)
  
}
interactions <- matrix(0, N, N_X)
for (i in 1:N_X){
  
  interactions[, i] <- dummies[, i] * Z
  
}
colnames(dummies) <- paste0("X_", 1:N_X)
colnames(interactions) <- paste0("Z_X_", 1:N_X)
data <- data.frame(Y, D, dummies, interactions)
iv_sat <- 
  ivreg(
    as.formula(
      paste0("Y ~ D + ",
             paste(colnames(dummies), collapse = " + "),
             "|",
             paste(colnames(dummies), collapse = " + "),
             " + ",
             paste(colnames(interactions), collapse = " + ")
             )
      ), 
    data = data)

gc()

iv_lin$coefficients[2]
iv_sat$coefficients[2]
sum(weights * LATEs) / sum(weights)
sum(weights * actual) / sum(weights)

# fitted v. true expectation of Z given X
L_Z <- lm(Z ~ X, data = data.frame(X = X, Z = Z))

# plot of approximated v true function
linear_approx <- L_Z$coefficients[1] + L_Z$coefficients[2] * (1:N_X)
true_expec <- pnorm(g)
data.frame(X = 1:N_X, LZ = linear_approx, EZ = true_expec) %>%
  ggplot() + 
  geom_line(aes(x = X, y = linear_approx, col = "linear approximation")) +
  geom_point(aes(x = X, y = EZ, col = "true expectation")) + 
  theme_bw() + 
  theme(legend.position = "null") + 
  labs(x = "X", y = "")

# E[Z] approximation error
sqrt(sum((linear_approx - true_expec) ^ 2)) / N_X

# E[D] approximation error
L_D <- lm(D ~ X, data = data.frame(X = X, D = D))
E_D <- pnorm(g) * pnorm(a[, 2]) + (1 - pnorm(g)) * pnorm(a[, 1])

e_D <- rep(0, N_X)
# check expectation
for (i in 1:N_X){
  
  e_D[i] <- mean(D[first[i]:last[i]])
  
}

cbind(E_D, e_D)
E_D - e_D

linear_approx_D <- L_D$coefficients[1] + L_D$coefficients[2] * (1:N_X)
data.frame(X = 1:N_X, LD = linear_approx_D, ED = E_D) %>%
  ggplot() + 
  geom_line(aes(x = X, y = LD, col = "linear approximation")) +
  geom_point(aes(x = X, y = ED, col = "true expectation")) + 
  theme_bw() + 
  theme(legend.position = "null") + 
  labs(x = "X", y = "")

# E[D] approximation error
sqrt(sum((linear_approx_D - E_D) ^ 2)) / N_X

# E[Y] approximation error
L_Y <- lm(Y ~ X, data = data.frame(X = X, Y = Y))
E_Y <- E_D * mu[, 2] + (1 - E_D) * mu[, 1]

e_Y <- rep(0, N_X)
# check expectation
for (i in 1:N_X){
  
  e_Y[i] <- mean(Y[first[i]:last[i]])
  
}

cbind(E_Y, e_Y)
E_Y - e_Y

linear_approx_Y <- L_Y$coefficients[1] + L_Y$coefficients[2] * (1:N_X)
data.frame(X = 1:N_X, LY = linear_approx_Y, EY = E_Y) %>%
  ggplot() + 
  geom_line(aes(x = X, y = LY, col = "linear approximation")) +
  geom_point(aes(x = X, y = EY, col = "true expectation")) + 
  theme_bw() + 
  theme(legend.position = "null") + 
  labs(x = "X", y = "")

# E[D] approximation error
sqrt(sum((linear_approx_Y - E_Y) ^ 2)) / N_X

# plot of errors against each other
plot(true_expec - linear_approx, E_Y - linear_approx_Y)
plot(true_expec - linear_approx, E_D - linear_approx_D)


#################### MANY X CASE (DEPENDENT Z SIMPLIFIED) ####################

# X is a discrete variable, and there is a separate model for each value
# of X.
set.seed(1234)

N <- 500000
N_X <- 100 # number of values X can take
X <- sample(1:N_X, N, replace = TRUE)

# to speed up calculations, sort everything 
X <- sort(X)
indices <- (2:N)[X[2:N] != X[1:(N - 1)]]
first <- c(1, indices)
last <- c(indices - 1, N)

# draw all parameters apart from g
# choose a0 and a1 for each X
a <- matrix(0, N_X, 2) #a0 in column 1, a1 in column 2
for (i in 1:N_X){
  
  draw <- runif(2, 0, 1)
  a[i, 1] <- qnorm(min(draw))
  a[i, 2] <- qnorm(max(draw))
  
}

# generate a sigma for each X
sigma <- diag(1, 3, 3)
sigma[2, 1] <- runif(1, -1, 1)
sigma[3, 1] <- runif(1, -1, 1)
sigma[3, 2] <- runif(1, -1, 1)
sigma[1, 2] <- sigma[2, 1]
sigma[1, 3] <- sigma[3, 1]
sigma[2, 3] <- sigma[3, 2]

mu <- matrix(runif(2 * N_X, 0, 1), N_X, 2) 
mu[, 2] <- mu[, 2] + 1
# matrix of average treatment effects, column 1 is D = 0, 
# column 2 is D = 1

# store 1: gbounds, true LATE, saturate and control LATE, linear estimate
results <- matrix(0, 40, 7)
g_lims <- matrix(0, 4, 2)
g_lims[1, ] <- c(0.5, 0.5)
g_lims[2, ] <- c(0.4, 0.6)
g_lims[3, ] <- c(0.3, 0.7)
g_lims[4, ] <- c(0.2, 0.8)
results[, 6] <- 1:nrow(results)
colnames(results) <- c("group", 
                       "true", 
                       "saturate", 
                       "linear", 
                       "percent error",
                       "index",
                       "L2P error")

for (k in 1:4){
  for (j in 1:10){
    
    results[(k - 1) * 10 + j, 1] <- k
    
    # Probability that Z = 1 depends on X
    g <- qnorm(runif(N_X, g_lims[k, 1], g_lims[k, 2]))
    # pnorm(g) gives the probability that Z = 1 at each value of X 
    
    Z <- rep(0, N)
    U_Z <- rnorm(N)
    
    for (i in 1:N_X){
      
      Z[first[i]:last[i]] <- as.numeric(g[i] - U_Z[first[i]:last[i]] > 0)
      
    }
    
    U <- rmvnorm(N, mean = rep(0, 3), sigma = sigma)
    
    D <- rep(0, N)
    
    for (i in 1:N_X){
      
      D[first[i]:last[i]] <- 
        as.numeric(
          a[i, 1] * (1 - Z[first[i]:last[i]]) + 
            a[i, 2] * Z[first[i]:last[i]] -
            U[first[i]:last[i], 3] > 0)
      
    }
    
    Y <- rep(0, N)
    
    for (i in 1:N_X){
      
      Y[first[i]:last[i]] <- 
        (1 - D[first[i]:last[i]]) * (mu[i, 1] + U[first[i]:last[i], 1]) +
        D[first[i]:last[i]] * (mu[i, 2] + U[first[i]:last[i], 2])
      
    }
    
    data <- data.frame(Z = Z, D = D, Y = Y, X = X)
    
    # individual LATES
    # estimated
    LATEs <- rep(0, length = N_X)
    # true
    actual <- rep(0, length = N_X)
    
    for (i in 1:N_X){
      
      temp <- data[first[i]:last[i], ]
      LATEs[i] <- ivreg(Y ~ D|Z, data = temp)$coefficients[2]
      
      o0 <- sqrt(sigma[1, 1])
      o1 <- sqrt(sigma[2, 2])
      oV <- sqrt(sigma[3, 3])
      rho_01 <- sigma[2, 1] / (o0 * o1)
      rho_0V <- sigma[3, 1] / (o0 * oV)
      rho_1V <- sigma[3, 2] / (o1 * oV)
      
      actual[i] <- mu[i, 2] - mu[i, 1] + (o1 * rho_1V - o0 * rho_0V) * 
        (dnorm(a[i, 1]) - dnorm(a[i, 2])) / 
        (pnorm(a[i, 2]) - pnorm(a[i, 1]))
      
    }
    
    # weighted average
    weights <- (pnorm(a[, 2]) - pnorm(a[, 1])) * (pnorm(g) * (1 - pnorm(g)))
    results[(k - 1) * 10 + j, 3] <- sum(weights * LATEs) / sum(weights)
    results[(k - 1) * 10 + j, 2] <- sum(weights * actual) / sum(weights)
    
    # linear specification
    iv_lin <- ivreg(Y ~ D + X|Z + X, data = data)
    
    results[(k - 1) * 10 + j, 4] <- iv_lin$coefficients[2]
    
    # fitted v. true expectation of Z given X
    L_Z <- lm(Z ~ X, data = data.frame(X = X, Z = Z))
    linear_approx <- L_Z$coefficients[1] + L_Z$coefficients[2] * (1:N_X)
    true_expec <- pnorm(g)
    
    # E[Z] approximation error
    results[(k - 1) * 10 + j, 7] <- 
      sqrt((sum((linear_approx - true_expec) ^ 2)) / N_X) 
   
  }
  
  
} 

results[, 5] <- ((results[, 2] - results[, 4]) / results[, 2]) * 100 

results <- data.frame(results)

results %>% ggplot() + 
  geom_point(aes(x = index, y = percent.error, col = as.character(group))) +
  theme_bw() +
  theme(legend.position = "null") +
  labs(x = "", y = "bias (%)") 
  
results %>% ggplot() + 
  geom_point(aes(x = index, y = L2P.error, col = as.character(group))) +
  theme_bw() +
  theme(legend.position = "null") +
  labs(x = "", y = "L2P error") 
