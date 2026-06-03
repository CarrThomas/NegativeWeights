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


############################## ROUGH MODEL ##############################

rm(list = ls())
gc()

# X takes three values
# Fix E[Z] and E[XZ], get covariance
n <- 3
p <- 1 / n
a <- 1
b <- 3
E_X <- (b + a) / 2
V_X <- sum(p * (a:b - E_X) ^ 2)
sd_X <- sqrt(V_X)

E_Z <- 1 / 2
V_Z <- E_Z * (1 - E_Z)
sd_Z <- sqrt(V_Z) 
corr_XZ <- 1 / 2  
cov_XZ <- corr_XZ * sd_X * sd_Z
E_XZ <- cov_XZ + E_Z * E_X

# back out probabilities from parameters
d_2 <- seq(0, 1, length.out = 101)

A <- E_Z - p * d_2
B <- E_XZ - p * 2 * d_2

e_2 <- d_2 - (g_0 + g_1 * 2)

d_a <- (b * A - B) / (p * (b - a))
d_b <- (B - a * A) / (p * (b - a))

data.frame(d_2 = d_2, d_a = d_a, d_b = d_b) %>% ggplot() + 
  geom_line(aes(x = d_2, y = d_a, col = "d_a")) +
  geom_line(aes(x = d_2, y = d_b, col = "d_b")) + 
  theme_bw()

# check
ind_1 <- as.logical((d_a >= 0) * (d_a <= 1))
ind_2 <- as.logical((d_b >= 0) * (d_b <= 1))
valid <- (1:101)[as.logical(ind_1 * ind_2)]

i <- valid[sample(1:length(valid), 1)]

d <- c(d_a[i], d_2[i], d_b[i])

N <- 100000
X <- sample(1:3, N, replace = T, prob = rep(1 / 3, 3))
X <- sort(X)
change <- X[1:(N-1)] != X[2:N]
first <- c(1, (2:N)[change])
last <- c((1:(N-1))[change], N)

Z <- rep(0, N)

for (j in 1:3){
  
  Z[first[j]:last[j]] <- rbinom(last[j] + 1 - first[j], 1, d[j])
  
}

data <- data.frame(Z = Z, X = X)
mean(X * Z)
E_XZ
mean(Z)
E_Z

# figure out the L2 approximation error in each case
g_1 <- (E_XZ - E_X * E_Z) / V_X
g_0 <- E_Z - g_1 * E_X

L2 <- rep(0, length(d_a))

for (i in 1:length(d_a)){
  
  e <- c(d_a[i] - g_0 - g_1 * a,
         d_2[i] - g_0 - g_1 * 2,
         d_b[i] - g_0 - g_1 * b)
  
  L2[i] <- sqrt(sum(p * e ^ 2))
  
}

data.frame(d_2 = d_2, L2 = L2) %>% ggplot() + 
  geom_line(aes(x = d_2, y = L2, col = "L2")) +
  theme_bw()

e_2 <- d_2 - (g_0 + g_1 * 2)
e_a <- d_a - (g_0 + g_1 * a)
e_b <- d_b - (g_0 + g_1 * b)

# Second stage
# need P(D = 1|X, Z) -> six values, all bounded between 0 and 1

############################# SECOND ATTEMPT #############################

rm(list = ls())
gc()

# x = (-1, 0, 1) with equal probability
X_vals <- c(-1, 0, 1)

n <- 100000
X <- sample(c(-1, 0, 1), n, replace = T)
X <- sort(X)
ind <- (1:(n - 1))[X[1:(n - 1)] != X[2:n]]
first <- c(1, ind + 1)
last <- c(ind, n)

#### E[Z|X] ####

# linear part
d_0 <- 1 / 2
d_1 <- 0

# non-linear part
e <- rep(0, 3)
e[2] <- -0.1
e[1] <- -e[2] / 2
e[3] <- e[1]

# generate Z
Z <- rep(0, n)

for (i in 1:3){
  
  p <- d_0 + d_1 * X_vals[i] + e[i]
  n_X <- last[i] + 1 - first[i]
  Z[first[i]:last[i]] <- rbinom(n_X, 1, p)
  
}

# check that everything is working properly
for (i in 1:3){
  
  p <- d_0 + d_1 * X_vals[i] + e[i]
  print(mean(Z[first[i]:last[i]]) - p)
  
}

lm(Z ~ X, data = data.frame(Z = Z, X = X))$coefficients
c(d_0, d_1)

# Construct P(D = 1|X) so that E[D(E[Z|X] - L[Z|X])] = 0
Delta <- rep(0, 3)
Delta[2] <- 0.1
Delta[1] <- -Delta[2] / 2
Delta[3] <- Delta[1]

a <- rep(0, 3)
a[2] <- - e[2] ^ 2 * Delta[2] / e[2]


############################# THIRD ATTEMPT #############################

rm(list = ls())
gc()

# x = (-1, 0, 1) with equal probability
X_vals <- c(-4, 0, 8)

n <- 100000
X <- sample(c(-1, 0, 1), n, replace = T)
X <- sort(X)
ind <- (1:(n - 1))[X[1:(n - 1)] != X[2:n]]
first <- c(1, ind + 1)
last <- c(ind, n)

#### E[Z|X] ####

# linear part
d_0 <- 1 / 2
d_1 <- 0

# non-linear part
w <- (X_vals[3] - X_vals[2]) / (X_vals[3] - X_vals[1])
e <- rep(0, 3)
e[2] <- -0.1
e[1] <- -w * e[2]
e[3] <- - (1 - w) * e[2]

# generate Z
Z <- rep(0, n)

for (i in 1:3){
  
  p <- d_0 + d_1 * X_vals[i] + e[i]
  n_X <- last[i] + 1 - first[i]
  Z[first[i]:last[i]] <- rbinom(n_X, 1, p)
  
}

# check that everything is working properly
for (i in 1:3){
  
  p <- d_0 + d_1 * X_vals[i] + e[i]
  print(mean(Z[first[i]:last[i]]) - p)
  
}

lm(Z ~ X, data = data.frame(Z = Z, X = X))$coefficients
c(d_0, d_1)

# Solve for Delta
A <- rbind(c(w, 1 - w), c(w ^ 2, (1 - w) ^ 2))
a <- rep(0, 3)
a[2] <- 0.1
a[c(1, 3)] <- as.vector(solve(A) %*% as.matrix(c(a[2], -a[2]), 2, 1))
a

# Numerator

############################# FOURTH ATTEMPT #############################

rm(list = ls())
gc()

# x = (-2, -1, 1, 1) with equal probability
X_vals <- c(-2, -1, 1, 2)
n_vals <- length(X_vals)
p <- 1 / n_vals

d_0 = 1 / 2
d_1 = 1 / 5

# combinations of approximation errors giving the same L2(P) error

# implicit weights
w <- (X_vals[2:(n_vals - 1)] - X_vals[1]) / (X_vals[n_vals] - X_vals[1])
w <- c(0, w, 0)

# Form the matrix
A <- p * ((1 - w[2]) ^ 2 + 1 + w[2] ^ 2)
B <- 2 * p * ((1 - w[2]) * (1 - w[3]) + w[2] * w[3])
C <- p * ((1 - w[3]) ^ 2 + 1 + w[3] ^ 2)

# Our equation is 
# A e_2 ^ 2 + B e_2 e_3 + C e_3 ^ 2 = xi ^ 2
# Want to transform this into an ellipse

mat <- rbind(c(A, B / 2),
             c(B / 2, C))
# Write this equation as
# e^T mat e = xi ^ 2
# where e = (e_2, e_3) is a column vector

e_decomp <- eigen(mat)

# the matrix is decomposed as
# V %*% D %*% t(V)
# V is e_decomp$vectors, the columns of V are the eigenvectors
# D is a diagonal matrix whose entries are the eigenvalues

# We then have
# t(e) %*% V %*% D %*% t(V) %*% e = xi ^ 2
# (t(e) %*% V) %*% D %*% (t(V) %*% e) = xi ^ 2
# where t(e) %*% V = (x y) are our transformed variables. Thus
# (x y) %*% D %*% t(x y) = D[1, 1] * x ^ 2 + D[2, 2] * y ^ 2 

# the equation for an ellipse is
# (x ^ 2 / a ^ 2) + (y ^ 2 / b ^ 2) = 1
# Hence, for us
# a = sqrt(xi ^ 2 / D[1, 1])
# b = sqrt(xi ^ 2 / D[2 , 2])
xi <- 0.05
a <- xi / sqrt(e_decomp$values[1])
b <- xi / sqrt(e_decomp$values[2])

# The ellipse can be traced out as
# (x, y) = (a cos(t), b sin (t))
t <- seq(0, 2 * pi, length.out = 101)
x <- a * cos(t)
y <- b * sin(t)

e <- cbind(x, y) %*% t(e_decomp$vectors)
plot(e[, 1], e[, 2], type = "l")

# verify the conditions hold
e1 <- -((1 - w[2]) * e[, 1] + (1 - w[3]) * e[, 2])
e4 <- -(w[2] * e[, 1] + w[3] * e[, 2])
e <- cbind(e1, e, e4)
rowSums(e ^ 2) * p 
xi ^ 2

par(mfrow = c(2, 2))

plot(e[, 2], e[, 1], type = "l")
plot(e[, 2], e[, 3], type = "l")
plot(e[, 2], e[, 4], type = "l")

par(mfrow = c(1, 1))

# add in constraints
E_Z <- d_0 + matrix(X_vals * d_1, nrow(e), 4, byrow = T) + e
l_cons <- apply(E_Z, 1, function(x){return(min(x) >= 0)})
u_cons <- apply(E_Z, 1, function(x){return(max(x) <= 1)})
cons <- u_cons * l_cons
cons

# plots with the feasible sections highlighted
plot(e[, 2], e[, 1], type = "l")
lines(e[as.logical(cons), 2], e[as.logical(cons), 1], col = "red")

plot(e[, 2], e[, 3], type = "l")
lines(e[as.logical(cons), 2], e[as.logical(cons), 3], col = "red")

plot(e[, 2], e[, 4], type = "l")
lines(e[as.logical(cons), 2], e[as.logical(cons), 4], col = "red")

# Calculate the C and D coefficients
C <- p * e
D <- p * e * E_Z

# caculate the norm
norm <- sqrt(rowSums(C ^ 2) + rowSums(D ^ 2))

plot(1:length(norm), norm)
points((1:length(norm))[as.logical(cons)], norm[as.logical(cons)], col = "red")

par(mfrow = c(2, 2))

plot(e[, 2], norm, type = "l")
lines(e[as.logical(cons), 2], norm[as.logical(cons)], col = "red")

plot(e[, 2], e[, 1], type = "l")
lines(e[as.logical(cons), 2], e[as.logical(cons), 1], col = "red")

plot(e[, 2], e[, 3], type = "l")
lines(e[as.logical(cons), 2], e[as.logical(cons), 3], col = "red")

plot(e[, 2], e[, 4], type = "l")
lines(e[as.logical(cons), 2], e[as.logical(cons), 4], col = "red")

par(mfrow = c(1, 1))

y_lim <- c(min(e) - 0.01, max(e) + 0.01)
plot(e[, 2], e[, 1], type = "l", ylim = y_lim)
lines(e[, 2], e[, 3], col = "red")
lines(e[, 2], e[, 4], col = "blue")

par(mfrow = c(2, 2))

plot(norm, e[, 1], type = "l")
plot(norm, e[, 2], type = "l")
plot(norm, e[, 3], type = "l")
plot(norm, e[, 4], type = "l")

par(mfrow = c(1, 1))

# Plot out points
# maximum and minimum second stage norms
min_ind <- which.min(norm)
max_ind <- which.max(norm)

# E[Z|X] function
plot(X_vals, d_0 + d_1 * X_vals, type = "l")
points(X_vals, d_0 + d_1 * X_vals + e[min_ind, ], col = "red")
points(X_vals, d_0 + d_1 * X_vals + e[max_ind, ], col = "blue")

data <- data.frame(X = X_vals, 
                   linear = d_0 + d_1 * X_vals, 
                   min = d_0 + d_1 * X_vals + e[min_ind, ],
                   max = d_0 + d_1 * X_vals + e[max_ind, ])

data %>% ggplot() + 
  geom_line(aes(x = X, y = linear, col = "linear approximation")) +
  geom_point(aes(x = X, y = min, col = "maximum width")) +
  geom_point(aes(x = X, y = max, col = "minimum width")) + 
  theme_bw() + 
  theme(legend.title = element_blank())


#### Implied polynomial ####
Vandermode <- cbind(rep(1, n_vals), X_vals, X_vals ^ 2, X_vals ^ 3)
coefs <- matrix(0, nrow(E_Z), ncol(E_Z))
for (i in 1:nrow(coefs)){
  
  coefs[i, ] <- solve(Vandermode, E_Z[i, ])
  
}

par(mfrow = c(2, 2))

plot(norm, coefs[, 1])
plot(norm, coefs[, 2])
plot(norm, coefs[, 3])
plot(norm, coefs[, 4])

par(mfrow = c(1, 1))

y_lim <- c(min(coefs) - 0.05, max(coefs) + 0.05)

plot(norm, coefs[, 1], type = "l", ylim = y_lim)
lines(norm, coefs[, 2], col = "blue")
lines(norm, coefs[, 3], col = "red")
lines(norm, coefs[, 4], col = "green")

# check on the calculations
N <- 100000
X <- sample(X_vals, N, replace = T)
X <- sort(X)

e_dist <- rep(0, n_vals)
for (i in 1:n_vals){
  
  e_dist[i] <- sum(X == X_vals[i])
  
}

Z <- rep(0, N)
bot <- 1
top <- e_dist[1]

j <- 1

for (i in 1:4){
  
 Z[bot:top] <- rbinom(e_dist[i], 1, E_Z[j, i])  
 bot <- bot + e_dist[i]
 top <- top + e_dist[i + 1]
  
}

data <- data.frame(X = X, Z = Z)
lm(Z ~ X, data = data)$coefficients
c(d_0, d_1)
