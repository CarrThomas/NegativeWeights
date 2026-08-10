rm(list=ls())
gc()

library(ggplot2)
library(dplyr)
library(tidyr)

n <- 50
p_1 <- runif(n, 0, 1)
p_0 <- runif(n, 0, 1)

data <- data.frame(Z = c(rep(1, n), rep(0, n)), 
                   p = c(p_1, p_0))
data <- data[order(data$p), ]

data$Z1_CDF <- cumsum(data$Z) / n
data$Z0_CDF <- cumsum(1 - data$Z) / n

data %>% ggplot() +
  geom_step(aes(x = p, y = Z1_CDF, col = "Z1")) +
  geom_step(aes(x = p, y = Z0_CDF, col = "Z0")) + 
  theme_bw()
data$D <- data$Z1_CDF - data$Z0_CDF
max(data$D)

data %>% ggplot() + 
  geom_step(aes(x = p, y = D, col = "violation")) +
  theme_bw() +
  theme(legend.position = "null")

# # tau_1 fixed, get tau_0
# tau_1 <- 0.24
# 
# data$S1 <- rep(1, 2 * n)
# data$S1 <- data$S1 - data$Z * (data$Z1_CDF <= tau_1)
# 
# temp <- (data$D - tau_1) * (1 - tau_1) / (data$Z1_CDF - tau_1)
# temp <- temp[(data$Z1_CDF - tau_1) > 0]
# tau_0 <- tau_1 + max(temp)
# 
# data$S1 <- data$S1 - (1 - data$Z) * (data$Z0_CDF > (1 - tau_0))
# 
# data$Z1_trimmed <- 
#   cumsum(data$Z * data$S1) / (sum(data$Z * data$S1))
# data$Z0_trimmed <- 
#   cumsum((1 - data$Z) * data$S1) / sum((1 - data$Z) * data$S1)
# 
# data %>% ggplot() +
#   geom_step(aes(x = p, y = Z1_trimmed, col = "Z1")) +
#   geom_step(aes(x = p, y = Z0_trimmed, col = "Z0")) + 
#   theme_bw()

# tau_0 fixed, get tau_1
tau_0 <- 0.24

data$S1 <- rep(1, 2 * n)
data$S1 <- data$S1 - (1 - data$Z) * (data$Z0_CDF > (1 - tau_0))

temp <- (data$D - tau_0 * data$Z1_CDF) / (1 - data$Z0_CDF - tau_0)
temp <- temp[(1 - data$Z0_CDF - tau_0) > 0]
tau_1 <- ceiling(max(temp) * n) / n

data$S1 <- data$S1 - data$Z * (data$Z1_CDF <= tau_1)

data$Z1_trimmed <- 
  cumsum(data$Z * data$S1) / (sum(data$Z * data$S1))
data$Z0_trimmed <- 
  cumsum((1 - data$Z) * data$S1) / sum((1 - data$Z) * data$S1)

data %>% ggplot() +
  geom_step(aes(x = p, y = Z1_trimmed, col = "Z1")) +
  geom_step(aes(x = p, y = Z0_trimmed, col = "Z0")) + 
  theme_bw()

# Loop over tau_0, get tau_1
N <- 101
tau_0 <- seq(0, 1, length.out = N)
tau_1 <- rep(0, N)

for (i in 1:(N - 1)){
  
  temp <- (data$D - tau_0[i] * data$Z1_CDF) / (1 - data$Z0_CDF - tau_0[i])
  temp <- temp[(1 - data$Z0_CDF - tau_0[i]) > 0]
  tau_1[i] <- min(max(temp), 1)
  
}

data.frame(tau_0 = tau_0[2:(N - 1)], tau_1 = tau_1[2:(N - 1)]) %>%
  ggplot() + 
  geom_step(aes(x = tau_0, y = tau_1, col = "frontier")) +
  theme_bw() +
  theme(legend.position = "null")

tau_0 + tau_1
which.min(tau_0 + tau_1)

data.frame(tau_0 = tau_0[2:(N - 1)], tau_1 = tau_1[2:(N - 1)]) %>%
  ggplot() + 
  geom_step(aes(x = 1 - tau_0, y = 1 - tau_1, col = "frontier")) +
  theme_bw() +
  theme(legend.position = "null")

#########################################################################

# Violation in the middle? 

n <- 50

# need to cook up very particular distributions
p_1 <- rbinom(n, 1, 2 / 3)
p_1 <- p_1 * runif(n, 1 / 4, 1 / 2) + (1 - p_1) * runif(n, 3 / 4, 1)

p_0 <- rbinom(n, 1, 1 / 3)
p_0 <- p_0 * runif(n, 0,1 / 4) + (1 - p_0) * runif(n, 1 / 2,  3 / 4)

data <- data.frame(Z = c(rep(1, n), rep(0, n)), 
                   p = c(p_1, p_0))
data <- data[order(data$p), ]

data$Z1_CDF <- cumsum(data$Z) / n
data$Z0_CDF <- cumsum(1 - data$Z) / n

data %>% ggplot() +
  geom_step(aes(x = p, y = Z1_CDF, col = "Z1")) +
  geom_step(aes(x = p, y = Z0_CDF, col = "Z0")) + 
  theme_bw()
data$D <- data$Z1_CDF - data$Z0_CDF
max(data$D)

data %>% ggplot() + 
  geom_step(aes(x = p, y = D, col = "violation")) +
  theme_bw() +
  theme(legend.position = "null")

# tau_0 fixed, get tau_1
tau_0 <- 0.24

data$S1 <- rep(1, 2 * n)
data$S1 <- data$S1 - (1 - data$Z) * (data$Z0_CDF > (1 - tau_0))

temp <- (data$D - tau_0 * data$Z1_CDF) / (1 - data$Z0_CDF - tau_0)
temp <- temp[(1 - data$Z0_CDF - tau_0) > 0]
tau_1 <- ceiling(max(temp) * n) / n

data$S1 <- data$S1 - data$Z * (data$Z1_CDF <= tau_1)

data$Z1_trimmed <- 
  cumsum(data$Z * data$S1) / (sum(data$Z * data$S1))
data$Z0_trimmed <- 
  cumsum((1 - data$Z) * data$S1) / sum((1 - data$Z) * data$S1)

data %>% ggplot() +
  geom_step(aes(x = p, y = Z1_trimmed, col = "Z1")) +
  geom_step(aes(x = p, y = Z0_trimmed, col = "Z0")) + 
  theme_bw()

# Loop over tau_0, get tau_1
N <- 101
tau_0 <- seq(0, 1, length.out = N)
tau_1 <- rep(0, N)

for (i in 1:(N - 1)){
  
  temp <- (data$D - tau_0[i] * data$Z1_CDF) / (1 - data$Z0_CDF - tau_0[i])
  temp <- temp[(1 - data$Z0_CDF - tau_0[i]) > 0]
  tau_1[i] <- max(0, min(max(temp), 1))
  
}

data.frame(tau_0 = tau_0[2:(N - 1)], tau_1 = tau_1[2:(N - 1)]) %>%
  ggplot() + 
  geom_step(aes(x = tau_0, y = tau_1, col = "frontier")) +
  theme_bw() +
  theme(legend.position = "null")

tau_0 + tau_1
which.min(tau_0 + tau_1)

# different preferences
k <- - 1 / 2 
((1 - tau_0) ^ (1 / k) + (1 - tau_1) ^ (1 / k)) ^ k
which.max((((1 - tau_0) ^ (1 / k) + (1 - tau_1) ^ (1 / k))) ^ k)

log(1 - tau_0) + log(1 - tau_1)
which.max(log(1 - tau_0) + log(1 - tau_1))

which.max(pmin(1 - tau_0, 1 - tau_1))

data.frame(tau_0 = tau_0[2:(N - 1)], tau_1 = tau_1[2:(N - 1)]) %>%
  ggplot() + 
  geom_step(aes(x = 1 - tau_0, y = 1 - tau_1, col = "frontier")) +
  theme_bw() +
  theme(legend.position = "null")

# quick calculations
w <- 1
A_1 <- 0.55
A_0 <- 0.5
p <- -1

# constraint set
N <- 101
max_0 <- (A_1 - A_0) / A_1
tau_0 <- seq(0, max_0, length.out = N)
tau_1 <- rep(0, N)

for (i in 1:(N - 1)){
  
  temp <- A_0 / (1 - tau_0[i])
  tau_1[i] <- (A_1 - temp) / (1 - temp)
  
}

data.frame(tau_0 = tau_0, tau_1 = tau_1) %>%
  ggplot() + 
  geom_line(aes(x = 1 - tau_0, y = 1 - tau_1, col = "frontier")) +
  theme_bw() +
  theme(legend.position = "null")

# utility
U <- (w * (1 - tau_1) ^ p + (1 - tau_0) ^ p) ^ (1 / p)
data.frame(tau_0 = tau_0, tau_1 = tau_1, U = U) %>%
  ggplot() + 
  geom_line(aes(x = tau_0, y = U, col = "utlity")) +
  theme_bw() +
  theme(legend.position = "null")


((w * (1 - A_1) ^ p + (1 - A_0) ^ p) ^ (1 / p)) / (1 - A_0)
t_1 <- (A_1 - A_0) / (1 - A_0)
t_1

((w * (A_1) ^ p + (A_0) ^ p) ^ (1 / p)) / A_1
t_0 <- (A_1 - A_0) / A_1
t_0

((1 - A_1) ^ (p / (1 + p)) * w ^ (1 / (1 + p)) + A_0 ^ (p / (1 + p))) ^ ((p + 1) / p)
C <- (A_0 * w / (1  - A_1)) ^ (1 / (p + 1))
t_1 <- A_1 - A_0 / C
t_0 <- 1 - C * (1 - t_1)
t_1 
t_0

#### Case with multiple violations

w <- 2
A_1 <- c(0.52, 0.58)
A_0 <- c(0.47, 0.53)
p <- - 2

# constraint set
N <- 101
max_0 <- max((A_1 - A_0) / A_1)
tau_0 <- seq(0, max_0, length.out = N)
tau_1 <- rep(0, N)
binding <- rep(0, N)

for (i in 1:N){
  
  temp <- A_0 / (1 - tau_0[i])
  temp_tau <- (A_1 - temp) / (1 - temp)
  temp_tau <- temp_tau[(1 - A_0) > tau_0[i]]
  tau_1[i] <- max(temp_tau)
  
  temp <- 1:length(A_1)
  temp <- temp[(1 - A_0) > tau_0[i]]
  binding[i] <- temp[which.max(temp_tau)]
  
}

data.frame(tau_0 = tau_0, tau_1 = tau_1) %>%
  ggplot() + 
  geom_line(aes(x = 1 - tau_0, y = 1 - tau_1, col = "frontier")) +
  theme_bw() +
  theme(legend.position = "null")
binding

# slope of the constraint
plot(1 - tau_0[1:(N - 1)],
     (-tau_1[2:N] + tau_1[1:(N - 1)]) / (tau_0[2:N] - tau_0[1:(N - 1)]),
     xlab = "1 - tau_0",
     ylab = "slope")

# utility
U <- (w * (1 - tau_1) ^ p + (1 - tau_0) ^ p) ^ (1 / p)
data.frame(tau_0 = tau_0, tau_1 = tau_1, U = U) %>%
  ggplot() + 
  geom_line(aes(x = tau_0, y = U, col = "utlity")) +
  theme_bw() +
  theme(legend.position = "null")

tau_1[which.max(U)]
tau_0[which.max(U)]

# Indifference curves
opt_1 <- 1 - tau_1[which.max(U)]
opt_0 <- 1 - tau_0[which.max(U)]

e <- 0.0005
N <- 100
x_0_plus <- opt_0 + 1:N * e
x_1_plus <- rep(0, N)

x_1_plus[1] <- opt_1 - (1 / w) * ((opt_0 / opt_1) ^ (p - 1)) * e
for (i in 2:N){
  
  x_1_plus[i] <- x_1_plus[i - 1] - 
    (1 / w) * ((x_0_plus[i - 1] / x_1_plus[i - 1]) ^ (p - 1)) * e
  
}

(w * x_1_plus ^ p + x_0_plus ^ p) ^ (1 / p)

x_0_minus <- opt_0 - 1:N * e
x_1_minus <- rep(0, N)

x_1_minus[1] <- opt_1 + (1 / w) * ((opt_0 / opt_1) ^ (p - 1)) * e
for (i in 2:N){
  
  x_1_minus[i] <- x_1_minus[i - 1] + 
    (1 / w) * ((x_0_minus[i - 1] / x_1_minus[i - 1]) ^ (p - 1)) * e
  
}

(w * x_1_plus ^ p + x_0_plus ^ p) ^ (1 / p)

constraint <- data.frame(tau_0 = tau_0, tau_1 = tau_1)
indiff <- data.frame(x_0 = c(x_0_minus, opt_0, x_0_plus),
                     x_1 = c(x_1_minus, opt_1, x_1_plus))

ggplot() + 
  geom_line(aes(x = 1 - tau_0, y = 1 - tau_1, col = "frontier"), data = constraint) +
  geom_line(aes(x = x_0, y = x_1, col = "indiff"), data = indiff) + 
  theme_bw() +
  theme(legend.position = "null")


# perfect substitutes
U <- w * (1 - tau_1) + (1 - tau_0) 
data.frame(tau_0 = tau_0, tau_1 = tau_1, U = U) %>%
  ggplot() + 
  geom_line(aes(x = tau_0, y = U, col = "utlity")) +
  theme_bw() +
  theme(legend.position = "null")
tau_1[which.max(U)]
tau_0[which.max(U)]

#### Leontief preferences

w <- 2
A_1 <- c(0.15, 0.9999)
A_0 <- c(0.1, 0.99985)
p <- -2

# constraint set
N <- 101
max_0 <- max((A_1 - A_0) / A_1)
tau_0 <- seq(0, max_0, length.out = N)
tau_1 <- rep(0, N)
binding <- rep(0, N)

for (i in 1:N){
  
  temp <- A_0 / (1 - tau_0[i])
  temp_tau <- (A_1 - temp) / (1 - temp)
  temp_tau <- temp_tau[(1 - A_0) > tau_0[i]]
  tau_1[i] <- max(temp_tau)
  
  temp <- 1:length(A_1)
  temp <- temp[(1 - A_0) > tau_0[i]]
  binding[i] <- temp[which.max(temp_tau)]
  
}

data.frame(tau_0 = tau_0, tau_1 = tau_1) %>%
  ggplot() + 
  geom_line(aes(x = 1 - tau_0, y = 1 - tau_1, col = "frontier")) +
  theme_bw() +
  theme(legend.position = "null")
binding

# utility
U <- pmin(w * (1 - tau_1), 1 - tau_0)
data.frame(tau_0 = tau_0, tau_1 = tau_1, U = U) %>%
  ggplot() + 
  geom_line(aes(x = tau_0, y = U, col = "utlity")) +
  theme_bw() +
  theme(legend.position = "null")

tau_1[which.max(U)]
tau_0[which.max(U)]

#### CROSSING POINTS

tau_1 <- 0.1
tau_0 <- 0.2
A_0 <- c(0.1, 0.105, 0.4, 0.695, 0.7)
A_1 <- tau_1 + A_0 * (1 - tau_1) / (1 - tau_0)

# constraint set
N <- 101
max_0 <- max((A_1 - A_0) / A_1)
tau_0 <- seq(0, max_0, length.out = N)
tau_1 <- rep(0, N)
binding <- rep(0, N)

for (i in 1:N){
  
  temp <- A_0 / (1 - tau_0[i])
  temp_tau <- (A_1 - temp) / (1 - temp)
  temp_tau <- temp_tau[(1 - A_0) > tau_0[i]]
  tau_1[i] <- max(temp_tau)
  
  temp <- 1:length(A_1)
  temp <- temp[(1 - A_0) > tau_0[i]]
  binding[i] <- temp[which.max(temp_tau)]
  
}

data.frame(tau_0 = tau_0, tau_1 = tau_1) %>%
  ggplot() + 
  geom_line(aes(x = 1 - tau_0, y = 1 - tau_1, col = "frontier")) +
  theme_bw() +
  theme(legend.position = "null")
binding

#### Continuous stuff ####

w <- 1

# just take draws? 

# Z = 1 distribution
Z_1 <- function(p){
  
  return(2 * 
           (min(1 / 6, max((p - 1 / 6), 0)) +
              min(1 / 6, max(p -  1 / 2, 0)) + 
              min(1 / 6, max(p - 5 / 6, 0)))
  )
  
}

Z_0 <- function(p){
  
  return((6 / 4) * min(1 / 6, p) +
           (6 / 2) * min(1 / 6, max(p - 1 / 3, 0)) +
           (6 / 4) * min(1 / 6, max(p - 2 / 3, 0)))
  
}

N <- 101
p <- seq(0, 1, length.out = N)

plot(p, sapply(p, Z_1), type = "l")
lines(p, sapply(p, Z_0), lty = 2)

data <- data.frame(p = p, Z_1 = sapply(p, Z_1), Z_0 = sapply(p, Z_0))
data$D <- data$Z_1 - data$Z_0

data %>% ggplot() +
  geom_line(aes(x = p, y = D, col = "violation")) +
  theme_bw() + 
  theme(legend.position = "null")


violation <- data$D > 0
A_1 <- data$Z_1[violation]
A_0 <- data$Z_0[violation]


# constraint set
N <- 101
max_0 <- max((A_1 - A_0) / A_1)
tau_0 <- seq(0, max_0, length.out = N)
tau_1 <- rep(0, N)
binding <- rep(0, N)

for (i in 1:N){
  
  temp <- A_0 / (1 - tau_0[i])
  temp_tau <- (A_1 - temp) / (1 - temp)
  temp_tau <- temp_tau[(1 - A_0) > tau_0[i]]
  tau_1[i] <- max(temp_tau)
  
  temp <- 1:length(A_1)
  temp <- temp[(1 - A_0) > tau_0[i]]
  binding[i] <- temp[which.max(temp_tau)]
  
}

data.frame(tau_0 = tau_0, tau_1 = tau_1) %>%
  ggplot() + 
  geom_line(aes(x = 1 - tau_0, y = 1 - tau_1, col = "frontier")) +
  theme_bw() +
  theme(legend.position = "null")
binding

# utility
U <- (w * (1 - tau_1) ^ p + (1 - tau_0) ^ p) ^ (1 / p)
data.frame(tau_0 = tau_0, tau_1 = tau_1, U = U) %>%
  ggplot() + 
  geom_line(aes(x = tau_0, y = U, col = "utlity")) +
  theme_bw() +
  theme(legend.position = "null")

tau_1[which.max(U)]
tau_0[which.max(U)]

# perfect substitutes
U <- w * (1 - tau_1) + (1 - tau_0) 
data.frame(tau_0 = tau_0, tau_1 = tau_1, U = U) %>%
  ggplot() + 
  geom_line(aes(x = tau_0, y = U, col = "utlity")) +
  theme_bw() +
  theme(legend.position = "null")
tau_1[which.max(U)]
tau_0[which.max(U)]

#### Case with multiple violations

w <- 1
p <- -2

# plot each individual constraint and then the combined constraint
A_1 <- c(0.9, 0.55, 0.1)
A_0 <- c(0.85, 0.45, 0.05)

N_p <- length(A_1)
N <- 101
constraints <- matrix(0, N * N_p, 3)
constraints[, 3] <- rep(1:N_p, each = N)
index <- 0

for (i in 1:N_p){
  
  max_0 <- (A_1[i] -A_0[i]) / A_1[i]
  tau_0 <- seq(0, max_0, length.out = N)
  tau_1 <- rep(0, N)
  
  for (j in 1:N){
    
    temp <- A_0[i] / (1 - tau_0[j])
    tau_1[j] <- (A_1[i] - temp) / (1 - temp)
    
  }
  
  constraints[index + 1:N, 1:2] <- cbind(tau_0, tau_1)
  index <- index + N
  
}

colnames(constraints) <- c("tau_0", "tau_1", "constraint")
constraints <- data.frame(constraints)

# plot individual constraints 
constraints %>%
  ggplot() + 
  geom_line(aes(x = 1 - tau_0, y = 1 - tau_1, col = as.factor(constraint))) +
  theme_bw() +
  theme(legend.position = "null")

constraints %>%
  ggplot() + 
  geom_line(aes(x = tau_0, y = tau_1, col = as.factor(constraint))) +
  theme_bw() +
  theme(legend.position = "null")

# utility
U <- (w * (1 - tau_1) ^ p + (1 - tau_0) ^ p) ^ (1 / p)
data.frame(tau_0 = tau_0, tau_1 = tau_1, U = U) %>%
  ggplot() + 
  geom_line(aes(x = tau_0, y = U, col = "utlity")) +
  theme_bw() +
  theme(legend.position = "null")

tau_1[which.max(U)]
tau_0[which.max(U)]

################################################################################

rm(list = ls())
gc()

x_0 <- runif(1, 1, 2)
N <- 1001
R <- seq(0.1, 9.9, length.out = N)

a <- runif(1, 0, 0.5)

U <- rep(0, N)

U <- (x_0 ^ (-1)) * (((a ^ p) * (R ^ (- p)) + ((1 - a) ^ p)) ^ (1 / p))

plot(R, U, type = "l")

################################## MAKE PLOTS ##################################

# Single plot showing constraint and indifference curves for different p
w <- 1
A_1 <- c(0.2, 0.55, 0.85)
A_0 <- c(0.15, 0.45, 0.8)
p <- - 2

# constraint set
N <- 1001
max_0 <- max((A_1 - A_0) / A_1)
tau_0 <- seq(0, max_0, length.out = N)
tau_1 <- rep(0, N)
binding <- rep(0, N)

for (i in 1:N){
  
  temp <- A_0 / (1 - tau_0[i])
  temp_tau <- (A_1 - temp) / (1 - temp)
  temp_tau <- temp_tau[(1 - A_0) > tau_0[i]]
  tau_1[i] <- max(temp_tau)
  
  temp <- 1:length(A_1)
  temp <- temp[(1 - A_0) > tau_0[i]]
  binding[i] <- temp[which.max(temp_tau)]
  
}

data.frame(tau_0 = tau_0, tau_1 = tau_1) %>%
  ggplot() + 
  geom_line(aes(x = 1 / (1 - tau_0), y = 1 / (1 - tau_1), col = "frontier")) +
  theme_bw() +
  theme(legend.position = "null")

# utility
p <- c(0.9, -1, -3)
N_p <- length(p)
opt <- rep(0, N_p)

for (i in 1:N_p){
  
  U <- (w * (1 - tau_1) ^ p[i] + (1 - tau_0) ^ p[i]) ^ (1 / p[i])
  opt[i] <- which.max(U)
  
}

# Indifference curves
x_0_range <- c(1 / (1 - tau_0[1]), 1 / (1 - tau_0[N]))
e <- 0.00005
N_i <- ceiling((x_0_range[2] - x_0_range[1]) / e)
x_0 <- x_0_range[1] + e * (1:N_i)
indiff <- matrix(0, N_i, N_p) # each column contains a single indifference curve 

for (i in 1:N_p){
  
  x_1 <- rep(0, N_i)
  
  x_1_star <- 1 / (1 - tau_1[opt[i]])
  x_0_star <- 1 / (1 - tau_0[opt[i]])
  
  index <- max((1:N_i)[x_0 < x_0_star])
  
  # area above optimal x_0
  x_1[index + 1] <- x_1_star - 
    (1 / w) * (x_0_star / x_1_star) ^ (-(p[i] + 1)) * (x_0[index + 1] - x_0_star)
  
  for (j in (index + 2):N_i){
   
    x_1[j] <- x_1[j - 1] - (1 / w) * (x_0[j - 1] / x_1[j - 1]) ^ (-(p[i] + 1)) * e
      
  }
  
  # area below optimal x_0
  x_1[index] <- x_1_star -
    (1 / w) * (x_0_star / x_1_star) ^ (-(p[i] + 1)) * (x_0[index] - x_0_star)
  
  for (j in (index - 1):1){
    
    x_1[j] <- x_1[j + 1] + (1 / w) * (x_0[j + 1] / x_1[j + 1]) ^ (-(p[i] + 1)) * e
    
  }
  
  indiff[, i] <- x_1
  
}

constraint <- data.frame(x_0 = 1 / (1 - tau_0), x_1 = 1 / (1 - tau_1))
curves <- data.frame(p = rep(p, each = N_i),
                     x_0 = rep(x_0, N_p),
                     x_1 = c(indiff))

# first plot
ggplot() +
  geom_line(aes(x = x_0, y = x_1, col = "constraint"), data = constraint) + 
  geom_ribbon(aes(x = x_0, ymin = x_1, ymax = max(x_1)), 
              data = constraint, fill = 'red', alpha = 0.2) + 
  geom_line(aes(x = x_0, y = x_1, col = "indifference curve"),
            data = curves[(2 * N_i) + 1:N_i, ]) +
  annotate(geom = "text", x = 1.2, y = 1.2, label = "feasible choices") + 
  geom_segment(aes(x = 1.1, y = 1.1, xend = 1.05, yend = 1.05), arrow = arrow()) + 
  annotate(geom = "text", x = 1.1, y = 1, label = "increasing payoff") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "null") +
  labs(x = expression(x[0]),
       y = expression(x[1]))

ggplot() +
  geom_line(aes(x = x_0, y = x_1, col = "constraint"), data = constraint) + 
  geom_line(aes(x = x_0, y = x_1, col = "indifference curve"),
            data = curves[(2 * N_i) + 1:N_i, ]) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "null") +
  labs(x = expression(x[0]),
       y = expression(x[1]))

# Convex indifference curve case 
ggplot() +
  geom_line(aes(x = x_0, y = x_1, col = "constraint"), data = constraint) + 
  geom_line(aes(x = x_0, y = x_1, col = "indifference curve"),
            data = curves[1:N_i, ]) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "null") +
  labs(x = expression(x[0]),
       y = expression(x[1]))

# line indifference curve case 
ggplot() +
  geom_line(aes(x = x_0, y = x_1, col = "constraint"), data = constraint) + 
  geom_line(aes(x = x_0, y = x_1, col = "indifference curve"),
            data = curves[N_i + 1:N_i, ]) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "null") +
  labs(x = expression(x[0]),
       y = expression(x[1]))
