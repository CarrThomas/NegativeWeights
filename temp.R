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

# tau_1 fixed, get tau_0
tau_1 <- 0.24

data$S1 <- rep(1, 2 *n)
data$S1 <- data$S1 - data$Z * (data$Z1_CDF <= tau_1)

temp <- (data$D - tau_1) * (1 - tau_1) / (data$Z1_CDF - tau_1)
temp <- temp[(data$Z1_CDF - tau_1) > 0]
tau_0 <- tau_1 + max(temp)

data$S1 <- data$S1 - (1 - data$Z) * (data$Z0_CDF > (1 - tau_0))

data$Z1_trimmed <- 
  cumsum(data$Z * data$S1) / (sum(data$Z * data$S1))
data$Z0_trimmed <- 
  cumsum((1 - data$Z) * data$S1) / sum((1 - data$Z) * data$S1)

data %>% ggplot() +
  geom_step(aes(x = p, y = Z1_trimmed, col = "Z1")) +
  geom_step(aes(x = p, y = Z0_trimmed, col = "Z0")) + 
  theme_bw()

# tau_0 fixed, get tau_1
tau_0 <- 0.24

data$S1 <- rep(1, 2 *n)
data$S1 <- data$S1 - (1 - data$Z) * (data$Z0_CDF > (1 - tau_0))

temp <- (data$D - tau_0) * (1 - tau_0) / (1 - data$Z0_CDF - tau_0)
temp <- temp[(1 - data$Z0_CDF - tau_0) > 0]
tau_1 <- tau_0 + max(temp)

data$S1 <- data$S1 - data$Z * (data$Z1_CDF < tau_1)

data$Z1_trimmed <- 
  cumsum(data$Z * data$S1) / (sum(data$Z * data$S1))
data$Z0_trimmed <- 
  cumsum((1 - data$Z) * data$S1) / sum((1 - data$Z) * data$S1)

data %>% ggplot() +
  geom_step(aes(x = p, y = Z1_trimmed, col = "Z1")) +
  geom_step(aes(x = p, y = Z0_trimmed, col = "Z0")) + 
  theme_bw()
