################################# INTRODUCTION #################################

# A script to work out the probability that we classify switchers correctly.

# The propensity score is estimated by the cell mean
# sum(D_i * 1(X_i = x, Z_i = z)) / sum(1(X_i = x, Z_i = z))

# The first stage is specified as D = 1{a + dZ >= V} with V normally distributed
# p(x, 1) = pnorm(a + d)
# p(x, 0) = pnorm(a)

# We fix n (number of observations with X = x), Pr(Z = 1|X = x), a and d. We 
# then work out the probability that switchers are correctly classified.

################################# FRONT MATTER #################################

rm(list=ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(xtable)

################################### FUNCTIONS ##################################

get_prob_correct <- function(k, n_Z0, n_Z1, a, d){
  
  # k is the number of D = 1 observations with Z = 0
  # n_Z0, n_Z1 are the numbers of observations with Z = 0 and Z = 1
  # a and d are propensity score parameters (p(x,z) = pnorm(a + dZ))
  
  # Given n_Z0 observations with Z = 0 and n_Z1 with Z = 1, the function returns
  # the joint probability that
  # i) there are k observations with D = 1 and Z = 0
  # ii) sum(D * Z) / n_Z1 >= k / n_Z0
  # (i.e. the sample average of D with Z = 1 weakly exceeds the Z = 0 average)
  # Rearranging, this inequality holds if
  # sum(D * Z) >= k * n_Z1 / n_Z0 -> d_Z1 >= floor(d_Z0 * n_Z1 / n_Z0)
  # as all numbers are integers.
  
  return(
    pbinom(floor((k - 0.00001) * n_Z1 / n_Z0), n_Z1, pnorm(a + d), 
           lower.tail = FALSE) * 
      dbinom(k, n_Z0, pnorm(a)))
  
}

################################### COMPLIERS ##################################

# The case where d > 0, so that switchers are compliers. 
delta <- c(0, 0.05, 0.1, 0.15, 0.2)
a <- 0 
d <- qnorm(pnorm(a) + delta) - a
p_Z1 <- 1 / 2 # Pr(Z = 1|X)

n <- seq(2, 100) # possible sample sizes
P_correct <- matrix(0, length(n), length(d))

for (l in 1:length(d)){
  
  for (i in 1:length(n)){
    
    # loop over sample sizes
    
    # enumerate possible combinations of Z = 1 and Z = 0 
    n_Z1 <- seq(0, n[i]) # observations with Z = 1
    n_Z0 <- n[i] - n_Z1 # observations with Z = 0
    P_temp <- rep(0, n[i] - 1)
    
    # convention: if we have 0 observations with Z = 0 or Z = 1, then the cell is
    # zeroed out. 
    
    for (j in 2:n[i]){
      
      # i fixes n_Z0 and n_Z1
      # loop over the number of D = 1 draws for Z = 0
      P_temp[j - 1] <- sum(sapply(0:n_Z0[j],
                                  get_prob_correct,
                                  n_Z0 = n_Z0[j],
                                  n_Z1 = n_Z1[j],
                                  a = a,
                                  d = d[l])) * dbinom(n_Z1[j], n[i], p_Z1)
      
    }
    
    P_correct[i, l] <- sum(P_temp)
    
  }
  
}

plot_data <- matrix(0, length(n) * length(d), 3)
plot_data[, 1] <- c(P_correct)
plot_data[, 2] <- rep(n, 5)
plot_data[, 3] <- rep(delta, each = length(n))
plot_data <- data.frame(plot_data)
names(plot_data) <- c("P", "n", "switchers")

plot_data %>% dplyr::filter(switchers > 0) %>% ggplot() +
  geom_point(aes(x = n, 
                 y = P, 
                 col = as.factor(switchers), 
                 shape = as.factor(switchers))) +
  theme_bw() +
  labs(color = expression(delta),
       shape = expression(delta))

############################ VALIDATE BY SIMULATION ############################

sims <- 1000
sim_correct <- matrix(0, length(n), length(d))

for (l in 1:length(d)){
  # loop over parameter values
  
  for (i in 1:length(n)){
    # loop over sample sizes
    
    for (s in 1:sims){
      
      Z <- rbinom(n[i], 1, p_Z1)
      V <- rnorm(n[i])
      D <- (a + d[l] * Z >= V) 
      
      if (sum(Z) > 0 & sum(1 - Z > 0)){
        
        sim_correct[i, l] <- sim_correct[i, l] + 
          as.numeric(sum(D * Z) / sum(Z) >= sum(D * (1 - Z)) / sum(1 - Z))
        
      }
      
      
    }
    
  }
  
}

sim_correct <- sim_correct / sims

sim_plot_data <- matrix(0, length(n) * length(d), 3)
sim_plot_data[, 1] <- c(sim_correct)
sim_plot_data[, 2] <- rep(n, 5)
sim_plot_data[, 3] <- rep(delta, each = length(n))
sim_plot_data <- data.frame(sim_plot_data)
names(sim_plot_data) <- c("P", "n", "switchers")

sim_plot_data %>% dplyr::filter(switchers > 0) %>% ggplot() +
  geom_point(aes(x = n, 
                 y = P, 
                 col = as.factor(switchers), 
                 shape = as.factor(switchers))) +
  theme_bw() +
  labs(color = expression(delta),
       shape = expression(delta))

#################################### DEFIERS ###################################

# The case where d < 0, so that switchers are defiers. 
delta <- - c(0, 0.05, 0.1, 0.15, 0.2)
a <- 0 
d <- qnorm(pnorm(a) + delta) - a
p_Z1 <- 1 / 2 # Pr(Z = 1|X)

n <- seq(2, 100) # possible sample sizes
P_correct <- matrix(0, length(n), length(d))

for (l in 1:length(d)){
  
  for (i in 1:length(n)){
    
    # loop over sample sizes
    
    # enumerate possible combinations of Z = 1 and Z = 0 
    n_Z1 <- seq(0, n[i]) # observations with Z = 1
    n_Z0 <- n[i] - n_Z1 # observations with Z = 0
    P_temp <- rep(0, n[i] - 1)
    
    # convention: if we have 0 observations with Z = 0 or Z = 1, then the cell is
    # zeroed out. 
    
    for (j in 2:n[i]){
      
      # i fixes n_Z0 and n_Z1
      # loop over the number of D = 1 draws for Z = 0
      P_temp[j - 1] <- sum(sapply(0:n_Z0[j],
                                  get_prob_correct,
                                  n_Z0 = n_Z0[j],
                                  n_Z1 = n_Z1[j],
                                  a = a,
                                  d = d[l])) * dbinom(n_Z1[j], n[i], p_Z1)
      
    }
    
    P_correct[i, l] <- 1 - sum(P_temp)
    
  }
  
}

plot_data <- matrix(0, length(n) * length(d), 3)
plot_data[, 1] <- c(P_correct)
plot_data[, 2] <- rep(n, 5)
plot_data[, 3] <- rep(delta, each = length(n))
plot_data <- data.frame(plot_data)
names(plot_data) <- c("P", "n", "switchers")

plot_data %>% dplyr::filter(switchers < 0) %>% ggplot() +
  geom_point(aes(x = n, 
                 y = P, 
                 col = as.factor(switchers), 
                 shape = as.factor(switchers))) +
  theme_bw() +
  labs(color = expression(delta),
       shape = expression(delta))

############################ VALIDATE BY SIMULATION ############################

sims <- 1000
sim_correct <- matrix(0, length(n), length(d))

for (l in 1:length(d)){
  # loop over parameter values
  
  for (i in 1:length(n)){
    # loop over sample sizes
    
    for (s in 1:sims){
      
      Z <- rbinom(n[i], 1, p_Z1)
      V <- rnorm(n[i])
      D <- (a + d[l] * Z >= V) 
      
      if (sum(Z) > 0 & sum(1 - Z > 0)){
        
        sim_correct[i, l] <- sim_correct[i, l] + 
          as.numeric(sum(D * Z) / sum(Z) < sum(D * (1 - Z)) / sum(1 - Z))
        
      }
      
      
    }
    
  }
  
}

sim_correct <- sim_correct / sims

sim_plot_data <- matrix(0, length(n) * length(d), 3)
sim_plot_data[, 1] <- c(sim_correct)
sim_plot_data[, 2] <- rep(n, 5)
sim_plot_data[, 3] <- rep(delta, each = length(n))
sim_plot_data <- data.frame(sim_plot_data)
names(sim_plot_data) <- c("P", "n", "switchers")

sim_plot_data %>% ggplot() +
  geom_point(aes(x = n, 
                 y = P, 
                 col = as.factor(switchers), 
                 shape = as.factor(switchers))) +
  theme_bw() +
  labs(color = expression(delta),
       shape = expression(delta))

########################### COMPLIERS FOR SPECIFIC N ###########################

# The case where d > 0, so that switchers are compliers. 
delta <- c(0, 0.05, 0.1, 0.15, 0.2)
a <- 0 
d <- qnorm(pnorm(a) + delta) - a
p_Z1 <- 1 / 2 # Pr(Z = 1|X)

n <- c(50, 100, 200, 400) # possible sample sizes
P_correct <- matrix(0, length(n), length(d))

for (l in 1:length(d)){
  
  for (i in 1:length(n)){
    
    # loop over sample sizes
    
    # enumerate possible combinations of Z = 1 and Z = 0 
    n_Z1 <- seq(0, n[i]) # observations with Z = 1
    n_Z0 <- n[i] - n_Z1 # observations with Z = 0
    P_temp <- rep(0, n[i] - 1)
    
    # convention: if we have 0 observations with Z = 0 or Z = 1, then the cell is
    # zeroed out. 
    
    for (j in 2:n[i]){
      
      # i fixes n_Z0 and n_Z1
      # loop over the number of D = 1 draws for Z = 0
      P_temp[j - 1] <- sum(sapply(0:n_Z0[j],
                                  get_prob_correct,
                                  n_Z0 = n_Z0[j],
                                  n_Z1 = n_Z1[j],
                                  a = a,
                                  d = d[l])) * dbinom(n_Z1[j], n[i], p_Z1)
      
    }
    
    P_correct[i, l] <- sum(P_temp)
    
  }
  
}

colnames(P_correct) <- delta
rownames(P_correct) <- n
P_correct
P_correct_compliers <- P_correct

############################ DEFIERS FOR SPECIFIC N ############################

# The case where d < 0, so that switchers are defiers. 
delta <- - c(0, 0.05, 0.1, 0.15, 0.2)
a <- 0 
d <- qnorm(pnorm(a) + delta) - a
p_Z1 <- 1 / 2 # Pr(Z = 1|X)

n <- c(50, 100, 200, 400) # possible sample sizes
P_correct <- matrix(0, length(n), length(d))

for (l in 1:length(d)){
  
  for (i in 1:length(n)){
    
    # loop over sample sizes
    
    # enumerate possible combinations of Z = 1 and Z = 0 
    n_Z1 <- seq(0, n[i]) # observations with Z = 1
    n_Z0 <- n[i] - n_Z1 # observations with Z = 0
    P_temp <- rep(0, n[i] - 1)
    
    # convention: if we have 0 observations with Z = 0 or Z = 1, then the cell is
    # zeroed out. 
    
    for (j in 2:n[i]){
      
      # i fixes n_Z0 and n_Z1
      # loop over the number of D = 1 draws for Z = 0
      P_temp[j - 1] <- sum(sapply(0:n_Z0[j],
                                  get_prob_correct,
                                  n_Z0 = n_Z0[j],
                                  n_Z1 = n_Z1[j],
                                  a = a,
                                  d = d[l])) * dbinom(n_Z1[j], n[i], p_Z1)
      
    }
    
    P_correct[i, l] <- 1 - sum(P_temp)
    
  }
  
}

colnames(P_correct) <- delta
rownames(P_correct) <- n
P_correct

P_correct_defiers <- P_correct

############################## MAKE A MASTER TABLE #############################

# combine both sets of results
K <- length(delta)
P_correct <- cbind(P_correct_defiers, P_correct_compliers)
P_correct <- P_correct[, c(2:K, (K + 2):(2 * K))]
P_correct <- P_correct[, c((K - 1):1,K:(2 * (K - 1)))]

# add a column for the probability that all are correct
# apply(P_correct, 1, prod)
P_correct <- cbind(P_correct, apply(P_correct, 1, prod))
xtable(P_correct, digits = 3)

# replace with averages
temp <- matrix(0, length(n), K - 1)
for (i in 1:length(n)){
  for (j in 1:(K - 1)){
    temp[i, j] <- (P_correct[i, j] + P_correct[i, 2 * (K - 1) + 1 - j]) / 2
  }
}
temp <- temp[, (K - 1):1]
rownames(temp) <- n
colnames(temp) <- -1 * delta[2:K]
temp <- cbind(temp, P_correct[, 2 * (K - 1) + 1])
xtable(temp, digits = 3)

xtable(1 - temp, digits = 3)
