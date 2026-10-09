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

##################################### SETUP ####################################



# fix parameters
delta <- 0
a <- 0 
d <- qnorm(pnorm(a) + delta) - a
p_Z1 <- 1 / 2 # Pr(Z = 1|X)


n <- seq(2, 50) # possible sample sizes
P_correct <- rep(0, length(n))

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
                                d = d)) * dbinom(n_Z1[j], n[i], p_Z1)
    
  }
  
  P_correct[i] <- sum(P_temp)
  
}

plot(n, P_correct)

# validate by simulation
sims <- 10000
sim_correct <- rep(0, length(n))

for (i in 1:length(n)){
  
  for (s in 1:sims){
    
    Z <- rbinom(n[i], 1, p_Z1)
    V <- rnorm(n[i])
    D <- (a + d * Z >= V) 
    
    if (sum(Z) > 0 & sum(1 - Z > 0)){
      
      sim_correct[i] <- sim_correct[i] + 
        as.numeric(sum(D * Z) / sum(Z) >= sum(D * (1 - Z)) / sum(1 - Z))
      
    }
    
  }
  
  sim_correct[i] <- sim_correct[i] / sims
  
}

sim_correct

plot(n, P_correct)
points(n, sim_correct, col = "red")


#### BUILD A PLOT ####

# Plot the probability that switchers are classified correctly as a function of
# n for different values of the switcher mass

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

plot(n, P_correct[, 5])

plot_data <- matrix(0, length(n) * length(d), 3)
plot_data[, 1] <- c(P_correct)
plot_data[, 2] <- rep(n, 5)
plot_data[, 3] <- rep(delta, each = length(n))
plot_data <- data.frame(plot_data)
names(plot_data) <- c("P", "n", "switchers")

plot_data %>% ggplot() +
  geom_point(aes(x = n, 
                 y = P, 
                 col = as.factor(switchers), 
                 shape = as.factor(switchers))) +
  theme_bw() +
  labs(color = expression(delta),
       shape = expression(delta))


