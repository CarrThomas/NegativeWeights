################################# FRONT MATTER #################################

# For binary Z and binary D, this script works out the probability that the 
# finite sample covariance between Z and D is positive.

rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)

################################# FRONT MATTER #################################

# set parameters
p_z <- 1 / 2 # Pr(Z = 1)
p_a <- 1 / 3 # always takers
p_n <- 1 / 3 # never takers
p_c <- 1 - p_a - p_n # compliers

p_11 <- p_z * (p_a + p_c) # Pr(D = 1, Z = 1)
p_10 <- (1 - p_z) * p_a # Pr(D = 1, Z = 0)
p_01 <- p_z * p_n # Pr(D = 0, Z = 1)
p_00 <- (1 - p_z) * (p_n + p_c) # Pr(D = 0, Z = 0)

probs <- c(p_11, p_10, p_01, p_00)

# brute force enumeration
n <- 6 # sample size
P <- 0
count <- 0
record <- matrix(0, 100, 4)

for (n_11 in 1:(n - 1)){
  for (n_00 in 1:(n - n_11)){
    for (n_10 in 0:max(0, (n - n_11 - n_00))){
      
      n_01 <- max(0, (n - n_11 - n_00 - n_10))
      
      if ((n_00 * n_11) > (n_10 * n_01)){
        
        P <- P + dmultinom(x = c(n_11, n_10, n_01, n_00),
                           prob = probs)
        count <- count + 1
        record[count, ] <- c(n_11, n_10, n_01, n_00)
        # print(c(n_11, n_10, n_01, n_00))
        
      }
    }
  }
}

P
count
record[1:count, ]

# loop over n
N <- 15
P <- rep(0, N)

for (i in 1:N){
  
  n <- 1 + i # sample size

  for (n_11 in 1:(n - 1)){
    for (n_00 in 1:(n - n_11)){
      for (n_10 in 0:max(0, (n - n_11 - n_00))){
        
        n_01 <- max(0, (n - n_11 - n_00 - n_10))
        
        if ((n_00 * n_11) > (n_10 * n_01)){
          
          P[i] <- P[i] + dmultinom(x = c(n_11, n_10, n_01, n_00),
                             prob = probs)
          
        }
      }
    }
  }
}

