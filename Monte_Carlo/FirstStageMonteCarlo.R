#################################### SUMMARY ################################### 

# This scipt performs Monte Carlo for the maximum rank correlation estimator. 

# The first stage is assumed to be
# p(x,z) = Phi(a_0 + x'a + z(d_0 + x'd))
# where Phi is an unknown function. In this formulation, the parameters d_0 and
# d control the sign of the relationship between Z and D. These are our primary
# object of interest. 

# We use the maximum rank correlation estimator. As the constant a_0 does not
# affect ranks, it is *not* estimated by this approach. In addition, the 
# remaining parameters are estimated only up to scale.

################################# FRONT MATTER #################################

rm(list = ls())
gc()

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(boot)

path <- getwd()

################################### FUNCTIONS ##################################

# Generate parameter vectors a and d which match sd_a, sd_d, and rho_ad
gen_ad <- function(k, sd_a, sd_d, rho_ad){
  
  # create orthonormal vectors e_1 and e_2
  e_1 <- runif(k, -1, 1)
  e_1 <- e_1 / sqrt(sum(e_1 ^ 2))
  u <- runif(k, -1, 1)
  e_2 <- u - (sum(e_1 * u)) * e_1 / (sum(e_1 ^ 2))
  e_2 <- e_2 / sqrt(sum(e_2 ^ 2))
  
  # a is a scaled version of e_1
  a <- sd_a * e_1
  # d is a combination of e_1 and e_2
  d <- sd_d * rho_ad * e_1 + sd_d * sqrt(1 - rho_ad ^ 2) * e_2
  
  return(list(a = a,
              d = d))
  
}

# Efficient calculation of the MRC objective function (as suggested by 
# Windmeijer/Abrevata) 
# * theta is a vector of estimated parameters (a, d_0, d)
# * data is a data frame containing X, Z, Z * X
# * varnames is a list of all variables used in the regression (X_1, X_2, etc)
S_theta <- function(theta, data, varnames){
  
  theta <- theta / sqrt(sum(theta ^ 2))
  index <- as.matrix(data[varnames]) %*% theta
  D_sort <- data$D[order(index)] 
  tally <- cumsum(1 - D_sort)
  val <- sum(tally[as.logical(D_sort)])
  
  return(val)
  
}

# old_S_theta <- function(candidates, data, varnames){
#   
#   index <- as.matrix(data[varnames]) %*% candidates
#   D_sort <- data$D[order(index)] 
#   tally <- cumsum(1 - D_sort)
#   val <- sum(tally[as.logical(D_sort)])
#   
#   return(val)
#   
# }

# Estimation function written so that it can be passed to the boot function
# * varnames is a list of variables (X_1,..,X_k,Z,Z * X_1,...,Z * X_k)
# * data is the sample
# * indices is included to allow the boot function to operate
# The fun
est_theta <- function(data, indices, varnames){
  
  boot_data <- data[indices, ] 
  model <- as.formula(paste0("D ~ ", paste(varnames, collapse = " + ")))
  
  probit_est <- glm(model,
                    family = binomial(link = "probit"), 
                    data = boot_data)$coefficients
  theta_probit <- probit_est[names(probit_est) %in% varnames]
  
  logit_est <- glm(model,
                   family = binomial(link = "logit"), 
                   data = boot_data)$coefficients
  theta_logit <- logit_est[names(logit_est) %in% varnames]
  
  lpm_est <- lm(model,
                data = boot_data)$coefficients
  theta_lpm <- lpm_est[names(lpm_est) %in% varnames]
  
  # Run the Nelder-Mead-Head simplex from each pilot, take the maximum
  candidates <- vector(mode = "list", length = 3)
  
  candidates[[1]] <- optim(par = theta_probit, 
                           fn = S_theta, 
                           method = "Nelder-Mead", 
                           data = boot_data, 
                           varnames = varnames, 
                           control = list("fnscale" = -1))
  
  candidates[[2]] <- optim(par = theta_logit, 
                           fn = S_theta, 
                           method = "Nelder-Mead", 
                           data = boot_data, 
                           varnames = varnames, 
                           control = list("fnscale" = -1))
  
  candidates[[3]] <- optim(par = theta_lpm, 
                           fn = S_theta, 
                           method = "Nelder-Mead", 
                           data = boot_data, 
                           varnames = varnames, 
                           control = list("fnscale" = -1))
  
  # Choose the best performing value
  i <- which.max(sapply(candidates, getElement, "value"))
  theta_hat <- candidates[[i]]$par
  theta_hat <- theta_hat / sqrt(sum(theta_hat ^ 2)) # scale to be length 1
  
  return(theta_hat)
  
}



##################################### SETUP ####################################

# Load in parameters
parameters <- readRDS(paste0(path, "/Monte_Carlo/Parameters1.rds"))

for (i in 1:length(parameters)){
  assign(names(parameters)[i], parameters[[i]])
}

# Settings 
k <- 2 # number of parameters
n <- 100 # number of observations
pZ <- 1 / 2 # probabilty Z = 1

# Generate parameters according to the constraints
ad <- gen_ad(k, sd_a, sd_d, rho_ad)
a <- ad$a
d <- ad$d

g <- a + d
g_0 <- a_0 + d_0

# Generate sample
Z <- rbinom(n, 1, pZ)
X <- matrix(rnorm(n * k), n, k)
V <- rnorm(n)

index <- a_0 + X %*% as.matrix(a, k, 1) + Z * (d_0 + X %*% as.matrix(d, k, 1))
P <- pnorm(index)
D <- as.numeric(index + V > 0) 

# Build a data frame
data <- data.frame(D = D, Z = Z)
data <- cbind(data, X)
names(data)[2 + 1:k] <- paste0("X_", 1:k)

# Generate interactions
interactions <- apply(X, 2, function(X, Z){return(Z * X)}, Z = Z)
colnames(interactions) <- paste("Z_X", 1:k, sep = "_")
data <- cbind(data, interactions)

##################### MAXIMUM RANK CORRELATION ESTIMATOR ######################

# Here theta = (a, d_0, d) (scaled to have length 1)
varnames <- c(paste("X", 1:k, sep = "_"), "Z", paste("Z_X", 1:k, sep = "_"))

# Get pilot estimates by (1) probit (2) logit (3) lpm
probit_est <- glm(as.formula(paste0("D ~ ", paste(varnames, collapse = " + "))),
                  family = binomial(link = "probit"), 
                  data = data)$coefficients
theta_probit <- probit_est[names(probit_est) %in% varnames]

logit_est <- glm(as.formula(paste0("D ~ ", paste(varnames, collapse = " + "))),
                 family = binomial(link = "logit"), 
                 data = data)$coefficients
theta_logit <- logit_est[names(logit_est) %in% varnames]

lpm_est <- lm(as.formula(paste0("D ~ ", paste(varnames, collapse = " + "))),
              data = data)$coefficients
theta_lpm <- lpm_est[names(lpm_est) %in% varnames]

# Run the Nelder-Mead-Head simplex from each pilot, take the maximum
candidates <- vector(mode = "list", length = 3)

candidates[[1]] <- optim(par = theta_probit, 
                        fn = S_theta, 
                        method = "Nelder-Mead", 
                        data = data, 
                        varnames = varnames, 
                        control = list("fnscale" = -1))

candidates[[2]] <- optim(par = theta_logit, 
                        fn = S_theta, 
                        method = "Nelder-Mead", 
                        data = data, 
                        varnames = varnames, 
                        control = list("fnscale" = -1))

candidates[[3]] <- optim(par = theta_lpm, 
                        fn = S_theta, 
                        method = "Nelder-Mead", 
                        data = data, 
                        varnames = varnames, 
                        control = list("fnscale" = -1))

# Choose the best performing value
i <- which.max(sapply(candidates, getElement, "value"))
theta_hat <- candidates[[i]]$par
theta_hat <- theta_hat / sqrt(sum(theta_hat ^ 2)) # scale to be length 1

########################### BOOTSTRAP STANDARD ERRORS ##########################

# Obtain standard errors by non-parametric bootstrap
temp <- boot(data = data, statistic = est_theta, varnames = varnames, R = 10)
