#################################### SUMMARY ################################### 

# This script calibrates the Monte Carlo process to match the Angrist and 
# Evans data

# The first stage is assumed to be
# D = 1(z * (d0 + x'd)  >=  V)
# where the marginal distribution of V is N(0, 1), so the propensity score is
# p(x,z) = Phi(z * (d0 + x'd))
# d0 + x'd controls the sign of the relationship between Z and D. 

# The second stage is 
# Y(0) = U(0)
# Y(1) = B_0 + x'B + U(1)
# Y = Y(0) + D(Y(1) - Y(0))
# where (V, U(1), U(0)) are jointly normal. 

################################# FRONT MATTER #################################

rm(list = ls())
gc()

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(boot)
library(xtable)

path <- getwd()

################################### SETTINGS ###################################

parameters <- vector(mode = "list", length = 1)
# parameters[[1]] <- list()
# parameters[[2]] <- list()

AE_moments <- readRDS(file = paste0(path, "/Monte_Carlo/AE_moments.rds"))
AE_moments$tsls_est <- AE_moments$tsls_est / 35

targets <- vector(mode = "list", length = 1)
# targets are
# 1) Pr(D = 1|Z = 0)
# 2) TSLS denominator (compliers - defiers)
# 3) TSLS hours estimate 
# 4) share of defiers
# 5) True tau
targets$PD0 <- AE_moments$PD0
targets$net_comp <- AE_moments$tsls_denom / AE_moments$Z_var
targets$tsls_est <- AE_moments$tsls_est
targets$defiers <- 0.05
targets$compliers <- targets$defiers + targets$net_comp
targets$tau_r <- targets$tsls_est + 0.1

# targets[[1]] <- list(compliers = 0.11,
#                      defiers = 0.09,
#                      tau_B = -0.05,
#                      tau_U = -0.05)
# targets[[2]] <- list(compliers = 0.06,
#                      defiers = 0.04,
#                      tau_B = -0.05,
#                      tau_U = -0.05)

e_tol <- 0.00000001 # error tolerance
# search bounds for use in calibration block 1 
bounds <- list(top =  c(0.2, 0.08), bottom = c(0.1, 0.01))

############################## CALIBRATION BLOCK 1 #############################

# Choose a0 to attain the desired Pr(D = 1|Z = 0)
a0 <- qnorm(targets$PD0)
parameters$a0 <- a0

# Choose d0 and sqrt(d'd) (sd_d) to attain the desired shares of compliers and 
# defiers

# integrate this kernel from -d0 to Inf to get the share of compliers
# integrate -1 times this kernel from -Inf to -d0 to get the share of defiers
kernel <- function(z, a0, d0, sd_d)(
  
  return((pnorm(a0 + d0 + z) - pnorm(a0)) * dnorm(z, mean = 0, sd = sd_d))
  
)

# Loop to solve for sd_d which meets the compliers target
get_sd_d <- function(a0, d0, e_tol, targets){
  
  # initial guess
  n_points <- 101
  sd_d <- seq(0.01, 1, length.out = n_points)
  compliers <- rep(0, n_points)
  
  for (i in 1:n_points){
    
    compliers[i] <- integrate(kernel, lower = -d0, upper = Inf,
                              a0 = a0, d0 = d0, sd_d = sd_d[i])$value
    
  }
  
  index <- which.min(abs(compliers - targets$compliers))
  upper <- sd_d[min(n_points, index + 1)]
  lower <- sd_d[max(1, index - 1)]
  error <- abs(compliers[index] - targets$compliers)
  
  # quick loop to improve
  while(error > e_tol){
    
    sd_d <- seq(lower, upper, length.out = n_points)
    for (i in 1:n_points){
      
      compliers[i] <- integrate(kernel, lower = -d0, upper = Inf, 
                                a0 = a0, d0 = d0, sd_d = sd_d[i])$value
      
    }
    
    index <- which.min(abs(compliers - targets$compliers))
    upper <- sd_d[min(n_points, index + 1)]
    lower <- sd_d[max(1, index - 1)]
    error <- abs(compliers[index] - targets$compliers)
    
  }
  
  # final answer
  return(list(sd_d = sd_d[index], compliers = compliers[index]))
  
}

# Split the difference type algorithm to find d0 and sd_d 
# Assume that desired compliers share > desired defiers share
# Loop is as follows:
# * Take d0 as fixed
# * Choose sd_d to attain the desired compliers shares
# * For this (d0, sd_d) pair, evaluate the defier shares
# * If the defiers share is too small, d0 is too large -> revise upper bound
# * If the defiers share is too large, d0 is too small -> revise lower bound

for (i in 1:1){
  
  # upper and lower bounds of the search region
  top <- bounds$top[i]
  bottom <- bounds$bottom[i]
  
  count <- 0
  limit <- 1000 # maximum number of iterations
  
  while (count < limit){
    
    mid <- (top + bottom) / 2
    
    # solve for sd_d that matches compliers
    sd_d <- get_sd_d(a0 = a0, d0 = mid, e_tol = e_tol, targets = targets)$sd_d
    
    # check defiers
    defiers <- - integrate(kernel, lower = - Inf, upper = -mid, 
                           a0 = a0, d0 = mid, sd_d = sd_d)$value
    
    if (abs(defiers - targets$defiers) < e_tol){
      print("converged")
      break
    }else{
      if (defiers > targets$defiers){
        bottom <- mid
      }else{
        top <- mid
      }
    }
    
    count <- count + 1
    
  }
  
  parameters$d0 <- mid
  parameters$sd_d <- sd_d
  
}

############################## VALIDATE PARAMETERS #############################

# Check if these parameters give the correct shares for complier and defiers

K <- 3 # number of X variables
n <- 100000

d <- runif(K, min = -1, max = 1)
d <- d * parameters$sd_d / sqrt(sum(d * d))

X <- matrix(rnorm(n * K), nrow = n, ncol = K)
Xd <- X %*% matrix(d, nrow = K, ncol = 1)
Z <- rbinom(n, 1, 1 / 2)
V <- rnorm(n)
D <- as.numeric((parameters$a0 + Z * (parameters$d0 + Xd)) > V)

comparisons <- matrix(0, 4, 2)

comparisons[, 1] <- c(targets$PD0,
                      targets$compliers + targets$defiers,
                      targets$compliers,
                      targets$defiers)
comparisons[1, 2] <- sum(D * (1 - Z)) / sum(1 - Z) 
comparisons[2, 2] <- mean(abs(pnorm(parameters$a0 + parameters$d0 + Xd) - pnorm(parameters$a0)))
comparisons[3, 2] <- mean(((parameters$d0 + Xd) > 0) * 
                            abs(pnorm(parameters$a0 + parameters$d0 + Xd) - 
                                  pnorm(parameters$a0)))
comparisons[4, 2] <- mean(((parameters$d0 + Xd) < 0) * 
                            abs(pnorm(parameters$a0 + parameters$d0 + Xd) - 
                                  pnorm(parameters$a0)))

rownames(comparisons) <- c("PD0", "switchers", "compliers", "defiers")

comparisons

############################## CALIBRATION BLOCK 2 #############################

# Given a0, d0, and sd_d, solve for the parameters 
# We have the formulas
# tau = B0 + (B'd / d'd) * term_one + (rho_1v - rho_0v) * term_two
# tau_tsls = B0 + (B'd / d'd) * term_three + (rho_1v - rho_0v) * term_four
# where term_one/three captures the effect of XB and term two/four the effect 
# of U. We can simply evaluate these integrals to back out values for the
# multiplying constants

# kernel_one is s(X) * D * (Z - E[Z|X])
kernel_one <- function(z, a0, d0, sd_d){
  
  return((2 * (d0 + z >= 0) - 1) * (pnorm(a0 + d0 + z) - pnorm(a0)) * 
           dnorm(z, mean = 0, sd = sd_d))
  
}

# kernel_two is s(X) * D * (Z - E[Z|X]) * Xd 
kernel_two <- function(z, a0, d0, sd_d){
  
  return((2 * (d0 + z >= 0) - 1) * (pnorm(a0 + d0 + z) - pnorm(a0)) * 
           z * dnorm(z, mean = 0, sd = sd_d))
  
}

# kernel_three is s(X) * ((1 - D) * U0 + D * U1) * (Z - E[Z|X])
kernel_three <- function(z, a0, d0, sd_d){
  
  return((2 * (d0 + z >= 0) - 1) * (dnorm(a0) - dnorm(a0 + d0 + z)) * 
           dnorm(z, mean = 0, sd = sd_d))
  
}

# kernel_four is D * (Z - E[Z|X])
kernel_four <- function(z, a0, d0, sd_d){
  
  return((pnorm(a0 + d0 + z) - pnorm(a0)) * dnorm(z, mean = 0, sd = sd_d))
  
}

# kernel_five is D * (Z - E[Z|X]) * Xd 
kernel_five <- function(z, a0, d0, sd_d){
  
  return((pnorm(a0 + d0 + z) - pnorm(a0)) * z * dnorm(z, mean = 0, sd = sd_d))
  
}

# kernel_six is ((1 - D) * U0 + D * U1) * (Z - E[Z|X])
kernel_six <- function(z, a0, d0, sd_d){
  
  return((dnorm(a0) - dnorm(a0 + d0 + z)) * dnorm(z, mean = 0, sd = sd_d))
  
}


reorder_denom <- integrate(kernel_one, lower = -Inf, upper = Inf, 
                           a0 = parameters$a0,
                           d0 = parameters$d0,
                           sd_d = parameters$sd_d)$value
  
term_one <- integrate(kernel_two, lower = -Inf, upper = Inf, 
                      a0 = parameters$a0,
                      d0 = parameters$d0,
                      sd_d = parameters$sd_d)$value / reorder_denom

term_two <- integrate(kernel_three, lower = -Inf, upper = Inf, 
                      a0 = parameters$a0,
                      d0 = parameters$d0,
                      sd_d = parameters$sd_d)$value / reorder_denom

tsls_denom <- integrate(kernel_four, lower = -Inf, upper = Inf, 
                        a0 = parameters$a0,
                        d0 = parameters$d0,
                        sd_d = parameters$sd_d)$value

term_three <- integrate(kernel_five, lower = -Inf, upper = Inf, 
                        a0 = parameters$a0,
                        d0 = parameters$d0,
                        sd_d = parameters$sd_d)$value / tsls_denom

term_four <- integrate(kernel_six, lower = -Inf, upper = Inf, 
                       a0 = parameters$a0,
                       d0 = parameters$d0,
                       sd_d = parameters$sd_d)$value / tsls_denom


# System of equations:
# tau_r = chi_two * term_one + chi_three * term_two
# tau_tsls = chi_two * term_three + chi_three * term_four
# with
# chi_two = ||B|| * B'd / ||d||
# chi_three = sd_1 rho_1V - sd_0 rho_0V

chi <- solve(matrix(c(term_one, term_two, term_three, term_four), 
                    ncol = 2,
                    nrow = 2,
                    byrow = T)) %*% 
  matrix(c(targets$tau_r, targets$tsls_est))

# dealing with chi_two
rho_Bd <- sign(chi[1]) * 0.5
sd_B <- chi[1] * parameters$sd_d / rho_Bd

# dealing with chi_three 
# fix correlations and adjust standard deviations to match the target
sigma <- diag(1, 3)
sigma[2, 1] <- -0.25
sigma[3, 1] <- -0.5
sigma[3, 2] <- 0.5
sigma[1, 2] <- -0.25
sigma[1, 3] <- -0.5
sigma[2, 3] <- 0.5

# assume sd_U is common for U(1) and U(0). Back it out.
sd_U <- abs(chi[2] / (sigma[3, 1] - sigma[3, 2]))
sigma[3, 1] <- - (sigma[3, 1] * sign(chi[2])) 
sigma[3, 2] <- - (sigma[3, 2] * sign(chi[2])) 

# multiply by sd_U so that sigma is a covariance matrix
sigma[, 1:2] <- sigma[, 1:2] * sd_U
sigma[1:2, ] <- sigma[1:2, ] * sd_U

# # Randomly draw Sigma
# A <- matrix(rnorm(9), nrow = 3, ncol = 3)
# sigma <- t(A) %*% A
# sd_V <- sqrt(sigma[3, 3])
# sigma[3, ] <- sigma[3, ] / sd_V
# sigma[, 3] <- sigma[, 3] / sd_V
# chi_2 <- sigma[3, 1] - sigma[3, 2]
parameters$sd_B <- sd_B
parameters$rho_Bd <- rho_Bd
parameters$sigma <- sigma
  
############################## VALIDATE PARAMETERS #############################

# Check if parameters give the correct asymptotic estimand
K <- 3
n <- 100000
runs <- 100

# other parameters
pZ <- 1 / 2
Var_Z <- pZ * (1 - pZ)
B0 <- 0

analytical <- rep(0, 6)
names(analytical) <- c("denom", "num_1", "num_2", "num_3", "num", "tau")

analytical[1] <- Var_Z * 
  integrate(kernel_one, lower = -Inf, upper = Inf,
            a0 = parameters$a0,
            d0 = parameters$d0, 
            sd_d = parameters$sd_d)$value
analytical[2] <- Var_Z * B0 * 
  integrate(kernel_one, lower = -Inf, upper = Inf, 
            a0 = parameters$a0,
            d0 = parameters$d0, 
            sd_d = parameters$sd_d)$value
analytical[3] <- Var_Z * 
  (parameters$sd_B * parameters$rho_Bd / parameters$sd_d) *
  integrate(kernel_two, lower = -Inf, upper = Inf, 
            a0 = parameters$a0,
            d0 = parameters$d0, 
            sd_d = parameters$sd_d)$value 
analytical[4] <- Var_Z * 
  (parameters$sigma[3, 1] - parameters$sigma[3, 2]) * 
  integrate(kernel_three, lower = -Inf, upper = Inf, 
            a0 = parameters$a0,
            d0 = parameters$d0, 
            sd_d = parameters$sd_d)$value    
analytical[5] <- sum(analytical[2:4])
analytical[6] <- analytical[5] / analytical[1]

tsls_analytical <- rep(0, 6)
names(tsls_analytical) <- c("denom", "num_1", "num_2", "num_3", "num", "tau")

tsls_analytical[1] <- Var_Z * 
  integrate(kernel_four, lower = -Inf, upper = Inf,
            a0 = parameters$a0,
            d0 = parameters$d0, 
            sd_d = parameters$sd_d)$value
tsls_analytical[2] <- Var_Z * B0 * 
  integrate(kernel_four, lower = -Inf, upper = Inf, 
            a0 = parameters$a0,
            d0 = parameters$d0, 
            sd_d = parameters$sd_d)$value
tsls_analytical[3] <- Var_Z * 
  (parameters$sd_B * parameters$rho_Bd / parameters$sd_d) *
  integrate(kernel_five, lower = -Inf, upper = Inf, 
            a0 = parameters$a0,
            d0 = parameters$d0, 
            sd_d = parameters$sd_d)$value 
tsls_analytical[4] <- Var_Z * 
  (parameters$sigma[3, 1] - parameters$sigma[3, 2]) * 
  integrate(kernel_six, lower = -Inf, upper = Inf, 
            a0 = parameters$a0,
            d0 = parameters$d0, 
            sd_d = parameters$sd_d)$value    
tsls_analytical[5] <- sum(tsls_analytical[2:4])
tsls_analytical[6] <- tsls_analytical[5] / tsls_analytical[1]

data <- matrix(0, runs, 6)
colnames(data) <- c("denom", "num_1", "num_2", "num_3", "num", "tau")

tsls_data <- matrix(0, runs, 6)
colnames(tsls_data) <- c("denom", "num_1", "num_2", "num_3", "num", "tau")

for (j in 1:runs){
  
  # Generate B and d subject to constraints
  # Construct orthogonal vectors 
  e_1 <- runif(K, min = - 1, max = 1)
  u <- runif(K, min = - 1, max = 1)
  e_2 <- u - sum(e_1 * u) * e_1 / (sum(e_1 ^ 2))
  e_1 <- e_1 / sqrt(sum(e_1 ^ 2))
  e_2 <- e_2 / sqrt(sum(e_2 ^ 2))
  
  # Use e_1 and e_2 to build d and B
  d <- e_1 * parameters$sd_d
  B <- parameters$sd_B * 
    (parameters$rho_Bd * e_1 + sqrt(1 - parameters$rho_Bd ^ 2) * e_2)
  
  X <- matrix(rnorm(n * K), nrow = n, ncol = K)
  Xd <- X %*% as.matrix(d, nrow = K, ncol = 1)
  XB <- X %*% as.matrix(B, nrow = K, ncol = 1)
  
  s <- 2 * as.numeric(parameters$d0 + Xd >= 0) - 1
  
  Z <- rbinom(n, 1, pZ)
  Z_tilde <- Z - mean(Z)
  
  U <- rmvnorm(n, mean = rep(0, K), sigma = parameters$sigma)
  colnames(U) <- c("U1", "U0", "V")
  U <- data.frame(U)
  
  D <- as.numeric((parameters$a0 + Z * (parameters$d0 + Xd)) >= U$V)
  
  Y <- (1 - D) * U$U0 + D * (B0 + XB + U$U1)
  
  #### Reordered IV ####
  
  # Denominator
  data[j, 1] <- mean(s * D * Z_tilde)
  
  # Numerator
  data[j, 2] <- B0 * mean(s * D * Z_tilde)
  data[j, 3] <- mean(s * D * XB * Z_tilde)
  data[j, 4] <- mean(s * ((1 - D) * U$U0 + D * U$U1) * Z_tilde)
  data[j, 5] <- sum(data[j, 2:4])
  
  # Estimand
  data[j, 6] <- data[j, 5] / data[j, 1]
  
  #### TSLS ####
  
  # Denominator
  tsls_data[j, 1] <- mean(D * Z_tilde)
  
  # Numerator
  tsls_data[j, 2] <- B0 * mean(D * Z_tilde)
  tsls_data[j, 3] <- mean(D * XB * Z_tilde)
  tsls_data[j, 4] <- mean(((1 - D) * U$U0 + D * U$U1) * Z_tilde)
  tsls_data[j, 5] <- sum(tsls_data[j, 2:4])
  
  # Estimand
  tsls_data[j, 6] <- tsls_data[j, 5] / tsls_data[j, 1]
  
}

results <- list(analytical = analytical, data = data,
                tsls_analytical = tsls_analytical,
                tsls_data = tsls_data)

#### Reordered IV Comparison ####

# denominator
plot(results$data[, 1], rep(results$analytical[1], runs))
mean(results$data[, 1])
results$analytical[1]

# numerator
plot(results$data[, 5], rep(results$analytical[5], runs))
mean(results$data[, 5])
results$analytical[5]

# numerator part 1
plot(results$data[, 2], rep(results$analytical[2], runs))
mean(results$data[, 2])
results$analytical[2]

# numerator part 2
plot(results$data[, 3], rep(results$analytical[3], runs))
mean(results$data[, 3])
results$analytical[3]

# numerator part 3
plot(results$data[, 4], rep(results$analytical[4], runs))
mean(results$data[, 4])
results$analytical[4]

# tau
plot(results$data[, 6], rep(results$analytical[6], runs))
mean(results$data[, 6])
results$analytical[6]
  
#### TSLS Comparison ####

# denominator
plot(results$tsls_data[, 1], rep(results$tsls_analytical[1], runs))
mean(results$tsls_data[, 1])
results$tsls_analytical[1]

# numerator
plot(results$tsls_data[, 5], rep(results$tsls_analytical[5], runs))
mean(results$tsls_data[, 5])
results$tsls_analytical[5]

# numerator part 1
plot(results$tsls_data[, 2], rep(results$tsls_analytical[2], runs))
mean(results$tsls_data[, 2])
results$tsls_analytical[2]

# numerator part 2
plot(results$tsls_data[, 3], rep(results$tsls_analytical[3], runs))
mean(results$tsls_data[, 3])
results$tsls_analytical[3]

# numerator part 3
plot(results$tsls_data[, 4], rep(results$tsls_analytical[4], runs))
mean(results$tsls_data[, 4])
results$tsls_analytical[4]

# tau
plot(results$tsls_data[, 6], rep(results$tsls_analytical[6], runs))
mean(results$tsls_data[, 6])
results$tsls_analytical[6]


############################### SWITCHER DENSITY ###############################

# Given these parameters, plot the density of switchers (compliers and defiers)
n_points <- 1001
C <- seq(-pnorm(parameters$a0), 1 - pnorm(parameters$a0), length.out = n_points)

filenames <- "switcher_dens_parameters_ae"

dens <- rep(0, n_points)
dist <- rep(0, n_points)

for (j in 1:n_points){
  
  dist[j] <- pnorm(qnorm(pnorm(parameters$a0) + C[j]) - parameters$d0 - parameters$a0, 
                   mean = 0, sd = parameters$sd_d) 
  dens[j] <-  dnorm(qnorm(pnorm(parameters$a0) + C[j]) - parameters$d0 - parameters$a0, 
                    mean = 0, sd = parameters$sd_d) / 
    dnorm(qnorm(pnorm(parameters$a0) + C[j]) - parameters$d0)
  
}

plot(C, dist)
plot(C, dens)

# verify by simulation
K <- 3
n <- 100000
d <- runif(K, min = - 1, max  = 1)
d <- d * parameters$sd_d / sqrt(sum(d ^ 2))
X <- matrix(rnorm(n * K), nrow = n, ncol = K)
Xd <- X %*% as.matrix(d, nrow = K, ncol = 1)
Delta <- pnorm(Xd + parameters$d0 + parameters$a0) - pnorm(parameters$a0)

emp_dist <- ecdf(Delta)
data_dist <- rep(0, n_points)

for (j in 1:n_points){
  
  data_dist[j] <- emp_dist(C[j])
  
}

plot_data <- data.frame(C = C, data = data_dist, analytical = dist)
plot_data %>% ggplot() + 
  geom_line(aes(x = C, y = data, col = "data")) + 
  geom_line(aes(x = C, y = analytical, col = "analytical")) + 
  theme_bw() + 
  theme(legend.title = element_blank())

data_dens <- density(Delta)
plot(data_dens)
lines(C, dens, lty = 2)

dens[1] <- 0
dens[n_points] <- 0

# Nice plot of the density of switchers
pdf(file = paste0(filenames[i], ".pdf"), width = 8.98, height = 5.60)
plot_data <- data.frame(C = C, dens = dens)
print(plot_data %>% ggplot() + 
        geom_line(aes(x = C, y = dens, col = "density")) + 
        geom_vline(aes(xintercept = pnorm(parameters$d0 + parameters$a0) - pnorm(parameters$a0), col = "d0"),
                   linetype = "dotted") + 
        theme_bw() + 
        scale_color_manual(values = c("density" = "#F8766D",
                                      "d0" = "black")) + 
        labs(x = "p(x,1) - p(x,0)", y = "") +
        theme(legend.position = "null"))
dev.off()
  
################################## SAVE OUTPUT #################################

calibration <- list(targets = targets, 
                    parameters = parameters)
saveRDS(calibration, 
        file = paste0(path, "/Monte_Carlo/ae_calibration.rds"))
