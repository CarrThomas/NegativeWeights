############################## CHOOSING PARAMETERS #############################

# Choose parameter values to attain certain targets 
# Specify the propensity score as
# p(X,Z) = \Phi(a_0 + X'a + Z(d_0 + X'd))
# With X ~ N(0, I), the parameters we need values for are a_0, d_0, and the 
# covariance matrix of X'a and X'd (sd_a, sd_d and rho_ad)

rm(list = ls())
gc()

#### Functions ####

kernel <- function(Xd, t_0, sd_d, sd_t, rho_td){
  
  # function to be integrated
  # consists of the density for Xd multiplied by either
  # 1) Pr(Xa + V < - a_0|Xd) -> set t_0 = a_0, sd_t = sd_a, rho_td = rho_ad
  # 2) Pr(Xg + V < - g_0|Xd) -> set t_0 = g_0, sd_t = sd_g, rho_td = rho_gd
  sd_temp <- sqrt(1 + (sd_t ^ 2) * (1 - rho_td ^ 2))
  bound <- t_0 + rho_td * sd_t * Xd / sd_d 
  
  return(pnorm(q = bound, mean = 0, sd = sd_temp) * dnorm(x = Xd, sd = sd_d))
  
}

# Calculates shares of always takers, never takers, defiers and compliers
get_shares <- function(a_0, d_0, g_0, sd_a, sd_d, sd_g, rho_ad, rho_dg){
  
  always_takers_1 <- integrate(kernel, lower = - d_0, upper = Inf,
                               t_0 = a_0,
                               sd_d = sd_d,
                               sd_t = sd_a,
                               rho_td = rho_ad)$value
  
  never_takers_1 <- 1 - pnorm(- d_0 / sd_d) -
    integrate(kernel, lower = - d_0, upper = Inf,
              t_0 = g_0,
              sd_d = sd_d,
              sd_t = sd_g,
              rho_td = rho_dg)$value
  
  compliers <- 1 - pnorm(q = -d_0, sd = sd_d) - always_takers_1 - never_takers_1
  
  always_takers_2 <- integrate(kernel, lower = - Inf, upper = - d_0,
                               t_0 = g_0,
                               sd_d = sd_d,
                               sd_t = sd_g,
                               rho_td = rho_dg)$value
  
  never_takers_2 <- pnorm(- d_0 / sd_d) -
    integrate(kernel, lower = - Inf, upper = - d_0,
              t_0 = a_0,
              sd_d = sd_d,
              sd_t = sd_a,
              rho_td = rho_ad)$value
  
  defiers <- pnorm(q = -d_0, sd = sd_d) - always_takers_2 - never_takers_2
  
  shares <- c(always_takers_1,
              always_takers_2,
              never_takers_1,
              never_takers_2,
              compliers,
              defiers)
  names(shares) <- c("always_takers_1",
                     "always_takers_2",
                     "never_takers_1",
                     "never_takers_2",
                     "compliers",
                     "defiers")
  
  return(shares)
  
}

#### Figuring Out Stuff ####

a_0 <- qnorm(0.45)
sd_a <- 3
sd_d <- 3
rho_ad <- -0.9

# # d_0 is chosen so that a_0 + d_0 = - a_0
# d_0 <- - 2 * a_0

# Calculate derived parameters 
g_0 <- a_0 + d_0
sd_g <- sqrt(sd_a ^ 2 + 2 * sd_a * sd_d * rho_ad + sd_d ^ 2)
rho_dg <- (sd_d + sd_a * rho_ad) / sd_g

shares <- get_shares(a_0, d_0, g_0, sd_a, sd_d, sd_g, rho_ad, rho_dg)
shares

# Generate a and d subject to these constraints
pZ <- 1 / 2
k <- 2

# construct orthonormal vectors
e_1 <- runif(k, -1, 1)
e_1 <- e_1 / sqrt(sum(e_1 ^ 2))
u <- runif(k, -1, 1)
# create a vector which is orthogonal to e_1
e_2 <- u - (sum(e_1 * u)) * e_1 / (sum(e_1 ^ 2))
e_2 <- e_2 / sqrt(sum(e_2 ^ 2))

# a is a scaled version of e_1
a <- sd_a * e_1

# d is a combination of e_1 and e_2
d <- sd_d * rho_ad * e_1 + sd_d * sqrt(1 - rho_ad ^ 2) * e_2

sd_a - sqrt(sum(a ^ 2))
sd_d - sqrt(sum(d ^ 2))
sd_a * sd_d * rho_ad - sum(a * d)

g <- a + d

# Generate a sample to check calculations
n <- 100000
X <- matrix(rnorm(n * k), n, k)

P_x0 <- pnorm(a_0 + X %*% as.matrix(a, k, 1))
P_x1 <- pnorm(g_0 + X %*% as.matrix(g, k, 1))
ind <- (P_x1 > P_x0)

shares[1]
sum(P_x0[ind]) / n

shares[2]
sum(P_x1[!ind]) / n

shares[3]
sum((1 - P_x1[ind])) / n

shares[4]
sum((1 - P_x0[!ind])) / n

shares[5]
sum((P_x1[ind] - P_x0[ind])) / n

shares[6]
sum((P_x0[!ind] - P_x1[!ind])) / n

# Treatment rates for Z = 0 and Z = 1 observations
pnorm(0, mean = - a_0, sd = sqrt(1 + sd_a ^ 2))
pnorm(a_0, mean = 0, sd = sqrt(1 + sd_a ^ 2))
pnorm(a_0 / sqrt(1 + sd_a ^ 2), mean = 0, sd = 1)
mean(P_x0)

pnorm(0, mean = - g_0, sd = sqrt(1 + sd_g ^ 2))
pnorm(g_0, mean = 0, sd = sqrt(1 + sd_g ^ 2))
pnorm(g_0 / sqrt(1 + sd_g ^ 2), mean = 0, sd = 1)
mean(P_x1)

# Expectation of the squared Z = 0 and Z = 1 treatment rates
pmvnorm(upper = c(0, 0), 
        mean = rep(-a_0, 2), 
        sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                      c(sd_a ^ 2, 1 + sd_a ^ 2)))
mean(P_x0 ^ 2)

pmvnorm(upper = c(0, 0), 
        mean = rep(-g_0, 2), 
        sigma = rbind(c(1 + sd_g ^ 2, sd_g ^ 2),
                      c(sd_g ^ 2, 1 + sd_g ^ 2)))
mean(P_x1 ^ 2)

# Variance of the Z = 0 and Z = 1 treatment rates
pmvnorm(upper = c(0, 0), 
        mean = rep(-a_0, 2), 
        sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                      c(sd_a ^ 2, 1 + sd_a ^ 2)),
        keepAttr = FALSE) - pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2)) ^ 2
var(P_x0)

pmvnorm(upper = c(0, 0), 
        mean = rep(-g_0, 2), 
        sigma = rbind(c(1 + sd_g ^ 2, sd_g ^ 2),
                      c(sd_g ^ 2, 1 + sd_g ^ 2)),
        keepAttr = FALSE) - pnorm(0, mean = -g_0, sd = sqrt(1 + sd_g ^ 2)) ^ 2
var(P_x1)

#### Tuning ####

# Targets are
# 1) The fraction of observations where the effect of the instrument is positive
# 2) The share of compliers
# 3) The share of defiers 
# 4) The overall treatment take-up rate (always-takers + compliers with Z = 1 + 
# defiers with Z = 0)
# 5) The standard deviation of the treatment rate with Z = 0
# The idea is to keep the total number of "switchers" (compliers + defiers) 
# fixed, but to vary the individual shares

# Calculate the targets
get_targets <- function(a_0, d_0, g_0, sd_a, sd_d, sd_g, rho_ad, rho_dg, pZ){
  
  # target 1: share of observations where X encourages treatment
  X_plus <- 1 - pnorm(0, mean = d_0, sd = sd_d)
  
  # target 2: share of compliers
  always_takers_1 <- integrate(kernel, lower = - d_0, upper = Inf,
                               t_0 = a_0,
                               sd_d = sd_d,
                               sd_t = sd_a,
                               rho_td = rho_ad)$value
  
  never_takers_1 <- 1 - pnorm(- d_0 / sd_d) -
    integrate(kernel, lower = - d_0, upper = Inf,
              t_0 = g_0,
              sd_d = sd_d,
              sd_t = sd_g,
              rho_td = rho_dg)$value
  
  compliers <- 1 - pnorm(q = -d_0, sd = sd_d) - always_takers_1 - never_takers_1
  
  # target 3: share of defiers
  always_takers_2 <- integrate(kernel, lower = - Inf, upper = - d_0,
                               t_0 = g_0,
                               sd_d = sd_d,
                               sd_t = sd_g,
                               rho_td = rho_dg)$value
  
  never_takers_2 <- pnorm(- d_0 / sd_d) -
    integrate(kernel, lower = - Inf, upper = - d_0,
              t_0 = a_0,
              sd_d = sd_d,
              sd_t = sd_a,
              rho_td = rho_ad)$value
  
  defiers <- pnorm(q = -d_0, sd = sd_d) - always_takers_2 - never_takers_2
  
  # target 4: overall treatment take-up rate
  pD <- always_takers_1 + always_takers_2 + pZ * compliers + (1 - pZ) * defiers
  
  # target 5: standard deviation of P(D = 1|X, Z = 0)
  sd_P_x0 <- sqrt(pmvnorm(upper = c(0, 0), 
                          mean = rep(-a_0, 2), 
                          sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                                        c(sd_a ^ 2, 1 + sd_a ^ 2)),
                          keepAttr = FALSE) - 
                    pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2)) ^ 2)
  
  targets <- c(X_plus,
               compliers, 
               defiers,
               pD,
               sd_P_x0)
  names(targets) <- c("X_plus", "compliers", "defiers", "pD", "sd_P_x0")
  
  return(targets)
  
}

get_targets(a_0, d_0, g_0, sd_a, sd_d, sd_g, rho_ad, rho_dg, pZ)

# Let theta = c(a_0, d_0, sd_a, sd_d, rho_ad) be our vector of parameter values
# to be calibrated. 
# Function to be optimised takes as inputs theta, pZ, and the target values
# Given theta, it calculates the derived parameters g_0, sd_g and rho_dg
# For these parameters, it calls the function `get_targets'
# Returns the distance between the output of the function and the target values

discrepancy <- function(theta, pZ, targets){
  
  a_0 <- theta[1]
  d_0 <- theta[2]
  sd_a <- theta[3]
  sd_d <- theta[4]
  rho_ad <- theta[5]
  
  g_0 <- a_0 + d_0
  sd_g <- sqrt(sd_a ^ 2 + 2 * sd_a * sd_d * rho_ad + sd_d ^ 2)
  rho_dg <- (sd_d + sd_a * rho_ad) / sd_g
  
  vals <- get_targets(a_0 = a_0,
                      d_0 = d_0,
                      g_0 = g_0,
                      sd_a = sd_a,
                      sd_d = sd_d,
                      sd_g = sd_g,
                      rho_ad = rho_ad,
                      rho_dg = rho_dg,
                      pZ = pZ)
  
  return(sqrt(sum((vals - targets) ^ 2)))
  
}

# initial targets
# 1) effect is positive 2 / 3 of the time
# 2) compliers 0.15, defiers 0.05
# 3) take up rate 0.55
# 4) standard deviation of the Z = 0 take up rate is 0.1

targets <- c(0.66, 0.15, 0.05, 0.55, 0.1)
theta <- c(a_0, d_0, sd_a, sd_d, rho_ad)

discrepancy(theta = theta, pZ = pZ, targets = targets)

# Initial guess: suppose there are no defiers. Then the only parameters to be 
# set are a_0, d_0 and sd_a. Let the targets be
# 1) compliers 0.2
# 2) take up rate 0.55
# 3) standard deviation of the Z = 0 take up rate is 0.1
# We proceed as follows. a_0 and sd_a are chosen to target
# 1) take up rate - P(Z = 1) * compliers
# 2) standard deviation of the Z = 0 take up rate
# Given a_0 and sd_a, d_0 can be backed out from the third condition

simple_discrepancy <- function(theta, targets){
  
  a_0 <- theta[1]
  sd_a <- theta[2]
  
  vals <- rep(0, 2)
  vals[1] <- pnorm(a_0 / sqrt(1 + sd_a ^ 2))
  
  vals[2] <- sqrt(pmvnorm(upper = rep(0, 2),
                          mean = rep(-a_0, 2),
                          sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                                        c(sd_a ^ 2, 1 + sd_a ^ 2)),
                          keepAttr = FALSE) -
                    vals[1] ^ 2)
  
  # return(sqrt(sum((targets - vals) ^ 2)))
  return(max(abs(targets - vals)))
  
}

theta <- c(a_0, sd_a)
targets <- rep(0, 2)
targets[1] <- 0.55 - pZ * 0.2 # share of always-takers
targets[2] <- 0.1 # standard deviation

simple_discrepancy(theta, targets)

# even simpler... given sd_a we can solve for a_0, and reduce this to a 
# univariate problem
# a_0 = sqrt(1 + sd_a ^ 2) * qnorm(targets[1])

AT <- 0.55 - pZ * 0.2
get_sd_Z0 <- function(AT, sd_a){
  
  mu <- qnorm(AT) * sqrt(1 + sd_a ^ 2)
  return(sqrt(pmvnorm(upper = rep(0, 2), 
                      mean = rep(-mu, 2),
                      sigma = rbind(c(1 + sd_a ^ 2, sd_a ^ 2),
                                    c(sd_a ^ 2, 1 + sd_a ^ 2)),
                      keepAttr = FALSE) - AT ^ 2))
  
}

temp <- data.frame(sd_a = seq(0, 2, length.out = 101))
temp$sd_Z0 <- sapply(temp$sd_a, get_sd_Z0, AT = AT)

temp %>% ggplot() + 
  geom_line(aes(x = sd_a, y = sd_Z0, col = "sd_Z0")) + 
  theme_bw() + 
  theme(legend.position = "null")

# choose the approximate minimiser
index <- which.min(abs(temp$sd_Z0 - targets[2]))

# refine
temp <- data.frame(sd_a = seq(temp$sd_a[index - 1], temp$sd_a[index + 1], length.out = 101))
temp$sd_Z0 <- sapply(temp$sd_a, get_sd_Z0, AT = AT)

index <- which.min(abs(temp$sd_Z0 - targets[2]))

sd_a <- temp$sd_a[index]
a_0 <- qnorm(AT) * sqrt(1 + sd_a ^ 2)
d_0 <- sqrt(1 + sd_a ^ 2) * qnorm(AT + 0.2) - a_0
g_0 <- a_0 + d_0

# Check
k <- 2
a <- runif(k, -1, 1)
a <- sd_a * a / sqrt(sum(a ^ 2))

# Generate a sample to check calculations
n <- 100000
X <- matrix(rnorm(n * k), n, k)

P_x0 <- pnorm(a_0 + X %*% as.matrix(a, k, 1))
P_x1 <- pnorm(g_0 + X %*% as.matrix(a, k, 1))

sd(P_x0)
targets[2]

mean(P_x0)
AT

mean(P_x0) + pZ * mean(P_x1 - P_x0)
pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2)) + 
  pZ * (pnorm(0, mean = -g_0, sd = sqrt(1 + sd_a ^ 2)) - 
          pnorm(0, mean = -a_0, sd = sqrt(1 + sd_a ^ 2))) 

theta_no_def <- c(a_0, d_0, sd_a)
names(theta_no_def) <- c("a_0", "d_0", "sd_a")

#### Main Calibration ####

# Start with the no defiers case, gradually increase the proportion of defiers
# With each new iteration, let the initial guess be the solution from the 
# previous iteration.

theta <- matrix(0, nrow = 5, ncol = 10)
row.names(theta) <- c("a_0", "d_0", "sd_a", "sd_d", "rho_ad")
colnames(theta) <- paste0("defiers = ", seq(0, 0.09, length.out = 10))
theta[1:3, 1] <- theta_no_def
theta[4:5, 1] <- 0.001

# Targets are
# 1) The fraction of observations where the effect of the instrument is positive
# (starts close to 1 and declines)
# 2) The share of compliers
# 3) The share of defiers (0.2 - compliers)
# 4) overall treatment take-up rate 0.55
# 5) The standard deviation of the treatment rate with Z = 0 is 0.1

defiers <- seq(0, 0.09, length.out = 10)
compliers <- 0.2 - defiers
X_plus <- 1 - seq(0, 0.45, length.out = 10)

runs <- 10
for (i in 2:runs){
  
  targets <- c(X_plus[i], compliers[i], defiers[i], 0.55, 0.1) 
  guess <- as.numeric(theta[, i - 1])
  discrepancy(guess, pZ, targets)
  
}

temp <- optim(par = guess, 
              fn = discrepancy, 
              method = "Nelder-Mead", 
              pZ = pZ, 
              targets = targets)$par

temp <- optim(par = temp, 
              fn = discrepancy, 
              method = "Nelder-Mead", 
              pZ = pZ, 
              targets = targets)$par

temp_guess <- guess
temp_guess[4:5] <- c(0.1, -0.8)
temp <- optim(par = temp_guess, 
              fn = discrepancy, 
              method = "Nelder-Mead", 
              pZ = pZ, 
              targets = targets)$par



a_0 <- temp[1]
d_0 <- temp[2]
sd_a <- temp[3]
sd_d <- temp[4]
rho_ad <- temp[5]
g_0 <- a_0 + d_0
sd_g <- sqrt(sd_a ^ 2 + 2 * sd_a * sd_d * rho_ad + sd_d ^ 2)
rho_dg <- (sd_d + sd_a * rho_ad) / sd_g

get_targets(a_0 = a_0,
            d_0 = d_0,
            g_0 = g_0,
            sd_a = sd_a,
            sd_d = sd_d,
            sd_g = sd_g,
            rho_ad = rho_ad,
            rho_dg = rho_dg,
            pZ = pZ)
