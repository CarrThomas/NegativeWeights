#################################### SUMMARY ################################### 

# This scipt performs Monte Carlo for corrected TSLS estimator.

# The first stage is assumed to be
# D = 1(a_0 + x'a + z(d_0 + x'd) + V > 0)
# where the marginal distribution of V is N(0, 1), so the propensity score is
# p(x,z) = Phi(a_0 + x'a + z(d_0 + x'd))
# In this formulation, the parameters d_0 and d control the sign of the 
# relationship between Z and D. 

# We use the maximum rank correlation estimator to obtain estimates of the first
# stage parameters. This estimator does not estimate a_0, and the remaining 
# parameters are only estimated up to scale (adding a constant to all terms 
# does not perturb ranks, and multiplying all non-constant terms by a positive
# constant would also not affect ranks)

# The second stage is assumed to be
# Y(0) = x'B + U(0)
# Y(1) = zeta_0 + x'B + U(1)
# Y = Y(0) + D(Y(1) - Y(0))
# where zeta_0 is the ATE and (V, U(1), U(0)) are jointly normal. 

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
k <- 3 # number of parameters
N <- c(100, 500, 1000, 2000) # number of observations
pZ <- 1 / 2 # probabilty Z = 1

# Monte Carlo settings
runs <- 1000

############################### MONTE CARLO LOOP ###############################

results <- vector(mode = "list", length = length(N))

for (j in 1:length(N)){
  
  tsls <- rep(0, runs)
  corrected <- rep(0, runs)
  oracle <- rep(0, runs)
  probit <- rep(0, runs)
  logit <- rep(0, runs)
  lpm <- rep(0, runs)
  diagnostics <- list(s_cor = rep(0, runs),
                      theta_hat = matrix(0, nrow = 1 + (2 * k), ncol = runs),
                      theta = matrix(0, nrow = 1 + (2 * k), ncol = runs))
  
  n <- N[j]
  
  for (r in 1:runs){
    
    #### GENERATE DATA ####
    
    # Generate parameters according to the constraints
    ad <- gen_ad(k, sd_a, sd_d, rho_ad)
    a <- ad$a
    d <- ad$d
    
    g <- a + d
    g_0 <- a_0 + d_0
    
    B <- runif(k, min = -1, max = 1)
    B <- B * sd_B / sqrt(sum(B ^ 2))
    
    # Generate sample
    Z <- rbinom(n, 1, pZ)
    X <- matrix(rnorm(n * k), n, k)
    V <- rnorm(n)
    
    index <- a_0 + X %*% as.matrix(a, k, 1) + Z * (d_0 + X %*% as.matrix(d, k, 1))
    P <- pnorm(index)
    D <- as.numeric(index + V > 0) 
    
    U <- rmvnorm(n, mean = rep(0, k), sigma = sigma)
    colnames(U) <- c("V", "U1", "U0")
    U <- data.frame(U)
    
    Y_0 <- B_0 + X %*% as.matrix(B, k, 1) + U$U0
    Y_1 <- B_0 + zeta_0 + X %*% as.matrix(B, k, 1) + U$U1
    Y <- Y_0 + D * (Y_1 - Y_0)
    
    # Build a data frame
    data <- data.frame(D = D, Z = Z, Y = Y)
    data <- cbind(data, X)
    names(data) <- c("D", "Z", "Y", paste0("X_", 1:k))
    
    # Generate interactions
    interactions <- apply(X, 2, function(X, Z){return(Z * X)}, Z = Z)
    colnames(interactions) <- paste("Z_X", 1:k, sep = "_")
    data <- cbind(data, interactions)
    
    #### MRC ESTIMATOR ####
    
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
    
    # Run the Nelder-Mead-Head simplex from each pilot (probit, logit, lpm), 
    # random perturbations of each pilot, and random points, take the maximum
    
    # 3 pilot points
    n_perturb <- 5 # n_perturb random perturbations of each pilot point
    n_random <- 20 # n_random randomly generated starting points
    candidates <- vector(mode = "list", length = 3 + 3 * n_perturb + n_random)
    
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
    
    # perturbations of theta probit
    theta_probit <- theta_probit / sqrt(sum(theta_probit ^ 2))
    for (p in 1:n_perturb){
      probit_perturbed <- theta_probit * 
        sample(c(0.9, 1.1), length(theta_probit), replace = TRUE)
      candidates[[3 + p]] <- optim(par = probit_perturbed, 
                                   fn = S_theta, 
                                   method = "Nelder-Mead", 
                                   data = data, 
                                   varnames = varnames, 
                                   control = list("fnscale" = -1))
    }
    
    # perturbations of theta logit
    theta_logit <- theta_logit / sqrt(sum(theta_logit ^ 2))
    for (p in 1:n_perturb){
      
      logit_perturbed <- theta_logit * 
        sample(c(0.9, 1.1), length(theta_logit), replace = TRUE)
      
      candidates[[3 + n_perturb + p]] <- optim(par = logit_perturbed, 
                                               fn = S_theta, 
                                               method = "Nelder-Mead", 
                                               data = data, 
                                               varnames = varnames, 
                                               control = list("fnscale" = -1))
    }
    
    # perturbations of theta lpm
    theta_lpm <- theta_lpm / sqrt(sum(theta_lpm ^ 2))
    for (p in 1:n_perturb){
      
      lpm_perturbed <- theta_lpm * 
        sample(c(0.9, 1.1), length(theta_lpm), replace = TRUE)
      
      candidates[[3 + (2 * n_perturb) + p]] <- 
        optim(par = lpm_perturbed, 
              fn = S_theta, 
              method = "Nelder-Mead", 
              data = data, 
              varnames = varnames, 
              control = list("fnscale" = -1))
    }
    
    # randomly generated starting points
    for (p in 1:n_random){
      
      theta_random <- runif(length(theta_probit), -1, 1) 
      theta_random <- theta_random / sqrt(sum(theta_random ^ 2))
      
      candidates[[3 + (3 * n_perturb) + p]] <- 
        optim(par = theta_random, 
              fn = S_theta, 
              method = "Nelder-Mead", 
              data = data, 
              varnames = varnames, 
              control = list("fnscale" = -1))
      
    }
    
    # Choose the best performing value
    i <- which.max(sapply(candidates, getElement, "value"))
    theta_hat <- candidates[[i]]$par
    names(theta_hat) <- names(candidates[[1]]$par)
    theta_hat <- theta_hat / sqrt(sum(theta_hat ^ 2)) # scale to be length 1
    
    s_hat <- sign(theta_hat["Z"] + X %*% 
                    matrix(theta_hat[colnames(interactions)], k, 1))
    s_star <- sign(d_0 + as.matrix(X) %*% as.matrix(d, k, 1)) # oracle
    s_probit <- sign(probit_est["Z"] + X %*% 
                       matrix(probit_est[colnames(interactions)], k ,1))
    s_logit <- sign(logit_est["Z"] + X %*% 
                      matrix(logit_est[colnames(interactions)], k ,1))
    s_lpm <- sign(lpm_est["Z"] + X %*% 
                    matrix(lpm_est[colnames(interactions)], k ,1))
    
    # Estimates
    Z_hat <- lm(as.formula(paste0("Z ~ ", paste(paste0("X_", 1:k), collapse = " + "))),
                data = data)$fitted.values
    Z_tilde <- Z - Z_hat
    # corrected TSLS
    oracle[r] <- sum(Y * s_star * Z_tilde) / sum(D * s_star * Z_tilde)
    corrected[r] <- sum(Y * s_hat * Z_tilde) / sum(D * s_hat * Z_tilde)
    probit[r] <- sum(Y * s_probit * Z_tilde) / sum(D * s_probit * Z_tilde)
    logit[r] <- sum(Y * s_logit * Z_tilde) / sum(D * s_logit * Z_tilde)
    lpm[r] <- sum(Y * s_lpm * Z_tilde) / sum(D * s_lpm * Z_tilde)
    tsls[r] <- cov(Y, Z_tilde) / cov(D, Z_tilde)
    
    # Diagnostics
    diagnostics$s_cor[r] <- cor(s_hat, s_star)
    diagnostics$theta[ , r] <- c(a, d_0, d) / sqrt(sum(c(a, d_0, d) ^ 2))
    diagnostics$theta_hat[, r] <- theta_hat 
    
  }
  
  run_results <- list(oracle = oracle, 
                      corrected = corrected,
                      tsls = tsls, 
                      probit = probit,
                      logit = logit,
                      lpm = lpm,
                      diagnostics = diagnostics)
  
  results[[j]] <- run_results
  
}

names(results) <- N

saveRDS(results, file = paste0(path, "/Monte_Carlo/results.rds"))

############################# VISUALISING RESULTS ############################## 

sims <- length(N)

for (i in 1:sims){
  
  x_lims <- rep(0, 2)
  x_lims[1] <- quantile(results[[i]]$tsls[!is.na(results[[1]]$corrected)], 0.05)
  x_lims[2] <- quantile(results[[i]]$tsls[!is.na(results[[1]]$corrected)], 0.95)
  plot_points <- list()
  
  plot_data <- data.frame(tsls = results[[i]]$tsls, 
                          oracle = results[[i]]$oracle, 
                          corrected = results[[i]]$corrected,
                          probit = results[[i]]$probit,
                          logit = results[[i]]$logit,
                          lpm = results[[i]]$lpm)
  
  # plot_data %>% ggplot() + 
  #   geom_density(aes(x = tsls, col = "TSLS")) + 
  #   geom_density(aes(x = corrected, col = "corrected")) + 
  #   geom_density(aes(x = oracle, col = "oracle"))
  
  plot_data %>% ggplot() + 
    geom_density(aes(x = corrected, col = "corrected")) + 
    geom_density(aes(x = oracle, col = "oracle")) +
    geom_density(aes(x = probit, col = "probit")) + 
    geom_density(aes(x = logit, col = "logit")) +
    geom_density(aes(x = lpm, col = "lpm")) + 
    theme_bw() +
    theme(legend.title = element_blank()) +
    labs(x = "Estimate", y = "")
  
  plot_data %>% ggplot() + 
    geom_density(aes(x = corrected, col = "corrected")) + 
    geom_density(aes(x = oracle, col = "oracle")) +
    geom_density(aes(x = probit, col = "probit")) + 
    theme_bw() +
    theme(legend.title = element_blank()) +
    labs(x = "Estimate", y = "")
  
}


# Summary statistics
table <- matrix(0, 6 * sims, 7)

get_stats <- function(results){
  
  stats <- rep(0, 7)
  results <- results[!is.infinite(results)]
  stats[1] <- mean(results, na.rm = TRUE)
  stats[2] <- sd(results, na.rm = TRUE)
  stats[3:7] <- quantile(results, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  return(stats)
  
}

rownames(table) <- rep(N, 6)
table[1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "tsls"))},
           results = results))
table[sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "corrected"))},
           results = results))
table[2 * sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "oracle"))},
           results = results))
table[3 * sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "probit"))},
           results = results))
table[4 * sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "logit"))},
           results = results))
table[5 * sims + 1:sims, ] <- 
  t(sapply(1:sims, 
           function(i, results){get_stats(getElement(results[[i]], "lpm"))},
           results = results))


xtable(table[1:sims, ], digits = 3) # tsls
xtable(table[sims + 1:sims, ], digits = 3) # corrected
xtable(table[2 * sims + 1:sims, ], digits = 3) # oracle
xtable(table[3 * sims + 1:sims, ], digits = 3) # probit
xtable(table[4 * sims + 1:sims, ], digits = 3) # logit
xtable(table[5 * sims + 1:sims, ], digits = 3) # lpm

# 
# ########################### BOOTSTRAP STANDARD ERRORS ##########################
# 
# # Obtain standard errors by non-parametric bootstrap
# temp <- boot(data = data, statistic = est_theta, varnames = varnames, R = 10)
