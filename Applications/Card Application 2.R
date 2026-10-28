############################ EMPIRICAL APPLICATIONS ############################

# Test applications to the data from Card (1993) 

rm(list = ls())
gc()

library(ggplot2)
library(dplyr)
library(tidyr)
library(wooldridge)
library(haven)
library(AER)
library(ivreg)
library(xtable)

################################### FUNCTIONS ################################## 

# Returns the p-value from a Wald test of joint significance of the Z * X 
# interaction terms in the propensity score estimation
# - est should be the output from an lm/glm estimation
# - interactions_names is a vector containing the names of all interaction terms
Wald_test <- function(est, interactions_names){
  
  theta_hat <- est$coefficients
  R <- names(theta_hat) %in% interactions_names
  r <- sum(R)
  Sigma <- vcov(est)
  
  return(1 - pchisq(matrix(theta_hat[R], 1, r) %*% 
                      solve(Sigma[R, R]) %*% 
                      matrix(theta_hat[R], r, 1),
                    r))
  
}

# Function to calculate the fitted values of the sign function
# - theta_hat is estimated parameters (have names)
# - data is the data (including interactions)
get_s_hat <- function(theta_hat, data){
  
  relevant <- names(theta_hat)[grep("Z", names(theta_hat))]
  theta_hat <- theta_hat[relevant]
  
  relevant <- relevant[2:length(relevant)]
  relevant <- gsub("Z_", "", relevant)
  k <- length(relevant)
  
  return(sign(theta_hat["Z"] + as.matrix(data[relevant]) %*% 
                matrix(theta_hat[1 + 1:k], k, 1)))
  
}

S_theta <- function(theta, data, varnames){
  
  theta <- theta / sqrt(sum(theta ^ 2))
  index <- as.matrix(data[varnames]) %*% theta
  D_sort <- data$D[order(index)] 
  tally <- cumsum(1 - D_sort)
  val <- sum(tally[as.logical(D_sort)])
  
  return(val)
  
}

calc_sd <- function(data, Z_tilde, tau_feasible, s_hat, X_names){
  
  # returns the estimate of the asymptotic SD
  
  # data: a data frame containing the data (including generated variables)
  # Z_tilde: Z - E[Z|X]
  # tau_feasible: the estimate (from the feasible estimator)
  # s_hat: the estimated sign function
  # X_names: list of covariate names
  
  n <- nrow(data)
  K <- length(X_names)
  
  # "errors"
  e_Y <- data$Y - data$D * tau_feasible
  e_Z <- Z_tilde
  
  X <- data[names(data) %in% X_names]
  
  # sigma_X = E_n[XX^T]
  sigma_X <- (t(as.matrix(X)) %*% as.matrix(X)) / n
  
  # C = E_n[e_Y * s_hat(X) * X]
  C <- apply(X, 2, function(X_col, e_Y, s_hat){
    return(X_col * e_Y * s_hat)
  }, e_Y = e_Y, s_hat = s_hat)
  C <- colMeans(C) 
  
  # denominator
  denom <- sum(data$D * s_hat * e_Z) / n
  
  # numerator term 1: E_n[e_Y ^ 2 * e_Z ^ 2 * s_hat]
  num_1 <- sum((e_Y ^ 2) * (e_Z ^ 2)) / n 
  
  # numerator term 2: C^T %*% sigma_X^{-1} %*% E_n[e_Y * s_hat * (e_Z ^ 2) * X]
  # 
  num_2 <- apply(X, 2, function(X_col, e_Y, s_hat, e_Z){
    return(e_Y * s_hat * (e_Z ^ 2) * X_col)
  }, e_Y = e_Y, s_hat = s_hat, e_Z = e_Z)
  num_2 <- matrix(colMeans(num_2), nrow = K, ncol = 1)
  num_2 <- matrix(C, nrow = 1, ncol = K) %*% solve(sigma_X) %*% num_2
  
  # numerator term 3: C^T %*% sigma_X^{-1} %*% E_n[e_Z^2 XX^T] %*% sigma_X^{-1}
  # %*% C
  X_temp <- apply(X, 2, function(X_col, e_Z){
    return(X_col * e_Z)
  }, e_Z = e_Z)
  X_temp <- (t(X_temp) %*% X_temp) / n
  num_3 <- matrix(C, nrow = 1, ncol = K) %*% solve(sigma_X) %*%
    X_temp %*% solve(sigma_X) %*% matrix(C, nrow = K, ncol = 1)
  
  # numerator = num_1 - 2 * num_2 + num_3
  num <- num_1 - 2 * num_2 + num_3
  
  return(sqrt(num) / denom)
  
}

################################# LOAD IN DATA #################################

set.seed(1993)
data('card')

# Tidy up
colnames(card)[3] <- "Z"
colnames(card)[33] <- "Y"
card$D <- as.numeric(card$educ >= 16)

# Construct all variables used in the original analysis
# dummies for missing fatheduc and motheduc
card$motheduc_na <- as.numeric(is.na(card$motheduc))
card$fatheduc_na <- as.numeric(is.na(card$fatheduc))
# Replace missing values with the mean
card$motheduc[is.na(card$motheduc)] <- mean(card$motheduc, na.rm = T)
card$fatheduc[is.na(card$fatheduc)] <- mean(card$fatheduc, na.rm = T)

# Create parental education class dummies by brute force
card$edclass1 <- as.numeric(card$fatheduc < 12 & card$motheduc < 12)
card$edclass2 <- as.numeric(card$fatheduc < 12 & card$motheduc == 12)
card$edclass3 <- as.numeric(card$fatheduc < 12 & card$motheduc > 12)
card$edclass4 <- as.numeric(card$fatheduc == 12 & card$motheduc < 12)
card$edclass5 <- as.numeric(card$fatheduc == 12 & card$motheduc == 12)
card$edclass6 <- as.numeric(card$fatheduc == 12 & card$motheduc > 12)
card$edclass7 <- as.numeric(card$fatheduc > 12 & card$motheduc < 12)
card$edclass8 <- as.numeric(card$fatheduc > 12 & card$motheduc == 12)

# Scale expersq as in the original paper
card$expersq <- card$expersq / 100

varlist <- c("exper", "expersq", "black", "smsa", "south", "smsa66", "reg661", 
             "reg662", "reg663","reg664", "reg665", "reg666", "reg667", "reg668", 
             "motheduc", "fatheduc", "motheduc_na","fatheduc_na", "edclass1", 
             "edclass2", "edclass3", "edclass4", "edclass5", "edclass6", 
             "edclass7", "edclass8", "momdad14", "sinmom14")

# COVARIATE GROUPS FOR ESTIMATION
# 1) "exper", "expersq", "black", "south", "smsa"
# 2) "smsa66", "reg661", "reg662", "reg663","reg664", "reg665", "reg666", 
# "reg667", "reg668"
# 3) "motheduc", "fatheduc", "motheduc_na","fatheduc_na"
# 4) "edclass1", "edclass2", "edclass3", "edclass4", "edclass5", "edclass6", 
# "edclass7", "edclass8"
# 5) "momdad14", "sinmom14"

######################### REPLICATE CARD (1993) RESULTS ########################

# ESTIMATES WITH DIFFERENT COVARIATE GROUPS
varnames <- list("group1" = c("exper", "expersq", "black", "south", "smsa"),
                 "group2" = c("smsa66","reg661", "reg662", "reg663","reg664", 
                              "reg665", "reg666", "reg667","reg668"),
                 "group3" = c("motheduc", "fatheduc", "motheduc_na","fatheduc_na"),
                 "group4" = c("edclass1", "edclass2", "edclass3", "edclass4", 
                              "edclass5", "edclass6", "edclass7", "edclass8"),
                 "group5" = c("momdad14", "sinmom14"))

# replicate the results from Card (1993)
ols_replication <- vector(mode = "list", length = 5)
ols_table <- matrix(0, 13, 5)

for (i in 1:5){
  
  ols_formula <- as.formula(
    paste0("Y ~ educ + ", 
           paste(unlist(varnames[1:i]), collapse = " + ")))
  ols_replication[[i]] <- lm(ols_formula, data = card)
  results <- summary(ols_replication[[i]])$coefficients
  # coefficients
  ols_table[seq(1, 11, by = 2), i] <- results[, 1][c("educ", varnames[[1]])]
  # standard errors
  ols_table[seq(2, 12, by = 2), i] <- results[, 2][c("educ", varnames[[1]])]
  # R^2
  ols_table[13, i] <- summary(ols_replication[[i]])$r.squared
  
}

rownames(ols_table) <- rep("", nrow(ols_table)) 
rownames(ols_table)[seq(1, 11, by = 2)] <- c("educ", varnames[[1]])
rownames(ols_table)[13] <- "R^2"
xtable(ols_table, digits = 3)

# TSLS estimates
tsls_replication <- vector(mode = "list", length = 5)
tsls_table <- matrix(0, 2, 5)

for (i in 1:5){
  
  tsls_formula <-  as.formula(
    paste0("Y ~ educ + ", 
           paste(c(unlist(varnames[1:i])), collapse = " + "),
           "| Z + ",
           paste(c(unlist(varnames[1:i])), collapse = " + ")))
  tsls_replication[[i]] <- ivreg(formula = tsls_formula, data = card)
  results <- summary(tsls_replication[[i]])$coefficients
  # coefficients
  tsls_table[1, i] <- results[, 1]["educ"]
  # standard errors
  tsls_table[2, i] <- results[, 2]["educ"]
  
}

xtable(tsls_table[, c(2,5)], digits = 3)

######################### CONVENTIONAL 2SLS ESTIMATION #########################

tsls <- vector(mode = "list", length = 5)

for (i in 1:5){
  
  tsls_formula <-  as.formula(
    paste0("Y ~ D + ", 
           paste(c(unlist(varnames[1:i])), collapse = " + "),
           "| Z + ",
           paste(c(unlist(varnames[1:i])), collapse = " + ")))
  tsls[[i]] <- ivreg(formula = tsls_formula, 
                     data = card)
  
}

# make a table
tsls_table <- matrix(0, 2, 5)

for (i in 1:5){
  
  temp <- summary(tsls[[i]])$coefficients
  temp <- temp[grep("D", rownames(temp)), ]
  tsls_table[1, i] <- temp["Estimate"]
  tsls_table[2, i] <- temp["Std. Error"]
  
}

rownames(tsls_table) <- c("Estimate", "Std. Error")
xtable(tsls_table, digits = 3)

################################ CORRECTED TSLS ################################ 

# We consider the following variations
# 1) correction with probit
# 2) correction with logit estimates
# 3) correction with lpm estimates
# 4) correction with MRC estimates

# Each version is estimated with each combination of control variables
# From the `first stage' (estimation of the sign function), we report
# 1) Joint test of the significance of the interaction terms
# 2) Proportions of observations where s(X) is positive/negative
# From the `second stage' (the corrected TSLS estimator), we report
# 1) The estimated LATE
# 2) Standard errors

# storage for results
probit <- vector(mode = "list", length = 5)
logit <- vector(mode = "list", length = 5)
lpm <- vector(mode = "list", length = 5)
mrc <- vector(mode = "list", length = 5)

for (i in 1:5){
  
  #### PREPARE DATA ####
  
  data <- card[c("Y", "D", "Z", unlist(varnames[1:i]))]
  control_vars <- unlist(varnames[1:i])
  
  # generate the orthogonalised instrument
  Z_hat <- 
    lm(as.formula(
      paste0("Z ~ ", 
             paste(control_vars, collapse = " + "))),
      data = data)$fitted.values
  Z_tilde <- data$Z - Z_hat
  
  # generate interactions
  interactions <- matrix(
    unlist(
      sapply(control_vars, 
             function(X_name, data){return(data["Z"] * data[X_name])},
             data = data)),
    nrow = nrow(data), 
    ncol = length(control_vars))
  colnames(interactions) <- paste0("Z_", control_vars)
  interactions <- data.frame(interactions)
  
  data <- cbind(data, interactions)
  
  #### PROBIT ####
  
  # coefficient estimates
  probit_est <- glm(as.formula(
    paste0("D ~ Z + ", 
           paste(control_vars, collapse = " + "),
           " + ",
           paste(names(interactions), collapse = " + "))),
    family = binomial(link = "probit"),
    data = data)
  
  # Wald test of the significance of the interaction terms
  probit_wald <- Wald_test(probit_est, colnames(interactions))
  
  # calculate S_hat 
  probit_s_hat <- get_s_hat(probit_est$coefficients, data)
  
  # sum of observations where we have defiers 
  probit_share <- sum(probit_s_hat == -1) / nrow(data)
  
  # get the estimator
  probit_tau <- sum(data$Y * probit_s_hat * Z_tilde) / 
    sum(data$D * probit_s_hat * Z_tilde)
  
  # standard error
  probit_se <- calc_sd(data, Z_tilde, probit_tau, probit_s_hat, control_vars) / 
    sqrt(nrow(data))
  
  probit[[i]] <- list(est = probit_est,
                      wald = probit_wald,
                      share = probit_share,
                      tau = probit_tau,
                      se = probit_se)
  
  #### LOGIT ####
  
  # coefficient estimates
  logit_est <- glm(as.formula(
    paste0("D ~ Z + ", 
           paste(control_vars, collapse = " + "),
           " + ",
           paste(names(interactions), collapse = " + "))),
    family = binomial(link = "logit"),
    data = data)
  
  # Wald test of the significance of the interaction terms
  logit_wald <- Wald_test(logit_est, colnames(interactions))
  
  # calculate S_hat 
  logit_s_hat <- get_s_hat(logit_est$coefficients, data)
  
  # sum of observations where we have defiers 
  logit_share <- sum(logit_s_hat == -1) / nrow(data)
  
  # get the estimator
  logit_tau <- sum(data$Y * logit_s_hat * Z_tilde) / 
    sum(data$D * logit_s_hat * Z_tilde)
  
  # standard error
  logit_se <- calc_sd(data, Z_tilde, logit_tau, logit_s_hat, control_vars) / 
    sqrt(nrow(data))
  
  logit[[i]] <- list(est = logit_est,
                     wald = logit_wald,
                     share = logit_share,
                     tau = logit_tau,
                     se = logit_se)
  
  #### LPM ####
  
  # coefficient estimates
  lpm_est <- lm(as.formula(
    paste0("D ~ Z + ", 
           paste(control_vars, collapse = " + "),
           " + ",
           paste(names(interactions), collapse = " + "))),
    data = data)
  
  # Wald test of the significance of the interaction terms
  lpm_wald <- Wald_test(lpm_est, colnames(interactions))
  
  # calculate S_hat 
  lpm_s_hat <- get_s_hat(lpm_est$coefficients, data)
  
  # sum of observations where we have defiers 
  lpm_share <- sum(lpm_s_hat == -1) / nrow(data)
  
  # get the estimator
  lpm_tau <- sum(data$Y * lpm_s_hat * Z_tilde) / 
    sum(data$D * lpm_s_hat * Z_tilde)
  
  # standard error
  lpm_se <- calc_sd(data, Z_tilde, lpm_tau, lpm_s_hat, control_vars) / 
    sqrt(nrow(data))
  
  lpm[[i]] <- list(est = lpm_est,
                   wald = lpm_wald,
                   share = lpm_share,
                   tau = lpm_tau,
                   se = lpm_se)
  
  #### MRC ####
  
  theta_probit <- 
    probit_est$coefficients[!names(probit_est$coefficients) %in% "(Intercept)"]
  
  theta_logit <-
    logit_est$coefficients[!names(logit_est$coefficients) %in% "(Intercept)"]
  
  theta_lpm <- 
    lpm_est$coefficients[!names(lpm_est$coefficients) %in% "(Intercept)"]
  
  
  # Run the Nelder-Mead-Head simplex from each pilot, take the maximum
  n_perturb <- 5 # n_perturb random perturbations of each pilot point
  n_random <- 10 # n_random randomly generated starting points
  
  temp_varlist <- names(theta_probit)
  
  candidates <- vector(mode = "list", length = 3 + 3 * n_perturb + n_random)
  
  candidates[[1]] <- optim(par = theta_probit, 
                           fn = S_theta, 
                           method = "Nelder-Mead", 
                           data = data, 
                           varnames = temp_varlist, 
                           control = list("fnscale" = -1))
  
  candidates[[2]] <- optim(par = theta_logit, 
                           fn = S_theta, 
                           method = "Nelder-Mead", 
                           data = data, 
                           varnames = temp_varlist, 
                           control = list("fnscale" = -1))
  
  candidates[[3]] <- optim(par = theta_lpm, 
                           fn = S_theta, 
                           method = "Nelder-Mead", 
                           data = data, 
                           varnames = temp_varlist, 
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
                                 varnames = temp_varlist, 
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
                                             varnames = temp_varlist, 
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
            varnames = temp_varlist, 
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
            varnames = temp_varlist, 
            control = list("fnscale" = -1))
    
  }
  
  # Choose the best performing value
  j <- which.max(sapply(candidates, getElement, "value"))
  mrc_est <- candidates[[j]]$par
  mrc_est <- mrc_est / sqrt(sum(mrc_est ^ 2)) # scale to be length 1
  
  # calculate S_hat 
  mrc_s_hat <- get_s_hat(mrc_est, data)
  
  # sum of observations where we have defiers 
  mrc_share <- sum(mrc_s_hat == -1) / nrow(data)
  
  # get the estimator
  mrc_tau <- sum(data$Y * mrc_s_hat * Z_tilde) / 
    sum(data$D * mrc_s_hat * Z_tilde)
  
  # standard error
  mrc_se <- calc_sd(data, Z_tilde, mrc_tau, mrc_s_hat, control_vars) / 
    sqrt(nrow(data))
  
  mrc[[i]] <- list(est = mrc_est,
                   share = mrc_share,
                   tau = mrc_tau,
                   se = mrc_se)
  
}

temp_table <- matrix(0, 12, 5)
stat_names <- c("tau", "se", "share")

closeness <- matrix(0, 3, 5)

for (i in 1:5){
  
  temp_table[1:3, i] <- sapply(stat_names, getElement, object = probit[[i]])
  temp_table[4:6, i] <- sapply(stat_names, getElement, object = logit[[i]])
  temp_table[7:9, i] <- sapply(stat_names, getElement, object = lpm[[i]])
  temp_table[10:12, i] <- sapply(stat_names, getElement, object = mrc[[i]])
  
  mrc_est <- mrc[[i]]$est
  
  probit_est <- probit[[i]]$est$coefficients
  probit_est <- probit_est[names(probit_est) %in% names(mrc_est)]
  closeness[1, i] <- cor(probit_est, mrc_est)
  
  logit_est <- logit[[i]]$est$coefficients
  logit_est <- logit_est[names(logit_est) %in% names(mrc_est)]
  closeness[2, i] <- cor(logit_est, mrc_est)
  
  lpm_est <- lpm[[i]]$est$coefficients
  lpm_est <- lpm_est[names(lpm_est) %in% names(mrc_est)]
  closeness[3, i] <- cor(lpm_est, mrc_est)
  
}

rownames(temp_table) <- rep(stat_names, 4)
closeness

xtable(temp_table, digits = 3)

path = getwd()

card_results <- list(tsls = tsls,
                     probit = probit,
                     logit = logit,
                     lpm = lpm,
                     mrc = mrc)

saveRDS(card_results, file = paste0(path, "/card_results.rds"))
