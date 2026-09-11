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


probit_version <- function(data, indices, varnames){
  
  # select data
  boot_data <- data[indices, ]
  
  # probit
  
  
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

# number of boostrap draws
B <- 10 

for (i in 1:5){
  
  #### PREPARE DATA ####
  
  data <- card[c("Y", "D", "Z", unlist(varnames[1:i]))]
  control_vars <- unlist(varnames[1:i])
  
  # # principal components
  # PC <- prcomp(data[control_vars], retx = T)$x
  
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
  
  probit[[i]] <- list(est = probit_est,
                      wald = probit_wald,
                      share = probit_share,
                      tau = probit_tau)
  
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
  
  logit[[i]] <- list(est = logit_est,
                      wald = logit_wald,
                      share = logit_share,
                      tau = logit_tau)
  
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
  
  lpm[[i]] <- list(est = lpm_est,
                     wald = lpm_wald,
                     share = lpm_share,
                     tau = lpm_tau)
  
  #### MRC ####
  
  theta_probit <- 
    probit_est$coefficients[!names(probit_est$coefficients) %in% "(Intercept)"]
  
  theta_logit <-
    logit_est$coefficients[!names(logit_est$coefficients) %in% "(Intercept)"]
  
  theta_lpm <- 
    lpm_est$coefficients[!names(lpm_est$coefficients) %in% "(Intercept)"]
  
  # Run the Nelder-Mead-Head simplex from each pilot, take the maximum
  candidates <- vector(mode = "list", length = 3)
  temp_varlist <- names(theta_probit)
  
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
  
  mrc[[i]] <- list(est = mrc_est,
                   share = mrc_share,
                   tau = mrc_tau)
  
  ################################# BOOTSTRAP ##################################
  
  boot_probit <- rep(0, B)
  boot_logit <- rep(0, B)
  boot_lpm <- rep(0, B)
  boot_mrc <- rep(0, B)
  boot_mrc_theta <- vector(mode = "list", length = B) 
  
  # obtain bootstrap standard errors
  for (b in 1:B){
    
    draw <- sample(1:nrow(data), replace = T)
    boot_data <- data[draw, ]
    
    # bootstrap orthogonalised instrument
    boot_Z_hat <- 
      lm(as.formula(
        paste0("Z ~ ", 
               paste(control_vars, collapse = " + "))),
        data = boot_data)$fitted.values
    boot_Z_tilde <- boot_data$Z - boot_Z_hat
    
    
    # Probit
    boot_probit_est <- glm(as.formula(
      paste0("D ~ Z + ", 
             paste(control_vars, collapse = " + "),
             " + ",
             paste(names(interactions), collapse = " + "))),
      family = binomial(link = "probit"),
      data = boot_data)
    
    boot_probit_s_hat <- get_s_hat(boot_probit_est$coefficients, boot_data)
    
    boot_probit[b] <- sum(boot_data$Y * boot_probit_s_hat * boot_Z_tilde) / 
      sum(boot_data$D * boot_probit_s_hat * boot_Z_tilde)
    
    # Logit
    boot_logit_est <- glm(as.formula(
      paste0("D ~ Z + ", 
             paste(control_vars, collapse = " + "),
             " + ",
             paste(names(interactions), collapse = " + "))),
      family = binomial(link = "logit"),
      data = boot_data)
    
    boot_logit_s_hat <- get_s_hat(boot_logit_est$coefficients, boot_data)
    
    boot_logit[b] <- sum(boot_data$Y * boot_logit_s_hat * boot_Z_tilde) / 
      sum(boot_data$D * boot_logit_s_hat * boot_Z_tilde)
    
    # LPM
    boot_lpm_est <- lm(as.formula(
      paste0("D ~ Z + ", 
             paste(control_vars, collapse = " + "),
             " + ",
             paste(names(interactions), collapse = " + "))),
      data = boot_data)
    
    boot_lpm_s_hat <- get_s_hat(boot_lpm_est$coefficients, boot_data)
    
    boot_lpm[b] <- sum(boot_data$Y * boot_lpm_s_hat * boot_Z_tilde) / 
      sum(boot_data$D * boot_lpm_s_hat * boot_Z_tilde)
    
    # MRC
    boot_theta_probit <- 
      boot_probit_est$coefficients[!names(boot_probit_est$coefficients) %in% "(Intercept)"]
    
    boot_theta_logit <-
      boot_logit_est$coefficients[!names(boot_logit_est$coefficients) %in% "(Intercept)"]
    
    boot_theta_lpm <- 
      boot_lpm_est$coefficients[!names(boot_lpm_est$coefficients) %in% "(Intercept)"]
    
    # Run the Nelder-Mead-Head simplex from each pilot, take the maximum
    candidates <- vector(mode = "list", length = 3)
    temp_varlist <- names(boot_theta_probit)
    
    candidates[[1]] <- optim(par = boot_theta_probit, 
                             fn = S_theta, 
                             method = "Nelder-Mead", 
                             data = boot_data, 
                             varnames = temp_varlist, 
                             control = list("fnscale" = -1))
    
    candidates[[2]] <- optim(par = boot_theta_logit, 
                             fn = S_theta, 
                             method = "Nelder-Mead", 
                             data = boot_data, 
                             varnames = temp_varlist, 
                             control = list("fnscale" = -1))
    
    candidates[[3]] <- optim(par = boot_theta_lpm, 
                             fn = S_theta, 
                             method = "Nelder-Mead", 
                             data = boot_data, 
                             varnames = temp_varlist, 
                             control = list("fnscale" = -1))
    
    # Choose the best performing value
    j <- which.max(sapply(candidates, getElement, "value"))
    boot_mrc_est <- candidates[[j]]$par
    boot_mrc_est <- boot_mrc_est / sqrt(sum(boot_mrc_est ^ 2)) # scale to be length 1
    boot_mrc_theta[[b]] <- boot_mrc_est
    
    # calculate S_hat 
    boot_mrc_s_hat <- get_s_hat(boot_mrc_est, boot_data)
    
    # get the estimator
    boot_mrc[b] <- sum(boot_data$Y * boot_mrc_s_hat * boot_Z_tilde) / 
      sum(boot_data$D * boot_mrc_s_hat * boot_Z_tilde)
    
  }
  
  probit[[i]]$boot <- boot_probit
  probit[[i]]$sd <- sd(boot_probit)
  
  logit[[i]]$boot <- boot_logit
  logit[[i]]$sd <- sd(boot_logit)
  
  lpm[[i]]$boot <- boot_lpm
  lpm[[i]]$sd <- sd(boot_lpm)
  
  mrc[[i]]$boot <- boot_mrc
  mrc[[i]]$sd <- sd(boot_mrc)
  mrc[[i]]$boot_theta <- boot_mrc_theta
  
  # Wald test for the MRC estimates
  n_vars <- length(boot_mrc_theta[[1]])
  store <- matrix(unlist(boot_mrc_theta), B, n_vars, byrow = T)
  colnames(store) <- names(boot_mrc_theta[[1]])
  store <- data.frame(store)
  
  R <- names(mrc_est) %in% colnames(interactions)
  r <- sum(R)
  Sigma <- cov(store[, R])
  
  mrc[[i]]$wald <- try(
    1 - pchisq(matrix(mrc_est[R], 1, r) %*%
                 solve(Sigma) %*%
                 matrix(mrc_est[R], r, 1),
               r), silent = TRUE)
  
}

temp_table <- matrix(0, 15, 5)
stat_names <- c("tau", "sd", "share", "wald")

closeness <- matrix(0, 3, 5)

for (i in 1:5){
  
  temp_table[1:4, i] <- sapply(stat_names, getElement, object = probit[[i]])
  temp_table[5:8, i] <- sapply(stat_names, getElement, object = logit[[i]])
  temp_table[9:12, i] <- sapply(stat_names, getElement, object = lpm[[i]])
  
  temp_table[13:15, i] <- sapply(c("tau", "sd", "share"), getElement,
                                 object = mrc[[i]])
  
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

rownames(temp_table) <- c(rep(stat_names, 3), c("tau", "sd", "share"))
temp_table[seq(1, 13, by = 4), ]
temp_table[seq(2, 14, by = 3), ]
closeness

temp_table

path = getwd()

saveRDS(card_results, file = paste0(path, "/card_results.rds"))
