################################# FRONT MATTER #################################

# Implements the estimator with the data from Angrist and Evans (1998)

library(ivreg)
library(xtable)
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list = ls())

path <- getwd()
path <- paste0(path, "/Applications/")

################################### FUNCTIONS ##################################

# Function to calculate the fitted values of the sign function
# - theta_hat is estimated parameters (have names)
# - data is the data (including interactions)
get_s_hat <- function(theta_hat, data){
  
  relevant <- names(theta_hat)[grep("Z", names(theta_hat))]
  theta_hat <- theta_hat[relevant]
  theta_hat <- theta_hat[!is.na(theta_hat)]
  
  relevant <- names(theta_hat)[2:length(theta_hat)]
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

calc_sd <- function(data, Z_tilde, tau_feasible, s_hat, X_names, Yname){
  
  # returns the estimate of the asymptotic SD
  
  # data: a data frame containing the data (including generated variables)
  # Z_tilde: Z - E[Z|X]
  # tau_feasible: the estimate (from the feasible estimator)
  # s_hat: the estimated sign function
  # X_names: list of covariate names

  # "errors"
  s_hat <- s_hat[!is.na(data[[Yname]])]
  
  e_Z <- Z_tilde
  e_Z <- e_Z[!is.na(data[[Yname]])]
  
  data <- data[!is.na(data[[Yname]]), ]
  e_Y <- data[[Yname]] - data$D * tau_feasible

  n <- nrow(data)
  K <- length(X_names)
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

################################## PREPARATION #################################

# Load in cleaned and prepared data
data <- readRDS(file = paste0(path, "1980_data.rds"))

# Drop observations with missing years of education
data <- dplyr::filter(data, is.na(yr_ed) == F)

data <- dplyr::rename(data, Z = same_sex, D = morekids)

Ylist <- c("wkswork1", "worked", "uhrswork", "inctot")
Xlist <- c("black", "oth", "hispan", "bpl", "sp_bpl", "age", "age_fst", "yr_ed")

data$Z <- as.numeric(data$Z)
data$D <- as.numeric(data$D)

############################## SUMMARY STATISTICS ##############################

# immigrant mother, US father
sum((data$bpl != 0) * (data$sp_bpl == 0))
# US mother, immigrant father
sum((data$bpl == 0) * (data$sp_bpl != 0))
# US couple
sum((data$bpl == 0) * (data$sp_bpl == 0))
# immigrant couple
sum((data$bpl != 0) * (data$sp_bpl != 0))

# counts for each birthplace
vals <- length(unique(data$bpl))
data <- data[order(data$bpl), ]

bpl_table <- matrix(0, nrow = vals, ncol = 2)
bpl_table[, 1] <- sort(unique(data$bpl))
bpl_table[, 2] <- tabulate(as.factor(data$bpl))

vals <- length(unique(data$sp_bpl))
data <- data[order(data$sp_bpl), ]

sp_bpl_table <- matrix(0, nrow = vals, ncol = 2)
sp_bpl_table[, 1] <- sort(unique(data$sp_bpl))
sp_bpl_table[, 2] <- tabulate(as.factor(data$sp_bpl))

# drop any bpl or spl_bpl that occurs less than 10 times (arbitrary)
bpl_drop <- bpl_table[, 1][bpl_table[, 2] <= 10]
sp_bpl_drop <- sp_bpl_table[, 1][sp_bpl_table[, 2] <= 10]

data <- data[!(data$bpl %in% bpl_drop), ]
data <- data[!(data$sp_bpl %in% sp_bpl_drop), ]


##################################### TSLS #####################################

store <- data
data$bpl <- as.factor(data$bpl)
data$sp_bpl <- as.factor(data$sp_bpl)

tsls <- vector(mode = "list", length = 4)

for (i in 1:4){
  
  tsls_formula <-  as.formula(
    paste0(Ylist[i], " ~ D + ", 
           paste(Xlist, collapse = " + "),
           "| Z + ",
           paste(Xlist, collapse = " + ")))
  tsls[[i]] <- ivreg(formula = tsls_formula, 
                     data = data)
  
}

# make a table
tsls_table <- matrix(0, 3, 4)

for (i in 1:4){
  
  temp <- summary(tsls[[i]])$coefficients
  temp <- temp[grep("D", rownames(temp)), ]
  tsls_table[1, i] <- temp["Estimate"]
  tsls_table[2, i] <- temp["Std. Error"]
  tsls_table[3, i] <- tsls[[i]]$df.residual
  
}

rownames(tsls_table) <- c("Estimate", "Std. Error", "n")
colnames(tsls_table) <- Ylist
xtable(tsls_table, digits = 3)

########################### REORDERED IV PREPARATIONS ##########################

# generate the orthogonalised instrument
Z_hat <- 
  lm(as.formula(
    paste0("Z ~ ", 
           paste(Xlist, collapse = " + "))),
    data = data)$fitted.values
Z_tilde <- data$Z - Z_hat
rm(Z_hat)

# go back to the non-factor bpl and sp_bpl
rm(data)
data <- store

# generate fixed effects
n <- nrow(data)
bpl_vals <- sort(unique(data$bpl))
n_vals <- length(bpl_vals)
bpl_fe <- matrix(0, n, n_vals)
colnames(bpl_fe) <- paste0("bpl_", bpl_vals)

for (i in 1:n_vals){
  
  bpl_fe[, i] <- as.numeric(data$bpl == bpl_vals[i])
  
}

bpl_fe <- data.frame(bpl_fe)

sp_bpl_vals <- sort(unique(data$sp_bpl))
n_vals <- length(sp_bpl_vals)
sp_bpl_fe <- matrix(0, n, n_vals)
colnames(sp_bpl_fe) <- paste0("sp_bpl_", sp_bpl_vals)

for (i in 1:n_vals){
  
  sp_bpl_fe[, i] <- as.numeric(data$sp_bpl == sp_bpl_vals[i])
  
}

sp_bpl_fe <- data.frame(sp_bpl_fe)

data <- cbind(data, bpl_fe, sp_bpl_fe)
Xlist <- Xlist[grep("bpl", Xlist, invert = T)]
Xlist <- c(Xlist, names(bpl_fe), names(sp_bpl_fe))

interactions <- matrix(
  unlist(
    sapply(Xlist, 
           function(X_name, data){return(data["Z"] * data[X_name])},
           data = data)),
  nrow = nrow(data), 
  ncol = length(Xlist))
colnames(interactions) <- paste0("Z_", Xlist)
interactions <- data.frame(interactions)
data <- cbind(data, interactions)

############################ REORDERED IV ESTIMATION ###########################

# We consider the following variations
# 1) correction with probit
# 2) correction with logit estimates
# 3) correction with lpm estimates
# 4) correction with MRC estimates

# Each version is estimated with each combination of control variables
# From the `first stage' (estimation of the sign function), we report
# 1) Proportions of observations where s(X) is positive
# From the `second stage' (the corrected TSLS estimator), we report
# 1) The estimated LATE
# 2) Standard errors

# storage for results
probit <- vector(mode = "list", length = 4)
logit <- vector(mode = "list", length = 4)
lpm <- vector(mode = "list", length = 4)
mrc <- vector(mode = "list", length = 4)

# As the conditioning covariates stay fixed, we need only run the first stage
# once!!!

#### PROBIT ####

# coefficient estimates
probit_est <- glm(as.formula(
  paste0("D ~ Z + ", 
         paste(Xlist, collapse = " + "),
         " + ",
         paste(names(interactions), collapse = " + "))),
  family = binomial(link = "probit"),
  data = data)

# calculate S_hat 
probit_s_hat <- get_s_hat(probit_est$coefficients, data)

#### LOGIT ####

# coefficient estimates
logit_est <- glm(as.formula(
  paste0("D ~ Z + ", 
         paste(Xlist, collapse = " + "),
         " + ",
         paste(names(interactions), collapse = " + "))),
  family = binomial(link = "logit"),
  data = data)

# calculate S_hat 
logit_s_hat <- get_s_hat(logit_est$coefficients, data)

#### LPM ####

# coefficient estimates
lpm_est <- lm(as.formula(
  paste0("D ~ Z + ", 
         paste(Xlist, collapse = " + "),
         " + ",
         paste(names(interactions), collapse = " + "))),
  data = data)

# calculate S_hat 
lpm_s_hat <- get_s_hat(lpm_est$coefficients, data)

#### MRC ####

theta_probit <- 
  probit_est$coefficients[!names(probit_est$coefficients) %in% "(Intercept)"]
theta_probit <- theta_probit[!is.na(theta_probit)]

theta_logit <-
  logit_est$coefficients[!names(logit_est$coefficients) %in% "(Intercept)"]
theta_logit <- theta_logit[!is.na(theta_logit)]

theta_lpm <- 
  lpm_est$coefficients[!names(lpm_est$coefficients) %in% "(Intercept)"]
theta_lpm <- theta_lpm[!is.na(theta_lpm)]

# Run the Nelder-Mead-Head simplex from each pilot, take the maximum
n_perturb <- 5 # n_perturb random perturbations of each pilot point
n_random <- 10 # n_random randomly generated starting points

temp_varlist <- names(theta_probit)
temp_varlist <- gsub("TRUE", "", temp_varlist)

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

# drop interactions from the data
data <- dplyr::select(data, - starts_with("Z_"))

#### SECOND STAGE ####

for (i in 1:4){
  
  #### PROBIT ####
  
  # sum of observations where we have defiers 
  probit_share <- sum((probit_s_hat == -1)[!is.na(data[Ylist[i]])]) / 
    sum(!is.na(data[Ylist[i]]))
  
  # get the estimator
  probit_tau <- sum(data[Ylist[i]] * probit_s_hat * Z_tilde, na.rm = T) / 
    sum((data$D * logit_s_hat * Z_tilde)[!is.na(data[[Ylist[i]]])])
  
  # standard error
  vars <- gsub("TRUE", "", names(probit_est$coefficients))
  vars <- Xlist[Xlist %in% vars[!is.na(probit_est$coefficients)]]
  probit_se <- calc_sd(data, Z_tilde, probit_tau, probit_s_hat, vars, Ylist[i]) / 
    sqrt(sum(!is.na(data[Ylist[i]])))
  
  probit[[i]] <- list(est = probit_est,
                      share = probit_share,
                      tau = probit_tau,
                      se = probit_se)
  
  #### LOGIT ####
  
  # sum of observations where we have defiers 
  logit_share <- sum((logit_s_hat == -1)[!is.na(data[Ylist[i]])]) / 
    sum(!is.na(data[Ylist[i]]))
  
  # get the estimator
  logit_tau <- sum(data[Ylist[i]] * logit_s_hat * Z_tilde, na.rm = T) / 
    sum((data$D * logit_s_hat * Z_tilde)[!is.na(data[[Ylist[i]]])])
  
  # standard error
  vars <- gsub("TRUE", "", names(logit_est$coefficients))
  vars <- Xlist[Xlist %in% vars[!is.na(logit_est$coefficients)]]
  logit_se <- calc_sd(data, Z_tilde, logit_tau, logit_s_hat, vars, Ylist[i]) / 
    sqrt(sum(!is.na(data[Ylist[i]])))
  
  logit[[i]] <- list(est = logit_est,
                     share = logit_share,
                     tau = logit_tau,
                     se = logit_se)
  
  #### LPM ####
  
  # sum of observations where we have defiers 
  lpm_share <- sum((lpm_s_hat == -1)[!is.na(data[Ylist[i]])]) / 
    sum(!is.na(data[Ylist[i]]))
  
  # get the estimator
  lpm_tau <- sum(data[Ylist[i]] * lpm_s_hat * Z_tilde, na.rm = T) / 
    sum((data$D * lpm_s_hat * Z_tilde)[!is.na(data[[Ylist[i]]])])
  
  # standard error
  vars <- gsub("TRUE", "", names(lpm_est$coefficients))
  vars <- Xlist[Xlist %in% vars[!is.na(lpm_est$coefficients)]]
  lpm_se <- calc_sd(data, Z_tilde, lpm_tau, lpm_s_hat, vars, Ylist[i]) / 
    sqrt(sum(!is.na(data[Ylist[i]])))
  
  lpm[[i]] <- list(est = lpm_est,
                   share = lpm_share,
                   tau = lpm_tau,
                   se = lpm_se)
  
  #### MRC ####
  
  # sum of observations where we have defiers 
  mrc_share <- sum((mrc_s_hat == -1)[!is.na(data[Ylist[i]])]) / 
    sum(!is.na(data[Ylist[i]]))
  
  # get the estimator
  mrc_tau <- sum(data[Ylist[i]] * mrc_s_hat * Z_tilde, na.rm = T) / 
    sum((data$D * mrc_s_hat * Z_tilde)[!is.na(data[[Ylist[i]]])])
  
  # standard error
  vars <- gsub("TRUE", "", names(mrc_est))
  vars <- Xlist[Xlist %in% vars[!is.na(mrc_est)]]
  mrc_se <- calc_sd(data, Z_tilde, mrc_tau, mrc_s_hat, vars, Ylist[i]) / 
    sqrt(sum(!is.na(data[Ylist[i]])))
  
  mrc[[i]] <- list(est = mrc_est,
                   share = mrc_share,
                   tau = mrc_tau,
                   se = mrc_se)
  
}

lpm_est <- lpm_est$coefficients
probit_est <- probit_est$coefficients
logit_est <- logit_est$coefficients

for (i in 1:4){
  
  probit[[i]]$est <- probit_est
  logit[[i]]$est <- logit_est
  lpm[[i]]$est <- lpm_est
  
}

for (i in 1:4){
  
  tsls[[i]] <- tsls[[i]]$coefficients
  
}

angrist_evans_results <- list(tsls = tsls,
                              probit = probit,
                              logit = logit,
                              lpm = lpm,
                              mrc = mrc)

saveRDS(angrist_evans_results, file = paste0(path, "angrist_evans_results.rds"))

reordered_IV_table <- matrix(0, 3, 4)

for (i in 1:4){
  
  reordered_IV_table[1, i] <- mrc[[i]]$tau
  reordered_IV_table[2, i] <- mrc[[i]]$se
  reordered_IV_table[3, i] <- mrc[[i]]$share
  
  
}

colnames(reordered_IV_table) <- Ylist
rownames(reordered_IV_table) <- c("Estimate", "Std. Error", "share")

xtable(tsls_table[, c(2, 1, 3, 4)], digits = 3)
xtable(reordered_IV_table[, c(2, 1, 3, 4)], digits = 3)
