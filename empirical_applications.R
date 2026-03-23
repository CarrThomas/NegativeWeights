########################### EMPIRICAL APPLICATIONS ###########################

# Test applications to the data from Card (1993) and Angrist and Evans (1998)

rm(list = ls())
gc()

library(ggplot2)
library(dplyr)
library(tidyr)
library(wooldridge)
library(haven)
library(AER)

################################ CARD (1993) ################################

set.seed(1993)

# Load in the data
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

varlist <- c("exper", "expersq", "black", "smsa", "south", "smsa66", "reg661", 
             "reg662", "reg663","reg664", "reg665", "reg666", "reg667", "reg668", 
             "motheduc", "fatheduc", "motheduc_na","fatheduc_na", "edclass1", 
             "edclass2", "edclass3", "edclass4", "edclass5", "edclass6", 
             "edclass7", "edclass8", "momdad14", "sinmom14")

# First stage with no interactions
no_interactions <- lm(as.formula(paste0("D ~ Z + ", 
                                        paste(varlist, collapse = " + "))),
                      data = card)

summary(no_interactions)

# First stage with interactions
interactions <- matrix(0, nrow(card), length(varlist))
for (i in 1:length(varlist)){
  
  interactions[, i] <- getElement(card, varlist[i]) * card$Z
  
}
colnames(interactions) <- paste0("Z_", varlist)
interactions <- data.frame(interactions)

with_interactions <- lm(
  as.formula(
    paste0("D ~ Z + ", 
           paste(varlist, collapse = " + "),
           " + ",
           paste(colnames(interactions), collapse = "+"))
    ),
  data = cbind(card, interactions))

summary(with_interactions)

# find data such that the overall effect of Z is negative
coefs <- with_interactions$coefficients
Z0_coefs <- coefs[grep("Z", names(coefs), invert = T)]
Z1_coefs <- coefs[grep("Z", names(coefs))]

# change when Z = 1
Z_effect <- Z1_coefs[1] + 
  as.matrix(dplyr::select(card, all_of(varlist))) %*% 
  matrix(Z1_coefs[2:length(Z1_coefs)],
         nrow = length(Z1_coefs) - 1,
         ncol = 1)

# summary states
# 1) count of positive and negative
sum(Z_effect >= 0)
sum(Z_effect < 0)

# 2) mean effect conditional on positive and negative
mean(Z_effect[Z_effect >= 0])
mean(Z_effect[Z_effect < 0])

# plot of Z_effect
data.frame(Z_effect = Z_effect) %>% ggplot() +
  geom_density(aes(x = Z_effect, col = "density")) + 
  geom_vline(aes(xintercept = 0, col = "line")) + 
  theme_bw() + 
  labs(y = "", x = expression(theta)) +
  theme(legend.position = "null")

# LATE for each group
Z_pos <- card[Z_effect >= 0, ]

# positive 
Z_tilde <- lm(
  as.formula(
    paste0("Z ~ ", 
           paste(varlist, collapse = " + "))), 
  data = Z_pos)$residuals

# Wald estimator numerator and denominator
num <- lm("Y ~ Z_tilde - 1", 
          data = data.frame(Y = Z_pos$Y, 
                            Z_tilde = Z_tilde))$coefficients
denom <- lm("D ~ Z_tilde - 1", 
            data = data.frame(D = Z_pos$D, 
                              Z_tilde = Z_tilde))$coefficients
LATE_pos <- num / denom
LATE_pos

# negative 
Z_neg <- card[Z_effect < 0, ]
Z_neg$Z <- 1 - Z_neg$Z

Z_tilde <- lm(
  as.formula(
    paste0("Z ~ ", 
           paste(varlist, collapse = " + "))), 
  data = Z_neg)$residuals

# Wald estimator numerator and denominator
num <- lm("Y ~ Z_tilde - 1", 
          data = data.frame(Y = Z_neg$Y, 
                            Z_tilde = Z_tilde))$coefficients
denom <- lm("D ~ Z_tilde - 1", 
            data = data.frame(D = Z_neg$D, 
                              Z_tilde = Z_tilde))$coefficients
LATE_neg <- num / denom
LATE_neg

# overall 
Z_tilde <- lm(
  as.formula(
    paste0("Z ~ ", 
           paste(varlist, collapse = " + "))), 
  data = card)$residuals

# Wald estimator numerator and denominator
num <- lm("Y ~ Z_tilde - 1", 
          data = data.frame(Y = card$Y, 
                            Z_tilde = Z_tilde))$coefficients
denom <- lm("D ~ Z_tilde - 1", 
            data = data.frame(D = card$D, 
                              Z_tilde = Z_tilde))$coefficients
LATE <- num / denom
LATE

# Implied weight
Z_tilde <- lm(
  as.formula(
    paste0("Z ~ ", paste(varlist, collapse = " + "))
  ), data = card)$residuals

sum(card$D * Z_tilde)
sum((card$D * Z_tilde)[Z_effect < 0])
sum((card$D * Z_tilde)[Z_effect >= 0])

sum((card$D * Z_tilde)[Z_effect < 0]) / sum(card$D * Z_tilde)


########################### ANGRIST AND EVANS (1998) ###########################

rm(list = ls())

path <- getwd()

####  Data wrangling ####

data <-  read_sas(paste(path, "/AngEv98/m_d_806.sas7bdat", sep = ""))
names(data) <- tolower(names(data))

# flag for first child born before marriage
data$qtrmar <- as.integer(data$qtrmar)
data$agemar <- as.integer(data$agemar)
data$qtrbthm <- as.integer(data$qtrbthm)
data$qtrbkid <- as.integer(data$qtrbkid)
data$illegit <- 0
data_0 <- subset(data, qtrmar == 0)
data <- subset(data, qtrmar != 0)
data$qtrmar <- data$qtrmar - 1
yom <- data$yobm + data$agemar + (data$qtrbthm > data$qtrmar)
dom_q <- yom + (data$qtrmar / 4)
do1b_q <- data$yobk + (data$qtrbkid / 4)
data$illegit <- dom_q > do1b_q
data <- rbind(data_0, data)
rm(data_0, yom, dom_q, do1b_q)

# sexes of the first two kids
data$Z <- data$sexk == data$sex2nd
data$D <- data$kidcount > 2

# race indicators for mom 
data$blackm <- data$racem == "02"
data$hispm <- data$racem == "12"
data$whitem <- data$racem == "01"
data$othracem <- 1 - data$blackm - data$hispm - data$whitem

# ages of mom when 1st kid was born
data$agem <- as.integer(data$agem)
data$ageqm <- 4 * (80 - data$yobm) - data$qtrbthm - 1
data$agefstm <- floor((data$ageqm - data$ageqk) / 4)

data$aged <- as.integer(data$aged)
data$yobd <- 79 - data$aged + (data$qtrbthd == "0")
data$qtrbthd <- as.integer(data$qtrbthd)
data$ageqd <- 4 * (80 - data$yobd) - data$qtrbthd
data$agefstd <- floor((data$ageqd - data$ageqk) / 4)

# outcome variable (log family income)
data$faminc1 <- as.numeric(data$faminc)
data$faminc1 <- data$faminc1 * 2.099173554
data$Y <- log(pmax(data$faminc1 , 1))

# indicator for whether the observation belongs in the married sample
data$msample <- data$aged != "" & data$timesmar == "1" & data$marital == "0" & 
  data$illegit == 0 & data$agefstd >= 15 & data$agefstm >= 15 

data <- subset(data, agem >= 21 & agem <= 35 & kidcount >= 2 & 
                 ageq2nd > 4 & agefstm >= 15 & asex == 0 & 
                 aage == 0 & aqtrbrth == 0 & asex2nd == 0 & 
                 aage2nd == 0 & msample == 1)

varlist <- c("blackm", "hispm", "othracem", "agem", "agefstm")

data <- data[names(data) %in% c("Y", "D", "Z", varlist)]

data$Z <- as.numeric(data$Z)
data$D <- as.numeric(data$D)

gc()

#### ESTIMATION ####

# First stage with no interactions
no_interactions <- lm(as.formula(paste0("D ~ Z + ", 
                                        paste(varlist, collapse = " + "))),
                      data = data)

summary(no_interactions)

# First stage with interactions
interactions <- matrix(0, nrow(data), length(varlist))
for (i in 1:length(varlist)){
  
  interactions[, i] <- getElement(data, varlist[i]) * data$Z
  
}
colnames(interactions) <- paste0("Z_", varlist)
interactions <- data.frame(interactions)

with_interactions <- lm(
  as.formula(
    paste0("D ~ Z + ", 
           paste(varlist, collapse = " + "),
           " + ",
           paste(colnames(interactions), collapse = "+"))
  ),
  data = cbind(data, interactions))

summary(with_interactions)

# find data such that the overall effect of Z is negative
coefs <- with_interactions$coefficients
Z0_coefs <- coefs[grep("Z", names(coefs), invert = T)]
Z1_coefs <- coefs[grep("Z", names(coefs))]

# change when Z = 1
Z_effect <- Z1_coefs[1] + 
  as.matrix(dplyr::select(data, all_of(varlist))) %*% 
  matrix(Z1_coefs[2:length(Z1_coefs)],
         nrow = length(Z1_coefs) - 1,
         ncol = 1)

# count of positive and negative
sum(Z_effect >= 0)
sum(Z_effect < 0)

#### WEIGHTS BY CELLS ####

# Add the referance group
data$whitem <- 1 - data$blackm - data$othracem - data$hispm

# get Z residuals 
Z_tilde <- lm(
  as.formula(
    paste0("Z ~ ", paste(varlist, collapse = " + "))
  ), data = data)$residuals
data$Z_tilde <- Z_tilde

# calculate the implied weight for each cell
agem_vals <- sort(unique(data$agem))
agem_dummies <- matrix(0, nrow(data), length(agem_vals))
for (i in 1:length(agem_vals)){
  
  agem_dummies[, i] <- as.numeric(data$agem == agem_vals[i])
  
}

agefstm_vals <- sort(unique(data$agefstm))
agefstm_dummies <- matrix(0, nrow(data), length(agefstm_vals))
for (i in 1:length(agefstm_vals)){
  
  agefstm_dummies[, i] <- as.numeric(data$agefstm == agefstm_vals[i])
}

counts <- rep(0, 4 * length(agem_vals) * length(agefstm_vals))
index <- 1

# loop over dummies and get the count
for (i in 1:3){
  for (j in 1:length(agem_vals)){
    for (k in 1:length(agefstm_vals)){
      
      counts[index] <- sum(data[varlist[i]] * 
                             agem_dummies[, j] * 
                             agefstm_dummies[, k]) 
      index <- index + 1
      
    }
  }
}

# add in the reference group
for (j in 1:length(agem_vals)){
  for (k in 1:length(agefstm_vals)){
    
    counts[index] <- sum(data$whitem * 
                           agem_dummies[, j] * 
                           agefstm_dummies[, k]) 
    index <- index + 1
    
  }
}

nrow(data) - sum(counts)
cells <- 4 * length(agem_vals) * length(agefstm_vals)

counts <- rep(0, cells)
numerator <- rep(0, cells)
denominator <- rep(0, cells)
index <- 1

data$negative <- rep(0, nrow(data))

for (i in 1:3){
  for (j in 1:length(agem_vals)){
    for (k in 1:length(agefstm_vals)){
      
      ind <- getElement(data, varlist[i]) * 
        agem_dummies[, j] * 
        agefstm_dummies[, k]
      count <- sum(ind) 
      
      if (count > 0){
        
        counts[index] <- count
        temp <- data[as.logical(ind), ]
        
        numerator[index] <- sum(temp$D * temp$Z_tilde)
        denominator[index] <- sum(temp$Z_tilde ^ 2)
        
        
      }
      
      if(numerator[index] < 0 ){
        data$negative <- data$negative + ind
      }
      
       index <- index + 1
    }
  }
}

for (j in 1:length(agem_vals)){
  for (k in 1:length(agefstm_vals)){
    
    ind <- data$whitem * 
      agem_dummies[, j] * 
      agefstm_dummies[, k]
    count <- sum(ind) 
    
    if (count > 0){
      
      counts[index] <- count
      temp <- data[as.logical(ind), ]
      
      numerator[index] <- sum(temp$D * temp$Z_tilde)
      denominator[index] <- sum(temp$Z_tilde ^ 2)
    
    }
    
    if(numerator[index] < 0 ){
      data$negative <- data$negative + ind
    }
    
    index <- index + 1
  }
}

sum(numerator) / sum(counts) - cov(data$D, data$Z_tilde)
sum(denominator) / sum(counts) - var(data$Z_tilde)

# Cell covariances
Z_effect <- numerator / counts
Z_effect <- Z_effect[counts > 0]
weight <- counts
weight <- weight[weight > 0]
weight <- weight / sum(weight)

data.frame(Z_effect = Z_effect, weight = weight) %>% ggplot() +
  geom_density(aes(x = Z_effect, weight =weight, col = "density")) + 
  geom_vline(aes(xintercept = 0, col = "line")) + 
  theme_bw() + 
  labs(y = "", x = expression(theta)) +
  theme(legend.position = "null") +
  xlim(c(-0.1, 0.1))

# weight: empirical frequency * covariance / variance of Z tilde 
lm(
  as.formula(
    paste0("D ~ Z + ", paste(varlist, collapse = " + "))
    ),
  data = data)

lm("D ~ Z_tilde - 1", data = data)

sum(numerator) / sum(denominator)

negative <- numerator < 0 

# implied mass of negative weights
sum(numerator[negative]) / sum(counts)
sum(numerator[!negative]) / sum(counts)

# normalised implied mass
sum(numerator[negative]) / sum(numerator)

# split out the two groups
data_neg <- dplyr::filter(data, negative == 1)

Z_tilde <- lm(
  as.formula(
    paste0("Z ~ ", 
           paste(varlist, collapse = " + ")
           )
    ), 
  data = data_neg)$residuals
  
num <- lm("Y ~ Z_tilde - 1", 
          data = data.frame(Y = data_neg$Y, 
                            Z_tilde = Z_tilde))$coefficients
denom <- lm("D ~ Z_tilde - 1", 
            data = data.frame(D = data_neg$D, 
                              Z_tilde = Z_tilde))$coefficients
LATE_neg <- num / denom
LATE_neg

data_pos <- dplyr::filter(data, negative == 0)

Z_tilde <- lm(
  as.formula(
    paste0("Z ~ ", 
           paste(varlist, collapse = " + ")
    )
  ), 
  data = data_pos)$residuals

num <- lm("Y ~ Z_tilde - 1", 
          data = data.frame(Y = data_pos$Y, 
                            Z_tilde = Z_tilde))$coefficients
denom <- lm("D ~ Z_tilde - 1", 
            data = data.frame(D = data_pos$D, 
                              Z_tilde = Z_tilde))$coefficients
LATE_pos <- num / denom
LATE_pos

# overall LATE
Z_tilde <- lm(
  as.formula(
    paste0("Z ~ ", 
           paste(varlist, collapse = " + ")
    )
  ), 
  data = data)$residuals

num <- lm("Y ~ Z_tilde - 1", 
          data = data.frame(Y = data$Y, 
                            Z_tilde = Z_tilde))$coefficients
denom <- lm("D ~ Z_tilde - 1", 
            data = data.frame(D = data$D, 
                              Z_tilde = Z_tilde))$coefficients
LATE_pos <- num / denom
LATE_pos

#### SATURATE AND CONTROL ####

# count the number of cells
cells <- 0

agem_vals <- sort(unique(data$agem))
agem_dummies <- matrix(0, nrow(data), length(agem_vals))
for (i in 1:length(agem_vals)){
  
  agem_dummies[, i] <- as.numeric(data$agem == agem_vals[i])
  
}

agefstm_vals <- sort(unique(data$agefstm))
agefstm_dummies <- matrix(0, nrow(data), length(agefstm_vals))
for (i in 1:length(agefstm_vals)){
  
  agefstm_dummies[, i] <- as.numeric(data$agefstm == agefstm_vals[i])
}

# loop over dummies and get the count
for (i in 1:3){
  for (j in 1:length(agem_vals)){
    for (k in 1:length(agefstm_vals)){
      
      cells <- cells + 
        as.numeric(
          sum(data[varlist[i]] * 
                agem_dummies[, j] * 
                agefstm_dummies[, k]) > 0)
      
    }
  }
}

counts <- rep(0, cells)
first_stage <- rep(0, cells)
second_stage <- rep(0, cells)
index <- 1

for (i in 1:3){
  for (j in 1:length(agem_vals)){
    for (k in 1:length(agefstm_vals)){
      
      ind <- getElement(data, varlist[i]) * 
        agem_dummies[, j] * 
        agefstm_dummies[, k]
      count <- sum(ind) 
      
      if (count > 0){
        
        counts[index] <- count
        temp <- data[as.logical(ind), ]
        
        first_stage[index] <- lm("D ~ Z", data = temp)$coefficients[2]
        second_stage[index] <- ivreg("Y ~ D|Z", data = temp)$coefficients[2]
        
        index <- index + 1
      }
    }
  }
}

# counts
check <- !is.na(first_stage)
counts <- counts[check]
first_stage <- first_stage[check]
second_stage <- second_stage[check]

check <- !is.na(second_stage)
counts <- counts[check]
first_stage <- first_stage[check]
second_stage <- second_stage[check]


sum(first_stage >= 0)
sum(first_stage < 0)

mean(second_stage[first_stage >= 0])
mean(second_stage[first_stage < 0])
