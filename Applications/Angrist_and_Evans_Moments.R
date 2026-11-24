########################### ANGRIST AND EVANS MOMENTS ##########################

# Get moments from Angrist and Evans (1998) to calibrate data

rm(list = ls())
gc()

library(ggplot2)
library(dplyr)
library(tidyr)
library(wooldridge)
library(haven)
library(AER)

################################ DATA WRANGLING ################################

path <- getwd()
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
data$log_faminc1 <- log(pmax(data$faminc1 , 1))

# indicator for whether the observation belongs in the married sample
data$msample <- data$aged != "" & data$timesmar == "1" & data$marital == "0" & 
  data$illegit == 0 & data$agefstd >= 15 & data$agefstm >= 15 

data <- subset(data, agem >= 21 & agem <= 35 & kidcount >= 2 & 
                 ageq2nd > 4 & agefstm >= 15 & asex == 0 & 
                 aage == 0 & aqtrbrth == 0 & asex2nd == 0 & 
                 aage2nd == 0 & msample == 1)

varlist <- c("blackm", "hispm", "othracem", "agem", "agefstm")
Ylist <- c("faminc1", "log_faminc1", "weeksd", "hoursd")

data <- data[names(data) %in% c("Y", "D", "Z", varlist)]

data$Z <- as.numeric(data$Z)
data$D <- as.numeric(data$D)

gc()

################################## GET MOMENTS #################################

moments <- list()

# Orthogonalise Z
Z_hat <- lm(as.formula(paste0("Z ~ ", paste(varlist, collapse = " + "))),
            data = data)$fitted.values
Z_tilde <- data$Z - Z_hat

moments$tsls_denom <- sum(data$D * Z_tilde) / nrow(data)
moments$Z_var <- var(data$Z)

sum(data$Y * Z_tilde) / sum(data$D * Z_tilde)
