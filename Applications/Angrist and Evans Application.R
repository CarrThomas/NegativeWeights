################################# FRONT MATTER #################################

# Implements the estimator with the data from Angrist and Evans (1998)
# Largely follows the replication files of Aaranson et al (2020)

library(ivreg)
library(xtable)

rm(list = ls())

path <- getwd()

################################ DATA WRANGLING ################################

#     ** inflation multiplier for post-1940 data (to 1990 dollars)
#     local infl1980 = 1.59
#     local infl1990 = 1
#     local infl2000 = 0.76
#     local infl2010 = 0.60
#     
#   # Data Cleaning
#
# 	** flag for group quarters
# 	gen fgq = gq!=1
# 
# 	** flags for age and gender
# 	gen fage = age==.|age<0|qage!=0
# 	gen fsex = (sex!=1&sex!=2)|qsex!=0
# 
#   ** Houehold flags (Any Imputation, GQ, and Inmate Relationships)
#   gen fhh_t = fage==1|fsex==1|fgq==1
# 
# 	** Reprocess race
# 	gen black = race==2
# 	gen other_race = race>2
# 
# 	** Process birth quarter/month. For US, quarter of birth is avilable for
# 		* 1960-1980 and ACS Samples (2010). For month, this is avilable for 1900.
# 		* The variables are available for other years, but these only provide
# 		* data for individuals with age code of 0 (<1 year old).
# 	gen bm = .
# 	if (`yr'>=1960 & `yr'<=1980 | `yr'>2000) {
# 		replace bm = birthqtr if birthqtr!=9
# 	}
# 
# ** labor force participation (all years)
# replace labforce = . if labforce!=1 & labforce!=2
# replace labforce = labforce==2 if labforce!=.
# 
# 
# 	** employed (standard employment status, only after 1960)
# 	cap gen empstat = . // generate if missing
# 	replace empstat = . if empstat!=1 & empstat!=2 & empstat !=3
# 	replace empstat = empstat==1 if empstat != .
# 
# 	** grandparents
# 	gen grandma = momloc!=0 if momloc!=.
# 	gen grandpa = poploc!=0 if poploc!=.
# 
# 	** children ever born
# 	cap gen chborn = .
# 	replace chborn = . if chborn >= 98
# 	replace chborn = . if chborn==0 // first code is missng, rest are off by 1
# 	replace chborn = chborn - 1
# 
# 	** flag for if number of children ever born not equal to number of children
# 	gen fchborn = chborn!=nchild if chborn!=.
# 	replace fchborn = 0 if fchborn==1 & chborn > 9 & nchild==9 // nchild topcoded at 9
# 
# 	** code the education variable
# 	gen education = .
# 	if `yr' >= 1940 & `yr' <= 1980 {
# 		replace education = 0 if educd == 2 | educd == 12
# 		replace education = 1 if educd == 14
# 		replace education = 2 if educd == 15
# 		replace education = 3 if educd == 16
# 		replace education = 4 if educd == 17
# 		replace education = 5 if educd == 22
# 		replace education = 6 if educd == 23
# 		replace education = 7 if educd == 25
# 		replace education = 8 if educd == 26
# 		replace education = 9 if educd == 30
# 		replace education = 10 if educd == 40
# 		replace education = 11 if educd == 50
# 		replace education = 12 if educd == 60
# 		replace education = 13 if educd == 70
# 		replace education = 14 if educd == 80
# 		replace education = 15 if educd == 90
# 		replace education = 16 if educd == 100
# 		replace education = 17 if educd == 110
# 		replace education = 18 if educd == 111
# 		replace education = 19 if educd == 112
# 		replace education = 21 if educd == 113
# 	}
# 	else if `yr' >= 1990 {
# 	  replace education = 0 if educd == 2 | educd == 12 | educd == 11
# 	  replace education = 2 if educd == 13
# 	  replace education = 2 if educd == 10
# 	  replace education = 6 if educd == 20
# 	  replace education = 5.5 if educd == 21
# 	  replace education = 7.5 if educd == 24
# 	  replace education = 9 if educd == 30
# 	  replace education = 10 if educd == 40
# 	  replace education = 11 if educd == 50
# 	  replace education = 11.5 if educd == 61
# 	  replace education = 12 if educd == 62
# 	  replace education = 12.5 if educd == 65
# 	  replace education = 13 if educd == 71
# 	  replace education = 14 if educd == 82 | educd == 83 | educd == 81
# 	  replace education = 16 if educd == 101
# 	  replace education = 18 if educd == 114 | educd == 115
# 	  replace education = 21 if educd == 116
# 	}
# 	
# 	** married
# 		gen married = marst==1 if marst>=1&marst<=6
# 		gen evermarried = married==0 & marst!=6 if marst!=. & marst!=9
# 		gen marr_divorce = marst==4 if marst!=9 & marst!=.
# 		gen marr_separate = marst==3 if marst!=9 & marst!=.
# 		gen marr_single = marst==6 if marst!=9 & marst!=.
# 		gen marr_widow = marst==5 if marst!=9 & marst!=.
# 
# 	** clean wage and income data
# 	  sum incwage
# 	  replace incwage = . if incwage<=0
# 	  replace incwage = . if incwage>=999998
# 	  replace incwage = incwage * `infl`yr''
# 	  
# 	  ** weeks worked
# 	  gen weeks = .
# 	  cap replace weeks = wkswork1
# 	  
# 	  ** hours
# 	  gen hours = .
# 	  cap {
# 	    if (`yr'==2010|`yr'==2000) {
# 				replace hours = uhrswork
# 			}
# 			else {
# 				replace hours = hrswork1
# 			}
# 		}
# 
# 		** hourly wages
# 		gen wage = incwage / (hours*weeks)
# 		replace wage = . if wage==0
# 		_pctile wage, p(1 99)
# 		replace wage = r(r1) if wage<=r(r1) & wage != .
# 		replace wage = r(r2) if wage>=r(r2) & wage != .
# 	}
# 	else {
# 		gen incwage = .
# 		gen weeks = .
# 		gen hours = .
# 		gen wage = .
# 	}
# 
# /* Weights */
#   ** modification to the weights
# replace perwt = . if perwt <=0
# replace hhwt = . if hhwt <=0
# 
# 
# /* Household Data */
#   sort sample statecode serial pernum
# 
# ** number of older women, 16+
#   gen t_var = sex==2 & age > 15
# by sample statecode serial: egen nfemale_p15 = sum(t_var)
# drop t_var
# 
# ** number of older women, 16+
#   gen t_var = sex==2 & age > 15 & labforce==1
# by sample statecode serial: egen nfemale_lab_p15 = sum(t_var)
# drop t_var
# 
# ** number of older women, 36+
#   gen t_var = sex==2 & age >= 36
# by sample statecode serial: egen nfemale_p36 = sum(t_var)
# drop t_var
# 
# ** number of older women in labor force, 36+
#   gen t_var = sex==2 & age >= 36 & labforce==1
# by sample statecode serial: egen nfemale_lab_p36 = sum(t_var)
# drop t_var
# 
# ** indicator for multigenerational household
# cap drop multgen
# gen	t_var = nchild>=1 & (grandpa==1|grandma==1)
# by sample statecode serial : egen multgen = max(t_var)
# drop t_var
# 
# ** flag for household data quality
# by sample statecode serial : egen fhh = max(fhh_t)
# 
# 
# /* Variables to Keep */
#   ** be sure this variable is available
# cap gen momrule = 0
# cap gen poprule = 0
# 
# ** these are variables we keep. this step forces them to exist in all samples
# keep source sample country year statecode serial pernum perwt hhwt age sex ///
#   labforce labforceb manuf boarder empstat occupied workedyr wage_worker occagg occprofess ///
#   hours weeks wage incwage education occscore attend_school ///
#   grandma grandpa multgen bpl moved urban_cat momrule poprule labforceb boarder ///
#   momloc poploc sploc marst married evermarried occupation ///
#   fgq fage fsex fhh fchborn chborn nchild bm ///
#   black other_race south deep_south nfemale_* ///
#   marr_divorce marr_separate marr_single marr_widow
# qui compress

# /* Extract Mothers */
# use if age>=16 & age<=45 & sex==2 & nchild!=0 ///
# 
# 
# 
# 
# /* Extact Children */
#   ** load the data
# clear
# use if momloc!=0 using data/temp_`f'_data
# 
# ** sort into family groups, then descending by age. This will put
# 	* the oldest kids at the top. Save this ordering
# gsort sample statecode serial momloc -age pernum
# by sample statecode serial momloc: gen kid_idx = _n
# 
# ** store how many kids there are, total, for a given mother, but only
# 	* keep around the first three for the reshape to save time
# keep sample statecode serial kid_idx momloc sex age fage fsex ///
#      attend_school bm momrule
# sort sample statecode serial momloc
# by sample statecode serial momloc: egen num_kids_found = max(kid_idx)
# by sample statecode serial momloc: egen youngest_age = min(age)
# drop if kid_idx>6
# 
# ** reshape so the age and gender of each kid is its own variable
# reshape wide ///
# sex age fage fsex youngest_age attend_school bm momrule, ///
# i(sample statecode serial momloc) j(kid_idx)
# rename youngest_age1 youngest_age
# drop youngest_age2 youngest_age3 youngest_age4 ///
# youngest_age5 youngest_age6
# 
# ** merge to moms
# rename momloc pernum
# sort sample statecode serial pernum
# merge 1:1 sample statecode serial pernum using data/moms_`f', ///
#   keep(2 3)
# replace num_kids_found = 0 if _merge==2
# drop _merge
# 
# ** Now we can compare the reported number of children with the number
# * of children  who were observed and drop if they don't match.
# 	* NAPP Top-codes nchild at 9.
# if ("`f'"=="napp" ) {
# 	keep if  ///
# 	nchild==num_kids_found | (num_kids_found > 9 & nchild==9)
# }
# else {
# 	keep if nchild==num_kids_found
# }
# 
# 
# /* Additional Variables */
# **Find age at first birth, spouse's age at first birth.
# 	* Will drop mothers that are unusually long
# 	* before analysis.
# gen age_first_birth = age - age1
# gen spouse_age_first_birth = spouse_age - age1
# 
# **Create indicators for whether a woman's first two children are of
# 	* the same sex
# gen same_sex = (sex1==sex2) if nchild>=2
# gen two_boys = (sex1==1 & sex2==1) if nchild>=2
# gen two_girls = (sex1==2 & sex2==2) if nchild>=2
# 
# ** Indicators for twins and "real" twins for each child. Only
# 	* populate month/quarter of birth twins if we
# 	* can observe this.
# forval i=1/5 {
# 	local j = `i' + 1
# 	gen twins`i' = age`i'==age`j' if nchild>=`i'
# 	gen rtwins`i' = twins`i'==1 & bm`i'==bm`j' ///
# 	if nchild>=`i' & bm`i'!=.
# }
# 
# ** Baseline Treatment
# gen fertility = nchild >=3 if nchild >= 2
# 
# ** Count number of twins
# foreach v in twins rtwins {
# 	gen n`v' = 0
# 	forval i=1/6 {
# 		local previ = `i' - 1
# 		if (`i'==1) {
# 			replace n`v'=n`v' + 1 if `v'1==1
# 		}
# 		else if (`i'==6) {
# 			replace n`v'=n`v' + 1 if `v'5==1
# 		}
# 		else {
# 			replace n`v'=n`v' + 1 if `v'`i'==1
# 			replace n`v'=n`v' + 1 if ///
# 			`v'`i'==0 & `v'`previ'==1
# 		}
# 	}
# }
# replace nrtwins = . if bm1==.
# 
# ** flag for family with triplets
# forval i = 1/4 {
# 	local i1 = `i' + 1
# 	local i2 = `i' + 2
# 	gen ftrip`i' = age`i'==age`i1' & age`i1'==age`i2'
# }
# gen fmult_twin = ntwins>2
# 
# ** indicator for baseline sample. Useful for testing
# gen _samp = age>=21 & age<=35 & nchild>=2 & fage1==0 & fage2==0 & ///
# fsex1==0 & fsex2==0 & fgq==0 & ftrip1==0 & fage==0 & fsex==0 & ///
# labforce!=. & hhwt!=.
# replace _samp = 0 if age_first_birth < 15 | age1 > 17
# replace _samp = 0 if nchild>=3 & (fage3!=0|fsex3!=0)
# 
# 
# 
# /* Finalize Data */
# ** keep these variables
# keep sample source statecode year serial pernum _samp ///
# perwt hhwt fertility nchild age age_first_birth youngest_age ///
# same_sex two_boys two_girls black other_race workedyr education ///
# empstat occupied labforce spouse spouse_* occagg chborn fchborn ///
# urban_cat twins* rtwins* bm* age1 age2 age3 age4 age5 age6 ///
# sex1 sex2 sex3 sex4 sex5 sex6 fage* fsex* marst married momrule* ///
# attend_school1 attend_school2 attend_school3 attend_school4 ///
# attend_school5 attend_school6 nchild fgq fhh fmult_twin ftrip* ///
# incwage hours weeks wage wage_worker nfemale_* ///
# labforceb manuf boarder occprofess occupation grandma grandpa ///
# multgen evermarried south deep_south occscore bpl moved ///
# oldersister12 oldersister15 ntwins nrtwins marr_divorce ///
# marr_separate marr_single marr_widow
# 
# ** merge with sample aggregates if we are working with full
# 	* samples. If we are using full country states, we will
# 	* merge in later to save time.
# merge m:1 sample statecode using ///
# 	data/samp_stat/`f'_summary, ///
#   nogen assert(2 3) keep(3)
# 



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

# outcome variables 
data$faminc1 <- as.numeric(data$faminc)
data$faminc1 <- data$faminc1 * 2.099173554
# data$Y <- log(pmax(data$faminc1 , 1))
data$weeksm1 <- as.numeric(data$weeksm)
data$workedm <- as.numeric(data$weeksm1 > 0)
data$hourswm <- as.numeric(data$hoursm)
data$incomem <- as.numeric(data$income1m) + pmax(0, as.numeric(data$income2m))
data$incomem <- data$incomem * 2.099173554

# indicator for whether the observation belongs in the married sample
data$msample <- data$aged != "" & data$timesmar == "1" & data$marital == "0" & 
  data$illegit == 0 & data$agefstd >= 15 & data$agefstm >= 15 

data <- subset(data, agem >= 21 & agem <= 35 & kidcount >= 2 & 
                 ageq2nd > 4 & agefstm >= 15 & asex == 0 & 
                 aage == 0 & aqtrbrth == 0 & asex2nd == 0 & 
                 aage2nd == 0 & msample == 1)

Ylist <- c("weeksm1", "workedm", "hourswm", "incomem")
varlist <- c("blackm", "hispm", "othracem", "agem", "agefstm")

data <- data[names(data) %in% c("D", "Z", Ylist, varlist)]

data$Z <- as.numeric(data$Z)
data$D <- as.numeric(data$D)

gc()

################################### FUNCTIONS ##################################

calc_sd <- function(data, Z_tilde, tau_feasible, s_hat, X_names, Y_name){
  
  # returns the estimate of the asymptotic SD
  
  # data: a data frame containing the data (including generated variables)
  # Z_tilde: Z - E[Z|X]
  # tau_feasible: the estimate (from the feasible estimator)
  # s_hat: the estimated sign function
  # X_names: list of covariate names
  
  n <- nrow(data)
  K <- length(X_names)
  
  # "errors"
  e_Y <- data[[Y_name]] - data$D * tau_feasible
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

get_s_hat <- function(theta_hat, data){
  
  relevant <- names(theta_hat)[grep("Z", names(theta_hat))]
  theta_hat <- theta_hat[relevant]
  
  relevant <- relevant[2:length(relevant)]
  relevant <- gsub("Z_", "", relevant)
  k <- length(relevant)
  
  return(sign(theta_hat["Z"] + as.matrix(data[relevant]) %*% 
                matrix(theta_hat[1 + 1:k], k, 1)))
  
}

##################################### TSLS #####################################

tsls <- vector(mode = "list", length = 5)

for (i in 1:4){
  
  tsls_formula <-  as.formula(
    paste0(Ylist[i], " ~ D + ", 
           paste(varlist, collapse = " + "),
           "| Z + ",
           paste(varlist, collapse = " + ")))
  tsls[[i]] <- ivreg(formula = tsls_formula, 
                     data = data)
  
}

# make a table
tsls_table <- matrix(0, 2, 4)

for (i in 1:4){
  
  temp <- summary(tsls[[i]])$coefficients
  temp <- temp[grep("D", rownames(temp)), ]
  tsls_table[1, i] <- temp["Estimate"]
  tsls_table[2, i] <- temp["Std. Error"]
  
}

rownames(tsls_table) <- c("Estimate", "Std. Error")
xtable(tsls_table, digits = 3)

################################### REORDERED ##################################

# generate the orthogonalised instrument
Z_hat <- 
  lm(as.formula(
    paste0("Z ~ ", 
           paste(varlist, collapse = " + "))),
    data = data)$fitted.values
Z_tilde <- data$Z - Z_hat

interactions <- matrix(
  unlist(
    sapply(varlist, 
           function(X_name, data){return(data["Z"] * data[X_name])},
           data = data)),
  nrow = nrow(data), 
  ncol = length(varlist))
colnames(interactions) <- paste0("Z_", varlist)
interactions <- data.frame(interactions)
data <- cbind(data, interactions)

probit_results <- matrix(0, 2, 4)
lpm_results <- matrix(0, 2, 4)

for (i in 1:4){
  
  #### PROBIT ####
  
  # coefficient estimates
  probit_est <- glm(as.formula(
    paste0("D ~ Z + ", 
           paste(varlist, collapse = " + "),
           " + ",
           paste(names(interactions), collapse = " + "))),
    family = binomial(link = "probit"),
    data = data)
  
  # calculate S_hat 
  probit_s_hat <- get_s_hat(probit_est$coefficients, data)
  
  # get the estimator
  probit_tau <- sum(data[[Ylist[i]]] * probit_s_hat * Z_tilde) / 
    sum(data$D * probit_s_hat * Z_tilde)
  
  # standard error
  probit_se <- calc_sd(data, Z_tilde, probit_tau, probit_s_hat, varlist, Ylist[i]) / 
    sqrt(nrow(data))
  
  probit_results[1, i] <- probit_tau
  probit_results[2, i] <- probit_se
  
  #### LPM ####
  
  # coefficient estimates
  lpm_est <- lm(as.formula(
    paste0("D ~ Z + ", 
           paste(varlist, collapse = " + "),
           " + ",
           paste(names(interactions), collapse = " + "))),
    data = data)
  
  # calculate S_hat 
  lpm_s_hat <- get_s_hat(lpm_est$coefficients, data)
  
  # get the estimator
  lpm_tau <- sum(data[[Ylist[i]]] * lpm_s_hat * Z_tilde) / 
    sum(data$D * lpm_s_hat * Z_tilde)
  
  # standard error
  lpm_se <- calc_sd(data, Z_tilde, lpm_tau, lpm_s_hat, varlist, Ylist[i]) / 
    sqrt(nrow(data))
  
  lpm_results[1, i] <- lpm_tau
  lpm_results[2, i] <- lpm_se
  
}

