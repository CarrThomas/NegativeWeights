################################# FRONT MATTER #################################

# Cleans 1980 census data for the Angrist and Evans application

# Input should be the 1980 5% sample.
# I use an extract that is cut down to households with women aged between 21 and
# 35 who have at least two children.

# Variables are
# "YEAR": 1980 for all observations, can be dropped
# "SAMPLE": 1980 5% for all observations, can be dropped
# "SERIAL": household unique identifier. 
# "HHWT": household weight
# "HHTYPE": household type, 1 for married couple households
# "CLUSTER": cluster for variance estimation
# "STRATA": strata for variance estimation
# "GQ": indicator for household type, retain observations with GQ = 1 (housholds)
# "PERNUM": person identifier, unique within household. 
# "PERWT": person weight
# "FAMUNIT": indexes family units (relatives) within the same household (e.g. 
# 1 for the first set of relatives, 2 for the second...)
# "MOMLOC": for children, contains the "PERNUM" of their mother, 0 for anyone
# not in the same household as their mother. 
# "MOMRULE": rule used to match mothers and children
# "SPLOC": contains the "PERNUM" of their spouse    
# "SPRULE": rule used to match spouse
# "NCHILD": number of own children in the household
# "RELATE": relation to head of household (2 digits)
# "RELATED": detailed version of RELATED (4 digits)
# "SEX": 1 for male, 2 for female, 9 for missing/blank
# "AGE": 999 for missing
# "BIRTHQTR": quarter of birth (1 for Jan - Mar, 2,...0 for N/A, 9 for missing)
# "MARST": marital status (1 for married spouse present, 2 for married absent,
# 9 for blank/missing)
# "MARRNO": times married, 0 for N/A, 7,8,9 are unknown, illegible, missing)
# "MARRQTR": quarter of first marriage (0 for N/A, 8 for missing)
# "AGEMARR": age at first marriage (0 for missing)
# "CHBORN": children ever born (0 for N/A, 1 for 0, 2 for 1 etc ..., 97,98,99 
# for unknown, illegible, missing)
# "RACE": 1 for white, 2 for black, 3-9 for other
# "RACED": detailed race code (3 digit version of RACE)
# "HISPAN": hispanic origin (0 for not hispanic, 1-4 for hispanic, 9 for not 
# reported)
# "HISPAND": detailed 3 digit version of HISPAN
# "BPL": 3 digit code for birthplace. 0 - 99 for US, 100 - 120 for US 
# territories, 150 - 800 for other countries, 900 for at sea, 950 for other, 
# 997 for unknown, 999 for missing/blank
# "BPLD": 5 digit version of BPL
# "YRIMMIG": year of immigration to the US, 0000 for NA
# "EDUC": highest grade completed, 00 for N/A or none, 99 for missing
# "EDUCD": 3 digit version
# "EMPSTAT": employment status, 1 for employed, 2 for unemployed, 3 for not in 
# labour force, 0 for N/A, 9 for unknown/illegible
# "EMPSTATD": 2 digit version of empstat
# "LABFORCE": 1 for not in labour force, 2 for in labour force, 0 for N/A, 9 
# for unknown
# "WKSWORK1": weeks worked last year, 00 for N/A
# "UHRSWORK": usual hours worked in the past year, 00 for N/A, top code 99
# "INCWAGE": wage and salary income, 999999 for N/A, 999998 for missing, top
# code 75000
# "INCBUS": business income, 999999 for N/A, bottom code -09995 (=-9900), top 
# code 75000
# "INCFARM": farm income, 999999 for N/A, bottom code -09995 (=-9900), top
# code 7500

# Quality indicators
# "QAGE": 0 for unaltered, != for issues
# "QAGEMARR": 0 for unaltered, != for issues
# "QCHBORN": 0 for unaltered, != for issues   
# "QMARRNO": 0 for unaltered, != for issues    
# "QMARST": 0 for unaltered, != for issues    
# "QSEX": 0 for unaltered, != for issues      
# "QBIRTHMO": 0 for unaltered, != for issues  
# "QMARRQTR": 0 for unaltered, != for issues 
# "QBPL": 0 for unaltered, != for issues     
# "QHISPAN": 0 for unaltered, != for issues   
# "QRACE": 0 for unaltered, != for issues     
# "QRACE2": 0 for unaltered, != for issues    
# "QYRIMM": 0 for unaltered, != for issues    
# "QEDUC": 0 for unaltered, != for issues     
# "QEMPSTAT": 0 for unaltered, != for issues  
# "QOCC": 0 for unaltered, != for issues     
# "QUHRSWOR": 0 for unaltered, != for issues  
# "QWKSWORK1": 0 for unaltered, != for issues
# "QINCBUS": 0 for unaltered, != for issues   
# "QINCFARM": 0 for unaltered, != for issues 
# "QINCWAGE": 0 for unaltered, != for issues

# Codes for RELATE:
# 1         Head/Householder
# 2                   Spouse
# 3                    Child
# 4             Child-in-law
# 5                   Parent
# 6            Parent-in-Law
# 7                  Sibling
# 8           Sibling-in-Law
# 9               Grandchild
# 10          Other relatives
# 11 Partner, friend, visitor
# 12      Other non-relatives
# 13    Institutional inmates

# Codes for RELATED:
# 101                                             Head/Householder
# 201                                                       Spouse
# 301                                                        Child
# 401                                                 Child-in-law
# 501                                                       Parent
# 601                                                Parent-in-Law
# 701                                                      Sibling
# 801                                               Sibling-in-Law
# 901                                                   Grandchild
# 1011                                                  Grandparent
# 1021                                                Aunt or Uncle
# 1031                                                Nephew, Niece
# 1041                                                       Cousin
# 1061                                         Other relatives, nec
# 1113                           Residual category for non relative
# 1240                 Roomers/boarders/lodgers and foster children
# 1250                                                    Employees
# 1260                Other non-relatives (1990 includes employees)
# 1291                                                     Military
# 1292                                          College dormitories
# 1293                                  Residents of rooming houses
# 1294 Other non-inmate 1980 (includes employees and non-inmates in
# 1301                                        Institutional inmates

# Possible values for YRIMMIG
# 1949  1949 or earlier
# 1950  1950 - 1959
# 1960  1960 - 1964
# 1965  1965 - 1969
# 1970  1970 - 1974
# 1975  1975 - 1979

# Codes for educd in 1980
# 0                            N/A or no schooling
# 2                         No schooling completed
# 11                      Nursery school, preschool
# 12                                   Kindergarten
# 14                                        Grade 1
# 15                                        Grade 2
# 16                                        Grade 3
# 17                                        Grade 4
# 22                                        Grade 5
# 23                                        Grade 6
# 25                                        Grade 7
# 26                                        Grade 8
# 30                                        Grade 9
# 40                                       Grade 10
# 50                                       Grade 11
# 60                                       Grade 12
# 65             Some college, but less than 1 year
# 70                              1 year of college
# 80                             2 years of college
# 90                             3 years of college
# 100                             4 years of college
# 110                            5+ years of college
# 111           6 years of college (6+ in 1960-1970)
# 112                             7 years of college
# 113                            8+ years of college

rm(list = ls())
gc()

library(ipumsr)
library(dplyr)
library(tidyr)

path <- getwd()
ddi_file <- "/Applications/usa_00007.xml"
data_file <- "/Applications/usa_00007.csv.gz"

################################# READ IN DATA #################################

ddi <- read_ipums_ddi(paste0(path, ddi_file))
data <- read_ipums_micro(ddi, data_file = paste0(path, data_file))

names(data) <- tolower(names(data))

# drop some unnecessary variables
# "year" and "sample" dropped as all data is census 1980 5% sample, "strata" and
# "cluster" dropped as we do not calculate those standard errors
data <- data[!(names(data) %in% c("cluster", "strata", "year", "sample"))]
gc()

################################### CLEANING ###################################

#### TOTAL EARNED INCOME ####

# Construct inctot = incwage + incbus + incfarm and a qinctot

# step 1) indicator for bottom coding 
# incbus and incfarm are bottom coded at -9995 (which means -9990)
data$inctot_bot_code <- ((data$incbus == -9995) + (data$incfarm == -9995)) > 0

# step 2) indicator for top coding
# incwage, incbus, incfarm all top coded at 75000
data$inctot_top_code <- 
  ((data$incwage == 75000) + 
     (data$incbus == 75000) + 
     (data$incfarm == 75000)) > 0

# step 3) replace codes for missing values with NA
# 999999 is N/A for all three. 999998 is missing for incwage
data$incwage[data$incwage >= 999998] <- NA
data$incbus[data$incbus >= 999998] <- NA
data$incfarm[data$incfarm >= 999998] <- NA

# step 4) replace bottom codes for incbus/incfarm 
#  -9995 means bottom code, which is -9990
data$incbus[data$incbus == -9995] <- -9990
data$incfarm[data$incfarm == -9995] <- -9990

# step 5) generate inctot
data$inctot <- data$incbus + data$incfarm + data$incwage

# step 6) generate qinctot (will take value 1 when there is a flag)
data$qinctot <- as.numeric((data$qincwage + data$qincbus + data$qincfarm) > 0)

# remove unnecessary variables
data <- data[!(names(data) %in% c("incwage", "incbus", "incfarm",
                                  "qincwage", "qincbus", "qincfarm"))]

gc()

#### LABOUR MARKET OUTCOMES ####

# Replace "uhrswork" (usual hours worked) values with quality issues
data$uhrswork[data$quhrswor > 0] <- NA

# Generate the "worked" variable (wkswork1 > 0).
data$worked <- data$wkswork1 > 0
data$worked[data$qwkswork1 > 0] <- NA

# Replace "wkswork1" values with quality issues
data$wkswork1[data$qwkswork1 > 0] <- NA

# remove unnecessary variables
data <- data[!(names(data) %in% c("quhrswor", "qwkswork1"))]

# DEMOGRAPHICS

# Reduce the number of quality flags to natural combinations.
# Each generated flag takes a value of 1 whenever there is an issue with either
# variable. 

# Combine age and sex flags
data$flag_age_sex <- as.numeric((data$qage + data$qsex) > 0)

# Combine flags for race and hispan
data$flag_race <- as.numeric((data$qrace + data$qrace2 + data$qhispan) > 0)

# combined flag for yrimm and bpl (1 whenever there are quality issues with either)
data$flag_imm <- as.numeric((data$qbpl + data$qyrimm) > 0)

# remove unnecessary variables
data <- data[!(names(data) %in% c("qage", 
                                  "qsex", 
                                  "qhispan", 
                                  "qrace", 
                                  "qrace2", 
                                  "qbpl", 
                                  "qyrimm"))]

gc()

################################# BUILDING DATA ################################

n <- nrow(data)

# N.B. data starts off sorted by serial and the pernum.

# id number for everyone in the data
data$id <- 1:n

# Construct a household index. This is the index of the first observation for 
# each household. 
house_id <- data$id * c(TRUE, (data$serial[1:(n - 1)] < data$serial[2:n])) 
# second expression is TRUE for the first observation for each value of serial
for (i in 2:n){
  house_id[i] <- max(house_id[i], house_id[i - 1])
}
data$house_id <- house_id

# Construct a pointer for each person's mom (if available)
data$momloc[data$momloc == 0] <- NA
data$mom_id <- data$house_id + (data$momloc - 1)

# Construct a pointer for each person's spouse (if available)
data$sploc[data$sploc == 0] <- NA
data$sp_id <- data$house_id + (data$sploc - 1)

# Count children for each mother
data$kid_count <- tabulate(data$mom_id, nbins = n)

# Indicator for mother
data$is_mum <- data$kid_count > 0

# Indicator for spouse
data$is_spouse <- rep(FALSE, n)
data$is_spouse[data$sp_id[data$is_mum]] <- TRUE

# Indicator for child
data$is_child <- !is.na(data$momloc)

# Split into three datasets
mothers <- data[data$is_mum, ]
spouses <- data[data$is_spouse, ]
kids <- data[data$is_child, ]

n_kids <- nrow(kids)
n_mums <- nrow(mothers)
n_spouses <- nrow(spouses)

saveRDS(object = mothers, file = paste0(path, "/Applications/mothers.rds"))
saveRDS(object = spouses, file = paste0(path, "/Applications/spouses.rds"))
saveRDS(object = kids, file = paste0(path, "/Applications/kids.rds"))

rm(data)
gc()

##################################### KIDS #####################################

# kids <- readRDS(file = paste0(path, "/Applications/kids.rds"))

# Cut down to only the relevant variables
kids <- kids[names(kids) %in% 
               c("sex", "age", "birthqtr", "mom_id", "flag_age_sex")]

# 1980 Census: conducted on April 1 1980 -> quarter 2
# Given two children of age y
# Youngest will have birth quarter 1 (they turned y in the previous quarter)
# Oldest will have birth quarter 2 (they will turn y + 1 this quarter)
# Birth quarter 3 is second oldest
# Birth quarter 4 is third oldest
# Record birthqtr so that it is increasing with age
kids$birthqtr <- recode_values(kids$birthqtr, 
                               from = c(1,2,3,4),
                               to = c(1,4,3,2))

# Arrange by mom_id and then age (each set of children is sorted from oldest to
# youngest)
kids <- arrange(kids, mom_id, desc(age), desc(birthqtr))

# Cut down data to at most 3 kids per mum
kid_number <- rep(0, n_kids)
kid_number[1] <- 1
for (i in 2:n_kids){
  
  kid_number[i] <- 1 + 
    as.numeric(kids$mom_id[i] == kids$mom_id[i - 1]) * kid_number[i - 1]
    
}
kids$kid_number <- kid_number
kids <- kids[kids$kid_number <= 3, ]
rm(kid_number)

# Reshape to wide with data for kid1, kid2, kid3
kids <- reshape(data = data.frame(kids[, ]), 
                direction = "wide",
                idvar = "mom_id", 
                timevar = "kid_number",
                v.names = c("sex", "age", "birthqtr", "flag_age_sex"))

# Tidy up some variables
kids$age.2[is.na(kids$age.2)] <- 0
kids$sex.2[is.na(kids$sex.2)] <- 0
kids$birthqtr.2[is.na(kids$birthqtr.2)] <- 0
kids$flag_age_sex.2[is.na(kids$flag_age_sex.2)] <- 0

kids$age.3[is.na(kids$age.3)] <- 0
kids$sex.3[is.na(kids$sex.3)] <- 0
kids$birthqtr.3[is.na(kids$birthqtr.3)] <- 0
kids$flag_age_sex.3[is.na(kids$flag_age_sex.3)] <- 0

# Generate new variables

# same sex instrument
kids$same_sex <- (kids$sex.1 == kids$sex.2) 
# oldest over 18
kids$flag_18_plus <- (kids$age.1 >= 18) 
# twins at second birth
kids$multi2nd <- (kids$age.2 == kids$age.3) * 
  (kids$birthqtr.2 == kids$birthqtr.3) 
# triplets at first birth
kids$triplets <- (kids$age.1 == kids$age.2) * 
  (kids$birthqtr.1 == kids$birthqtr.2) * kids$multi2nd
# kid variables combined quality flags
kids$flag_qkids <- as.numeric(
  (kids$flag_age_sex.1 + kids$flag_age_sex.2 + kids$flag_age_sex.3) > 0)


# Cut down to what we need (NEED TO RETAIN AGE OF THE OLDEST CHILD!)
kids <- kids[names(kids) %in% c("mom_id", 
                                "same_sex", 
                                "flag_18_plus", 
                                "multi2nd",
                                "triplets",
                                "flag_qkids",
                                "age.1",
                                "birthqtr.1")]
names(kids)[1] <- "id"

saveRDS(object = kids, file = paste0(path, "/Applications/kid_data.rds"))

#################################### SPOUSES ###################################

# spouses <- readRDS(file = paste0(path, "/Applications/spouses.rds"))

# Want bpl, years of immigration and quality flags
spouses <- 
  spouses[names(spouses) %in% c("sp_id", "bpl", "yrimmig", "flag_imm")]

spouses <- spouses[c("sp_id", "bpl", "yrimmig", "flag_imm")]
names(spouses) <- c("id", "sp_bpl", "sp_yrimmg", "sp_flag")

##################################### MERGE ####################################

# mothers <- readRDS(file = paste0(path, "/Applications/mothers.rds"))

# Every "mother" record has a "kids" record, so these too can just be combined
mothers <- cbind(mothers, kids)

# Drop mothers without spouses records
mothers <- mothers[!is.na(mothers$sploc), ]

# Sort mothers and spouses (just in case)
mothers <- mothers[order(mothers$id), ]
spouses <- spouses[order(spouses$id), ]

mothers <- cbind(mothers, spouses)

rm(kids, spouses)

saveRDS(object = mothers, file = paste0(path, "/Applications/mom_data.rds"))

########################### ADD VARIABLES AND TIDY UP ##########################

# mothers <- readRDS(file = paste0(path, "/Applications/mom_data.rds"))

# drop unused columns
keep <- c("nchild", "age", "birthqtr", "marst", "marrno", "agemarr", "chborn", 
          "race", "hispan", "bpl", "yrimmig", "educ", "educd", "wkswork1", 
          "uhrswork", "qagemarr", "qchborn", "qmarrno", "qmarst","qbirthmo", 
          "qmarrqtr", "qeduc", "qempstat", "qocc", "inctot_bot_code", 
          "inctot_top_code", "inctot", "qinctot", "worked", "flag_age_sex", 
          "flag_race", "flag_imm", "kid_count", "age.1", "birthqtr.1", 
          "same_sex", "flag_18_plus", "multi2nd", "triplets", "flag_qkids", 
          "sp_bpl", "sp_yrimmg", "sp_flag")         

mothers <- mothers[names(mothers) %in% keep]

gc()

# drop iffy observations
mothers <- mothers[mothers$nchild == (mothers$chborn - 1), ]
# number of children in household should equal children born 
# (N.B. chborn is 0 for NA, 1 for 0, 2 for 1 etc...)
mothers <- mothers[mothers$nchild == mothers$kid_count, ]
# reported number of children should equal number actually in sample

# cut down to mothers aged 21 - 35
mothers <- dplyr::filter(mothers, age >= 21, age <= 35)

# oldest child must < 18
mothers <- dplyr::filter(mothers, flag_18_plus == 0)

# Aged > 16 when married and only married once
mothers <- dplyr::filter(mothers, agemarr > 15, marrno == 1)

# drop any observations with quality issues
mothers <- mutate(mothers, flag = qagemarr +
                    qchborn +
                    qmarst +
                    flag_age_sex +
                    flag_race +
                    flag_imm +
                    flag_qkids +
                    sp_flag)

mothers <- dplyr::filter(mothers, flag == 0)

# drop all the quality variables
mothers <- 
  mothers[!(names(mothers) %in% c("flag", 
                                  "qagemarr", 
                                  "qchborn", 
                                  "qmarst",
                                  "flag_age_sex", 
                                  "flag_race", 
                                  "flag_imm", 
                                  "flag_qkids", 
                                  "sp_flag"))]

# drop all triplet first births and all twin second births
mothers <- dplyr::filter(mothers, triplets == 0, multi2nd == 0)

# define the treatment
mothers$morekids <- mothers$nchild > 2

mothers <- 
  mothers[!(names(mothers) %in% c("marst", "marrno", "multi2nd",
                                  "triplets", "flag_18_plus", "qmarrno",
                                  "qbirthmo", "qmarrqtr", "nchild", "chborn",
                                  "kid_count"))]

# work out age at first birth

# recode birthqtr.1 (redefined it earlier so that kids would be sorted properly)
mothers$birthqtr.1 <- recode_values(mothers$birthqtr.1, 
                                    from = c(1, 2, 3, 4),
                                    to = c(1, 4, 3, 2))

# mothers and childs age in quarters:
# 4 * age in years plus 
# 0 if born in quarter 1
# 3 if born in quarter 2
# 2 if born in quarter 3
# 1 if born in quarter 4
# (as the census is April 1 1980)

# recode birthqtr and birthqtr.1
mothers$birthqtr <- recode_values(mothers$birthqtr,
                                  from = c(1, 2, 3, 4),
                                  to = c(5, 2, 3, 4))
mothers$birthqtr.1 <- recode_values(mothers$birthqtr.1,
                                    from = c(1, 2, 3, 4),
                                    to = c(5, 2, 3, 4))

mothers <- mutate(mothers, age_quarters =  age * 4 + 5 - birthqtr)
mothers <- mutate(mothers, age_quarters.1 = age.1 * 4 + 5 - birthqtr.1)
mothers <- mutate(mothers, age_fst = floor((age_quarters - age_quarters.1) / 4))

# View(mothers[c("age", "birthqtr", "age_quarters", 
#                "age.1", "birthqtr.1", "age_quarters.1",
#                "age_fst")])

# Drop any mothers whose first child was born before marriage
mothers <- dplyr::filter(mothers, age_fst >= agemarr)

# drop extra variables
mothers <- dplyr::select(mothers, 
                         -starts_with("birthqtr"), 
                         -starts_with("age_quarters"),
                         -age.1,
                         -agemarr)

# form race variables (as in Angrist and Evans)
mothers <- mutate(mothers, black = (race == 2), oth = (race > 2))
mothers$hispan <- (mothers$hispan > 0)
mothers <- select(mothers, - race)

# coarsen bpl and sp_bpl: codes < 100 are mainland US, count these as natives
mothers$bpl <- as.numeric(mothers$bpl) * (as.numeric(mothers$bpl) >= 100)
mothers$sp_bpl <- as.numeric(mothers$sp_bpl) * (as.numeric(mothers$sp_bpl) >= 100)
# zero out any where the year of immigration is 0
mothers$bpl <- mothers$bpl * (mothers$yrimmig != 0)
mothers$sp_bpl <- mothers$sp_bpl * (mothers$sp_yrimmg != 0)

mothers <- dplyr::select(mothers, - yrimmig, - sp_yrimmg)

# construct years of education from educd (detailed education)
mothers <- mutate(mothers, yr_ed = educd)
mothers$yr_ed <- recode_values(mothers$yr_ed,
                               from = c(2, 11, 12, 14, 15, 16, 17, 22, 23, 25, 
                                        26, 30, 40, 50, 60, 65, 70, 80, 90, 
                                        100, 110, 111, 112, 113),
                               to = c(0,0,0, 1:12, 12, 13:20))
# replace with NA if qeduc > 0
mothers$yr_ed[mothers$qeduc > 0] <- NA
mothers <- dplyr::select(mothers, - educ, - educd, - qeduc)
mothers <- dplyr::select(mothers, - qempstat, -qocc)

# inflate incomes to 1990 dollars to be comparable to Angrist Evans results
mothers$inctot <- mothers$inctot * 2.099173554
mothers$inctot[mothers$qinctot > 0] <- NA
mothers <- dplyr::select(mothers, -qinctot)

save(mothers, file =  paste0(path, "/Applications/1980_data.rds"))
