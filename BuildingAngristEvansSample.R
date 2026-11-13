################################# BUILDING DATA ################################

# Using Census data to build a sample similar to Angrist and Evans

# Census data is structured as follows
# * SERIAL is the household ID
# * NFAMS is the number of families in the household
# * FAMUNIT codes families within a household (1 for family 1, 2 for family 2,...)
# * FAMSIZE codes the number of family members living in that household
# * RELATE codes how members of the household are related: (1 is head of the 
# house/householder, 2 is spouse of the householder, 3 is child of the 
# householder...)
# * NCHILD number of own children residing with each individual (includes step
# children and adopted children)
# * SEX  (1 for male, 2 for female)
# * AGE (999 for missing)
# * BIRTHYR year of birth 
# * BIRTHQTR quarter of birth
# * MARST marital status (1: married spouse present, 2: married spouse absent,
#   3 separated, 4 divorced, 5 widowed, 6 never married, single)
# * MARRNO times married (1980 only)
# * MARRQTR quarter of first marriage (1980 only)
# * AGEMARR age at first marriage (1980 only)
# * EDUC educational attainment (level COMPLETED, incomplete doesn't count!)
#   - 00 NA or no schooling 
#   - 01 nursery to grade 4
#   - 02 grade 5 to 8
#   - 03, 04, 05, 06 are grade 9, 10, 11, 12
#   - 07, 08, 09, 10 are 1, 2, 3, 4 years of college
#   - 11 is 5+ years of college

# * MOMLOC indicator of whether mom is in the household, gives her PERNUM if so
# 

# Outcome variables
# * EMPSTAT employment status (0 = N/A, 1 = employed, 2 = unemployed, 3 = not in
# labour force)
# * CLASSWKR type of worker (0 = N/A, 1 = self-employed, 2 = works for wages)
# * WKSWORK1 weeks worked in the last year (00 = N/A) [WKSWORK2 is the same 
# thing coded in intervals]
# * HRSWORK1 hours worked in the last week (00 = N/A, top code = 99), available
# in 1980 and 1990, but not 2000 [HRSWORK2  is the same thing coded in intervals]
# * UHRSWORK usual hours worked in the last year (00 = N/A, top code = 99)
# * FTOTINC total family income, available all years, rounded, may be top-coded
# * INCTOT total pre-tax personal income, available all years, rounded, may be 
#   top-coded

# To match mothers and children 

library(dplyr)
library(tidyr)
library(ggplot2)

path <- getwd()
data <- read.csv(file = paste0(path, "/usa_00002.csv"))
# SERIAL in the main sample ranges from 1 to 4711341

# sex of first child and sex of second child
# 1) loop over SERIAL (household indicator)
# 2) loop over FAMUNIT (within serial)
# 3) for each FAMUNIT, sexk is SEX of the oldest child, sex2nd is SEX for the 
# second oldest. 

# construct the child variables (age of first child, age of second child)
# drop all households without at least two children?
data <- dplyr::filter(data, FAMSIZE >= 4)
gc()

names(data) <- tolower(names(data))
data$sexk <- rep(0, nrow(data))
data$sex2nd <- rep(0, nrow(data))

hholds <- unique(data$serial)

for (i in 1:length(hholds)){
  
  temp <- dplyr::filter(data, serial == hholds[i])
  temp <- temp[order(temp$age), ]
  nkids <- 
  oldest <- which.max(temp$age * (temp$relate == 3))
  
  
}

################################# SHOPPING LIST ################################ 

# PERSON NUMBER (PERNUM) UNIQUE IDENTIFIER
# 

# OUTCOMES
# EMPLOYMENT STATUS (EMPSTAT) 1 FOR EMPLOYED
# HOURS WORKED (HRSWORK1)
# USUAL HOURS WORKED (UHRSWORK)
# WEEKS WORKED (WKSWORK1)
# INCOME 
# After 1990, use INCEARN which is INCWAGE + INCBUS + INCFARM in 1990 and 
# INCWAGE + INCBUS00 from 2000 on. 
# In 1980, use INCWAGE + INCBUS + INCFARM

# TREATMENT
# NUMBER OF CHILDREN

# INSTRUMENT
# SEX (ALL)

# OTHERS
# AGE
# AGE IN QUARTERS

# RACE DUMMIES (BLACK, WHITE, HISPANIC, OTHER)
# COUNTRY OF BIRTH
# YEARS OF EDUCATION? 

# MOMLOC (PERSON NUMBER OF PROBABLE MOTHER)
# THEIR CONSTRUCTION OF `MOTHER'


#### SAMPLE DEFINITIONS ####

# MOTHERS AGED 21 - 35
# OLDEST CHILD < 18
# MARRIED SAMPLE AND UNMARRIED SAMPLE
# 1980 MARRIED: MARRIED AT THE TIME OF THE CENSUS, MARRIED ONLY ONCE, MARRIED
# AT THE TIME OF THE FIRST BIRTH
# 1990 MARRIED: ALL WOMEN MARRIED AT THE TIME OF THE CENSUS
# ONLY 1980 HAS DATE OF FIRST MARRIAGE! 


# VARIABLES TO DOWNLOAD
# serial, pernum

# gq
# age, qage
# sex, qsex
# race
# hispan
# labforce
# empstat
# uhrswork
# wskwork1
# incearn (for 1980, incwage + incbus + incfarm)

# chborn
# nchild
# educ

# marst
# momrule
# momloc

# perwt
# hhwt

# ADDITIONAL
# nativity
# bpl, qbpl (joint flag for nativity and birthplace quality)
# yrimmig

