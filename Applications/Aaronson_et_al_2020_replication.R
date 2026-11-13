# Aaronson et al (2020) replication files


# /* ----- IPUMS US (Random Samples) ----- */
#   if ("`f'"=="1860"|"`f'"=="1870"|"`f'"=="1881"|"`f'"=="1901"|"`f'"=="1911"| ///
#       "`f'"=="1921"|"`f'"=="1931"|"`f'"=="1941"|"`f'"=="1950"|"`f'"=="1960"| ///
#       "`f'"=="1970"|"`f'"=="1980"|"`f'"=="1990"|"`f'"=="2000"|"`f'"=="2010") {
#     
#     ** inflation multiplier for post-1940 data (to 1990 dollars)
#     local infl1940 = 9.34
#     local infl1950 = 5.42
#     local infl1960 = 4.42
#     local infl1970 = 3.37
#     local infl1980 = 1.59
#     local infl1990 = 1
#     local infl2000 = 0.76
#     local infl2010 = 0.60
#     
#     ** source
#     gen source = 0
#     
#     ** in the good ol' US of A
# 	gen country = "USA"
# 
# 	** state we are in
# 	gen statecode = statefip
# 
# 	**Identify some notions of Southern states
# 	gen south = (statefip==10 | statefip==24 | statefip==11 | statefip==51 | ///
# 	statefip==54 |statefip==37 |statefip==45 |statefip==13 |statefip==12 | ///
# 	statefip==21 |statefip==47 |statefip==1 |statefip==28 |statefip==5 | ///
# 	statefip==22 |statefip==40 |statefip==48)
# 	gen deep_south = (statefip==45 |statefip==13 |statefip==12 |statefip==1 | ///
# 	statefip==28 |statefip==22 |statefip==48)
# 
# 	** flag for group quarters
# 	gen fgq = gq!=1
# 
# 	** flags for age and gender
# 	gen fage = age==.|age<0|qage!=0
# 	gen fsex = (sex!=1&sex!=2)|qsex!=0
# 
# 	** indicator for having boarders
# 	if (`yr'<1940) {
# 	  gen boarder_t = related>=1201 & related<=1205
# 	  bys statecode serial : egen boarder = max(boarder_t)
# 	}
# else {
#   gen boarder = 0
# }
# 
# ** Houehold flags (Any Imputation, GQ, and Inmate Relationships)
# * There are no missing/imputed relationships in these samples,
# * so we do not have to worry about missing relationship variables.
# gen fhh_t = fage==1|fsex==1|fgq==1
# 
# ** rename
# if (`yr'==1860) {
# 		gen occ=occ1950
# 	}
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
# 	if (`yr'==1900) {
#   replace bm = birthmo if birthmo!=99
# }
# 
# ** labor force participation (all years)
# replace labforce = . if labforce!=1 & labforce!=2
# replace labforce = labforce==2 if labforce!=.
# 
# ** occupation string derived
# cap {
#   replace occstr=lower(occstr)
#   replace occstr = "" if regexm(occstr,"school")==1|regexm(occstr,"student")==1|regexm(occstr,"none")==1|occstr=="n"|occstr=="*"
#   replace occstr = "" if regexm(occstr,"wife")
#   replace occstr = "" if regexm(occstr,"house")
#   gen occupieds = occstr!=""
# }
# if _rc != 0 {
#   gen occstr=""
#   gen occupieds = .
# }
# 
# ** occupation code derived
# gen occupiedc = occ1950<980 if occ1950!=.
# 
# ** combine occupation codes
# gen occupied = occupiedc
# if (`yr'<1940) {
# 		replace occupied = labforce
# 	}
# 
# 	** employed (standard employment status, only after 1960)
# 	cap gen empstat = . // generate if missing
# 	replace empstat = . if empstat!=1 & empstat!=2 & empstat !=3
# 	replace empstat = empstat==1 if empstat != .
# 
# 	** wage worker
# 	cap gen classwkr=.
# 	if _rc==0 {
# 		gen wage_worker=.
# 	}
# 	else {
# 		gen wage_worker = labforce==1 & classwkr==2 if labforce!=.
# 	}
# 
# 	** worked last year (2 is code specific to ACS)
# 	cap gen workedyr=.
# 	replace workedyr = . if workedyr != 3 & workedyr!= 1 & workedyr!= 2
# 	replace workedyr = workedyr == 3 if workedyr != .
# 
# 	** boarder labor force
# 	gen labforceb = labforce
# 	replace labforceb = 1 if boarder==1
# 
# 	** indicator for manufacturing worker
# 	gen manuf = ind1950>=306 & ind1950<=499
# 
# 	** urban/rural
# 	cap gen urban=.
# 	gen urban_cat = urban==2 if urban==1|urban==2
# 
# 	** attend school
# 	cap gen school=.
# 	gen attend_school = school==2 if school==1|school==2
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
# 	if (`yr'==1860|`yr'==1870) {
# 		gen married=.
# 		gen marst=.
# 		gen evermarried = .
# 		gen marr_divorce = .
# 		gen marr_separate = .
# 		gen marr_single = .
# 		gen marr_widow = .
# 	}
# 	else {
# 		gen married = marst==1 if marst>=1&marst<=6
# 		gen evermarried = married==0 & marst!=6 if marst!=. & marst!=9
# 		gen marr_divorce = marst==4 if marst!=9 & marst!=.
# 		gen marr_separate = marst==3 if marst!=9 & marst!=.
# 		gen marr_single = marst==6 if marst!=9 & marst!=.
# 		gen marr_widow = marst==5 if marst!=9 & marst!=.
# 	}
# 
# 	** clean wage and income data
# 	if (`yr'>=1940) {
# 	  sum incwage
# 	  replace incwage = . if incwage<=0
# 	  replace incwage = . if incwage>=999998
# 	  replace incwage = incwage * `infl`yr''
# 	  
# 	  ** weeks worked
# 	  gen weeks = .
# 	  cap replace weeks = wkswork1
# 	  if _rc!=0 {
# 	    replace weeks = 7  if wkswork2==1
# 	    replace weeks = 20 if wkswork2==2
# 	    replace weeks = 33 if wkswork2==3
# 	    replace weeks = 43.5 if wkswork2==4
# 	    replace weeks = 48.5 if wkswork2==5
# 	    replace weeks = 51 if wkswork2==6
# 	  }
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
# 		if _rc!=0 {
# 			replace hours = 7.5 if hrswork2==1
# 			replace hours = 22 if hrswork2==2
# 			replace hours = 32 if hrswork2==3
# 			replace hours = 37 if hrswork2==4
# 			replace hours = 40 if hrswork2==5
# 			replace hours = 44.5 if hrswork2==6
# 			replace hours = 54 if hrswork2==7
# 			replace hours = 65 if hrswork2==8
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
# 	** occ is professional
# 	gen occprofess = labforce==1 & ((occ1950 >= 0 & occ1950 < 100)|(occ1950 >= 200 & occ1950 < 300))
# 	replace occprofess = . if labforce==. | occ1950==.
# 
# 	** occ is agricultural
# 	gen occagg = labforce==1 & ((occ1950>=100&occ1950<200)|(occ1950>=800&occ1950<900))
# 	replace occagg = . if labforce==.|occ1950==.
# 	rename occ1950 occupation
# 
# 	** move from birthplace
# 	gen moved = bpl!=statefip if bpl<=56 & bpl!=999 & bpl!=.
# 	replace moved = 1 if bpl>56 & bpl!=999 & bpl!=.
# 
# 	** grandparents
# 	gen grandma = momloc!=0 if momloc!=.
# 	gen grandpa = poploc!=0 if poploc!=.
# }


# /* ----- Final Steps ----- */
#   /* Needed Variables */
#   ** This will error out if the variable already exists. Create missing variable
# * if it is not missing
# foreach v in `needvars' {
# 	cap gen `v' = .
# }
# 
# 
# /* Samples and Labels */
#   ** indicator for samples
# cap drop sample
# gen sample = country + " " + string(year) + ": " + string(source)
# 
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

# /*
#   * Author:  Karl Schulze (kschulze@frbchi.org)
# * Created: 10/2/2017
# * Purpose: This file processes the mothers for the data.
# */
#   
#   /* Read File */
#   local f `1'
# 
# 
# /* Extract Mothers */
# use if age>=16 & age<=45 & sex==2 & nchild!=0 ///
# using data/temp_`f'_data
# save data/moms_`f', replace
# 
# 
# /* Extract Spouses */
# ** load the data
# clear
# use if sex==1 & sploc!=0 using data/temp_`f'_data
# 
# ** rename these variables
# rename workedyr spouse_workedyr
# rename labforce spouse_labforce
# rename empstat spouse_empstat
# rename education spouse_education
# rename age spouse_age
# rename black spouse_black
# rename other_race spouse_other_race
# rename perwt spouse_perwt
# rename occupied spouse_occupied
# rename wage_worker spouse_wage_worker
# rename weeks spouse_weeks
# rename hours spouse_hours
# rename wage spouse_wage
# rename incwage spouse_incwage
# rename occprofess spouse_occprofess
# rename occagg spouse_occagg
# rename occupation spouse_occupation
# rename grandma spouse_grandma
# rename grandpa spouse_grandpa
# rename occscore spouse_occscore
# rename bpl spouse_bpl
# rename moved spouse_moved
# 
# ** by renaming pernum as spouse location, merging by
# * pernum will match the spouse of the husbands
# * with the person code of the wives.
# drop pernum
# rename sploc pernum
# 
# ** do the merge
# keep sample statecode serial pernum spouse*
#   sort sample statecode serial pernum
# merge 1:1 sample statecode serial pernum using data/moms_`f'
# drop if _merge==1
# gen spouse = (_merge==3)
# drop _merge
# 
# ** save the results
# save data/moms_`f', replace
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
# ** indicator for older sister
# gen oldersister12 = 0
# gen oldersister15 = 0
# forval i = 1/6 {
# 	replace oldersister12 = 1 if age`i'>=12 & sex`i'==2
# 	replace oldersister15 = 1 if age`i'>=15 & sex`i'==2
# }
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
# ** merge GDP
# gen country_orig = country
# replace country = lower(subinstr(country," ","",.))
# gen year_true = year
# replace year = 2010 if year >= 2010
# replace year = 1820 if country=="denmark" & year == 1787
# replace year = 1820 if country=="denmark" & year == 1801
# replace year = 1820 if country=="norway" & year == 1801
# replace year = 1800 if country=="germany" & year==1819
# replace year = 2008 if country=="puertorico" & year==2010
# replace year = 2008 if country=="panama" & year==2010
# cap drop gdp gdp7
# merge m:1 country year using data/maddison, keep(1 3) nogen
# drop year country
# rename year_true year
# rename country_orig country
# 
# ** save
# qui compress
# save exports/data_`f'_mom, replace
# 
# ** remove extras
# rm data/moms_`f'.dta
# 

# /*==============================================================================
#   Author:  Karl Schulze (kschulze@princeton.edu)
# Date:	 	 5/29/2019
# Purpose: This code analyzes each sample used in the fertility-labor projct
# seperately. It reads in the same options as analyze_bin.do, except
# for options specific to binning (sample, specification, instrument).
# The final product is a stata dataset that contains all estimates and
# summary statistics.
# ==============================================================================*/
#   
#   ** Housekeeping
# clear all
# cap log close
# set more off
# set matsize 1000
# global newdate = substr(c(current_date),-4,.) + "_" + ///
#   substr(c(current_date),-8,3) + "_" + trim(substr(c(current_date),1,2))
# cap mkdir "log/$newdate"
# log using "log/${newdate}/analyze_sing_samp_`1'_reg_`2'_est_`3'.log", ///
#   text replace
# 
# ** read directories
# do read_dir.do
# 
# * Load necessary programs for these specifications
# qui do afv_progs.do
# 
# /* ----------------------------- SAMPLE SELECTION ----------------------------*/
#   // See analyze_bin.do for descriptions.
# local samp_sel = `1'
# 
# /* -------------------------- SPECIFICATION SELECTION ------------------------*/
# local reg_sel = `2'
# local reg_sel1 = substr("`reg_sel'",1,1)
# local reg_sel2 = substr("`reg_sel'",2,2)
# 
# ** if in first set (<10), change
# if (strlen("`reg_sel'")==1) {
#   local reg_sel2 = `reg_sel1'
# 	local reg_sel1 = 0
# }
# 
# ** if in special, set primary and secondary to 0
# if (strlen("`reg_sel'")==3) {
# 	local reg_sel2 = 0
# 	local reg_sel1 = 0
# }
# 
# /* ---------------------------- ESTIMATOR SELECTION --------------------------*/
# local est_sel = `3'
# 
# 
# /*=========================== PREPARE RECENTERING ============================*/
# ** Here we get the complier means for the 1980-US sample. We need to do this here
# 	* since sample_select.do is called after the data to analyze has been loaded
# 	* and we only want to do this once.
# if (`reg_sel2'==7|`reg_sel2'==8) {
# 
# 	* Error out if we have incompatible options specified
# 	if `est_sel'==1 | `est_sel'==5 {
# 		di "Error: do not use these instruments with AFV Centering"
# 		exit
# 	}
# 
# 	* load the 1980 data
# 	use if sample=="USA 1980: 0" using ///
# 	  "${dirData}/data_mom_main_usa${trimmed_sample}", clear
# 
# 	* Do the same sample preparation we do with the other data
# 	qui do sample_select `samp_sel' `reg_sel' `est_sel'
# 
# 	* Store the complier means for this sample
# 	di "====== Get United States 1980 Means ======"
# 	local means80
# 	if (`reg_sel2'==7) {
# 		* complier means
# 		cmean $covars [aw=hhwt] , z(twins2 ) d(fertility)
# 		foreach xv in $covars {
# 			local rm = r(`xv')
# 			local means80 `means80' `rm'
# 		}
# 	}
# 	if (`reg_sel2'==8) {
# 		* population means
# 		foreach xv in $covars {
# 			sum `xv' [aw=hhwt]
# 			local rm = r(mean)
# 			local means80 `means80' `rm'
# 		}
# 	}
# }
# 
# 
# /*========================== ITERATE OVER SAMPLES ============================*/
# ** get list of samples from summary file
# use data/samples_to_keep, clear
# qui count
# local sN = r(N)
# expand 7 in `sN'
# replace sample = "USA 1880: 0" if _n==_N
# replace sample = "USA 1900: 0" if _n==_N-1
# replace sample = "USA 1910: 0" if _n==_N-2
# replace sample = "USA 1920: 0" if _n==_N-3
# replace sample = "USA 1930: 0" if _n==_N-4
# replace sample = "USA 1940: 0" if _n==_N-5
# levelsof sample, local(sample_list)
# 
# ** create matrix of results
# local NS : word count `sample_list'
# mat sample_mat = J(`NS',90,.)
# 
# ** iterate over each file and then iterate over each sample in the data set
# forval row = 1/`NS' {
# 
# 	/* ---------------- Preparation ---------------- */
# 	** get sample we are looking at
# 	local samp : word `row' of `sample_list'
# 
# 	** identify file for main samples and random sample U.S.
# 	if ("`samp'"=="USA 1880: 0"|"`samp'"=="USA 1910: 0"|"`samp'"=="USA 1920: 0"| ///
# 	    "`samp'"=="USA 1930: 0"|"`samp'"=="USA 1940: 0"|"`samp'"=="USA 1900: 0") {
# 		local dtsub alt_usa
# 	}
# 	else if ("`samp'"=="USA 1860: 0"|"`samp'"=="USA 1870: 0"|"`samp'"=="USA 1880: 1"| ///
# 	  "`samp'"=="USA 1900: 1"|"`samp'"=="USA 1910: 1"|"`samp'"=="USA 1920: 1"| ///
# 	  "`samp'"=="USA 1930: 1"|"`samp'"=="USA 1940: 1"|"`samp'"=="USA 1950: 0"| ///
# 	  "`samp'"=="USA 1960: 0"|"`samp'"=="USA 1970: 0"|"`samp'"=="USA 1980: 0"| ///
# 	  "`samp'"=="USA 1990: 0"|"`samp'"=="USA 2000: 0"|"`samp'"=="USA 2010: 0") {
# 		local dtsub main_usa
# 	}
# 	else {
# 		local dtsub main_int
# 	}
# 
# 	** Load the data set
# 	clear
# 	use if sample=="`samp'" using "${dirData}/data_mom_`dtsub'${trimmed_sample}"
# 	if (`reg_sel'==102) {
# 	  assert "${trimmed_sample}"=="" // cannot do this analysis in the trimmed replication files	
# 	  save "${dirTemp}/temp_sing_sing_samp`samp_sel'_reg`reg_sel'_est`est_sel'", replace
# 	  use if sample=="`samp'" using ///
# 	    "${dirData}/data_sing_`dtsub'"
# 	  append using "${dirTemp}/temp_sing_sing_samp`samp_sel'_reg`reg_sel'_est`est_sel'.dta"
# 	  rm "${dirTemp}/temp_sing_sing_samp`samp_sel'_reg`reg_sel'_est`est_sel'.dta"
# 	}
# assert sample=="`samp'"
# 
# * print what we are doing
# di "====== [`row'] Analyzing: `samp' ======"
# 
# * get number of observations in total sample
# if (_N==0) {
#   di "SKIPPING"
#   continue
# }
# 
# ** prep sample fixed effects
# encode sample, gen(s)
# 
# ** Sample Selection and Covariates
# qui do sample_select `samp_sel' `reg_sel' `est_sel'
# 
# ** This is a weird restriction. Basically, we drop samples where the outcome
# * is constant within samples. The idea behind this (and the next restriction)
# * is that we don't want to include samples in the bin where we could not
# 		* estimate alone. This really has no discernable effect on the results and
# 		* only quantitativley affects estimates where we chop the data within
# 		* country-years or change the outcome variable
# 	if (`reg_sel1'!=6 & `reg_sel'!=117) {
# 
# 		* drop if labor force participation is constant within a sample
# 		di "Outcome: $yvar"
# 		local NS = _N
# 		qui count if labforce==1
# 		local NL = r(N)
# 		if (`NS'==`NL' | `NL'==0) {
# 			di "SKIPPING"
# 			continue
# 		}	
# 
# 		* if we are not using labor force (e.g. employed), we want to be sure
# 			* there are enough of all types within labor force
# 		if ("$yvar"!="labforce") {
# 			qui sum $yvar
# 			local NY = r(sum)
# 			if ((`reg_sel'>=111&`reg_sel'<=116)| ///
# 			    (`reg_sel'>=121 & `reg_sel'<=126)| ///
# 			    (`reg_sel'>=131 & `reg_sel'<=136)| ///
# 			    (`reg_sel'>=161 & `reg_sel'<=166)) {
# 				if (`NS'==`NY' | `NL'==`NY' | `NY'==0) {
# 					di "SKIPPING"  // Occupational heterogeneity
# 					continue
# 				}	  	
# 			}
# 			else {
# 				if (`NS'==`NY' | `NY'==0) {
# 					di "SKIPPING" // _NL==_NY okay for some (e.g. employed)
# 					continue
# 				}	  				
# 			}
# 		}
# 	}
# 
# 	** drop samples with no twins (or ss children but this rarely happens)
# 		* unless looking at the extensive margin
# 	if (`reg_sel'!=102 & (`reg_sel'<161|`reg_sel'>166)) {
# 		qui count if twins2==1
# 		local NTW = r(N)
# 		qui count if same_sex==1
# 		local NSS = r(N)
# 		di "`NTW'"
# 		if ("$ivvar"!="" & "$ivvar"!="twins2" & "$ivvar"!="same_sex") {
# 			local NIV=sum($ivvar)
# 			if (`NIV'==0) {
# 				di "SKIPPING"
# 				continue
# 			}	
# 		}
# 		if (`NTW'==0|`NSS'==0) {
# 			di "SKIPPING"
# 			continue
# 		}		
# 	}
# 	if (_N==0) {
# 		di "SKIPPING"
# 		continue
# 	}
# 
# 	* Get the year of the sample
# 	qui sum year if sample == "`samp'"
# 	mat sample_mat[`row',1] = round(r(mean),0)
# 	
# 	* Get N in sample
# 	mat sample_mat[`row',2] = _N
# 
# 	** store id in matrix
# 	mat sample_mat[`row',70] = `row'
# 
# 
# 	/* --------------- First-Stage --------------- */
# 	if (`est_sel'!=1 & `reg_sel2'!=7 & `reg_sel2'!=8) {
# 
# 		* first-stage regression
# 		reg $treat $ivvar $ivvar2 $ivvar3 $covars $sexvar if sample == "`samp'" [aw=hhwt], robust
#     local N_FS = e(N)
# 
# 		* if we are using two or three instruments
# 		if ("$ivvar2"!="" & "$ivvar3"!="") {
# 			mat sample_mat[`row',46] = _b[ $ivvar ]
# 			mat sample_mat[`row',47] = _se[ $ivvar ]
# 			mat sample_mat[`row',48] = _b[ $ivvar ] + invnormal(0.975) * _se[ $ivvar ]
# 			mat sample_mat[`row',49] = _b[ $ivvar ] - invnormal(0.975) * _se[ $ivvar ]
# 			mat sample_mat[`row',50] = _b[ $ivvar2 ]
# 			mat sample_mat[`row',51] = _se[ $ivvar2 ]
# 			mat sample_mat[`row',52] = _b[ $ivvar2 ] + invnormal(0.975) * _se[ $ivvar2 ]
# 			mat sample_mat[`row',53] = _b[ $ivvar2 ] - invnormal(0.975) * _se[ $ivvar2 ]
# 			mat sample_mat[`row',54] = _b[ $ivvar3 ]
# 			mat sample_mat[`row',55] = _se[ $ivvar3 ]
# 			mat sample_mat[`row',56] = _b[ $ivvar3 ] + invnormal(0.975) * _se[ $ivvar3 ]
# 			mat sample_mat[`row',57] = _b[ $ivvar3 ] - invnormal(0.975) * _se[ $ivvar3 ]
# 
# 			* combined coefficient
# 			qui lincom $ivvar + $ivvar2 + $ivvar3
# 			mat sample_mat[`row',4] = r(estimate)
# 			mat sample_mat[`row',5] = r(se)
# 			mat sample_mat[`row',6] = r(estimate) + invnormal(0.975) * r(se)
# 			mat sample_mat[`row',7] = r(estimate) - invnormal(0.975) * r(se)
# 
# 			* F-statistic
# 			test ($ivvar=0) ($ivvar2=0) ($ivvar3=0)
# 			mat sample_mat[`row',8] = r(F)
# 		}
# 		else if ("$ivvar2"!="") {
# 
# 			* store other coefficients seperately
# 			mat sample_mat[`row',46] = _b[ $ivvar ]
# 			mat sample_mat[`row',47] = _se[ $ivvar ]
# 			mat sample_mat[`row',48] = _b[ $ivvar ] + invnormal(0.975) * _se[ $ivvar ]
# 			mat sample_mat[`row',49] = _b[ $ivvar ] - invnormal(0.975) * _se[ $ivvar ]
# 			mat sample_mat[`row',50] = _b[ $ivvar2 ]
# 			mat sample_mat[`row',51] = _se[ $ivvar2 ]
# 			mat sample_mat[`row',52] = _b[ $ivvar2 ] + invnormal(0.975) * _se[ $ivvar2 ]
# 			mat sample_mat[`row',53] = _b[ $ivvar2 ] - invnormal(0.975) * _se[ $ivvar2 ]
# 
# 			* combined coefficient
# 			qui lincom $ivvar + $ivvar2 = 0
# 			mat sample_mat[`row',4] = r(estimate)
# 			mat sample_mat[`row',5] = r(se)
# 			mat sample_mat[`row',6] = r(estimate) + invnormal(0.975) * r(se)
# 			mat sample_mat[`row',7] = r(estimate) - invnormal(0.975) * r(se)
# 
# 			* F-statistic
# 			test ($ivvar=0) ($ivvar2=0)
# 			mat sample_mat[`row',8] = r(F)
# 		}
# 		else {
# 
# 			* store single coefficient
# 			mat sample_mat[`row',4] = _b[ $ivvar ]
# 			mat sample_mat[`row',5] = _se[ $ivvar ]
# 			mat sample_mat[`row',6] = _b[ $ivvar ] + invnormal(0.975) * _se[ $ivvar ]
# 			mat sample_mat[`row',7] = _b[ $ivvar ] - invnormal(0.975) * _se[ $ivvar ]
# 
# 			* F-statistic
# 			test $ivvar
# 			mat sample_mat[`row',8] = r(F)
# 		}
# 	}
# 
# 
# 	/* ---------- Unweighted Second-Stage ---------- */
# 	if (`est_sel' == 1) {
# 		* OLS
# 		reg $yvar $treat $covars $sexvar if sample == "`samp'" [aw=hhwt], robust
# 		mat sample_mat[`row',9] = _b[$treat]
# 		mat sample_mat[`row',10] = _se[$treat]
# 		mat sample_mat[`row',11] = _b[$treat] + invnormal(0.975) * _se[$treat]
# 		mat sample_mat[`row',12] = _b[$treat] - invnormal(0.975) * _se[$treat]
# 		mat sample_mat[`row',3] = e(N)
# }
# else if (`reg_sel2'!=7 & `reg_sel2'!=8 & `reg_sel2'!=9 & (`reg_sel'<121|`reg_sel'>129)) {
# 		ivregress 2sls $yvar $covars $sexvar ($treat=$ivvar $ivvar2 $ivvar3 ) if sample == "`samp'" [aw=hhwt], robust
# 		mat sample_mat[`row',9] = _b[$treat]
# 		mat sample_mat[`row',10] = _se[$treat]
# 		mat sample_mat[`row',11] = _b[$treat] + invnormal(0.975) * _se[$treat]
# 		mat sample_mat[`row',12] = _b[$treat] - invnormal(0.975) * _se[$treat]
# 		mat sample_mat[`row',3] = e(N)
#     local N_2S = e(N)
# 	}
# 
# 	** check that the number of observations is the same in the first and second
# 		* stage regression if we are running both
#   if ("`N_FS'"!="" & "`N_2S'"!="") {
#     assert round(`N_2S')==round(`N_FS')
#   }
# 
# 
# 	/* ---------- Centered Regressions ---------- */
# 	** here is a separate command written for the re-weighted/centered regression
# 	** AFV: regressions are centered on means for groups. Selection 7 center
# 		* on the means of compliers, 8 centers on population means in the 1980 U.S.
# 		* Those means are calculated at the beginning of this file.
# 	if (`reg_sel2'==7|`reg_sel2'==8) {
# 		assert `est_sel'!=1 & "$treat"=="fertility"
#                                  rlate $yvar fertility $ivvar $covars $sexvar [aw=hhwt] ///
#                                    if sample=="`samp'" , xmeans(`means80')
# 		mat sample_mat[`row',3] = r(N)
# 		mat sample_mat[`row',4] = r(fs_b)
# 		mat sample_mat[`row',5] = r(fs_se)
# 		mat sample_mat[`row',6] = r(fs_b) + invnormal(0.975) * r(fs_se)
# 		mat sample_mat[`row',7] = r(fs_b) - invnormal(0.975) * r(fs_se)
# 		mat sample_mat[`row',8] = r(fs_F)
# 		mat sample_mat[`row',9] = r(b)
# 		mat sample_mat[`row',10] = r(se)
# 		mat sample_mat[`row',11] = r(b) + invnormal(0.975) * r(se)
# 		mat sample_mat[`row',12] = r(b) - invnormal(0.975) * r(se)
# 	}
# 
# 
# 	/* ---------- Rescaled LATE ---------- */
# 	if (`reg_sel2'==9|(`reg_sel'>=120 & `reg_sel'<=129)) {
# 		assert `est_sel'!=1 & "$treat"=="fertility"
# 	                   rescale_late $yvar fertility $ivvar $covars $sexvar [aw=hhwt], sample(s)
# 	                   mat sample_mat[`i',3] = r(N)
# 		mat sample_mat[`i',9] = r(b_r)
# 	                   mat sample_mat[`i',10] = r(se_r)
# 		mat sample_mat[`i',11] = r(ul_r)
# 	                   mat sample_mat[`i',12] = r(ll_r)
# 		mat sample_mat[`i',58] = r(b_late)
# 	                   mat sample_mat[`i',59] = r(se_late)
# 		mat sample_mat[`i',60] = r(ul_late)
# 	                   mat sample_mat[`i',61] = r(ll_late)
# 		mat sample_mat[`i',62] = r(b_m)
# 	                   mat sample_mat[`i',63] = r(se_m)
# 		mat sample_mat[`i',64] = r(ul_m)
# 	                   mat sample_mat[`i',65] = r(ll_m)
# 	}
# 
# 
# 	/* ---------- Summary Statistics ---------- */
# 	** sample means of whole population
# 	gen child2 = nchild==2
# 	gen child3 = nchild==3
# 	gen child4 = nchild==4
# 	gen child5 = nchild>=5
# 	local c = 13
# 	foreach var in $yvar labforce nchild child2 child3 child4 child5 fertility ///
# 	spouse age1 age2 sex1 sex2 same_sex two_boys two_girls twins2 age age_first_birth gdp7 {
# 		sum `var' [aw=hhwt] if sample=="`samp'"
# 	                                   mat sample_mat[`row',`c'] = r(mean)
# 		local c = `c' + 1
# 	                                                  }
# 
# 
# /* ---------- Complier Summary Statistics ---------- */
#   ** complier means within population
# if (inlist(`est_sel',0,2,3,4,8,9,10,14)==1) {
# 		qui cmean $yvar sex1 sex2 age age_first_birth [aw=hhwt] if sample=="`samp'", z($ivvar ) d(fertility) ex($yvar )
# 		foreach v in $yvar sex1 sex2 age age_first_birth {
# 			local rm = r(`v')
# 			mat sample_mat[`row',`c'] = `rm'
# 			local c = `c' + 1
# 		}
# 	}
# 
# 	** control means within population
# 	if (inlist(`est_sel',0,2,3,4,8,9,10,14)==1) {
# 		cap drop z0
# 		gen z0 = $ivvar == 0
# 		qui cmean $yvar sex1 sex2 age age_first_birth [aw=hhwt] if sample=="`samp'", ///
# 			z($ivvar ) d(fertility) ex($yvar ) treat(0)
# 		foreach v in $yvar sex1 sex2 age age_first_birth {
# 			local rm = r(`v')
#     mat sample_mat[`row',`c'] = `rm'
#                    local c = `c' + 1
# 		}
# 	}
# 
# 	** shares of complier, always taker, never taker
# 	if (inlist(`est_sel',0,2,3,4,8,9,10,14)==1) {
# 	  qui sum $treat if $ivvar==0
# 	  local shA = r(mean)
# 	  qui sum $treat if $ivvar==1
# 	  local shAC = r(mean)
# 	  mat sample_mat[`row',66] = `shAC' - `shA'
# 	                 mat sample_mat[`row',67] = `shA'
# 		mat sample_mat[`row',68] = 1 - `shAC'
# 	}
# 
# 
# 	/*---------------- Dummies for Variable Availability -----------------*/
# 	mat sample_mat[`row',43] = 0
# 	  mat sample_mat[`row',44] = 0
# 
# 	* education
# 	cap count if education!=. & sample=="`samp'"
# 	if (_rc==0) {
# 		if (r(N)>0) {
# 			mat sample_mat[`row',43] = 1
# 		}
# 	}
# 
# 	* month of birth
# 	cap count if bm1!=. & sample=="`samp'"
# 	if (_rc==0) {
# 		if (r(N)>0) {
# 			mat sample_mat[`row',44] = 1
# 	}
# }
# 
# 
# /* ---------- Edits and Output ---------- */
#   ** Assign the source of the sample
# qui sum source
# local source = round(r(mean))
# mat sample_mat[`row',45] = `source'
# }
# 
# 
# 
# /*============================= OUTPUT RESULTS ===============================*/
# clear
# svmat sample_mat
# rename sample_mat* resmat*
# 
# ** Extract sample titles from the stored local macros
# gen sample = ""
# forval i = 1/`NS' {
#   local samp : word `i' of `sample_list'
# 	replace sample = "`samp'" if resmat70==`i'
# }
# 
# ** Get Countries from Sample Names
# gen country = substr(sample,1,strpos(sample,": ")-6)
# replace country = lower(subinstr(country," ","",.))
# rename resmat1 year
# rename resmat2 N_sample
# rename resmat3 N
# rename resmat4 FS_b
# rename resmat5 FS_se
# rename resmat6 FS_up
# rename resmat7 FS_dw
# rename resmat8 F
# rename resmat9 iv_b
# rename resmat10 iv_se
# rename resmat11 iv_up
# rename resmat12 iv_dw
# rename resmat13 m_yvar
# rename resmat14 m_labforce
# rename resmat15 m_nchild
# rename resmat16 m_child2
# rename resmat17 m_child3
# rename resmat18 m_child4
# rename resmat19 m_child5
# rename resmat20 m_fertility
# rename resmat21 m_spouse
# rename resmat22 m_age1
# rename resmat23 m_age2
# rename resmat24 m_sex1
# rename resmat25 m_sex2
# rename resmat26 m_same_sex
# rename resmat27 m_two_boys
# rename resmat28 m_two_girls
# rename resmat29 m_twins2
# rename resmat30 m_age
# rename resmat31 m_age_first_birth
# rename resmat32 gdp
# rename resmat33 cm_yvar
# rename resmat34 cm_sex1
# rename resmat35 cm_sex2
# rename resmat36 cm_age
# rename resmat37 cm_age_first_birth
# rename resmat38 cl_yvar
# rename resmat39 cl_sex1
# rename resmat40 cl_sex2
# rename resmat41 cl_age
# rename resmat42 cl_age_first_birth
# rename resmat43 d_educ
# rename resmat44 d_mob
# rename resmat45 source
# rename resmat46 fs1_b
# rename resmat47 fs1_se
# rename resmat48 fs1_up
# rename resmat49 fs1_dw
# rename resmat50 fs2_b
# rename resmat51 fs2_se
# rename resmat52 fs2_up
# rename resmat53 fs2_dw
# rename resmat54 fs3_b
# rename resmat55 fs3_se
# rename resmat56 fs3_up
# rename resmat57 fs3_dw
# rename resmat58 iv2_b
# rename resmat59 iv2_se
# rename resmat60 iv2_up
# rename resmat61 iv2_dw
# ** complier control mean we rescaled by
# if (`reg_sel2'==9) {
# 	rename resmat62 ccm_yvar
# 	rename resmat63 ccm_yvar_se
# 	rename resmat64 ccm_yvar_up
# 	rename resmat65 ccm_yvar_dw
# }
# rename resmat66 share_complier
# rename resmat67 share_always
# rename resmat68 share_never
# rename resmat70 id
# 
# ** Final Cleaning of Titles
# order sample country year source gdp
# sort gdp year
# 
# ** labelling for source varaible
# label define src_lab 0 "IPUMS-US" 1 "US 100%" 2 "IPUMS-PR" 3 "IPUMS-I: Euorpe" ///
# 4 "IPUMS-I: Asia" 5 "IPUMS-I: Africa" 6 "IPUMS-I: Americas" 7 "NAPP" 8 "DHS"
# label values source src_lab
# 
# ** save the data
# save "${dirEst}/res_samp_`samp_sel'_reg_`reg_sel'_est_`est_sel'", replace
# cap log close
# 
