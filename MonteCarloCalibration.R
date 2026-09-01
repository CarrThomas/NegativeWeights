#################################### SUMMARY ###################################

# A script used to calibrate the Monte Carlo processes

# The first stage is assumed to be
# p(x,z) = Phi(a_0 + x'a + z(d_0 + x'd))
# x is assumed to be a vector of independent normally distributed variables, so 
# that x'a and x'd will be jointly normally distributed variables. 



################################# FRONT MATTER #################################

rm(list = ls())
gc()

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)