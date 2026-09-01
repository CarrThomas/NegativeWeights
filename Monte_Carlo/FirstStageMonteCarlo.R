#################################### SUMMARY ################################### 

# This scipt performs Monte Carlo for the maximum rank correlation estimator. 

# The first stage is assumed to be
# p(x,z) = Phi(a_0 + x'a + z(d_0 + x'd))
# where Phi is an unknown function. In this formulation, the parameters d_0 and
# d control the sign of the relationship between Z and D. These are our primary
# object of interest. 

# We use the maximum rank correlation estimator. As the constant a_0 does not
# affect ranks, it is *not* estimated by this approach. In addition, the 
# remaining parameters are estimated only up to scale.

################################# FRONT MATTER #################################

rm(list = ls())
gc()

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)

################################### FUNCTIONS ##################################
