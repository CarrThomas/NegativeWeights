# simulations where we have discrete covariates only

rm(list = ls())
gc()


############################## TWO BINARY X ##############################

# Simplest case! 

set.seed(1234)

#### Specification for Y

# True model
# Y = B00 * (1 - X1) * (1 - X2) + B01 * (1 - X1) * X2 
#     + B10 * X1 * (1 - X2) + B11 * X1 * X2
# Alternative form
# Y = B00 + (B10 - B00) * X1 + (B01 - B00) * X2 
# + (B00 - B01 - B10 - B11) * X1 * X2

B <- rnorm(4, 0, 1)
B00 <- B[1]
B01 <- B[2]
B10 <- B[3]
B11 <- B[4]

#### Specification for X

# 4 joint probabilities
# P00 = P(X1 = 0, X2 = 0)
# P10 = P(X1 = 1, X2 = 0)
# P01 = P(X1 = 0, X2 = 1)
# P11 = P(X1 = 1, X2 = 1)
# Specify such that the variables are not independent 

P00 <- 0.1
P10 <- 0.4
P01 <- 0.4
P11 <- 0.1

#### Generate data

N <- 100

# generate X
X <- matrix(0, N, 2)
U <- runif(N, 0, 1)

X[, 2] <- as.numeric(U > (P00 + P10))
X[, 1] <- as.numeric((U > P00) * (U < (P00 + P10)))
X[, 1] <- X[, 1] + as.numeric(U > P00 + P10 + P01)
