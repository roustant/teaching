# -------------------
# Author: O. Roustant
# -------------------

# The aim of this code is to illustrate the main properties
# of the conditional distribution of Gaussian vectors.

library(MASS)

# draw realizations of a Gaussian vector (X1, X2) 
rho <- 0.95
Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
mu <- c(1, 2)
X <- mvrnorm(n = 10000, mu = mu, Sigma = Sigma)
plot(X)
abline(v = mu[1], h = mu[2], lty = 'dotted')

# let us condition on X1 = x1
# in practice: consider the points such that x1-h <= X1 <= x1+h (with a small h)
x1 <- 1.3
h <- 0.1

indices <- which(abs(X[,1] - x1) < h)   
X2 <- X[indices, 2]
points(X[indices, 1], X[indices, 2], col="blue")
abline(v = c(x1-h, x1+h))

# study the distribution of X2 conditioned on X1 = x1
# Is it a normal one ?
# What can you say of its mean as a function of x1 ? 
# (consider a reasonable domain for x1, such as [mu - sigma, mu + sigma],
# in which there are enough data points)
# Same question with the variance ?

hist(X2, breaks = 30)
qqnorm(X2); qqline(X2)

