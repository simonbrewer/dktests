## Example from https://medium.com/@natsunoyuki/radial-basis-function-interpolation-e1152476758e

library(tidyverse)

## Function to draw samples from
phi <- function(x) {
  return(sin(2 * pi * x) * cos(5 * pi * x))
}

## Function to create Gaussian weights
## Gives weight as a function of distance betwee2 x1 and x2
gaussian <- function(x1, x2, l = 1) {
  return(exp(-(x1 - x2)^2 / (2 * l **2)))
}

## Function to create matrix of kernel weights
kernel_matrix <- function(X, l = 1) {
  nx = length(X)
  G = matrix(0, nrow = nx, ncol = nx)
  for (i in 1:nx) {
    for (j in 1:nx) {
      G[i, j] = gaussian(X[i], X[j], l)
    }
  }
  return(G)
}

## Basis function
## Dot product of G^-1 and the data (d)
rbf_model <- function(G, d) {
  m = solve(G) %*% d
  return(m)
}

## Prediction function 
## Args: interpolation points, data points, weights, etc
rbf_predict <- function(x, X, m, l = 1) {
  S = rep(0, length(x))
  for (i in 1:length(m)) {
    S = S + m[i] * gaussian(x, X[i], l)
  }
  return(S)
}

## Variables
x = seq(0, 1, length.out = 100) ## Dimension for prediction
N = 10 ## Number of samples
X = seq(0, N, 1) / N
d = phi(X)

## Plot original function and samples
plot(x, phi(x), type = 'l')
points(X, d, pch = 16)

# Gaussian length scale.
L = 1 / 9

## Calculate the kernel matrix for these poitns
G = kernel_matrix(X, L)

matplot(G, type = 'l', col = 'grey')

# Invert d = Gm for the model parameters.
m = rbf_model(G, d)
m
## Note that these are symmetrical (as the obs are too)

## Make prediction with m
S = rbf_predict(x, X, m, L)


## Prediction plot
plot(x, phi(x), type = 'l', lwd = 1.5)
lines(x, S, lty = 2, col = 2, lwd = 1.5)
points(X, d, pch = 16)

## Plot the non-weighted gaussians
nG <- matrix(NA, nrow = length(x), ncol = length(m))
for (i in 1:length(m)) {
  nG[ ,i] = gaussian(x, X[i], L)
}
matplot(x, nG, type = 'l', col = 'gray')

## Plot the weighted gaussians
wG <- matrix(NA, nrow = length(x), ncol = length(m))
for (i in 1:length(m)) {
  wG[ ,i] = m[i] * gaussian(x, X[i], L)
}
matplot(x, wG, type = 'l', col = 'gray')

## And just to show the effect - this is the sum of these weighted kernels 
## i.e. the prediction
plot(x, apply(wG, 1, sum), type = 'l')
