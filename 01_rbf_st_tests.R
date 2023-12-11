library(sf)
library(tidyverse)

mynorm = function(x) {
  return(sqrt(sum(x^2)))
}

## ------------------------------------
## Code to make csv file
# df <- read.table("data/LOC_50000_univariate_spacetime_matern_stationary_1",
#                   sep = ',')
# 
# names(df) <- c("x", "y", "t")
# 
# y = read.table("data/Z1_50000_univariate_spacetime_matern_stationary_1")
# df$z = y$V1
# 
# write.csv(df, "./data/synthetic_50000.csv", row.names = FALSE)

## Read data
## Note that coordinates are pre-normalized
df <- read.csv("./data/synthetic_50000.csv")

## Output grid (for testing)
a = seq(0, 1, length.out = 100)
s_test = expand.grid(data.frame(x = a, y = a))

## Normalized time
s = df$t / max(df$t)
N_data = nrow(df)

## -------------------------
## Time slices (?) - these appear to be for the quantile interps
t = c(250, 350, 450) / 500
s_t = rep(t, each = nrow(s_test))

s = c(s,s_t)
N = N_data + length(s_t)
print(N)

## time basis (number of functions for three levels)
num_basis = c(70,250,410)
std_arr = c(0.2,0.09,0.009)
#std_arr = [0.3,0.15,0.05]

## Make up knots for time
mu_knots = NULL
for (i in 1:length(num_basis)) {
  mu_knots[[i]] = seq(0,1, length.out = num_basis[i])
  print(i)
}

## -------------------------
## Temporal basis functions (Gaussian)

## Begin test code
i = 1
K = 0
k = 1
# d = (abs(s - mu_knots[[k]][i]))^2
# plot(d, type = 'l')
# 
# phi_t = matrix(0, nrow = N, ncol = sum(num_basis))
# 
# std = std_arr[k]
# for (i in 1:num_basis[k]) {
#   d = (abs(s - mu_knots[[k]][i]))^2
#   for (j in 1:length(d)) {
#     if (d[j] >= 0 & d[j] <= 1) {
#       phi_t[j, i + K] = exp(-0.5 * d[j] / std^2)
#     } else {
#       phi_t[j, i + K] = 0
#     }
#   }
# }
# 
# plot(phi_t[, 1], type = 'l')
## End test code

## Full matrix starts here
phi_t = matrix(0, nrow = N, ncol = sum(num_basis))
K = 0

for (k in 1:length(num_basis)) { ## Basis level
  print(paste(k, num_basis[k]))
  std = std_arr[k]
  for (i in 1:num_basis[k]) { ## Build BF for each knot
    d = (abs(s - mu_knots[[k]][i]))^2 ## Distance from locations to knots
    for (j in 1:length(d)) {
      if (d[j] >= 0 & d[j] <= 1) {
        phi_t[j, i + K] = exp(-0.5 * d[j] / std^2) ## Gaussian weight
      } else {
        phi_t[j, i + K] = 0
      }
    }
  }
  K = K + num_basis[k]
}

## -----------------------------------------------
## Spatial basis functions (Wendland)
## Locations for basis (bases?)
s = df[, c("x", "y")]
## Add test locs
s = rbind(s, s_test, s_test, s_test)

## These are the weights for the Wendland kernels
num_basis = c(5**2,9**2,11**2)

## Make up knot locations
knots_1d = list()
for (i in 1:length(num_basis)) {
  print(num_basis[i])
  print(seq(0, 1, length.out = sqrt(num_basis[i])))
  knots_1d[[i]] = seq(0, 1, length.out = sqrt(num_basis[i]))
}

K = 0
phi = matrix(0, nrow = N, ncol = sum(num_basis))

for (k in 1:length(num_basis)) {
  print(paste(k, num_basis[k]))
  theta = 1/sqrt(num_basis[k])*2.5 ## Wendland weight
  print(theta)
  ## Grid of knots
  knots = expand.grid(data.frame(x = knots_1d[[k]],
                                 y = knots_1d[[k]]))
  # print(knots)
  
  for (i in 1:num_basis[k]) { ## Basis levels
    ## Knot distances
    d = cbind(s[,1] - knots[i, 1], s[,2] - knots[i, 2]) ## Distance
    d = apply(d, 1, mynorm)
    d = d / theta
    for (j in 1:length(d)) {
      if (d[j] >= 0 & d[j] <= 1) {
        phi[j,i + K] = (1-d[j])**6 * (35 * d[j]**2 + 18 * d[j] + 3)/3 ## Wendland eq
      } else {
        phi[j,i + K] = 0
      }
    }
  }
  K = K + num_basis[k]
}

## -----------------------------------------------
## Basis function viz
## Level 1
k = 1
knots = expand.grid(data.frame(x = knots_1d[[k]],
                               y = knots_1d[[k]]))

myloc = 100
knots$phi = phi[myloc, 1:25]
tmp <- df %>% slice(myloc)
ggplot(knots, aes(x = x, y = y)) +
  geom_raster(aes(fill = phi)) +
  geom_point(data = tmp, aes(x = x, y = y))

## Level 2
k = 2
knots = expand.grid(data.frame(x = knots_1d[[k]],
                               y = knots_1d[[k]]))
myloc = 2
knots$phi = phi[myloc, 26:106]
tmp <- df %>% slice(myloc)
ggplot(knots, aes(x = x, y = y)) +
  geom_raster(aes(fill = phi)) +
  geom_point(data = tmp, aes(x = x, y = y))

## Level 3
k = 3
knots = expand.grid(data.frame(x = knots_1d[[k]],
                               y = knots_1d[[k]]))
myloc = 2
knots$phi = phi[myloc, 107:227]
tmp <- df %>% slice(myloc)
ggplot(knots, aes(x = x, y = y)) +
  geom_raster(aes(fill = phi)) +
  geom_point(data = tmp, aes(x = x, y = y))
