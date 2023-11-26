library(sf)

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

## Output grid
a = seq(0, 1, length.out = 100)
s_test = expand.grid(data.frame(x = a, y = a))

## Normalized time
s = df$t / max(df$t)
N_data = nrow(df)

## Time slices (?)
t = c(250, 350, 450) / 500
s_t = rep(t, each = nrow(s_test))

s = c(s,s_t)
N = N_data + length(s_t)
print(N)

## time basis 
num_basis = c(70,250,410)
std_arr = c(0.2,0.09,0.009)
#std_arr = [0.3,0.15,0.05]

## Make up knots (same res in space and time?)
mu_knots = NULL
for (i in 1:length(num_basis)) {
  mu_knots[[i]] = seq(0,1, length.out = num_basis[i])
  print(i)
}

## Temporal basis functions
phi_t = matrix(0, nrow = N, ncol = sum(num_basis))
K = 0

for (k in 1:length(num_basis)) {
  print(paste(k, num_basis[k]))
  std = std_arr[k]
  for (i in 1:num_basis[k]) {
    print(i)
    d = (abs(s - mu_knots[[k]][i]))^2
    for (j in 1:length(d)) {
      if (d[j] >= 0 & d[j] <= 1) {
        phi_t[j, i + K] = exp(-0.5 * d[j] / std^2)
      } else {
        phi_t[j, i + K] = 0
      }
    }
  }
  stop()
}
