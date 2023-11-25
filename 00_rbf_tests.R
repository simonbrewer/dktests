library(sf)

grid <- st_read("./data/coarse_grid_pts/grid_pts_coarse.shp")

st_crs(grid) <- 4326

grid_utm <- st_transform(grid, 32612)

plot(st_geometry(grid_utm))


lon = st_coordinates(grid)[ ,1]
lat = st_coordinates(grid)[ ,2]
normalized_lon = (lon-min(lon))/(max(lon)-min(lon))
normalized_lat = (lat-min(lat))/(max(lat)-min(lat))
N = length(lon)

## Basis functions
num_basis = c(10**2,19**2,37**2)

xs = list()
knots_1dx = list()
knots_1dy = list()

for (i in 1:length(num_basis)) {
  print(num_basis[i])
  print(seq(0, 1, length.out = sqrt(num_basis[i])))
  xs[[i]] = rep(i, sqrt(num_basis[i]))
  knots_1dx[[i]] = seq(0, 1, length.out = sqrt(num_basis[i]))
  knots_1dy[[i]] = seq(0, 1, length.out = sqrt(num_basis[i]))
}

plot(knots_1dx[[1]], xs[[1]], ylim = c(0.5, 3.5))
points(knots_1dx[[2]], xs[[2]])
points(knots_1dx[[3]], xs[[3]])

knots_2d_1 = expand.grid(data.frame(x = knots_1dx[[1]],
                                    y = knots_1dy[[1]]))
knots_2d_2 = expand.grid(data.frame(x = knots_1dx[[2]],
                                    y = knots_1dy[[2]]))
knots_2d_3 = expand.grid(data.frame(x = knots_1dx[[3]],
                                    y = knots_1dy[[3]]))

plot(knots_2d_1, cex = 2)
points(knots_2d_2, cex = 1.5)
points(knots_2d_3, cex = 1)
points(normalized_lon, normalized_lat, pch = '.', col = 'blue', cex = 5)

## Make up the Wendland kernels
basis_size = 0
phi = matrix(0, nrow = N, ncol = sum(num_basis))

mynorm = function(x) {
  return(sqrt(sum(x^2)))
}

for (i in 1:length(num_basis)) {
  theta = 1/sqrt(num_basis[i])*2.5
  knots = expand.grid(data.frame(x = knots_1dx[[i]],
                                 y = knots_1dy[[i]]))
  
  for (j in 1:num_basis[i]) {
    ## Knot distances
    d = cbind(normalized_lon - knots[j, 1], normalized_lat - knots[j, 2]) / theta
    d = apply(d, 1, mynorm)
    for (k in 1:length(d)) {
      if (d[k] >= 0 & d[k] <= 1) {
        phi[k,j+basis_size] = (1-d[k])**6 * (35 * d[k]**2 + 18 * d[k] + 3)/3
      } else {
        phi[k,j+basis_size] = 0
      }
    }
  }
  basis_size = basis_size + num_basis[i]
  
}

## Find zero weights
phi_sum = apply(phi, 2, sum)
phi_reduce = phi[, phi_sum > 0]

image(matrix(phi[20, 1:100], 10, 10))
points(normalized_lon, normalized_lat)
points(normalized_lon[20], normalized_lat[20], pch = 16)

library(animation)
saveGIF({
  for (i in 1:N) {
    image(matrix(phi[i, 1:100], 10, 10))
    points(normalized_lon, normalized_lat)
    points(normalized_lon[i], normalized_lat[i], pch = 16)
  }
}, movie.name = "slc_rbf_1.gif", interval = 0.025)

saveGIF({
  for (i in 1:N) {
    image(matrix(phi[i, 1:100], 10, 10))
    points(normalized_lon, normalized_lat)
    points(normalized_lon[i], normalized_lat[i], pch = 16)
  }
}, movie.name = "slc_rbf_1.gif", interval = 0.01)

saveGIF({
  for (i in 1:N) {
    image(matrix(phi[i, 101:461], 19, 19))
    points(normalized_lon, normalized_lat)
    points(normalized_lon[i], normalized_lat[i], pch = 16)
  }
}, movie.name = "slc_rbf_2.gif", interval = 0.01)

saveGIF({
  for (i in 1:N) {
    image(matrix(phi[i, 462:1830], 37, 37))
    points(normalized_lon, normalized_lat)
    points(normalized_lon[i], normalized_lat[i], pch = 16)
  }
}, movie.name = "slc_rbf_3.gif", interval = 0.01)

