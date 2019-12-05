library("testthat")
source("halfspacemass-sol.R")

generate_data <- function(num_points, dimensionality) {
  data <- matrix(ncol = dimensionality, nrow = num_points)
  apply(data, 2, rnorm)
}

#create 3D data and compute its mass
num_points <- 100
data <- generate_data(num_points, 3)
trained_depth <- train_depth(data,
                             n_halfspace = 100,
                             scope = 1,
                             seed = 1234)
mass <- evaluate_depth(data, trained_depth, metric = "mass")
cols <- rev(heat.colors(num_points))
ord <- order(mass)
data <- data[ord,]
mass <- mass[ord]

#plot data and colour points by their depth
#points towards the center of the cloud should have higher mass than
#those further out
car::scatter3d(
  x = data[, 1],
  y = data[, 2],
  z = data[, 3],
  surface = FALSE,
  groups = as.factor(seq_len(nrow(data))),
  surface.col = cols
)

context("test halfspace mass computation")

#data is sampled from gaussian centered at 0, hence
#mass of the point at (0,0,0) should be higher than the average mass or mass
#at (1,1,1)

testthat::test_that("check that 3D points towards the center have more mass", {
  data <- generate_data(100, 3)
  trained_depth <- train_depth(data,
                               n_halfspace = 100,
                               scope = 1,
                               seed = 1234)
  mass_z <-
    evaluate_depth(data.frame(0, 0, 0), trained_depth, metric = "mass")
  mass_one <-
    evaluate_depth(data.frame(1, 1, 1), trained_depth, metric = "mass")
  testthat::expect_less_than(mean(mass), mass_z)
  testthat::expect_less_than(mass_one, mass_z)
})

testthat::test_that("check that 4D points towards the center have more mass", {
  data <- generate_data(100, 4)
  trained_depth <- train_depth(data,
                               n_halfspace = 100,
                               scope = 1,
                               seed = 1234)
  mass <- evaluate_depth(data, trained_depth, metric = "mass")
  mass_z <-
    evaluate_depth(data.frame(0, 0, 0, 0), trained_depth, metric = "mass")
  mass_one <-
    evaluate_depth(data.frame(1, 1, 1, 1), trained_depth, metric = "mass")
  testthat::expect_less_than(mean(mass), mass_z)
  testthat::expect_less_than(mass_one, mass_z)
})
