normalize <- function(x) {
  x / sqrt(sum(x ^ 2))
}

norm_vec <- function(x) sqrt(sum(x^2))

get_direction <- function(dimension) {
  coords <- runif(dimension)
  while (all(coords == 0)) {
    coords <- runif(dimension)
  }
  #runif draws from [0,1] but the directions should point to any quadrant
  signs <- as.numeric(runif(dimension) < 0.5)
  signs <- replace(signs, signs == 0,-1)
  normalize(coords * signs)
}

get_subsample <- function(data, subsample) {
  n <- subsample * dim(data)[[1]]
  #sample(data, n, replace = FALSE)
  data[sample(nrow(data), n, replace = FALSE), ]
}

project_df <- function(vector, points) {
  rowSums(points * vector)
}

get_split_point <- function(projections, scope) {
  max_projection <- max(projections)
  min_projection <- min(projections)
  mid_projection <- (max_projection + min_projection) / 2
  diff = max_projection - min_projection
  interval <- scope * 0.5 * diff
  left_border <- mid_projection - interval
  right_border <- mid_projection + interval
  runif(1, min = left_border, max = right_border)
}

train_depth <-
  function(data,
           n_halfspace,
           subsample = 1,
           scope = 1,
           seed) {
    dim <- dim(data)[[2]]
    directions <- matrix(nrow = n_halfspace, ncol = dim)
    split_points <- vector(mode = "numberic", length = n_halfspace)
    means <- matrix(nrow = n_halfspace, ncol = 2, dimnames = c("left_of_split", "right_of_split"))
    
    for (iteration in seq_len(n_halfspace)) {
      direction <- get_direction(dimension = (dim))# - 1))
      sample_points <- get_subsample(data, subsample)
      projections <- project_df(direction, sample_points)
      split_point <- get_split_point(projections, scope)
      mean_left <- sum(projections < split_point) / dim(sample_points)[[1]]
      mean_right <- sum(projections >= split_point) / dim(sample_points)[[1]]
      
      directions[iteration, ] <- direction
      split_points[iteration] <- split_point
      means[iteration, ] <- c(mean_left, mean_right)
    }
    list("directions" = directions, "split_points" = split_points, "means" = means, "number" = n_halfspace)
  }

evaluate_depth <- function(data, halfspaces) {
  masses <- replicate(dim(data)[[1]], 0)
  for (iteration in seq_len(halfspaces$number)) {
    projections <- project_df(halfspaces$directions[iteration, ], data)
    projection_smaller_split <- projections < halfspaces$split_points[iteration]
    masses[projection_smaller_split] <- masses[projection_smaller_split] + halfspaces$means[iteration, 1]
    masses[!projection_smaller_split] <- masses[!projection_smaller_split] + halfspaces$means[iteration, 2]
  }
  masses / halfspaces$number
}