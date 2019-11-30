normalize <- function(x) {
  x / sqrt(sum(x ^ 2))
}

norm_vec <- function(x)
  sqrt(sum(x ^ 2))

sample_direction <- function(dimension) {
  coords <- rnorm(dimension)#runif(dimension)
  while (all(coords == 0)) {
    coords <- runif(dimension)
  }
  normalize(coords)
}

get_subsample <- function(data, subsample) {
  n <- subsample * dim(data)[[1]]
  #sample(data, n, replace = FALSE)
  data[sample(nrow(data), n, replace = FALSE),]
}

project_df <- function(vector, points) {
  #rowSums(points * vector)
  rowSums(t(t(as.matrix(points)) * vector))
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
    split_points <- vector(mode = "numeric", length = n_halfspace)
    means <- matrix(nrow = n_halfspace, ncol = 2)
    
    for (iteration in seq_len(n_halfspace)) {
      direction <- sample_direction(dimension = (dim))
      sample_points <- get_subsample(data, subsample)
      projections <- project_df(direction, sample_points)
      split_point <- get_split_point(projections, scope)
      mean_left <-
        sum(projections < split_point) / dim(sample_points)[[1]]
      mean_right <-
        sum(projections >= split_point) / dim(sample_points)[[1]]
      
      directions[iteration,] <- direction
      split_points[iteration] <- split_point
      means[iteration,] <- c(mean_left, mean_right)
    }
    list(
      "data" = data,
      "directions" = directions,
      "split_points" = split_points,
      "means" = means,
      "number" = n_halfspace
    )
  }

update_mass <-
  function(masses,
           projections,
           split_point,
           mean_smaller,
           mean_greater) {
    projection_smaller_split <- projections < split_point
    masses[projection_smaller_split] <-
      masses[projection_smaller_split] + mean_smaller
    masses[!projection_smaller_split] <-
      masses[!projection_smaller_split] + mean_greater
    masses
  }

update_depth <- function(depths, projections, projections_z) {
  for (i in seq_len(length(projections))) {
    num_greater <- sum(projections_z > projections[i])
    depths[i] <- min(depths[i], num_greater)
  }
  depths
}

evaluate_depth <- function(data, halfspaces, metric = "mass") {
  if (metric == "mass") {
    measures <- replicate(dim(data)[[1]], 0)
    for (iteration in seq_len(halfspaces$number)) {
      projections <- project_df(halfspaces$directions[iteration,], data)
      measures <-
        update_mass(
          measures,
          projections,
          halfspaces$split_points[iteration],
          halfspaces$means[iteration, 1],
          halfspaces$means[iteration, 2]
        )
    }
    measures <- measures / halfspaces$number
  } else {
    measures <- replicate(dim(data)[[1]], .Machine$double.xmax)
    for (iteration in seq_len(halfspaces$number)) {
      projections <- project_df(halfspaces$directions[iteration,], data)
      projections_z <-
        project_df(halfspaces$directions[iteration,], halfspaces$data)
      measures <- update_depth(measures, projections, projections_z)
    }
  }
  measures
}