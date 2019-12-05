library(checkmate)

normalize <- function(x) {
  assert_vector(x)
  x / sqrt(sum(x ^ 2))
}

#sample a direction that is not 0
sample_direction <- function(dimension) {
  assert_count(dimension)
  coords <- rnorm(dimension)#runif(dimension)
  while (all(coords == 0)) {
    coords <- rnorm(dimension)
  }
  normalize(coords)
}

#returns a subsample of the given data. The subsample parameter gives the
#fraction of the data that should be sampled (without replacement)
get_subsample <- function(data, subsample) {
  assert_numeric(subsample)
  assert(check_data_frame(data), check_matrix(data), combine = "or")
  n <- subsample * dim(data)[[1]]
  #sample(data, n, replace = FALSE)
  data[sample(nrow(data), n, replace = FALSE), ]
}

#projects a dataframe/matrix of points onto a unit vector
project_df <- function(vector, points) {
  assert_vector(vector)
  assert(check_data_frame(points), check_matrix(points), combine = "or")
  rowSums(t(t(as.matrix(points)) * vector))
}

#computes split points as defined by Chen et al.
#for a given numeric vector of projections and a numeric scope
get_split_point <- function(projections, scope) {
  assert_vector(projections)
  assert_numeric(scope)
  max_projection <- max(projections)
  min_projection <- min(projections)
  mid_projection <- (max_projection + min_projection) / 2
  difference = max_projection - min_projection
  interval <- scope * 0.5 * difference
  left_border <- mid_projection - interval
  right_border <- mid_projection + interval
  runif(1, min = left_border, max = right_border)
}

#Algorithm 1 from Chen et al.
#Trains the halfspace mass estimate for a given numer of halfspaces on the data.
#Outputs a list of input data, used directions, split points, means, and the
#number of halfspaces
train_depth <-
  function(data,
           n_halfspace,
           subsample = 1,
           scope = 1,
           seed) {
    assert(check_data_frame(data), check_matrix(data), combine = "or")
    assert_count(n_halfspace)
    assert_numeric(subsample)
    assert_numeric(scope)
    assert_count(seed)
    
    dimensionality <- dim(data)[[2]]
    directions <- matrix(nrow = n_halfspace, ncol = dimensionality)
    split_points <- vector(mode = "numeric", length = n_halfspace)
    means <- matrix(nrow = n_halfspace, ncol = 2)
    proj <- matrix(nrow = n_halfspace, ncol = subsample * dim(data)[[1]])
    
    for (iteration in seq_len(n_halfspace)) {
      direction <- sample_direction(dimension = dimensionality)
      sample_points <- get_subsample(data, subsample)
      n_subsample <- dim(sample_points)[[1]]
      projections <- project_df(direction, sample_points)
      split_point <- get_split_point(projections, scope)
      mean_left <- sum(projections < split_point) / n_subsample
      mean_right <- sum(projections >= split_point) / n_subsample
      
      directions[iteration, ] <- direction
      split_points[iteration] <- split_point
      means[iteration, ] <- c(mean_left, mean_right)
      proj[iteration, ] <- projections
    }
    list(
      "data" = data,
      "directions" = directions,
      "split_points" = split_points,
      "means" = means,
      "number" = n_halfspace,
      "projections" = proj
    )
  }

#updates the vector of masses for each query point when evaluating a new
#halfspace. Takes the current vector of masses, the projections of the query
#points, the split point and m_l and m_r from Chen et al.
update_mass <-
  function(masses,
           projections,
           split_point,
           mean_smaller,
           mean_greater) {
    assert_vector(masses)
    assert_vector(projections)
    assert_numeric(split_point)
    assert_numeric(mean_smaller)
    assert_numeric(mean_greater)
    
    projection_smaller_split <- projections < split_point
    masses[projection_smaller_split] <-
      masses[projection_smaller_split] + mean_smaller
    masses[!projection_smaller_split] <-
      masses[!projection_smaller_split] + mean_greater
    masses
  }

#updates the depths of the query points when evaluating a new halfspace
#takes the current depth values, the projected query points and the projected
#input points
update_depth <- function(depths, projections, projections_z) {
  assert_vector(depths)
  assert_vector(projections)
  assert_vector(projections_z)
  for (i in seq_len(length(projections))) {
    num_greater <- sum(projections_z > projections[i])
    depths[i] <- min(depths[i], num_greater)
  }
  depths
}

#evaluate Tukey halfspace depth or halfspace mass (Chen et al.) on query points.
#data contains the query points, halfspaces is the list obtained from
#train_depth() and the metric has to be either mass or depth.
evaluate_depth <- function(data, halfspaces, metric = c("mass", "depth")) {
  assert(check_data_frame(data), check_matrix(data), combine = "or")
  assert_list(halfspaces)
  assert_character(metric)
  metric <- match.arg(metric)
  
  if (metric == "mass") {
    measures <- replicate(dim(data)[[1]], 0)
    for (iteration in seq_len(halfspaces$number)) {
      projections <- project_df(halfspaces$directions[iteration, ], data)
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
  }
  if (metric == "depth") {
    measures <- replicate(dim(data)[[1]], .Machine$double.xmax)
    for (iteration in seq_len(halfspaces$number)) {
      projections <- project_df(halfspaces$directions[iteration, ], data)
      #projections_z <-
      #  project_df(halfspaces$directions[iteration,], halfspaces$data)
      projections_z <- halfspaces$projections[iteration, ]
      measures <- update_depth(measures, projections, projections_z)
    }
  } 
  measures
}