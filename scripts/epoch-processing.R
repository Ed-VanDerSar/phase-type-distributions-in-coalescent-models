#' Calculate each epoch duration and related quantities.
#' @param epochs_time_boundaries Numeric vector of epoch boundaries. This 
#'              should start from t_1>0, since we are assuming t_0=0
#'              and t_H=infinity. We expect it to have the form 
#'              (t_1, t_2, ... ,t_{H-1})

#' Determines if the duration of given epoch is finite and calculates
#' its duration.
#' @param h epoch index (1 ≤ h ≤ H)
#' @param epochs_time_boundaries Numeric vector of epoch boundaries.
#' @return Duration of epoch h (d_h = t_h - t_{h-1})
epoch_duration <- function(h, epochs_time_boundaries) {
  H <- length(epochs_time_boundaries)+ 1
  if (h < 1 || h > H) {
    stop("Input is not in the epochs range.")
  } else if (h == H) {
    stop("Last epoch is infinite.")
  } else if (h == 1) {
    return(epochs_time_boundaries[1])
  }
  return(epochs_time_boundaries[h] - epochs_time_boundaries[h-1])
}

#' Find which epoch contains time t
#' @param t Numeric, time point of interest
#' @param epochs_time_boundaries Numeric vector of epochs boundaries
#' @return Index of epoch containing time t (l(t) = min{h: t_{h-1} ≤ t < t_h})
epoch_index <- function(t, epochs_time_boundaries) {
  H <- length(epochs_time_boundaries) + 1
  
  if (t >= epochs_time_boundaries[length(epochs_time_boundaries)]) {
    return(H)
  } else if (0 <= t && t < epochs_time_boundaries[1]) {
    return(1)
  }
  
  for (h in 1:H) {
    if (epochs_time_boundaries[h-1] <= t && t < epochs_time_boundaries[h]) {
      return(h)
    }
  }
  
  stop("Time is not within any defined epoch")
}

#' Calculate elapsed time in current epoch at time t
#' @param t Numeric, time point of interest  
#' @param epochs_time_boundaries Numeric vector of epoch boundaries
#' @return Time elapsed in current epoch (d_{l(t)} = t - t_{l(t)-1})
epoch_elapsed_time <- function(t, epochs_time_boundaries) {
  h <- epoch_index(t, epochs_time_boundaries)
  if (h == 1) {
    return(t)
  } 
  return(t - epochs_time_boundaries[h-1])
}
