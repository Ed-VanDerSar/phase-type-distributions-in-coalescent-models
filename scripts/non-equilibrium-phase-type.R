library("expm")
library(here)

source(here("scripts", "epoch-processing.R"))

#' Calculate each epoch duration and related quantities.
#' @param epochs_time_boundaries Numeric vector of epoch boundaries. This 
#'              should start from t_1>0, since we are assuming t_0=0
#'              and t_H=infinity. We expect it to have the form 
#'              (t_1, t_2, ... ,t_{H-1})
#' @param rate_matrices List of rate matrices, the i-th matrix corresponds 
#'              to the rate matrix of the i-th epoch.             


#' Computes the core probability matrix S(t) at time t for a 
#' non-equilibrium phase-type model with given epochs and 
#' rate matrices.
#' 
#' This is the block of the probability matrix associated to the 
#' non-absorbing states.
#' 
#' We are assuming that H matrix are provided and that the the 
#' vector epochs_time_boundaries has length equals H-1
#'
#' @param t the time evaluated
#' @param rate_matrices List of rate matrices.
#' @param epochs_time_boundaries Numeric vector of epoch boundaries.                                
#' @return the density function.
compute_S <- function(t, rate_matrices, epochs_time_boundaries) {
  # Initialize with identity matrix 
  result <- diag(nrow(rate_matrices[[1]])-1)
  current_epoch_index <- epoch_index(t, epochs_time_boundaries)
  times_passed_current_epoch <- epoch_elapsed_time(t, epochs_time_boundaries)
  
  # Compute the product from h=1 to current_epoch_index
  for (h in 1:(current_epoch_index-1)) {
    this_epoch_rate_matrix = rate_matrices[[h]][1:(ncol(rate_matrices[[h]]) - 1),
                                                1:(ncol(rate_matrices[[h]]) - 1)]
    result <- result %*% expm::expm(this_epoch_rate_matrix * epoch_duration(h, epochs_time_boundaries))
  }
  
  # Multiply by the final term, corresponding to the current epoch
  result <- result %*% expm::expm(rate_matrices[[current_epoch_index]] * times_passed_current_epoch)
  
  return(result)
}