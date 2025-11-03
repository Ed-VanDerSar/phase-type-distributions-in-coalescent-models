library("expm")
library(here)

source(here("scripts", "epoch-processing.R"))

#' Calculates TMRCA density function for the non-equilibrium model.  
#' @param epochs_time_boundaries Numeric vector of epoch boundaries. This 
#'              should start from t_1>0, since we are assuming t_0=0
#'              and t_H=infinity. We expect it to have the form 
#'              (t_1, t_2, ... ,t_{H-1})
#' @param rate_matrices List of rate matrices, the i-th matrix corresponds 
#'              to the rate matrix of the i-th epoch.             


#' Computes the transient states probability matrix S(t) at time t for a 
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
transient_states_submatrix_S <- function(t, rate_matrices, epochs_time_boundaries) {
  result <- diag(nrow(rate_matrices[[1]])-1) # Initialize with identity matrix 
  current_epoch_index <- epoch_index(t, epochs_time_boundaries)
  times_passed_current_epoch <- epoch_elapsed_time(t, epochs_time_boundaries)
  # Compute the product from h=1 to current_epoch_index
  if (current_epoch_index > 1) {
    for (h in 1:(current_epoch_index-1)) {
      this_epoch_rate_matrix <- rate_matrices[[h]][1:(ncol(rate_matrices[[h]]) - 1),
                                                  1:(ncol(rate_matrices[[h]]) - 1)]
      result <- result %*% expm::expm(this_epoch_rate_matrix * epoch_duration(h, epochs_time_boundaries))
    }
  }
  # Multiply by the final term, corresponding to the current epoch
  current_epoch_rate_matrix <- rate_matrices[[current_epoch_index]][1:(ncol(rate_matrices[[current_epoch_index]]) - 1),
                                                                    1:(ncol(rate_matrices[[current_epoch_index]]) - 1)]
  result <- result %*% expm::expm(current_epoch_rate_matrix * times_passed_current_epoch)
  return(result)
}

#' Computes the probability density function at time t for a 
#' non-equilibrium phase-type model with given epochs and 
#' rate matrices.
#' 
#' We are assuming that H rate matrix are provided and that the the 
#' vector epochs_time_boundaries has length equals H-1
#'
#' @param t the time evaluated
#' @param rate_matrices List of rate matrices.
#' @param epochs_time_boundaries Numeric vector of epoch boundaries.  
#' @param initial_probability_vector The initial probability vector corresponding 
#'                                   to the transient states. Its length is the total number 
#'                                   of states minus one (the only absorbing state).
#'                                   Usually we start form the first one, so that the vector 
#'                                   becomes (1,0,0,...,0).                              
#' @return the density function.
non_equilibrium_density <- function(t,
                                    rate_matrices,
                                    epochs_time_boundaries,
                                    initial_probability_vector) {
  S_t <- transient_states_submatrix_S(t, rate_matrices, epochs_time_boundaries)
  K <- nrow(S_t)
  # Check if dimensions are compatible
  if (K != length(initial_probability_vector)) {
    stop("The lenght of the initial probability vector is incorrect.")
  }
  ones_vector <- rep(1, K)  # Create vector of ones (1,...,1) of length K
  exit_rate_vector <- ones_vector - (S_t %*% ones_vector)
  # Calculate the quantity:  initial_probability_vector * S(t) * s(t), where s(t) is the exit rate vector.
  result <- initial_probability_vector %*% S_t %*% exit_rate_vector
  return(as.numeric(result))  
}