library(parallel)
library(here)

source(here("scripts/non-equilibrium-phase-type.R"))
source(here("scripts/structured-coalescent-state-space-rate-matrix.R"))
source(here("scripts/phase-type-density-and-moments.R"))

#' List of rate matrices for the two island model. Each matrix corresponds to 
#' and epoch where the population increases exponentially in the past. This corresponds
#' to H=5 epochs: Population doubles in each epoch.
island_rate_matrices_exponential_demographics <- function() {
  h_values <- c(1, 2, 4, 8, 16)
  matrices <- lapply(h_values, function(h) {
    two_islands_rate_matrix(20, 20, h, h, 1, 1)
  })
  return(matrices)
}

#The density for the model with the exponential demographics
non_equilibrium_islands_density <- function() {
  epochs_time_boundaries <- c(4,8,12,16) # We change epoch each four generations.
  rate_matrices <- island_rate_matrices_exponential_demographics()
  initial_probability_vector <- diag(1, (nrow(two_islands_state_space(20,20))-2))[1, ] #Starts in the first state. 
  function(x) {
    non_equilibrium_density(x,
                            rate_matrices,
                            epochs_time_boundaries,
                            initial_probability_vector)
  }
}

#The density for the model in equilibrium 
equilibrium_islands_density <- function() {
  rate_matrix <- two_islands_rate_matrix(20, 20, 1, 1, 1, 1)
  tmrca_density(rate_matrix)
}

plot_densities <- function(limit_interval) {
  f_1 <- non_equilibrium_islands_density()
  f_2 <- equilibrium_islands_density()
  x_vals <- seq(0, limit_interval, length.out = 250)
  
  # Compute in parallel
  y_vals1 <- unlist(mclapply(x_vals, f_1, mc.cores = detectCores() - 1))
  y_vals2 <- unlist(mclapply(x_vals, f_2, mc.cores = detectCores() - 1))
  
  plot(x_vals, y_vals1, type = "l", col = "darkred", lwd = 2,
       xlab = "time", ylab = "density", 
       main = "TMRCA probability density",
       ylim = range(c(y_vals1, y_vals2)))
  
  lines(x_vals, y_vals2, col = "darkblue", lwd = 2)
  
  legend("topright", 
         legend = c("Exponential population growth", "Demograpich equilibrium"),
         col = c("darkred", "darkblue"), 
         lwd = 2)  
}
 