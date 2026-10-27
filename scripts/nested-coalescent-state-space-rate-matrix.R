library(partitions)

#' Given an integer n, returns all possible arithmetic partitions of n
#' in a vector representation of the Young diagram.
#'
#' @param n the size of the initial sample.
#' @return all the possible partitions
all_possible_integer_partitions <- function(n) {
  if (n == 1) {
    return(matrix(1, nrow = 1, ncol = 1))
  }
  
  # parts() returns a matrix where each column is a partition
  x <- partitions::parts(n)
  
  # Tabulate counts occurrences of integers 1 to n natively (Huge speedup)
  # We apply it to each column, then transpose to match original row structure
  r_matrix <- t(apply(x, 2, function(y) tabulate(y, nbins = n)))
  
  # Reverse row order to match your original right-to-left processing,
  # then reorder by the first column decreasingly
  r_matrix <- r_matrix[nrow(r_matrix):1, , drop = FALSE]
  r_matrix[order(r_matrix[, 1], decreasing = TRUE), , drop = FALSE]
}

#' Given integers n and b, return a matrix encoding all possible truncated states
#' of the nested Kingman model.
#'
#' @param n the size of the initial gene sample in each species.
#' @param b the size of the initial species sample.
#' @return the state space.
nested_state_space <- function(n, b) {
  max_val <- n * b
  valid_states_list <- vector("list", max_val)
  
  # Process states with total gene mass less than n*b
  for (k in seq_len(max_val - 1)) {
    truncated_states <- all_possible_integer_partitions(k)
    zero_cols <- matrix(0, nrow = nrow(truncated_states), ncol = max_val - k)
    valid_states_list[[k]] <- cbind(truncated_states, zero_cols)
  }
  
  initial_states <- all_possible_integer_partitions(max_val)
  
  # Find the index of the row with the maximum value in column n
  max_index <- which.max(initial_states[, n])
  
  # Reorder matrix
  reordered_initial_states <- rbind(
    initial_states[-max_index, , drop = FALSE],
    initial_states[max_index, , drop = FALSE]
  )
  valid_states_list[[max_val]] <- reordered_initial_states
  
  # Bind everything efficiently and reverse row order
  valid_states <- do.call(rbind, valid_states_list)
  valid_states[nrow(valid_states):1, , drop = FALSE]
}

#' Populates the rate matrix for the state space associated to the nested coalescent model.
#'
#' @param n the size of the initial gene sample in each species.
#' @param b the size of the initial species sample.
#' @return the rate matrix.
nested_rate_matrix <- function(n, b) {
  e <- nested_state_space(n, b)
  dim_e <- nrow(e)
  total_gene_sample <- ncol(e)
  
  rate <- matrix(0, ncol = dim_e, nrow = dim_e)
  gene_index <- seq_len(total_gene_sample)
  
  for (i in seq_len(dim_e)) {
    for (j in seq_len(dim_e)) {
      if (i == j) next # Skip diagonals for now
      
      # Difference between two states
      state_diff <- e[i, ] - e[j, ]
      gene_mass <- sum(state_diff * gene_index)
      sum_diff <- sum(state_diff)
      
      if (gene_mass == 0 && sum_diff == -1) {
        # SPECIES MERGING EVENT
        neg_merged_species <- ifelse(state_diff < 0, -state_diff, 0)
        # Vectorized replacement for the inner loop
        rate[j, i] <- prod(choose(e[j, ], neg_merged_species))
        
      } else if (gene_mass == -1 && sum_diff == 0) {
        # GENE MERGING EVENT
        w2 <- ifelse(state_diff < 0, 1, 0)
        # Vectorized replacement for the inner loop
        rate[j, i] <- sum(choose(gene_index * w2, 2))
      }
    }
  }
  
  # Fill the diagonal efficiently
  diag(rate) <- -rowSums(rate)
  
  rate
}