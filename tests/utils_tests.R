library(partitions)

#' Returns the size of the sate space
#' associated to the nested Kingman starting with s species with
#'l genes each.
#'
#' @param l the size of the initial gene sample.
#' @param s the size of the initial species sample.
cardinality_state_space <- function(l, s) {
  e <- 0
  for (k in 1:(l * s)) {
    e <- e + P(k)
  }
  return(e)
}

#' Checks that the rows of the
#' state space matrix are in fact admissible states for the
#' nested Kingman starting with s species with
#'l genes each.
#'
#' @param l the size of the initial gene sample.
#' @param s the size of the initial species sample.
#' @param e the state space matrix.
consistency_state_space <- function(l, s, e) {
  for (i in 1:nrow(e)) {
    state <- e[i, ]
    total_gene_mass <- sum((1:length(state)) * state)
    if (total_gene_mass > (l * s)) {
      print(paste("State indexed with", i, "is not admissible"))
      break
    }
  }
}
