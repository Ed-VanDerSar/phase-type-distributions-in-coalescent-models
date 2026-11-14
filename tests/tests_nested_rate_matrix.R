library(testthat)
library(partitions)
library(here)

source(here("scripts", "nested-coalescent-state-space-rate-matrix.R"))

test_that("all_possible_integer_partitions returns correct structures", {
  # For n=1, should be a 1x1 matrix with value 1
  res_1 <- all_possible_integer_partitions(1)
  expect_equal(dim(res_1), c(1, 1))
  expect_equal(res_1[1,1], 1)
  
  # For n=3, partitions of 3 are (3), (2,1), (1,1,1). 
  # partitions::P(3) gives exactly 3 partitions.
  res_3 <- all_possible_integer_partitions(3)
  expect_equal(nrow(res_3), partitions::P(3))
  expect_equal(ncol(res_3), 3)
  
  # The row sums multiplied by indices must always equal n
  row_masses <- rowSums(res_3 * rep(1:3, each = nrow(res_3)))
  expect_true(all(row_masses == 3))
})

test_that("nested_state_space formats correctly", {
  n <- 2
  b <- 2
  max_mass <- n * b
  states <- nested_state_space(n, b)
  
  # Should return a matrix
  expect_true(is.matrix(states))
  
  # Number of columns should equal n * b
  expect_equal(ncol(states), max_mass)
  
  # States shouldn't have negative gene counts
  expect_true(all(states >= 0))
})

test_that("nested_rate_matrix adheres to phase-type and markov properties", {
  n <- 2
  b <- 2
  rate_mat <- nested_rate_matrix(n, b)
  
  # Must be a square matrix
  expect_true(nrow(rate_mat) == ncol(rate_mat))
  
  # The row sums of a continuous-time Markov chain transition matrix must equal 0
  row_sums <- rowSums(rate_mat)
  expect_true(all(abs(row_sums) < 1e-10))
  
  # Off-diagonals must be non-negative
  off_diags <- rate_mat[col(rate_mat) != row(rate_mat)]
  expect_true(all(off_diags >= 0))
  
  # Diagonals must be non-positive
  expect_true(all(diag(rate_mat) <= 0))
})
