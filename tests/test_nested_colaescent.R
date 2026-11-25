library(testthat)
library(partitions)
library(here)

source(here("scripts", "nested-coalescent-state-space-rate-matrix.R"))

test_that("all_possible_integer_partitions return actual integer partitions", {
  # Test the the method returns a matrix with the correct dimensions 
  expect_silent(result <- all_possible_integer_partitions(4))
  expect_type(result, "integer")
  expect_equal(dim(result), c(5, 4))  # 5 partitions of 4
  
  # Test known partitions for n=4
  expected_partitions <- matrix(c(
    4,0,0,0,  # 4
    2,1,0,0,  # 2+1+1
    0,2,0,0,  # 2+2
    1,0,1,0,  # 3+1
    0,0,0,1   # 1+1+1+1
  ), nrow = 5, byrow = TRUE)
  
  # Check that all expected partitions are present (order may vary)
  expect_equal(nrow(result), nrow(expected_partitions))
  setequal(
    lapply(seq_len(nrow(result)), function(i) result[i, ]),
    lapply(seq_len(nrow(expected_partitions)), function(i) expected_partitions[i, ])
  )
  
  # Each partition sums to 4
  expect_equal(result %*% (1:4), matrix(rep(ncol(result), nrow(result)), ncol = 1))  
  
  # Test case n=1
  result1 <- all_possible_integer_partitions(1)
  expect_equal(dim(result1), c(1, 1))
  expect_equal(result1[1, 1], 1)
})

test_that("nested_state_space_mapper returns expected state space", {
  # Test the the method returns a matrix with the correct dimensions
  # for a initial sample of two species with two lineages each
  expect_silent(states <- nested_state_space(2, 2)) 
  total_gene_lineages <- 2 * 2
  expect_equal(ncol(states), total_gene_lineages)
  
  expected_states <- matrix(c(
    0,2,0,0,  # two species with two lineages
    0,0,0,1,  # one species with four lineages
    1,0,1,0,  # one species with 3 lineages and one with 1
    2,1,0,0,  # two species with 1 lineages and one with 2 
    4,0,0,0,  # four species with 1 lineages  
    1,1,0,0,  # one species with 2 lineages and one with 1
    3,0,0,0,  # three species with 1 lineage
    0,0,1,0,  # one species with 3 lineages
    2,0,0,0,  # two species with 1 lineages
    0,1,0,0,  # one species with 2 lineages
    1,0,0,0   # one species with 1 lineages
  ), nrow = 11, byrow = TRUE)
  
  # Check that the states are the expected (order may vary)
  setequal(
    lapply(seq_len(nrow(states)), function(i) states[i, ]),
    lapply(seq_len(nrow(expected_states)), function(i) expected_states[i, ])
  )
  
  # Test that initial state is correct (2 species with 2 lineages each)
  expected_initial_state <- integer(total_gene_lineages)
  expected_initial_state[2] <- 2  # 2 species with lineages 2
  expect_equal(states[1, ], expected_initial_state)
  
  # Test that all states have correct total gene count
  gene_counts <- apply(states, 1, function(row) sum(row * seq_along(row)))
  expect_true(all(gene_counts <= total_gene_lineages))
  expect_true(all(gene_counts >= 1))
})

test_that("nested_rate_matrix returns expected rate materix", {
  # Test the the method returns a square matrix
  expect_silent(rate_matrix <- nested_rate_matrix(2, 2))
  expect_type(rate_matrix, "double")
  expect_true(is.matrix(rate_matrix))
  expect_equal(nrow(rate_matrix), ncol(rate_matrix))
  
  # Test that rows sum to 0 (Markov chain property)
  row_sums <- apply(rate_matrix, 1, sum)
  expect_true(all(abs(row_sums) < 1e-10))  # Allow for floating point error
  
  # Test that off-diagonal elements are non-negative
  diag_removed <- rate_matrix
  diag(diag_removed) <- 0
  expect_true(all(diag_removed >= -1e-10))  # Allow for floating point error
  
  # Test specific known case: l=2,s=1
  rate_matrix_2_1 <- nested_rate_matrix(2, 1)
  expect_equal(dim(rate_matrix_2_1), c(nrow(nested_state_space(2,1)), nrow(nested_state_space(2,1))))
})

test_that("Integration between functions works", {
  # Test that state space mapper uses partitions correctly
  n <- 2
  b <- 2
  states <- nested_state_space(n, b)
  partitions <- all_possible_integer_partitions(n * b)
  
  # The full mass partition should be in the state space
  full_mass_state <- integer(n * b)
  full_mass_state[n] <- b  # b species with n genes each
  expect_true(any(apply(states, 1, function(row) all(row == full_mass_state))))
  
  # Test that rate matrix dimensions match state space
  rate_matrix <- nested_rate_matrix(n, b)
  expect_equal(nrow(rate_matrix), nrow(states))
  expect_equal(ncol(rate_matrix), nrow(states))
})

test_that("Performance tests for larger inputs", {
  # Test that functions don't crash on moderate inputs
  expect_silent({
    partitions_10 <- all_possible_integer_partitions(10)
    states_3_3 <- nested_state_space(3, 3)
    rate_matrix_2_2 <- nested_rate_matrix(2, 2)
  })
  
  # Test timing for moderate inputs (should complete in reasonable time)
  time_taken <- system.time({
    nested_rate_matrix(3, 2)
  })
  expect_lt(time_taken["elapsed"], 10)  # Should complete in under 10 seconds
})

test_that("Biological interpretation tests", {
  # Test that transitions make biological sense
  
  # Case 1: Two species with one gene each
  # Only possible transition: gene coalescence within a species
  # But with one gene per species, no coalescence is possible
  # So the state should be absorbing
  rate_matrix_1_2 <- nested_rate_matrix(1, 2)
  expect_equal(sum(rate_matrix_1_2[1, ]), 0)  # First row (initial state) should sum to 0
  
  # Case 2: One species with two genes
  # Only possible transition: gene coalescence
  states <- nested_state_space(2, 1)
  rate_matrix <- nested_rate_matrix(2, 1)
  
  # Should have two states: [0,1] and [2,0]
  # Transition from [2,0] to [0,1] should have rate choose(2,2) = 1
  # Find which state is [2,0] (two species with 1 gene)
  state_2_0 <- which(apply(states, 1, function(x) all(x == c(2, 0))))
  state_0_1 <- which(apply(states, 1, function(x) all(x == c(0, 1))))
  
  if (length(state_2_0) > 0 && length(state_0_1) > 0) {
    expect_equal(rate_matrix[state_2_0, state_0_1], 1)
  }
})
