library(testthat)
library(here)

source(here("scripts", "epoch-processing.R"))

# Test for epoch_duration function
test_that("epoch_duration calculates correct durations", {
  t_boundaries <- c(5, 12, 20)  # H = 4 epochs
  
  # Test normal cases
  expect_equal(epoch_duration(1, t_boundaries), 5)   # 5 - 0 = 5
  expect_equal(epoch_duration(2, t_boundaries), 7)   # 12 - 5 = 7
  expect_equal(epoch_duration(3, t_boundaries), 8)   # 20 - 12 = 8
  expect_error(epoch_duration(4, t_boundaries), "Last epoch is infinite.")
  
  # Test error cases
  expect_error(epoch_duration(0, t_boundaries), "Input is not in the epochs range")
  expect_error(epoch_duration(5, t_boundaries), "Input is not in the epochs range")
  expect_error(epoch_duration(-1, t_boundaries), "Input is not in the epochs range")
})

# Test for epoch_index function
test_that("epoch_index correctly identifies epoch containing time t", {
  t_boundaries <- c(5, 12, 20)  # H = 4 epochs
  
  # Test times within each epoch
  expect_equal(epoch_index(0, t_boundaries), 1)    # t = t_0 (start of epoch 1)
  expect_equal(epoch_index(2.5, t_boundaries), 1)  # within epoch 1
  expect_equal(epoch_index(4.999, t_boundaries), 1) # still epoch 1
  
  expect_equal(epoch_index(5, t_boundaries), 2)    # t = t_1 (start of epoch 2)
  expect_equal(epoch_index(8, t_boundaries), 2)    # within epoch 2
  expect_equal(epoch_index(11.999, t_boundaries), 2) # still epoch 2
  
  expect_equal(epoch_index(12, t_boundaries), 3)   # t = t_2 (start of epoch 3)
  expect_equal(epoch_index(15, t_boundaries), 3)   # within epoch 3
  expect_equal(epoch_index(19.999, t_boundaries), 3) # still epoch 3
  
  # Test times in the infinite final epoch
  expect_equal(epoch_index(20, t_boundaries), 4)   # t = t_3 (start of infinite epoch)
  expect_equal(epoch_index(100, t_boundaries), 4)  # far in the future
  expect_equal(epoch_index(1000, t_boundaries), 4) # very far in the future
  
  # Test edge cases with different boundary sets
  expect_equal(epoch_index(0, c(10)), 1)
  expect_equal(epoch_index(5, c(10)), 1)
  expect_equal(epoch_index(10, c(10)), 2) # H=2
})

# Test for epoch_elapsed_time function
test_that("epoch_elapsed_time calculates correct elapsed time", {
  t_boundaries <- c(5, 12, 20)  # H = 4 epochs
  
  # Test times at epoch boundaries
  expect_equal(epoch_elapsed_time(0, t_boundaries), 0)    # just started epoch 1
  expect_equal(epoch_elapsed_time(5, t_boundaries), 0)    # just started epoch 2
  expect_equal(epoch_elapsed_time(12, t_boundaries), 0)   # just started epoch 3
  expect_equal(epoch_elapsed_time(20, t_boundaries), 0)   # just started infinite epoch
  
  # Test times within epochs
  expect_equal(epoch_elapsed_time(2.5, t_boundaries), 2.5)  # 2.5 - 0 = 2.5
  expect_equal(epoch_elapsed_time(8, t_boundaries), 3)      # 8 - 5 = 3
  expect_equal(epoch_elapsed_time(15, t_boundaries), 3)     # 15 - 12 = 3
  expect_equal(epoch_elapsed_time(25, t_boundaries), 5)     # 25 - 20 = 5
  
  # Test with different boundary sets
  expect_equal(epoch_elapsed_time(7, c(10)), 7)
  expect_equal(epoch_elapsed_time(15, c(10, 20)), 5) # 15 - 10 = 5
})
