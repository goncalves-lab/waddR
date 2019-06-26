library("testthat")
library("diffexpR")


#### TEST PERMUTATIONS OF ARMADILLO DATA STRUCTURES
test_that("permutations on cpp stl datatypes",{
  a <- rnorm(100)
  n <- 100
  df <- sapply(1:n, function(e) sample(a, replace = FALSE))
  
  permutations_matrix = permutations_internal_test_export(a, n)
  expect_equal(length(df), length(permutations_matrix))
  expect_true(all(apply(permutations_matrix, 2, function(perm) {a %in% perm})))
})


#### SANITY CHECK
test_that("Bogus tests", {
  
  x <- c(1, 2, 3)
  expect_false( length(x) == 2.7 )
  expect_false( typeof(x) == "data.frame") 
  expect_true(all(x %in% sample(x)))
})
