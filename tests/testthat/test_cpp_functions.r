library("testthat")
library("diffexpR")

##########################################################################
##                        CPP FUNCTIONS EXPOSED TO R                    ##
##########################################################################

#### NUMERIC VECTOR ABSOLUTE VALUE
test_that("abs_test_export", {
  v1 <- c(0,2.3,2,3,4,-3,-134)
  expect_equal(diffexpR::abs_test_export(v1), base::abs(v1))
})

#### Vector sum
test_that("sum_test_export", {
  v1 <- c(12,23,42.5, 22.5)
  expected <- 100
  expect_equal(diffexpR::sum_test_export(v1), expected)
})

#### NUMERIC VECTOR MEAN
test_that("mean_test_export", {
  set.seed(42)
  v <- rnorm(100)
  v2 <- c(-1,3.5,100)
  expect_equal(diffexpR::mean_test_export(v), base::mean(v))
  expect_equal(diffexpR::mean_test_export(v2), base::mean(v2))
})

#### NUMERIC VECTOR SD
test_that("sd_test_export", {
  set.seed(42)
  v <- rnorm(100)
  v2 <- c(-1,3.5,100)
  expect_equal(diffexpR::sd_test_export(v), sd(v))
  expect_equal(diffexpR::sd_test_export(v2), sd(v2))
})

#### Vector subtract
test_that("subtract_test_export", {
  v1 <- c(1,2,3,4)
  v2 <- c(2,2,2,2)
  expected <- c(-1,0,1,2)
  expect_equal(diffexpR::subtract_test_export(v1,v2), expected) # with vector
})

#### Vector divide
test_that("divide_test_export", {
  v1 <- c(1,2,3,4)
  v2 <- c(2,2,2,2)
  expected <- c(0.5,1.0,1.5,2)
  expect_equal(diffexpR::divide_test_export_vectors(v1,v2), expected) # with vector
  expect_equal(diffexpR::divide_test_export_sv(v1,2), expected) # with scalar
})

#### Vector multiply
test_that("multiply_test_export", {
  v1 <- c(1,2,3,4)
  v2 <- c(2,2,2,2)
  expected <- c(2,4,6,8)
  expect_equal(diffexpR::multiply_test_export(v1,v2), expected) # with vector
  expect_warning(diffexpR::multiply_test_export(v1,2)) 
  expect_equal(diffexpR::multiply_test_export_sv(v1,2), expected) # with scalar
})

#### Vector Pow
test_that("pow_test_export",{
  v1 <- c(1,2,3,4)
  expected <- c(1,4,9,16)
  expect_equal(diffexpR::pow_test_export(v1,2), expected) # with scalar
})

#### PERMUTATIONS OF NUMERIC VECTOR
test_that("permutations", {
  set.seed(123)
  a <- rnorm(100)
  n <- 1000
  df <- sapply(1:n, function(e) sample(a, replace = FALSE))
  
  permutations_matrix = permutations(a, n)
  expect_equal(length(df), length(permutations_matrix))
  expect_true(all( 
    apply(permutations_matrix, 2, function(perm) {a %in% perm})
  ))
})

#### CUMULATIVE SUM OF NUMERIC VECTOR
test_that("cumSum_test_export", {
  expect_equal(diffexpR::cumSum_test_export(c(1,2,3,4,5)), c(1,3,6,10,15))
  expect_equal(diffexpR::cumSum_test_export(c(1,2,3,4,5), 3), c(1,3,6))
  expect_equal(diffexpR::cumSum_test_export(c(1,2,3,4), 0), c(1,3,6,10))
})

#### INTERVAL TABLE
test_that("Interval_table_test_export", {
  # setup
  set.seed(42)
  wa<- rnorm(200, 20, 1)
  set.seed(43)
  wb <- rnorm(210, 21, 1) 
  ua <- (wa/sum(wa))[-2000]
  ub <- (wb/sum(wb))[-2010]
  cua <- c(cumsum(ua))
  cub <- c(cumsum(ub))
  
  wa2 <- c(1, 1)
  wb2 <- c(2) 
  ua2 <- (wa/sum(wa))[-2]
  ub2 <- (wb/sum(wb))[-1]
  cua2 <- c(cumsum(ua))
  cub2 <- c(cumsum(ub))
  
  set.seed(123)
  wa3 <- runif(c(1:60))
  set.seed(234)
  wb3 <- runif(c(1:100))
  ua3 <- (wa/sum(wa))[-2]
  ub3 <- (wb/sum(wb))[-1]
  cua3 <- c(cumsum(ua))
  cub3 <- c(cumsum(ub))
  
  # cpp implementation
  cppresult1 <- interval_table_test_export(cub, cua,1)
  cppresult2 <- interval_table_test_export(cua2, cub2,1)
  cppresult3 <- interval_table_test_export(cua3, cub3,1)
  
  # pure R
  temp <- cut(cub, breaks = c(-Inf, cua, Inf))
  rresult1 <- table(temp) + 1
  temp <- cut(cua2, breaks = c(-Inf, cub2, Inf))
  rresult2 <- table(temp) + 1
  temp <- cut(cua3, breaks = c(-Inf, cub3, Inf))
  rresult3 <- table(temp) + 1
  
  expect_true(all(rresult1 == cppresult1))
  expect_true(length(rresult1) == length(cppresult1))
  expect_true(all(rresult2 == cppresult2))
  expect_true(all(rresult3 == cppresult3))
  expect_type(interval_table_test_export(rnorm(1000, 10,0.3), rnorm(834, 20, 10)), "integer")
  
})

#### Repeat Weighted
test_that("rep_weighted_test_export", {
  expect_equal(rep_weighted_test_export(c(1,2,3),c(1,2,2)), rep(c(1,2,3), times=c(1,2,2)))
  expect_equal(rep_weighted_test_export(c(1,2,3,4), c(1,2,2,2)), c(1,2,2,3,3,4,4))
  expect_equal(rep_weighted_test_export(c(1,2), c(0,1)), c(2))
#  expect_error(rep_weighted_test_export(c(1,2), c()))
})

#### Concat NumericVectors
test_that("concat_test_export", {
  expect_equal(concat_test_export(c(1,2,3,4), c(5,6,7,8)), c(1,2,3,4,5,6,7,8))
  expect_equal(concat_test_export(c(1,2,3), c(1)), c(1,2,3,1))
  expect_equal(concat_test_export(c(1,2,3),c(1,2,3)), c(1,2,3,1,2,3))
#  expect_equal(concat_test_export(c(), c(1,2,3)), c(1,2,3))
})

#### correlation of vectors
test_that("cor_test_export", {
  skip("Likely to throw errors while emp_equi_quant... not fixed")
  expect_equal(cor_test_export(c(1,2,3,4,5), c(2,3,4,5,6)), 1)
  expect_equal(cor_test_export(c(1,2,3,4,5), c(2,3,4,5,6)*-1), -1)
  expect_equal(cor_test_export(c(1,34,134,13,50,5,1,2), c(2,20,55, 18, 55, 6, 2, 2)), cor(c(1,34,134,13,50,5,1,2), c(2,20,55, 18, 55, 6, 2, 2)))
})

#### empirical equidistant quantiles
test_that("emp_equi_quantiles_test_export",{
  
  expect_equal(emp_equi_quantiles_test_export(c(1:10),10), c(1:10))
  expect_equal(emp_equi_quantiles_test_export(c(1:20),10), c(2,4,6,8,10,12,14,16,18,20))
  expect_equal(emp_equi_quantiles_test_export(c(1:5),10), c(1,1,2,2,3,3,4,4,5,5))
})

