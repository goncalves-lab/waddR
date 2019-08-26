library("testthat")
library("waddR")

##########################################################################
##                        CPP FUNCTIONS EXPOSED TO R                    ##
##########################################################################

#### NUMERIC VECTOR ABSOLUTE VALUE
test_that("abs_test_export", {
  skip_if_not_exported()
  v1 <- c(0,2.3,2,3,4,-3,-134)
  expect_equal(abs_test_export(v1), base::abs(v1))
})

#### Vector sum
test_that("sum_test_export", {
  skip_if_not_exported()
  v1 <- c(12,23,42.5, 22.5)
  expected <- 100
  expect_equal(sum_test_export(v1), expected)
})

#### NUMERIC VECTOR MEAN
test_that("mean_test_export", {
  skip_if_not_exported()
  set.seed(42)
  v <- rnorm(100)
  v2 <- c(-1,3.5,100)
  expect_equal(mean_test_export(v), base::mean(v))
  expect_equal(mean_test_export(v2), base::mean(v2))
  expect_equal(mean_test_export(c(23)), base::mean(c(23)))
})

#### NUMERIC VECTOR SD
test_that("sd_test_export", {
  skip_if_not_exported()
  set.seed(42)
  v <- rnorm(100)
  v2 <- c(-1,3.5,100)
  expect_equal(sd_test_export(v), sd(v))
  expect_equal(sd_test_export(v2), sd(v2))
  # special case: Input vectors of length 1 should have sd=0
  expect_equal(sd_test_export(c(1)), 0)
})

#### Vector add
#### Vector subtract
test_that("add_test_export", {
  skip_if_not_exported()
  v1 <- c(1,2,3,4)
  v2 <- c(2,2,2,2)
  expected <- c(3,4,5,6)
  expect_equal(add_test_export(v1,v2), expected)
  expect_error(add_test_export(v1,2)) 
  expect_equal(add_test_export_sv(v1,2), expected)
})

#### Vector subtract
test_that("subtract_test_export", {
  skip_if_not_exported()
  v1 <- c(1,2,3,4)
  v2 <- c(2,2,2,2)
  expected <- c(-1,0,1,2)
  expect_equal(subtract_test_export(v1,v2), expected) # with vector
})

#### Vector divide
test_that("divide_test_export", {
  skip_if_not_exported()
  v1 <- c(1,2,3,4)
  v2 <- c(2,2,2,2)
  expected <- c(0.5,1.0,1.5,2)
  expect_equal(divide_test_export_vectors(v1,v2), expected) # with vector
  expect_equal(divide_test_export_sv(v1,2), expected) # with scalar
})

#### Vector multiply
test_that("multiply_test_export", {
  skip_if_not_exported()
  v1 <- c(1,2,3,4)
  v2 <- c(2,2,2,2)
  expected <- c(2,4,6,8)
  expect_equal(multiply_test_export(v1,v2), expected) # with vector
  expect_warning(multiply_test_export(v1,2)) 
  expect_equal(multiply_test_export_sv(v1,2), expected) # with scalar
})

#### Vector Pow
test_that("pow_test_export", {
  skip_if_not_exported()
  v1 <- c(1,2,3,4)
  expected <- c(1,4,9,16)
  expect_equal(pow_test_export(v1,2), expected) # with scalar
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
  skip_if_not_exported()
  expect_equal(cumSum_test_export(c(1,2,3,4,5)), c(1,3,6,10,15))
  expect_equal(cumSum_test_export(c(1,2,3,4,5), 3), c(1,3,6))
  expect_equal(cumSum_test_export(c(1,2,3,4), 0), c(1,3,6,10))
})

#### INTERVAL TABLE
test_that("Interval_table_test_export", {
  skip_if_not_exported()
  wa <- c(1)
  wb <- c(0,1,0,0,0,0)
  ua <- (wa/sum(wa))[-1]
  ub <- (wb/sum(wb))[-1]
  cua <- cumsum(ua)
  cub <- cumsum(ub)
  
  wa2 <- c(0, 0, 0, 10)
  wb2 <- c(13)
  ua2 <- (wa2/sum(wa2))[-1]
  ub2 <- (wb2/sum(wb2))[-1]
  cua2 <- cumsum(ua2)
  cub2 <- cumsum(ub2)
  
  set.seed(123)
  wa3 <- runif(c(1:60))
  set.seed(234)
  wb3 <- runif(c(1:100))
  ua3 <- (wa2/sum(wa2))[-1]
  ub3 <- (wb2/sum(wb2))[-1]
  cua3 <- cumsum(ua3)
  cub3 <- cumsum(ub3)
  
  empty_vector <- c(1)[-1]
  
  # cpp implementation
  cppresult1 <- interval_table_test_export(cub, cua,1)
  cppresult2 <- interval_table_test_export(cua2, cub2,1)
  cppresult3 <- interval_table_test_export(cua3, cub3,1)
  cppresult4 <- interval_table_test_export(c(1,2,3), empty_vector, 1)
  cppresult5 <- interval_table_test_export(empty_vector, c(1,2,3), 1)
  
  # pure R
  temp <- cut(cub, breaks = c(-Inf, cua, Inf))
  rresult1 <- table(temp) + 1
  temp <- cut(cua2, breaks = c(-Inf, cub2, Inf))
  rresult2 <- table(temp) + 1
  temp <- cut(cua3, breaks = c(-Inf, cub3, Inf))
  rresult3 <- table(temp) + 1
  rresult4 <- table(cut(c(1,2,3), breaks=c(-Inf, empty_vector, Inf))) + 1
  rresult5 <- table(cut(empty_vector, breaks=c(-Inf, c(1,2,3), Inf))) + 1
  
  expect_true(all(rresult1 == cppresult1))
  expect_true(length(rresult1) == length(cppresult1))
  expect_true(all(rresult2 == cppresult2))
  expect_true(all(rresult3 == cppresult3))
  expect_true(all(rresult4 == cppresult4))
  expect_true(all(rresult5 == cppresult5))
  expect_type(interval_table_test_export(rnorm(1000, 10,0.3),
                                         rnorm(834, 20, 10)),
              "integer")
  
})

#### Repeat Weighted
test_that("rep_weighted_test_export", {
  skip_if_not_exported()
  expect_equal(rep_weighted_test_export(c(1,2,3), c(1,2,2)),
               rep(c(1,2,3), times=c(1,2,2)))
  expect_equal(rep_weighted_test_export(c(1,2,3,4), c(1,2,2,2)),
               c(1,2,2,3,3,4,4))
  expect_equal(rep_weighted_test_export(c(1,2), c(0,1)), c(2))
#  expect_error(rep_weighted_test_export(c(1,2), c()))
})

#### Concat NumericVectors
test_that("concat_test_export", {
  skip_if_not_exported()
  expect_equal(concat_test_export(c(1,2,3,4), c(5,6,7,8)),
               c(1,2,3,4,5,6,7,8))
  expect_equal(concat_test_export(c(1,2,3), c(1)), c(1,2,3,1))
  expect_equal(concat_test_export(c(1,2,3),c(1,2,3)), c(1,2,3,1,2,3))
#  expect_equal(concat_test_export(c(), c(1,2,3)), c(1,2,3))
})

#### correlation of vectors
test_that("cor_test_export", {
  skip_if_not_exported()
  expect_equal(cor_test_export(c(1), c(11)), 1)
  expect_equal(cor_test_export(c(1.0,2.0,3,4,5), c(2,3,4,5,6)), 1)
  expect_equal(cor_test_export(c(1,2,3,4,5), c(2,3,4,5,6)*-1), -1)
  expect_equal(cor_test_export(c(1,34,134,13,50,5,1,2),
                               c(2,20,55, 18, 55, 6, 2, 2)),
               cor(c(1,34,134,13,50,5,1,2), c(2,20,55, 18, 55, 6, 2, 2)))
})

#### empirical equidistant quantile
test_that("equidist_quantile_test_export",{
  skip_if_not_exported()
  probs = (seq(1:10)) / 10
  probs2 = (seq(1:5)) / 5
  expect_true(all(equidist_quantile_test_export(c(1:10), 10) ==
                  quantile(seq(1:10), probs=probs, type=1)))
  expect_true(all(equidist_quantile_test_export(c(1:10), 5) ==
                  quantile(c(1:10), probs=probs2, type=1)))
  expect_true(all(equidist_quantile_test_export(c(1:5), 5) ==
                  quantile(c(1:5), probs=probs2, type=1)))
})

####quantile
test_that("quantile_test_export",{
  skip_if_not_exported()
  expect_true(all(quantile_test_export(c(1:10),seq(1:10)/10) ==
                  quantile(c(1:10), probs=seq(1:10)/10, type=1)))
  expect_true(all(quantile_test_export(c(1:10),seq(1:5)/5) ==
                  quantile(c(1:10), probs=seq(1:5)/5, type=1)))
  expect_true(all(quantile_test_export(c(1:5), seq(1:10)/10) ==
                  quantile(c(1:5), probs=seq(1:10)/10, type=1)))
})

