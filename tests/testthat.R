# setting R_TESTS to empty string because of
# https://github.com/hadley/testthat/issues/144
# revert this when that issue in R is fixed.
Sys.setenv("R_TESTS" = "")

library(testthat)
library(waddR)
library(devtools)
library(rprojroot)

# Workaround for issue on build systems, where package directory can't be found or attached
# and non-exported functions throw errors.
tryCatch({

  dir.start <- "."
  if ("waddR" %in% dir(".."))
    dir.start <- "../waddR/"
  crit <- has_dir("DESCRIPTION")
  abspath <- find_root(dir.start, criterion=crit)
  #load the non-exported functions defined in waddR
  load_all(abspath)
      
}, error = function(err){

  # Workaround of finding the package dir has not worked => catch rprojroot::find_root() error
  # Define dummy functions, assign to the non-exported functions and skip them
  dummy <- function(...) {return(1)}
  abs_test_export <- dummy
  sum_test_export <- dummy
  mean_test_export <- dummy
  sd_test_export <- dummy
  subtract_test_export <- dummy
  add_test_export <- dummy
  add_test_export_sv <- dummy
  divide_test_export_vectors <- dummy
  divide_test_export_sv <- dummy
  multiply_test_export_sv <- dummy
  multiply_test_export <- dummy
  pow_test_export <- dummy
  cumSum_test_export <- dummy
  interval_table_test_export <- dummy
  rep_weighted_test_export <- dummy
  concat_test_export <- dummy
  cor_test_export <- dummy
  emp_equi_quantiles_test_export <- dummy

}, finally = {

  # run tests:
  #   A) on all functions, both exported and non-exported
  #   B) just on exported, skip all non-exported
  test_check("waddR")

})


