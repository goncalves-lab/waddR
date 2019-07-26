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
if ("waddR" %in% dir("..")) {
  abspath <- find_root("../waddR", criterion=has_dir("DESCRIPTION"))
} else {
  abspath <- find_root(".", criterion=has_dir("DESCRIPTION"))
}
load_all(abspath)

test_check("waddR")
