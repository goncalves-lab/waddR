library(testthat)
library(diffexpR)
library(devtools)


if(!"DESCRIPTION" %in% dir()){
  package_base <- "DESCRIPTION" %in% dir()
} else {
  package_base <- "."
}
load_all(package_base)

test_check("diffexpR")
