# load asymptotic reference distribution
# contains:
#   value.integral  : a distribution
#   empcdf.ref      : an empirical cumulative distribution function
#                     based on the values in value.integral
load(system.file("data/ref_distr.dat", package="diffexpR"))


# on load
.First.lib <-function(lib,pkg){
  ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"),"Version")
  ver <- as.character(ver)
  # compiled cpp
  library.dynam("diffexpR", pkg, lib)
  
  cat("diffexpR", ver, "loaded\n")
}


# cleanup after our cpp libraries
#'@export
.onUnload <- function (libpath) {
  library.dynam.unload("diffexpR", libpath)
}
