
# load asymptotic reference distribution
# contains:
#   value.integral  : a distribution
#   empcdf.ref      : an empirical cumulative distribution function
#                     based on the values in value.integral
load(system.file("data/ref_distr.dat", package="waddR"))


# on load
.First.lib <-function(lib,pkg){
  ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"),"Version")
  ver <- as.character(ver)
  # compiled cpp
  library.dynam("waddR", pkg, lib)
  
  cat("waddR", ver, "loaded\n")
}


# cleanup after our cpp libraries
.onUnload <- function (libpath) {
  library.dynam.unload("waddR", libpath)
}
