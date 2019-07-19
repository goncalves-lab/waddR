
expect_known <- function(type=c("value","failure", "output", "hash"), function_call, ...) {
  KNOWN_VALUES_DIR = "known.values"
  if (!dir.exists(KNOWN_VALUES_DIR)){
    dir.create(KNOWN_VALUES_DIR)
  }
  kwargs = list(...)
  if (type == "value") 
    expect_known_value(function_call, unlist(kwargs))
  else if (type == "failure")
    expect_known_failure(function_call, unlist(kwargs))
  else if (type == "output")
    expect_known_output(function_call, unlist(kwargs))
  else # value must be "hash"
    expect_known_hash(function_call, unlist(kwargs))
}


skip_temporarily <- function() {
  skip("Hardskip: This has manually been placed in a test to skip it")
}
