if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Package 'testthat' is required but not installed.")
}
if (!requireNamespace("GenesPackage", quietly = TRUE)) {
  stop("Package 'GenesPackage' is required but not installed.")
}

testthat <- asNamespace("testthat")

testthat::test_check("GenesPackage")
