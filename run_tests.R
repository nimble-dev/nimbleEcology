
testfiles <- list.files("nimbleEcology/tests/testthat", full.names = T)

lapply(testfiles, source)
