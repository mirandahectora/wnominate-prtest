## define the test suite:
library(RUnit)
testsuite.cf <- defineTestSuite("wnominate", dirs="C:/Program Files/R/rw2001",testFileRegexp="runitwnominate.R")

## run test suite:
testResult <- runTestSuite(testsuite.cf)

## print text protocol to console:
printTextProtocol(testResult)
