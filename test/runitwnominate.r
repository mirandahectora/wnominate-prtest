##################################
## Unit Testing of wnominate()  ##
##################################

## How to run the tests (do not uncomment in this file,
## but execute the commands at the R prompt):
## All you have to do is to adapt the directory locations.
## ------------------------------------------------

## define the test suite:
#library(RUnit)
#testsuite.cf <- defineTestSuite("wnominate", dirs="C:/Program Files/R/rw2001",testFileRegexp="runitwnominate.R")

## run test suite:
#testResult <- runTestSuite(testsuite.cf)

## print text protocol to console:
#printTextProtocol(testResult)

## print HTML version to a file:
#printHTMLProtocol(testResult, fileName="someFileName.html")

warpData<-function(rcObject) {
    rcObject$votes[5,]<-9
    rcObject$votes[6,]<-9
    rcObject$votes[6,1]<-6
    rcObject$votes[,4]<-9
    rcObject$votes[7,]<-9
    rcObject$votes[,2]<-1
    rcObject$votes[,3]<-6
    return(rcObject)
}

## test functions:
## ---------------------

.setUp <- function() {  ## called before each test case, see also .tearDown()
  library(wnominate)
  cat("\n\nBeginning testing...")
}

test.class <- function() {
    set.seed(2)
    test<-generateTestData()
    class(test)<-c("other")
    checkException(wnominate(test,polarity=c(1,7,30)), msg="\n\t\t...Not a NomObject error failed")
    cat("\n\t\t====== NOT NOMOBJECT ERROR ======\n\n\n")
}

test.largedim <- function() {
    set.seed(2)
    test<-generateTestData()
    checkException(wnominate(test,polarity=c(1,7,30),dims=12), msg="\n\t\t...Dimensions >10 error failed")
    cat("\n\t\t====== DIMS >10 ERROR ======\n\n\n")
}

test.smalldim <- function() {
    set.seed(2)
    test<-generateTestData()
    checkException(wnominate(test,polarity=c(1,7,30),dims=-3), msg="\n\t\t...Dimensions <1 error failed")
    cat("\n\t\t====== DIMS <1 ERROR ======\n\n\n")
}


test.trials <- function() {
    set.seed(2)
    test<-generateTestData()
    checkException(wnominate(test,polarity=c(1,7,30),trials=-3), msg="\n\t\t...Trials <1 error failed")
    cat("\n\t\t====== TRIALS <1 ERROR ======\n\n\n")
}


test.nomprob <- function() {
    yp <- matrix(rep(0,10),nrow=10)
    np <- matrix(rep(0.1,10),nrow=10)
    ideal <- matrix(rep(0,10),nrow=10)
    checkTrue(is.matrix(nomprob(yp,np,ideal,15,0.5)),msg="nomprob failed")
    cat("\n\t\t======NOMPROB IS WORKING======\n\n\n")
}    

test.outOfRange <- function() {
    set.seed(2)
    test<-generateTestData()
    checkException(wnominate(test,polarity=c(1,30)), msg="\n\t\t...Out of Range error failed") 
    cat("\n\t\t====== OUT OF RANGE ERROR ======\n\n\n")
}

test.manyElements <- function() {
    set.seed(2)
    test<-generateTestData()
    checkException(wnominate(test,polarity=c(1,7,30)), msg="\n\t\t...Too Many Elements error failed")
    cat("\n\t\t====== TOO MANY ELEMENTS ERROR ======\n\n\n")
}

test.badName <- function() {
    set.seed(2)
    test<-generateTestData()
    checkException(wnominate(test,polarity=c("Legislator10","Legislator1","junk")), msg="\n\t\t...Bad Name error failed")  
    cat("\n\t\t====== BAD NAME ERROR ======\n\n\n")
}

test.LOPspecification <- function() {
    set.seed(2)
    test<-generateTestData()
    checkException(wnominate(test,lop=3, polarity=c("Legislator10","Legislator1")), msg="\n\t\t...'Lop' Specification error failed")
    cat("\n\t\t====== 'LOP' SPECIFICATION ERROR ======\n\n\n")
}

test.polarityString <- function() {
    set.seed(2)
    test<-generateTestData()
    checkTrue(class(wnominate(test, polarity=c("Legislator10","Legislator1")))=="nomObject", msg="\n\t\t...Polarity Specification by String failed")
    cat("\n\t\t====== THIS WILL WORK: POLARITY SPECIFICATION BY STRING ======\n\n\n")
}

test.polarityList <- function() {
    set.seed(2)
    test<-generateTestData()
    test$legis.data$cd[10]<-4
    test$legis.data$cd[9]<-5
    checkTrue(class(wnominate(test, polarity=list("cd",c(4,5))))=="nomObject", msg="\n\t\t...Polarity Specification by List failed") 
    cat("\n\t\t====== THIS WILL WORK: POLARITY SPECIFICATION BY LIST ======\n\n\n")
}

test.missingItem <- function() {
    set.seed(2)
    test<-generateTestData()
    checkException(wnominate(test, polarity=list("Junk",c(4,5))), msg="\n\t\t...Missing Item error failed")
    cat("\n\t\t====== POLARITY LIST ERROR: MISSING ITEM ======\n\n\n")
}

test.missingLegis <- function() {
    set.seed(2)
    test<-generateTestData()
    checkException(wnominate(test, polarity=list("cd",c(444,47))), msg="\n\t\t...Missing Legislator error failed")
    cat("\n\t\t====== POLARITY LIST ERROR: MISSING LEGISLATOR 444======\n\n\n")
}

test.duplicateLegis <- function() {
    set.seed(2)
    test<-generateTestData()
    test$legis.data$cd[8]<-4
    test$legis.data$cd[6]<-4
    checkException(wnominate(test, polarity=list("cd",c(4,47))), msg="\n\t\t...Duplicate Legislator error failed")
    cat("\n\t\t====== POLARITY LIST ERROR: DUPLICATE LEGISLATOR 4 ======\n\n\n")
}

test.purgedBills <- function() {
    set.seed(2)
    test<-warpData(generateTestData())
    checkTrue(class(wnominate(test,polarity=c("Legislator10","Legislator1")))=="nomObject", msg="\n\t\t...Purged Bills/Legislators processing failed")
    cat("\n\t\t====== THIS WILL WORK WITH PURGED BILLS/LEGISLATORS ======\n\n\n")
}

test.polarityDeleted <- function() {
    set.seed(2)
    test<-warpData(generateTestData())
    checkException(wnominate(test,polarity=c(1,7)), msg="\n\t\t...Deleted Polarity error failed") 
    cat("\n\t\t====== DELETED POLARITY ERROR ======\n\n\n")
}

test.invalidValue <- function() {
    set.seed(2)
    test<-generateTestData()
    test$votes[1,1]<-13
    checkException(wnominate(test,polarity=c("Legislator10","Legislator1")), msg="\n\t\t...Non 0/1/NA error failed")
    cat("\n\t\t====== NON 0/1/NA VALUE ERROR ======\n\n\n")
}

test.duplicateName <- function() {
    set.seed(2)
    test<-generateTestData()
    row.names(test$legis.data)<-c(paste("Legislator",1:19),"Legislator1")
    checkException(wnominate(test,polarity=c("Legislator10","Legislator1")), msg="\n\t\t...Duplicate Name error failed")
    cat("\n\t\t====== DUPLICATE NAME ERROR ======\n\n\n")
}
