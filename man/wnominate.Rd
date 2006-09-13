\name{wnominate}
\alias{wnominate}
\title{ W-NOMINATE Roll Call Scaling }
\description{ 
  \code{wnominate} is the function that takes a \code{rollcall} object and estimates Poole
  and Rosenthal W-NOMINATE scores with them.
}

\usage{
wnominate(rcObject, ubeta=15, uweights=0.5, dims=2, minvotes=20,
        lop=0.025,trials=3, polarity, verbose=FALSE)                        
}

\arguments{
  \item{rcObject}{ An object of class \code{rollcall}, from Simon Jackman's \code{pscl} package. }
  \item{ubeta}{ integer, beta parameter for NOMINATE. It is strongly recommended that
                you do not change the default.}
  \item{uweights}{ integer, weight parameter for NOMINATE. It is strongly recommended
                that you do not change the default.}
  \item{dims}{ integer, number of dimensions to estimate. Must be nonnegative and cannot
                exceed 10 dimensions.}
  \item{minvotes}{ minimum number of votes a legislator must vote in for them to be analyzed. }
  \item{lop}{ A proportion between 0 and 1, the cut-off used for excluding lopsided
                votes, expressed as the proportion of non-missing votes on the minority side.
                The default, \code{lop=0.025}, eliminates votes where the minority is smaller
                than 2.5%.  Inputting a valid proportion will generate a new vector that
                overwrites the \code{lopsided} attribute in the RC object inputted. }
  \item{trials}{ integer, number of bootstrap trials for standard errors. Any number set
                    below 4 here will not return any standard errors.  Setting this number
                    to be large will slow execution of W-NOMINATE considerably.}
  \item{polarity}{ a vector specifying the legislator in the data set who is conservative on
                    each dimension. For example, \code{c(3,5)} indicates legislator 3 is conservative
                    on dimension 1, and legislator 5 is conservative on dimension 2.
                    Alternatively, polarity can be specified as a string for legislator names
                    found in \code{legis.names} (ie. \code{c("Bush", "Gore")}) if every legislative name in
                    the data set is unique.  Finally, polarity can be specified as a list (ie.
                    \code{list("cd",c(4,5))}) where the first list item is a variable from the roll
                    call object's \code{legis.data}, and the second list item is a conservative
                    legislator on each dimension as specified by the first list item.
                    \code{list("cd",c(4,5))} thus specifies the legislators with congressional
                    district numbers of 4 and 5. }
  \item{verbose}{ logical, indicates whether bills and legislators to be deleted should be
                    printed while data is being checked before ideal points are estimated.}
}

\value{
    An object of class \code{nomObject}, which in this documentation is also referred to
    as a W-NOMINATE object.

  \item{legislators}{ data frame, containing all data from the old \code{nom31.dat} file about
        legislators. For a typical W-NOMINATE object run with an ORD file read using
        \code{readKH}, it will contain the following:
    \itemize{
    \item{\code{state}}{State name of legislator.}
    \item{\code{icpsrState}}{ICPSR state code of legislator.}
    \item{\code{cd}}{Congressional District number.}
    \item{\code{icpsrLegis}}{ICPSR code of legislator.}
    \item{\code{party}}{Party of legislator.}
    \item{\code{partyCode}}{ICPSR party code of legislator.}
    \item{\code{correctYea}}{Predicted Yeas and Actual Yeas.}
    \item{\code{wrongYea}}{Predicted Yeas and Actual Nays.}
    \item{\code{wrongNay}}{Predicted Nays and Actual Yeas.}
    \item{\code{correctNay}}{Predicted Nays and Actual Nays.}
    \item{\code{GMP}}{Geometric Mean Probability.}
    \item{\code{PRE}}{Proportional Reduction In Error.}
    \item{\code{coord1D}}{First dimension W-NOMINATE score, with all subsequent dimensions
                    numbered similarly.}
    \item{\code{se1D}}{Bootstrapped standard error of first dimension W-NOMINATE score,
                    with all subsequent dimensions numbered similarly.  This will be 
                    empty if trials is set below 4.}
    \item{\code{corr.1}}{Covariance between first and second dimension W-NOMINATE score,
                    with all subsequent dimensions numbered similarly.}
    }
  }

  \item{rollcalls}{ data frame, containing all data from the old \code{nom33.dat} file about
                    bills. For a typical W-NOMINATE object run with an ORD file read
                    using \code{readKH}, it will contain the following:
    \itemize{
    \item{\code{correctYea}}{Predicted Yeas and Actual Yeas.}
    \item{\code{wrongYea}}{Predicted Yeas and Actual Nays.}
    \item{\code{wrongNay}}{Predicted Nays and Actual Yeas.}
    \item{\code{correctNay}}{Predicted Nays and Actual Nays.}
    \item{\code{GMP}}{Geometric Mean Probability.}
    \item{\code{PRE}}{Proportional Reduction In Error.}
    \item{\code{spread1D}}{First dimension W-NOMINATE spread, with all subsequent dimensions
                    numbered similarly.}
    \item{\code{midpoint1D}}{First dimension W-NOMINATE midpoint, with all subsequent dimensions
                    numbered similarly.}
    }                 
  }
  \item{dimensions}{ integer, number of dimensions estimated. }
  \item{eigenvalues}{ A vector of roll call eigenvalues. }
  \item{beta}{ The beta value used in the final iteration. }
  \item{weights}{ A vector of weights used in each iteration. }
  \item{fits}{ A vector of length 3*dimensions with the classic measures of fit.
                In order, it contains the correct classifications for each dimension,
                the APREs for each dimension, and the overall GMPs for each dimension.} 
}

\author{ 
Keith Poole \email{kpoole@ucsd.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@ucla.edu}

Royce Carroll \email{rcarroll@ucsd.edu}
}

\references{
Jeffrey Lewis. \url{ http://adric.sscnet.ucla.edu/rollcall/}

Keith Poole and Howard Rosenthal. 1997. 'Congress: A Political-Economic History of
Roll Call Voting.' New York: Oxford University Press.

Keith Poole. \url{ http://voteview.ucsd.edu/ }

}

\seealso{
  '\link{generateTestData}','\link{plot.nomObject}','\link{summary.nomObject}'.
}

\examples{
    #This data file is the same as reading file using:
    #sen90 <- readKH("ftp://voteview.com/sen90kh.ord")
    #All ORD files can be found on www.voteview.com
    data(sen90)
    
    summary(sen90)
    result<-wnominate(sen90,polarity=c(2,5))
    #'result' is the same nomObject as found in 
    #data(sen90nomObject)
    summary(result)
    plot(result)
}
\keyword{ multivariate }
