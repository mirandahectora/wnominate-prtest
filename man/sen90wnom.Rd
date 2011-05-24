\name{sen90wnom}
\alias{sen90wnom}

\title{
  90th U.S. Senate Ideal Points
}

\description{
  This dataframe contains the estimated ideal points of the 90th U.S Senate
  using \code{wnominate}.  Although it can easily be obtained from calling
  the example in \code{wnominate}, it is included here to facilitate illustration
  of the examples for the plot and summary functions.
}

\usage{
data(sen90wnom)
}

\value{
    An object of class \code{nomObject}, which in this documentation is also referred to
    as a W-NOMINATE object.

  \item{legislators}{ data frame, containing all data from the old \code{nom33.dat} file about
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
   } }                     

  \item{dimensions}{ integer, number of dimensions estimated. }
  \item{eigenvalues}{ A vector of roll call eigenvalues. }
  \item{beta}{ The beta value used in the final iteration. }
  \item{weights}{ A vector of weights used in each iteration. }
  \item{fits}{ A vector of length 3*dimensions with the classic measures of fit.
                In order, it contains the correct classifications for each dimension,
                the APREs for each dimension, and the overall GMPs for each dimension.} 
}


\source{
  Keith Poole. 2005. \emph{90th Senate Roll Call Vote Data}.
  \url{http://www.voteview.com/}.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lo@uni-mannheim.de}

Royce Carroll \email{rcarroll@rice.edu}
}

\seealso{
  '\link{wnominate}'.
}

\examples{
    #This data file is the same as that obtained using:
    #data(sen90)
    #sen90wnom<-wnominate(sen90,polarity=c(2,5))
    data(sen90wnom)
    
    summary(sen90wnom)
    plot(sen90wnom)
}

\keyword{datasets}
