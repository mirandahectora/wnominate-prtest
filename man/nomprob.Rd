\name{nomprob}
\alias{nomprob}
\title{ NOMINATE Probability Matrix Generator}
\description{ 
  \code{nomprob} takes estimates from the W-NOMINATE model and returns a
  matrix of yeah choice probabilities.  It is used to generate a test
  \code{rollcall} object using \code{generateTestData}.

}

\usage{
nomprob(yea, nay, ideal, Beta, dimweight, normal=1)                       
}

\arguments{

For items below, m is the number of roll calls, n the number of legislators,
and d the number of dimensions.

\item{yea}{ m x d matrix of yeah locations. }
\item{nay}{ m x d matrix of no locations. }
\item{ideal}{ n x d matrix of legislator ideal points. }
\item{Beta}{ scalar giving beta parameter from W-NOMINATE. Usually set to 15.}
\item{dimweight}{ d x 1 vector of dimension weights. Usually set to 0.5.}
\item{normal}{ integer, '1' generates data using normal probabilities,
                    any other value generates data using logistic probabilities. }

}

\value{

An n x m matrix of probabilities giving the probability of yeah for each of n
legislators on each of m votes   

}

\author{ 
Keith Poole \email{kpoole@ucsd.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@ucla.edu}

Royce Carroll \email{rcarroll@ucsd.edu}
}

\seealso{
  '\link{generateTestData}' and '\link{wnominate}'.
}

\examples{

    yp <- matrix(rep(0,10),nrow=10)
    np <- matrix(rep(0.1,10),nrow=10)
    ideal <- matrix(rep(0,10),nrow=10)
    nomprob(yp,np,ideal,15,0.5)    #a matrix of yea probabilities

}
\keyword{ multivariate }
