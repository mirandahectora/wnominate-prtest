\name{qnprob}
\alias{qnprob}
\title{ Qudaratic Normal Probability Matrix Generator}
\description{ 
  \code{qnprob} takes estimates from the Quadratic Normal model and returns a
  matrix of yeah choice probabilities.  It is used to generate a test
  \code{rollcall} object using \code{generateTestData}. The function is set
  up to take identical arguments to \code{nomprob}, which explains why many
  of the arguments do not do anything.

}

\usage{
qnprob(yea,nay,ideal,Beta,dimweight,normal=1)                      
}

\arguments{

For items below, m is the number of roll calls, n the number of legislators,
and d the number of dimensions.

\item{yea}{ m x d matrix of yeah locations. }
\item{nay}{ m x d matrix of no locations. }
\item{ideal}{ n x d matrix of legislator ideal points. }
\item{Beta}{ Ignored.}
\item{dimweight}{ Ignored.}
\item{normal}{ integer, '1' generates data using normal probabilities,
                    any other value generates data using logistic probabilities. }

}

\value{

An n x m matrix of probabilities giving the probability of yea for each of n
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
    qnprob(yp,np,ideal,15,0.5)    #a matrix of yea probabilities

}
\keyword{ multivariate }
