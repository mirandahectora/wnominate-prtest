\name{generateTestData}
\alias{generateTestData}
\title{ Test Data Generator for W-NOMINATE }
\description{ 
  \code{generateTestData} is the function that generates a \code{rollcall} object used to
  test \code{wnominate}.  The description of the result below is copied from the
  documentation of package \code{pscl}, written by Simon Jackman.
}

\usage{
generateTestData(legislators=20, rcVotes=100,
            yea=matrix(runif(rcVotes,min=-0.2,max=0.7),nrow=rcVotes),
            nay=matrix(runif(rcVotes,min=-0.7,max=0.2),nrow=rcVotes),
            ideal=matrix(rnorm(legislators),nrow=legislators),
            Beta=15, dimweight=0.5,normal=1, seed = NULL, utility='nominate')

}

\arguments{
  \item{legislators}{ integer, number of Legislators ('n'). }
  \item{rcVotes}{ integer, number of roll calls ('m').}
  \item{yea}{ an m x d matrix of yea locations, where 'd' are the number of dimensions. }
  \item{nay}{ an m x d matrix of no locations, where 'd' are the number of dimensions. }
  \item{ideal}{ an n x d matrix of legislator ideal points. }
  \item{Beta}{ scalar giving beta parameter from W-NOMINATE.}
  \item{dimweight}{ d x 1 vector of dimension weights. }
  \item{normal}{ integer, '1' generates data using normal probabilities,
                    any other value generates data using logistic probabilities. }
  \item{seed}{ a single value, interpreted as an integer, used to set the seed.
		If seed is NULL, current seed is used. }
  \item{utility}{ String set to either ``nominate'' or ``qn''. `nominate' allows NOMINATE
		  logit or probit utilities, while ``qn'' allows for quadratic normal utilities
		  to be used when generating the roll call matrix.}
}

\value{
  An object of class \code{rollcall}
  
  \item{votes}{ n x m vote matrix in 0/1/NA format. }
  \item{n}{ integer, number of legislators. }
  \item{m}{ integer, number of roll call votes. }
  \item{lopsided}{ logical vector of length \code{m} indicating dropped vote. This is
                    recomputed in \code{wnominate} and is never used. }
  \item{legis.data}{ matrix, user-supplied data on legislators, containing data from
                    an ORD file.  Legislator names are rownames to this matrix.}
  \item{vote.data}{ user-supplied data on rollcall votes, set to NULL. }
  \item{desc}{ user-supplied description, set to NULL.}
}

\author{ 
Keith Poole \email{kpoole@ucsd.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@ucla.edu}

Royce Carroll \email{rcarroll@ucsd.edu}
}

\seealso{
  '\link{wnominate}', '\link{nomprob}'.
}

\examples{
    dat<-generateTestData()
    result<-wnominate(dat,polarity=c(1,2))
    summary(result)
    plot(result)
}
\keyword{ multivariate }
