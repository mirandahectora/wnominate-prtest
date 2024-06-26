\name{summary.nomObject}
\alias{summary.nomObject}
\title{ W-NOMINATE Summary }
\description{ 
  \code{summary.nomObject} reads a W-NOMINATE object and prints a summary.
}

\usage{
\method{summary}{nomObject}(object,verbose=FALSE, \dots)              
}

\arguments{
  \item{object}{ a \code{wnominate} output object. }
  \item{verbose}{ logical, includes all ideal points if TRUE, otherwise only returns
  the first 10 legislators. }
  \item{...}{ other arguments do nothing and are not passed to any functions. }
}

\value{
    A summary of a \code{wnominate} object. Correct classification, APRE, and GMP
    are reported separately for each dimension.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lojames@usc.edu}

Royce Carroll \email{rcarroll@rice.edu}
}


\seealso{
  '\link{wnominate}', '\link{plot.coords}', '\link{plot.scree}', '\link{plot.angles}',
  '\link{plot.cutlines}', '\link{plot.nomObject}'
}

\examples{
    #This data file is the same as that obtained using:
    #data(sen90)
    #sen90wnom<-wnominate(sen90,polarity=c(2,5))
    data(sen90wnom)
    summary(sen90wnom)
    plot(sen90wnom)
}
\keyword{ multivariate }
