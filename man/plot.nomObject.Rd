\name{plot.nomObject}
\alias{plot.nomObject}
\title{ W-NOMINATE Summary Plot }
\description{ 
  \code{plot.nomObject} reads a W-NOMINATE object in two user-specified dimensions
  and plots the coordinates, cutting lines, a Coombs mesh, and a Skree plot. For
  1-dimensional W-NOMINATE objects, it plots the coordinates against the ranks
  along with a Skree plot.
}

\usage{
plot.nomObject(x, dims=c(1,2), ...)              
}

\arguments{
  \item{x}{ a \code{wnominate} output object. }
  \item{dims}{ a vector of length 2, specifying the two dimensions to be plotted. }
  \item{...}{ other arguments do nothing and are not passed to any plot functions. }
}

\value{
    A summary plot of a \code{wnominate} object.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lo@uni-mannheim.de}

Royce Carroll \email{rcarroll@rice.edu}
}

\seealso{
  '\link{wnominate}', '\link{plot.coords}', '\link{plot.skree}', '\link{plot.angles}',
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
