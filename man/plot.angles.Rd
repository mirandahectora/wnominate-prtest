\name{plot.angles}
\alias{plot.angles}
\title{ W-NOMINATE Cutting Line Angles Plot }
\description{ 
  \code{plot.angles} reads a W-NOMINATE object and plots a histogram of
  the angles of the cutlines for two dimensions. \code{plot.angles} does not
  work for one-dimensional W-NOMINATE objects.
}

\usage{
\method{plot}{angles}(x, main.title="Cutting Line Angles",
        x.title="Angle in Degrees", y.title="Count",
        dims=c(1,2), \dots)
}

\arguments{
  \item{x}{ a \code{wnominate} output object. }
  \item{main.title}{ string, coordinate plot title. }
  \item{x.title}{ string, x-axis label.}
  \item{y.title}{ string, y-axis label.}
  \item{dims}{ vector of length 2, specifying the dimensions to be plotted.}
  \item{...}{ other arguments to \code{hist}.}
}

\value{
    A cutting line angle plot.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lojames@usc.edu}

Royce Carroll \email{rcarroll@rice.edu}
}

\seealso{
  '\link{wnominate}', '\link{plot.coords}', '\link{plot.scree}', '\link{plot.cutlines}',
  '\link{plot.nomObject}'
}

\examples{
    #This data file is the same as that obtained using:
    #data(sen90)
    #sen90wnom<-wnominate(sen90,polarity=c(2,5))
    data(sen90wnom)

    summary(sen90wnom)
    plot.angles(sen90wnom)
    plot(sen90wnom)
}
\keyword{ multivariate }
