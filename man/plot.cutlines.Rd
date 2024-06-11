\name{plot.cutlines}
\alias{plot.cutlines}
\title{ W-NOMINATE Cutline Plot }
\description{ 
  \code{plot.cutlines} reads a W-NOMINATE object and plots the cutting line of a
  specified proportion of all votes along two user-specified dimensions.  The
  default is to plot 50 cutting lines.  This is also known as a Coombs
  mesh. A unit circle is included to emphasize the constraints on the W-NOMINATE
  coordinates.  Only cutlines that are constrained to have midpoints lying in a
  unit circle are included. \code{plot.cutlines} does not work for 1D W-NOMINATE objects.
}

\usage{         
\method{plot}{cutlines}(x, main.title="Cutting Lines",
        d1.title="First Dimension", d2.title="Second Dimension",
        lines=50,dims=c(1,2),lwd=2, \dots)
}

\arguments{
  \item{x}{ a \code{wnominate} output object. }
  \item{main.title}{ string, coordinate plot title. }
  \item{d1.title}{ string, x-axis label.}
  \item{d2.title}{ string, y-axis label.}
  \item{lines}{ numeric, number of non-constrained cutlines to be plotted. If
  this number exceeds to total number of cutlines, then all cutlines are plotted.}
  \item{dims}{ numeric vector of length 2, specifying dimensions to be plotted.}
  \item{lwd}{ numeric, line width.}
  \item{...}{ other arguments to \code{symbols}.}
}

\value{
    A Coombs mesh.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lojames@usc.edu}

Royce Carroll \email{rcarroll@rice.edu}
}

\seealso{
  '\link{wnominate}', '\link{plot.coords}', '\link{plot.scree}', '\link{plot.angles}',
  '\link{plot.nomObject}'
}

\examples{
    #This data file is the same as that obtained using:
    #data(sen90)
    #sen90wnom<-wnominate(sen90,polarity=c(2,5))
    data(sen90wnom)

    summary(sen90wnom)
    plot.cutlines(sen90wnom)
    plot(sen90wnom)
}
\keyword{ multivariate }
