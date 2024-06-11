\name{plot.coords}
\alias{plot.coords}
\title{ W-NOMINATE Coordinate Plot }
\description{ 
  \code{plot.coords} reads a W-NOMINATE object in 2 user-specified dimensions
  and plots the coordinates of each member, applying separate colors and shapes
  to each party by default.  A unit circle is included to emphasize the
  constraints on the W-NOMINATE coordinates, and options to select non-party
  attributes of legislators are included.  For a 1D W-WNOMINATE object,
  W-NOMINATE scores are plotted against their ranks.
}

\usage{
\method{plot}{coords}(x, main.title="W-NOMINATE Coordinates",
    d1.title="First Dimension", d2.title="Second Dimension", dims=c(1,2),
    plotBy="party", color=TRUE, shape=TRUE, cutline=NULL, Legend=TRUE,
    legend.x=0.8, legend.y=1, \dots)            
}

\arguments{
  \item{x}{ a \code{wnominate} output object. }
  \item{main.title}{ string, coordinate plot title. }
  \item{d1.title}{ string, x-axis label.}
  \item{d2.title}{ string, y-axis label.}
  \item{dims}{ vector of length 2, specifying the dimensions to be plotted.}
  \item{plotBy}{ string, name of a variable in \code{nomObject\$data}.
                \code{plot.coords} will plot coordinates using this variable as a selector.}
  \item{color}{ logical, marks different groups specified by \code{plotBy} using
                different colors if TRUE. }
  \item{shape}{ logical, marks different groups specified by \code{plotBy} using different
                shapes if TRUE.}
  \item{cutline}{ vector, selects roll calls by row number for which a cutting line is desired.}
  \item{Legend}{ logical, include a generic legend.}
  \item{legend.x}{ numeric, corresponds to the 'x' argument of legend().}
  \item{legend.y}{ numeric, corresponds to the 'y' argument of legend().}
  \item{...}{ other arguments to \code{symbols}.}
}

\value{
    A coordinate plot.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lojames@usc.edu}

Royce Carroll \email{rcarroll@rice.edu}
}

\seealso{
  '\link{wnominate}', '\link{plot.scree}', '\link{plot.cutlines}', '\link{plot.angles}',
  '\link{plot.nomObject}'
}

\examples{
    #This data file is the same as that obtained using:
    #data(sen90)
    #sen90wnom<-wnominate(sen90,polarity=c(2,5))
    data(sen90wnom)

    summary(sen90wnom)
    plot.coords(sen90wnom)
    plot(sen90wnom)
}
\keyword{ multivariate }
