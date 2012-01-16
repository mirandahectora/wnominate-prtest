\name{plot.scree}
\alias{plot.scree}
\title{ W-NOMINATE Scree Plot }
\description{ 
  \code{plot.scree} is the function that takes a W-NOMINATE object and plots a Scree plot.
  Scree plots show the dimensionality of the voting by showing the sizes of the eigenvalues.
}

\usage{
\method{plot}{scree}(x, main.title="Scree Plot", x.title="Dimension",
        y.title="Eigenvalue",\ldots)
}

\arguments{
  \item{x}{ a \code{wnominate} output object. }
  \item{main.title}{ string, Skree plot title. }
  \item{x.title}{ string, x-axis label.}
  \item{y.title}{ string, y-axis label.}
  \item{...}{ other arguments to \code{plot}.}
}

\value{
    A Scree plot, showing the first 20 eigenvalues.
}

\author{ 
Keith Poole \email{ktpoole@uga.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lo@uni-mannheim.de}

Royce Carroll \email{rcarroll@rice.edu}
}

\seealso{
  '\link{wnominate}', '\link{plot.coords}', '\link{plot.cutlines}', '\link{plot.angles}',
  '\link{plot.nomObject}'
}

\examples{
    #This data file is the same as that obtained using:
    #data(sen90)
    #sen90wnom<-wnominate(sen90,polarity=c(2,5))
    data(sen90wnom)

    summary(sen90wnom)
    plot.scree(sen90wnom)
    plot(sen90wnom)
}
\keyword{ multivariate }
