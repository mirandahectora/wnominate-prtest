\name{plot.skree}
\alias{plot.skree}
\title{ W-NOMINATE Skree Plot }
\description{ 
  \code{plot.skree} is the function that takes a W-NOMINATE object and plots a Skree plot.
  Skree plots show the dimensionality of the voting by showing the sizes of the eigenvalues.
}

\usage{
plot.skree(x, main.title="Skree Plot", x.title="Dimension",
        y.title="Eigenvalue",...)
}

\arguments{
  \item{x}{ a \code{wnominate} output object. }
  \item{main.title}{ string, Skree plot title. }
  \item{x.title}{ string, x-axis label.}
  \item{y.title}{ string, y-axis label.}
  \item{...}{ other arguments to \code{plot}.}
}

\value{
    A Skree plot, showing the first 20 eigenvalues.
}

\author{ 
Keith Poole \email{kpoole@ucsd.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{jameslo@ucla.edu}

Royce Carroll \email{rcarroll@ucsd.edu}
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
    plot.skree(sen90wnom)
    plot(sen90wnom)
}
\keyword{ multivariate }
