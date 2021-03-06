\name{Bpest-package}
\alias{Bpest-package}
\alias{Bpest}
\docType{package}
\title{
  Bpest
}

\description{
  Estimate distribution of methylation patterns from a table of counts from a
  bisulphite sequencing experiment given a non-conversion rate and sequencing
  error rate.
}

\details{
  \tabular{ll}{
    Package: \tab Bpest\cr
    Type: \tab Package\cr
    Version: \tab 0.99\cr
    Date: \tab 2014-06-12\cr
    License: \tab GPL(>=3)\cr
  }
  The main component of this package is the function \code{\link{estimatePatterns}}, 
  which reads a table of read counts of bisulphite sequencing data for a given
  amplicon and generates a table and plot of the estimated distribution over
  methylation patterns.
}

\author{
  Peijie Lin, Sylvain Foret, Conrad Burden

  Maintainer: conrad.burden@anu.edu.au
}

\references{
}

\keyword{
  amplicon,
  bisulphite sequencing,
  bisufite sequencing,
  methylation
}

\seealso{
}

\examples{
  data(patternsExample)
  estimates <- estimatePatterns(patternsExample, epsilon=0.02, eta=0.01)
  estimates
  plotPatterns(estimates[[2]])
}

% vim:ts=2:sw=2:sts=2:expandtab:
