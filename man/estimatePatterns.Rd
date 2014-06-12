\name{estimatePatterns}
\alias{estimatePatterns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimate distribution of methylation patterns
}
\description{
  Estimate distribution of methylation patterns from a table of counts from a bisulphite sequencing experiment given a non-conversion rate and a sequencing error rate
}
\usage{
  estimatePatterns(patternCounts,
                   epsilon,
                   eta=0,
                   column=1,
                   fast=TRUE,
                   steps=20000,
                   plot=TRUE,
                   yLimit1=-1,
                   yLimit2=-1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{patternCounts}{
    data frame with methylation patterns in first column and pattern counts in subsequent columns
  }
  \item{epsilon}{
    non-converson rate, a value between 0 and 1
  }
  \item{eta}{
    error rate, either a vector of numbers between 0 and 1 of length equal to the number of CpG sites or a single value between 0 and 1 for a single error rate across all sites
  }
  \item{column}{
    index of the column of `patternCounts' to process, an integer value from 1 to the number of pattern counts columns in `patternCounts'
  }
  \item{fast}{
    logical, if TRUE, fast version implemented (default)
  }
  \item{steps}{
    number of steps for the optimiser, passed to \code{\link{constrOptim}}
  }
  \item{plot}{
    logical, if TRUE, plot the estimated distribution
  }
  \item{yLimit1}{
    upper limit of y-axis in the first graph
  }
  \item{yLimit2}{
    upper limit of y-axis in the second graph
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  data(patternsExample)
  estimatePatterns(patternsExample, epsilon=0.02, eta=0.01, column=2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

% vim:ts=2:sw=2:sts=2:expandtab:
