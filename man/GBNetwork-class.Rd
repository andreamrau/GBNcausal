\name{GBNetwork-class}
\Rdversion{1.1}
\docType{class}
\alias{GBNetwork-class}

\title{Class \code{"GBNetwork"}}
\description{
This class models gaussian bayesian networks. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GBNetwork", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{WeightMatrix}:}{Object of class \code{"matrix"} - Matrix of direct causal effects. }
    \item{\code{resMean}:}{Object of class \code{"vector"} - Means for the distribution law of the nodes. }
    \item{\code{resSigma}:}{Object of class \code{"vector"} - Variance for the distribution law of nodes. }
  }
}
\section{Methods}{
No methods defined with class "GBNetwork" in the signature.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
showClass("GBNetwork")
}
\keyword{classes}
