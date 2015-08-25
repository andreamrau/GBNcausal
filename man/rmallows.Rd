\name{rmallows}
\alias{rmallows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mallows proposal distribution function
}
\description{
The rmallows function will propose a new order for the elements of the parameter ref.
}
\usage{
rmallows(ref, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{
vector - elements of ref must be of the form "N#" where # is a number.
}
  \item{beta}{
positive - Mallows temperature. This parameter influences the way the elements of ref are mixed. 
}
}
\value{
The function returns a list of 2 items :

  \item{order}{order is a vector of integers of the size of ref. It is the order proposed by the rmallowsBlocks function.}
 \item{new}{ref[order] - Vector of ref elements mixed according to the new order.}
}

\seealso{
 \code{\link{MCMC.GBN}}, \code{\link{rmallowsBlocks}}
}
\examples{
 ref <- c("N1","N2","N3","N4","N5")
 r <- rmallows(ref, 0.7)
 r$order
 r$new

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
