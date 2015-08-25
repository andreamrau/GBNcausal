\name{validGBNetworkObject}
\alias{validGBNetworkObject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Validation of a GB network
}
\description{
This function verify if the GBN in argument is one that can be used by the main functions of this package. 
}
\usage{
validGBNetworkObject(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object - If it's a GBN, the function returns TRUE. 
}
}
\details{
The function tests the dimnames of all the elements of GBN (weightmatrix, mean and variance), and its acyclicity. 
}
\value{
It returns FALSE if the object in argument is not a GBN, and TRUE if it is.
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
