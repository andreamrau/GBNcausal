\name{uniformBlocks}
\alias{uniformBlocks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Uniform proposal law
}
\description{
The uniformBlocks function will propose a new order for the elements of the parameter ref with an uniform sampling.
 The list of blocks provided by the user allows to put a constraint on the new order : the elements put inside a block remain together. 
This is the proposal distribuion used in the MCMC.GBN algorithm when the data are observational only.  
}
\usage{
uniformBlocks(ref, listblocks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{
vector - elements of ref must be of the form "N#" where # is a number.
}
  \item{listblocks}{
list - A list of blocks. Listblocks is of the form (c("N1","N2"),c("N3","N4")), where "N1","N2","N3" and "N4" are elements of ref. 
 An element of a block can't be in another block.
}
}

\value{
The function returns a list of 2 items :

  \item{order}{order is a vector of integers of the size of ref. It is the order proposed by the UniformBlocks function.}
 \item{new}{ref[order] - Vector of ref elements mixed according to the new order.}
}

\seealso{
\code{\link{rmallowsBlocks}}, \code{\link{rmallows}}
}
\examples{
ref <- c("N1","N2","N3","N4","N5")
listBlocks = list(c("N1","N3"))
uniformBlocks(ref, listblocks = listBlocks)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
