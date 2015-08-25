\name{detriangularisation}
\alias{detriangularisation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ordering matrix columns and rows
}
\description{
Assume the names of the matrix are : "N3","N4","N2","N1". This function will order the matrix rows and columns so the names are in the order : "N1","N2","N3","N4".

}
\usage{
detriangularisation(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
matrix - rownames and colnames must be of the form "N#" where # is an integer. It's a square matrix.
} 
}

\value{
The initial matrix ordered in a new way. 
}

\examples{
W <- matrix(0,4,4)
W[1,2] = 1
W[2,3] = 1
W[3,4] = 1
colnames(W) = rownames(W) = c("N3","N4","N1","N2")
detriangularisation(W)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
