\name{Sfct}
\alias{Sfct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Covariance matrix
}
\description{
Computes the covariance matrix.
}
\usage{
Sfct(s, L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
vector - Variance. The function will take the diagonal matrix where the diagonals elements equals to the elements of s. 
}
  \item{L}{
matrix - Matrix of total causal effects of the GBN. 
}
}

\value{
Gives the covariance matrix between L and a diagonal matrix : the elements of this matrix are equal to the elements of s. 
}
\examples{
s <- rep(0.5,5)
L <- 1*upper.tri(matrix(0,5,5))+diag(rep(1,5))

Sfct(s,L)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
