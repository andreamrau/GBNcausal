\name{Lfct}
\alias{Lfct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Total causal effects 
}
\description{
Computes the inverse of (I-W) when W is upper triangular. 
}
\usage{
Lfct(W)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
matrix - W must be upper triangular.
}
}

\value{
L = I+W+W^2+...+W^(p-1) where p is the number of rows and columns of W. 
}

\seealso{
\code{\link{simulGBN}}
}
\examples{
W <- 1*upper.tri(matrix(0,5,5))
Lfct(W)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
