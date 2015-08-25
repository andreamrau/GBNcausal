\name{betaFromW}
\alias{betaFromW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes total causal effects
}
\description{
In the case of a gaussian bayesian network where W is the weightmatrix (direct causal effects), beta is the matrix of total causal effects. 

The function computes the inverse of (I-W)^p even if W is not upper triangular. In the case of a GBN, the graph is acyclic, it exists at least one topological order. W is upper triangular and nilpotent under these orderings. 

 The function is ordering rows and columns of W to make the matrix upper triangular, then computes I+W+W^2+...+W^(p-1) where p is the number of rows and columns of W.  

}
\usage{
betaFromW(Wgt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Wgt}{
matrix - Weightmatrix of the graph (direct causal effects). Must have dimnames of the form "N#" where # is a number. 
}
}

\value{
The function returns the matrix of total causal effects. 
}

\seealso{
\code{\link{Lfct}}
}
\examples{
p <- 10
W <- 1*upper.tri(matrix(0,10,10))
ord <- sample(1:p,10,replace = FALSE)
colnames(W) = rownames(W) = paste("N",1:p,sep = "")
W <- W[ord,ord]

betaFromW(W)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
