\name{swap}
\alias{swap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Switching rows and columns
}
\description{
This function is used to exchange two columns and rows of a matrix. 
}
\usage{
swap(mat, i, j)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
matrix - The matrix you want to modify
}
  \item{i}{
integer - The first column and row to be exchanged. 
}
  \item{j}{
integer - The column and row to switch with the ith column and ith row. 
}
}

\value{
The matrix mat where the columns i and j are switched and the rows i and j are switched. 
}

\examples{

p<-10
W1=1*upper.tri(matrix(0,p,p))
swap(W1,9,5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
