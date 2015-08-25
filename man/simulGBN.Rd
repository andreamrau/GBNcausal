\name{simulGBN}
\alias{simulGBN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate data from Gaussian Bayesian Networks}

\description{
simulGBN simulates data from a GBN. 
}
\usage{
simulGBN(N, m, s, W, int = numeric(0), int_data = matrix(runif(N * length(int)), nrow = N), seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
integer - Number of observations to simulate.
}
  \item{m}{
vector - Mean of the simulation law. Must be the of the same dimension as s and columns of W.
}
  \item{s}{
vector - Variance of the simulation law. Must be the of the same dimension as m and columns of W.
}
  \item{W}{
matrix - Weightmatrix of the GBNetwork. 
}
  \item{int}{
integer or vector of integers - Columns of the simulated data mentionned in int will have the corresponding value of int_data. 
}
  \item{int_data}{
vector - Values of the columns mentionned in int. 
}
  \item{seed}{
integer - The function simulate the data according to the seed. 
}
}

\value{
The function return a matrix of the simulated data : the number of rows is N and the number of columns is equal to the dimension of W. 
}

\examples{

# Simulate data from a GBN with intervention on nodes 1 and 9. 

p<-10

W1=1*upper.tri(matrix(0,p,p))
m1=rep(0,10)
sigma1=rep(0.1,10)

X=simulGBN(2*p,m1,sigma1,W1,seed = 1990)
X=rbind(X,simulGBN(p,m1,sigma1,W1,int=1,int_data=0,seed = 1990)) 
X=rbind(X,simulGBN(p,m1,sigma1,W1,int=9,int_data=0,seed = 1990)) 


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
