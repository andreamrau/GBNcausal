\name{dataOrder}
\alias{dataOrder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ordering data 
}
\description{
This function is used by the MCMC.GBN function to ordered the data by the new order given by the rmallows or rmallowsBlocks functions.

All elements in data are ordered according to the new ord given by the user.  
}
\usage{
dataOrder(data, ord)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data - Data is a format obtained by the function dataFormat. It's a list of observations, knock-out and means of the observatations. 
}
  \item{ord}{
vector of integers - The new order of the data. 
}
}
\value{
The elements of data, x, int.nodes and int.means are all ordered according to the paramter ord. 
}

\examples{

seed = 1990

p <- 10
m<-rep(0,10)
sigma<-rep(0.1,10) 

W <- 1*upper.tri(matrix(0,p,p))

data <- dataCreate(nbData = 2*p, p = 10,KO = list(1,9), nbKO = c(p,p), W = W , m = m,sigma = sigma, seed = seed)$data
ord <- sample(1:10,10,replace = FALSE)

dataOrder(data,ord)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
