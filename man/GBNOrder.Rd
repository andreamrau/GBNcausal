\name{GBNOrder}
\alias{GBNOrder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ordering GBN object
}
\description{
 This function is used in the MCMC algorithm to order the GBN according to the new order found by the rmallows or rmallowsBlocks functions. 
 
}
\usage{
GBNOrder(GBN, ord)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GBN}{
GBN - An object of class GBN. 
}
  \item{ord}{
vector of integers - The new order of the GBN. 
}
}

\value{
The function returns an object of type GBN : the slots Weightmatrix, ResMean and ResSigma of the new GBN are ordered according to the order given by the user.  
}

\examples{

seed = 1990
p <- 10
m<-rep(0,10)
sigma<-rep(0.1,10) 

W <- 1*upper.tri(matrix(0,p,p))

data <- dataCreate(nbData = 2*p, p = 10,KO = list(1,9), nbKO = c(p,p), W = W , m = m,sigma = sigma, seed = seed)$data

W1=1*upper.tri(matrix(0,p,p)) 
m1=rep(0,p)
s1=rep(10e-4,p)
colnames(W1)=names(m1)=names(s1)=rownames(W1)=paste("N",1:p,sep="")

GBnetwork = new("GBNetwork",WeightMatrix=W1,resMean=m1,resSigma=s1)

ord <- sample(1:p,10, replace = FALSE)
GBNOrder(GBnetwork,ord)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
