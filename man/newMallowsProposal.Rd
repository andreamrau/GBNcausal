\name{newMallowsProposal}
\alias{newMallowsProposal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum likelihood estimation of GBN with new node order. 
}
\description{
This function is used by the MCMC.GBN algorithm every step : it propose a new order, order nodes of the GBN according to the new order, computes the maximum likelihood estimator based on the data.
}
\usage{
newMallowsProposal(GBN, data, alpha, alpha2, lambda = 0, listblocks = list(), str = matrix(1, length(GBN@resSigma), length(GBN@resSigma)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GBN}{
GBN - An object of class GBN. If this is a step of the MCMC.GBN function, this argument is the maximum likelihood estimation of the previous iteration. 
}
  \item{data}{
data - can be obtained by the function dataFormat. 
}
  \item{alpha}{
double - First Mallows temperature. 
}
  \item{alpha2}{
double - Second Mallows temperature. If listblocks is empty, this argument won't be used. 
}
  \item{lambda}{
logarithmic - Coefficient of the penalty Ridge. 
}
  \item{listblocks}{
list - A list of nodes of the form (c("N1","N2"),c("N3","N4","N5")) where "N1","N2","N3","N4" and "N5" are elements of rownames and colnames of firstGBN elements and data elements. 
}
  \item{str}{
matrix - To improve the efficiency of the algorithm, a structure can be add. Colnames and rownames are not needed. 
}
}

\value{
 \item{GBN }{The maximum likelihood estimation of the first GBN calculated according to the data and the new node order. }
 \item{data }{The data are ordered according to the new node order proposed by the rmallows or rmallowsBlocks function using dataOrder.}
%% ...
}

\examples{
# Data creation

seed = 1990
n = 3000
p <- 10
m<-rep(0,10)
sigma<-rep(0.1,10) 

W <- 1*upper.tri(matrix(0,p,p))

data <- dataCreate(nbData = 2*p, p = 10,KO = list(1,9), nbKO = c(p,p), W = W , m = m,sigma = sigma, seed = seed)$data

# Initial Value

W1=1*upper.tri(matrix(0,p,p)) 
m1=rep(0,p)
s1=rep(10e-4,p)
colnames(W1)=names(m1)=names(s1)=rownames(W1)=paste("N",1:p,sep="")

firstGBN = new("GBNetwork",WeightMatrix=W1,resMean=m1,resSigma=s1)  

newMallowsProposal(firstGBN, data, alpha = 0.5, alpha2 = 0.01, lambda=0) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
