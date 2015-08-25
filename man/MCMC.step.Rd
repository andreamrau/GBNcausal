\name{MCMC.step}
\alias{MCMC.step}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main step of the MCMC.GBN algorithm
}
\description{
This function is used by the MCMC.GBN algorithm every step : it uses the newMallowsProposal function and accepts or rejects the proposed GBN. 

}
\usage{
MCMC.step(GBN, data, alpha = 0.05, alpha2 = 0.05, lambda, listblocks = list(), str = matrix(1, length(GBN@resSigma), length(GBN@resSigma)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GBN}{
GBN - An object of class GBN. If this is a step of the MCMC.GBN function, this argument is the maximum likelihood estimation of the previous iteration. 
}
  \item{data}{
data - can be obtained by the function dataFormat or dataCreate. 
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
 \item{newGBN }{If accept equals 1, new GBN is the maximum likelihood estimation of the previous GBN based on a new order of the nodes and the data. If accept equals 0, newGBN is the GBN of the previous step.  }
 \item{accept }{If 1, the new GBN is the new estimation. If 0, the new GBN is the old GBN. }

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
firstGBN = GBNmle(firstGBN,data,lambda= 0,sigmapre=s1)$GBN	  

# Algorithm

MCMC.step(firstGBN, data, alpha=0.05, lambda = 0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
