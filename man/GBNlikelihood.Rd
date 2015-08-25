\name{GBNlikelihood}
\alias{GBNlikelihood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Loglikelihhood
}
\description{
Compute the loglikelihood between a GBN and data. 
}
\usage{
GBNlikelihood(GBN, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GBN}{
GBN - an object of type GBN. colnames and rownames of the Weightmatrix must be the same as the colnames and rownamesof elements from data. 
}
  \item{data}{
data - Data is a format obtained by the function dataFormat.
}
}

\value{
It returns the loglikelihood of GBN and data. 
}

\examples{
# Data creation

seed = 1990

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

GBnetwork = new("GBNetwork",WeightMatrix=W1,resMean=m1,resSigma=s1)
GBNlikelihood(GBnetwork,data)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
