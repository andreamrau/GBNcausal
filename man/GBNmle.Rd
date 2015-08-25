\name{GBNmle}
\alias{GBNmle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum likelihood estimator
}
\description{
This function computes the maximum likelihood estimator of the object GBN calculated according to the data. 

The system to solve can be rewritten as : AW=b, where W is the weightmatrix of the GBN. 
}
\usage{
GBNmle(GBN, data, lambda = 0, sigmapre = rep(0, dim(data$x)[2]))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GBN}{
An object of the class GBN. 
}
  \item{data}{
data - Data is a format obtained by the function dataFormat.
}
  \item{lambda}{
logarithmic - coefficient of the penalty Ridge. 
}
  \item{sigmapre}{
vector - In the algorithm MCMC.GBN, sigmapre is the variance of the GBN of the previous iteration. If there is no previous iteration, ie if the function is used alone, sigmapre is a null vector.  
}
}

\value{

\item{GBN }{The estimator of the initial GBN}
\item{A }{A member of the equation AW = b solved by the function.}
\item{b }{b member of the equation AW = b solved by the function.}
\item{y }{centered observation (data$x).}
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
firstGBN = GBNmle(firstGBN,data,lambda=10e10,sigmapre=rep(0.1,p))$GBN

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
