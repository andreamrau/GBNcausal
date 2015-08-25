\name{plotingPostOrder}
\alias{plotingPostOrder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ploting post order
}
\description{
This function plot the probability a node has to be on the ith position in the graph, according to the results of a MCMC.GBN run or a ObsOnly run.

The more dark the color is, the higher is the probability.  
}
\usage{
plotingPostOrder(W, refNames, obs = "", GBN.all)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
matrix - The true weightmatrix (direct causal effects) of the GBN. 
}
  \item{refNames}{
vector - A vector of character, names of the nodes N# where # is a number. 
}
  \item{obs}{
If obs = "obsOnly", it means the type of results is produced by observational data only. 
}
  \item{GBN.all}{
list - A list of full.run (MCMC.GBN or ObsOnly results).
}
}

\value{
It returns as many graphs as there is full.run in GBN.all list. 
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

results=MCMC.GBN(data, firstGBN, nbSimulation=2000, burnIn=20, seq=1, verbose=TRUE,verbose.index=100, 
              alpha=1,lambda=0) 
			  

refNames <- paste("N",1:10,sep = "")
plotingPostOrder(W1,refNames,obs = "",GBN.all = list(results))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
