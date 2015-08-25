\name{causalEffects}
\alias{causalEffects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Direct and total causal effects trajectories
}
\description{
Computes direct and total causal effects of a list of GBN objects.

This function enables to get back the trajectories of direct and undirect causal effects through a MCMC.GBN run. 
}
\usage{
causalEffects(full.run)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{full.run}{
list - full.run can be a result of MCMC.GBN function. It's a list of gaussian bayesian networks. They must have the same number of nodes and the same names. 
}
}

\value{
This function returns a list of two arguments : 
\item{alphaRes}{
matrix - direct causal effects. The ith line of alphaRes corresponds to the direct causal effects matrix of the ith GBN of full.run.
}

\item{betaRes}{
matrix - total causal effects. The ith line of alphaRes corresponds to the direct causal effects matrix of the ith GBN of full.run.
}
}

\examples{

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

alphaRes <- causalEffects(results$full.run)$alphaRes
betaRes <- causalEffects(results$full.run)$betaRes

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
