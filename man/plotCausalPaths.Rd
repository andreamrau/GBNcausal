\name{plotCausalPaths}
\alias{plotCausalPaths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot causal effects trajectory 
}
\description{

Plot the trajectory of direct or total causal effects values
}
\usage{
plotCausalPaths(alphaRes, trueW = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alphaRes}{
matrix - Direct or total causal effects matricies, one per line. They are all estimators of trueW. 
}
  \item{trueW}{
matrix -  True causal effects of the GBNetwork. 
}
}

\value{
Plot as many graphs as there is causal effects. 
}

\seealso{
\code{\link{plotResPaths}}
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
			  
alphaRes <- causalEffects(results$full.run)$alphaRes
plotCausalPaths(alphaRes,W)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
