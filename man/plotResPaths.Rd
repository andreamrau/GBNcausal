\name{plotResPaths}
\alias{plotResPaths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot mean and variance trajectory
}
\description{
Plot the trajectory of the estimations of mean and variance of the GBNetwork estimated. 
}
\usage{
plotResPaths(full.run, trueM = NULL, trueS = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{full.run}{
list - full.run can be a result of MCMC.GBN function. It's a list of gaussian bayesian networks. They must have the same number of nodes and the same names.
}
  \item{trueM}{
vector - Mean of the true GBNetwork. 
}
  \item{trueS}{
vector - Variance of the true GBNetwork.
}
}

\value{
Plot as many graphs as there are elements in variance and mean vectors. 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plotCausalPaths}}
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

plotResPaths(results$full.run, m1, s1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
