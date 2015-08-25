\name{multiplot}
\alias{multiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ploting multiple graphs on the same page
}
\description{
Plot multiple graphs on the same page. 
}
\usage{
multiplot(..., plotlist = NULL, cols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
ggplot objects
}
  \item{plotlist}{
list of ggplot objects, if ... is not used.
}
  \item{cols}{
integer - Number of columns in the graph
}
}

\value{
The function plots the graph in the configuration asked by the user.
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
			  
results2=MCMC.GBN(data, firstGBN, nbSimulation=6000, burnIn=2000, seq=1, verbose=TRUE,verbose.index=100, 
              alpha=1,lambda=0) 		  
GBN.all <- list(results, results2)

refNames <- paste("N",1:10,sep = "")
h <- plotingPostOrder(W1,refNames,obs = "",GBN.all)

multiplot(h[[1]],h[[2]],cols = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
