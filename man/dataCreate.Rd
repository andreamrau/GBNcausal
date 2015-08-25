\name{dataCreate}
\alias{dataCreate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating data
}
\description{
This function does the same as simulGBN and dataFormat combined. 
}
\usage{
dataCreate(nbData, p, KO = list(), nbKO = c(), W = 1 * upper.tri(matrix(0, p, p)), m, sigma, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nbData}{
integer - Number of observational replicats.
}
  \item{p}{
integer - Number of nodes in the network. 
}
  \item{KO}{
list of vectors - Indicates on which nodes there is an intervention. 
}
  \item{nbKO}{
vector of integer - Number of replicat for each KO.
}
  \item{W}{
matrix - Matrix of direct causal effects of the graph. 
}
  \item{m}{
vector - Mean of the law that generate thes data.
}
  \item{sigma}{
vector - Variance of the law that generates the data
}
  \item{seed}{
integer - 
}
}

\value{
The function returns a list :

	\item{data }{data that can be used by the main functions of this package.}
	\item{X}{Simulated values of the GBNetwork nodes.}
	\item{W}{Matrix of direct causal effects of the GBNetwork}
	\item{intnode}{matrix - 0 indicates that there is no intervention on the node. 1 indicates that the node is knock-out.}

}

\examples{
seed = 1990
p <- 10
m<-rep(0,10)
sigma<-rep(0.1,10) 

W <- 1*upper.tri(matrix(0,p,p))

# Simple knock-out on two nodes

dataCreate(nbData = 2*p, p = 10,KO = list(1,9), nbKO = c(p,p), W = W , m = m,sigma = sigma, seed = seed)

# Multiple knock-out on 2 nodes

dataCreate(nbData = 2*p, p = 10,KO = list(c(6,10)), nbKO = c(p), W = W , m = m,sigma = sigma, seed = seed)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
