\name{dataFormat}
\alias{dataFormat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Put data into a format usable by the main functions
}
\description{
Put the data into a format usable by the main functions of the package. 
}
\usage{
dataFormat(x, int.nodes = c(), int.means = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
matrix - Observations or simulated observations from a GBN. 
}
  \item{int.nodes}{
matrix - 0 indicates that there is no intervention on the node. 1 indicates that the node is knock-out. Only those two values are supported by the MCMC.GBN function.
}
  \item{int.means}{
matrix - Values of the knock-out : if int.nodes = 1 then x = int.means. Only 0 is currently supported by the MCMC.GBN function.
}
}

\value{
data is a list of x, int.nodes and int.means with the same colnames/rownames. If there is no colnames/rownames in x, names given by the function are of the form "N#" where # is a number. 
}

\seealso{
 \code{\link{dataCreate}}
}
\examples{

p<-10

W1=1*upper.tri(matrix(0,p,p))
m1=rep(0,10)
sigma1=rep(0.1,10)

X=simulGBN(2*p,m1,sigma1,W1,seed = 1990)
X=rbind(X,simulGBN(p,m1,sigma1,W1,int=1,int_data=0,seed = 1990)) 
X=rbind(X,simulGBN(p,m1,sigma1,W1,int=9,int_data=0,seed = 1990)) 

intnode=matrix(0,4*p,p)
intnode[(2*p+1):(3*p),p-9]=1
intnode[(3*p+1):(4*p),p-1]=1

colnames(X)=colnames(intnode)=paste("N",1:p,sep="")

dataFormat(X,intnode) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
