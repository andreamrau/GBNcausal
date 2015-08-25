\name{MCMC.GBN}
\alias{MCMC.GBN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC algorithm to estimate causal effetcts of a GBN object
}
\description{

MCMC.GBN is used to infer causal effects (direct or undirect) of a GB network. It is a random walk MCMC Metropolis Hasting algorithm with a Mallows proposal distribution. 
The algorithm explores the space of nodes order and find the best estimation of the causal effects. In each iteration, it computes the maximum likelihood estimation of the GBN law, based on the data. 


%% La fonction MCMC.GBN est utilisée pour inférer sur les effets causaux (directs et indirects) d'un GBN. C'est un algorithme MCMC Metropolis Hasting à marche aléatoire, dont la loi de proposition est la loi Mallows.
%% Le principe est de parcourir l'espace des ordres des nodes du réseau et trouver la meilleure estimation des effets causaux. Pour cela, on calcule à chaque itération l'estimation du maximum de vraisemblance
%% de la loi du GBN, basé sur les observations (data). 


}
\usage{
MCMC.GBN(data, firstGBN, nbSimulation = 20000, burnIn = 5000, seq = 25, verbose = FALSE, verbose.index = 500, 
 alpha = 0.05, alpha2 = 0.05, lambda = 0, listblocks = list(), str = matrix(1, length(firstGBN@resSigma), length(firstGBN@resSigma)), type = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data - Can be obtained by the function \code{\link{dataFormat}}.
}
  \item{firstGBN}{
GBN - An object of type GBN which is the initial value of the algorithm. data and firstGBN must have the same rownames and colnames for each elements. 
}
  \item{nbSimulation}{
integer - The number of iterations the function has to do. 
}
  \item{burnIn}{
integer - The number of iterations to do before starting to save the trajectory of the algorithm.
}
  \item{seq}{
integer - Every seq, if the number of iterations done is superior to the burnIn, the GBN found is saved into full.run. 
}
  \item{verbose}{
logical - If TRUE, every verbose.index interation a message is written to give the number of iterations made so far by the algorithm.  
}
  \item{verbose.index}{
integer - Indicates the number of iteration to do before writting a message if verbose is TRUE. 
}
  \item{alpha}{
double - First temperature of the mallows function. 
}
  \item{alpha2}{
double - Second temperature of the mallows function. If listblocks is empty, this parameter will not be used. 
}
  \item{lambda}{
logarithmic - Coefficient of the penalty Ridge. Can be 0.
}
  \item{listblocks}{
list - A list of nodes of the form (c("N1","N2"),c("N3","N4","N5")) where "N1","N2","N3","N4" and "N5" are elements of rownames and colnames of firstGBN elements and data elements. 
}
  \item{str}{
matrix - To improve the efficiency of the algorithm, a structure can be add. Colnames and rownames are not needed. 
Once the structure is set, the algorithm can't add any interaction that's not in the structure. It can remove some of them.
}
	\item{type}{
	string - Type of observations. If type = obsOnly, it means the data are observational only, and the function will use a uniform sampling each step instead of a Mallows distribution.}
}
\details{

This function only works with acyclic directed networks (DAG).  

It allows to insert an a priori information, which can be a structure or a list of blocks. 

A structure is in general a skeleton of a graph, or a CPDAG. It can be obtained by a PC algorithm (PCAlg) or a Lasso (glasso) method. 

The structure is a higher constraint than the list of blocks, in the way that the algorithm can't add any direct causal effect that is not in the structure. It can remove some of them. However, this is a good way of improving the
efficiency of the algorithm : it's faster and give better results if the structure corresponds to the network to estimate. 

The list of blocks is a list composed of vector of nodes that we know they are close in the graph.
 If the list is not empty, the algorithm uses the rmallowsBlocks function to propose a new order of nodes. The elements of a block can't be separated. 
 
This method can take observational and interventional data (only knock-out are implemented). Infering on observational data only can't give any good results.  

%% Cette méthode permet d'introduire une information a priori dans l'algorithme : sous la forme d'une structure ou sous la forme d'une liste de voitures. 
%% L'incorporation d'une structure est plus contraignante pour l'algorithme que celle de voitures, puisqu'une fois la structure fixée, l'algorithme ne peut ajouter un lien entre deux nodes si celui-ci n'est pas présent dans 
%% la structure. 

%%  La liste de voitures est une liste composée de vecteurs de nodes, dont on sait qu'ils se situent à côté dans le réseau. Si elle n'est pas vide, l'algorithme utilise la loi de proposition rmallowsBlocks. 

%% Cette méthode a été conçue dans le cadre de réseaux acycliques seulement. 

%% Elle permet de prendre en compte des données observationnelles et interventionnelles (seulement des knock-out).  L'inférence sur des données observationnelles seulement n'aboutira pas. Pour réellement avoir une estimation
%% des effets causaux, il faut que les data comportent des internventions. 
}
\value{
  \item{full.run }{The full trajectory of the algorithm : after the burnIn iteration, every seq iteration the GBN is saved into the full.run. Inference is done on full.run. }
  \item{accept }{The acceptance rate of the algorithm.}
}

\seealso{
 \code{\link{MCMC.step}}, \code{\link{newMallowsProposal}}, \code{\link{rmallows}}, \code{\link{rmallowsBlocks}}}
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
# Without constraint
results=MCMC.GBN(data, firstGBN, nbSimulation=2*n, burnIn=n, seq=1, verbose=TRUE,
             verbose.index=n/10, alpha=0.7,lambda=0.1)

# With constraint
listBlocks <- list(c("N1","N2","N3"),c("N4","N5"))

results_Blocks=MCMC.GBN(data, firstGBN, nbSimulation=2*n, burnIn=n, seq=1, verbose=TRUE,
                 verbose.index=n/10, alpha=0.7, alpha2 = 0.5, lambda=0.1, listblocks = listBlocks)
				 
# Observational data only :
data <- dataCreate(nbData = 2*p, p = 10,KO = list(), nbKO = c(), W = W , m = m,sigma = sigma, seed = seed)$data

firstGBN = new("GBNetwork",WeightMatrix=W1,resMean=m1,resSigma=s1)
firstGBN = GBNmle(firstGBN,data,lambda= 0,sigmapre=s1)$GBN

results_Blocks=MCMC.GBN(data, firstGBN, nbSimulation=100, burnIn=0, seq=1, verbose=TRUE,
                 verbose.index=n/10, alpha=0.7, alpha2 = 0.5, lambda=0.1, listblocks = listBlocks, type = "obsOnly")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
