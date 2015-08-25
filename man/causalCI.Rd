\name{causalCI}
\alias{causalCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence interval for direct causal effects 
}
\description{
% Chaque colonne de alphaRes représente une valeur de la matrice des effets causaux totaux ou directs. 
% La matrice que retourne cette fonction est la matrice des moyennes de ces colonnes. 
% De plus, pour chaque moyenne, on vérifie qu'elle est entre la borne inférieure et la borne supérieure de l'intervalle, calculé comme : X[0.05*n], X[0.95*n] 
}
\usage{
causalCI(alphaRes, CIlb = 0.05, CIub = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alphaRes}{
matrix - contains a matrix of direct causal effects per line. 
}
  \item{CIlb}{
double in [0,1] - lower bound of the interval.
}
  \item{CIub}{
double in [0,1] - upper bound of the interval.
}
}

\value{
It returns a matrix of the values that are in the confidence interval. Those found out of it are equals to 0 in the matrix. 
}

\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
