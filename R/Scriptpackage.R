#########################################################
#########################################################
## Script pour la création du package GBN_causal : fonctions
#########################################################
#########################################################


########################################
# Class GBNetwork:
# W (WeightMatrix), m (resMean), sigma (resSigma)
########################################
setClass("GBNetwork",representation = list(WeightMatrix="matrix",resMean="vector",resSigma="vector"))

########################################
# Validation GBNetwork object (dimensions and acyclicity)
########################################
validGBNetworkObject <- function(object) {
  valid = 1
  if(length(object@resMean) != length(object@resSigma) ||
       length(object@resMean) != length(object@WeightMatrix[,1]) ||
       length(object@resSigma) != length(object@WeightMatrix[1,]))
  {
    valid = 0; paste("Unequal length of arguments")
  }
  if(length(names(object@resMean)) == 0 || length(names(object@resSigma)) == 0 ||
       length(colnames(object@WeightMatrix)) == 0 || length(rownames(object@WeightMatrix)) == 0)
  {
    valid = 0; paste("resMean, resSigma, and WeightMatrix must all have colnames", sep="")
  }
  if(sum(names(object@resMean)!=names(object@resSigma))>0 ||
       sum(names(object@resMean)!=rownames(object@WeightMatrix))>0 ||
       sum(names(object@resMean)!=colnames(object@WeightMatrix))>0 ||
       sum(names(object@resSigma)!=rownames(object@WeightMatrix))>0 ||
       sum(names(object@resSigma)!=colnames(object@WeightMatrix))>0)
  {
    valid = 0; paste("resMean, resSigma, and WeightMatrix must have same dimnames", sep="")
  }
  valid = valid*isAcyclic(abs(sign(object@WeightMatrix)))
  return(valid==1)
}
setValidity("GBNetwork", validGBNetworkObject)


########################################
# Format data
# data is an (n x p) matrix
# int.nodes is an (n x p) matrix with 1's indicating a knock-out for a given gene in a given sample
# int.means is an (n x p) matrix with corresponding values for knock-out mean
########################################
dataFormat <- function(x, int.nodes=c(), int.means=c()) {
  n <- nrow(x); p <- ncol(x)
  if(is.null(int.nodes[[1]])) {
    int.nodes=matrix(0, nrow=n, ncol=p)
    int.means=matrix(0, nrow=n, ncol=p)
  }
  if(is.null(int.means[[1]])) {
    int.means=matrix(0, nrow=n, ncol=p)
  }
  if(is.null(colnames(x)[1])==TRUE) {
    colnames(x)=colnames(int.nodes)=colnames(int.means)=paste("N", 1:ncol(x), sep="")
  }
  if(is.null(colnames(x)[1])!=TRUE) {
    colnames(int.nodes)=colnames(int.means)=colnames(x)
  }
  data=list(x=as.matrix(x), int.nodes=int.nodes, int.means=int.means)
  return(data)
}

########################################
# Permute node order in data and GBNnetwork class
########################################
dataOrder <- function(data, ord) {
  x0=data$x[,ord]
  int.nodesO=data$int.nodes[,ord]
  int.meansO=data$int.means[,ord]
  dataOrder=list(x=x0, int.nodes=int.nodesO, int.means=int.meansO)
  return(dataOrder)
}
GBNOrder <- function(GBN, ord) {
  WeightMatrixO=GBN@WeightMatrix[ord,ord]
  resMeanO=GBN@resMean[ord]
  resSigmaO=GBN@resSigma[ord]
  GBNO=new("GBNetwork", WeightMatrix=WeightMatrixO, resMean=resMeanO, resSigma=resSigmaO)
  return(GBNO)
}

#######################################################
# Metropolis Hastings algorithm: Max log-likelihood for a given graph
#######################################################

MCMC.GBN = function(data, firstGBN, nbSimulation=20000, burnIn=5000, seq=25, verbose=FALSE,
                     verbose.index=500, alpha=0.05, alpha2 = 0.05, lambda = 0, 
                     listblocks = list(),
                     str=matrix(1,length(firstGBN@resSigma),length(firstGBN@resSigma)),
                    type = "")
{
  if(is.null(colnames(data$x)[1])==TRUE) stop("data matrix must have column names");
       if(is.null(colnames(data$int.nodes)[1])==TRUE) stop("int.nodes must have column names");
       if(is.null(colnames(data$int.means)[1])==TRUE) stop("int.means must have column names");
       if(sum(colnames(data$x)!=colnames(data$int.nodes)) > 0 ||
            sum(colnames(data$x)!=colnames(data$int.means)) > 0 ||
            sum(colnames(data$int.means)!=colnames(data$int.nodes)) > 0 ||
            sum(colnames(data$x)!=colnames(firstGBN@WeightMatrix)) > 0)
         stop("dimnames must match");
		 
  if(type == "obsOnly"){
	p <- length(colnames(data$int.nodes))
    Trajectoire <- list()
    mu <- firstGBN@resMean
    sigma <- firstGBN@resSigma
    
    for(i in 1:nbSimulation){
      W = matrix(0, nrow=p, ncol=p)
      W[upper.tri(W)] = 1
      ref <- colnames(firstGBN@WeightMatrix)
      newOrder <- uniformBlocks(ref,listblocks)$order
      rownames(W) = colnames(W) = names(mu) = names(sigma) = paste("N",1:p,sep="")[newOrder]
      dataO = dataOrder(data, match(colnames(W), colnames(data$x))) 
      firstGBN = new("GBNetwork",WeightMatrix=W,resMean=mu,resSigma=sigma)
      mle = GBNmle(firstGBN, dataO, lambda = lambda)$GBN
      ord = order(as.numeric(substr(colnames(dataO$x), 2, 10))) ## Assume names are "N#"
      mle = GBNOrder(mle, ord)
      Trajectoire[[i]] <- mle
    }
    return(full.run = Trajectoire)
    
  }
  else{
      Trajectoire = list()
       count = 0; u = 1; accept = 0;
       if(verbose==TRUE) print("Running MCMC...")
       message("Function currently assumes node names are of the form N###\n")
       
       ## Start out with ordered names for GBN & data to facilitate matching up later
       ord = order(as.numeric(substr(colnames(data$x), 2, 10))) ## Assume names are "N#"
       firstGBN = GBNOrder(firstGBN, ord)
       str=str[ord,ord]
       data = dataOrder(data, ord)
       nextGBN = firstGBN
       for(i in 1:nbSimulation)
       {
         count = count + 1
         if(verbose==TRUE & floor(count/verbose.index)==count/verbose.index)
           print(paste("Iteration:", count));
         MCstep = MCMC.step(nextGBN, data, alpha=alpha, alpha2 =alpha2,lambda=lambda,
                            listblocks = listblocks, str=str)
         nextGBN = MCstep$newGBN
         accept = accept + MCstep$accept
         if(verbose==TRUE & count==burnIn+1) print("Burn-in complete.")
         if(count%%seq==0 && count>burnIn)
         {
           Trajectoire[[u]] = nextGBN       
           u = u+1
         }
       }
       return(list(full.run=Trajectoire, accept.rate=accept/nbSimulation))
       
       
    
  }
  
}

########################################
# Iteration in Metropolis Hastings algorithm
# type = c("perEdgeProposal", "nodeOrderProposal")
########################################
MCMC.step<- function(GBN, data, alpha=0.05, alpha2 = 0.05, lambda,listblocks = list(),
                       str=matrix(1,length(GBN@resSigma),
                                  length(GBN@resSigma)))
{ 
  ## Sanity check that all dimnames match:
  if(sum(colnames(data$x)!=colnames(data$int.nodes))>0 ||
       sum(colnames(data$x)!=colnames(GBN@WeightMatrix))>0 ||
       sum(colnames(data$x)!=names(GBN@resMean))>0 ||
       sum(colnames(data$x)!=names(GBN@resSigma))>0) stop("dimnames do not match")
  
  ord = topOrder(abs(sign(GBN@WeightMatrix))) # Put previous network back into upper triangual
  GBN.ord = GBNOrder(GBN, ord); data.ord = dataOrder(data, ord)
  new = newMallowsProposal(GBN.ord, data.ord, alpha,alpha2,lambda,
                            listblocks = listblocks,str[ord,ord])
  hatnewGBN = new$GBN; newdata = new$data
  
  ## Permute nodes in former GBN using topOrder()
  ord = topOrder(abs(sign(GBN@WeightMatrix)))
  oldGBN = GBNOrder(GBN, ord); olddata = dataOrder(data, ord)
  
  ## Calculate Metropolis-Hastings acceptance rate
  ## Watch out for how nodes are ordered
  acceptanceRate = min(1,exp(GBNlikelihood(hatnewGBN,newdata)-GBNlikelihood(oldGBN,olddata)))
  accept = 1
  if(runif(1)>acceptanceRate) {
    accept = 0
    hatnewGBN = oldGBN
  }
  ## Put nodes back into correct order
  ord = order(as.numeric(substr(colnames(hatnewGBN@WeightMatrix), 2, 10))) ## Assume names are "N#"
  newGBN = GBNOrder(hatnewGBN, ord)
  return(list(newGBN=newGBN, accept=accept))
}


########################################
# Mallows order proposal for new order (full network estimation)
# Full estimation of graph....
#
# sampling through RIM = Repeated Insertion Model (Doignon et al. 2004)
# found in Lu and Boutilier (2011)
# beta>=0 temperature: beta=0 dirac ref, beta=Inf uniform
# ref order and dispersion phi in [0,1]
# -log phi = lambda = 1/beta
## listblocks = liste des paquets de noeuds dont on sait qu'ils sont ensemble
## l = list(c("N1","N2","N3"),c("N5","N6","N7","N8","N9")) pour potentiellement plus de genes.
########################################

rmallows <- function(ref,beta) {
  m=length(ref);
  phi=exp(-1/beta);
  order=1;
  for (i in 2:m) {
    if (phi<1) {
      p=(1-phi)*phi^(i-(1:i))/(1-phi^i);
      j=sample(1:i,size=1,prob=p)
    } else {
      j=sample(1:i,size=1)
    }
    aux=rep(i,i); aux[-j]=order;
    order=aux;
  }
  return(list(order=order, new=ref[order]))
}

rmallowsBlocks <- function(ref,beta,beta2,listblocks){
  l <- length(listblocks)
  if(l>0 && l != length(ref)){
    finalOrder <- numeric() 
    
    listecomplete <- list()
    for(i in 1:length(listblocks)){
      listecomplete[[i]] = listblocks[[i]]
    }
    
    for(i in 1:length(ref)){
      if(length(grep(ref[i],listblocks))==0){
        listecomplete=c(listecomplete,ref[i])
      }
    } 
    
    blocksOrder <- rmallows(listecomplete,beta)$order
    
    listecomplete <- listecomplete[blocksOrder]
    
    vecteurcomplet <- c()
    
    for(k in 1:length(listecomplete)){
      m <-length(listecomplete[[k]])
      if(m==1){
        vecteurcomplet <- c(vecteurcomplet,listecomplete[[k]]) 
      }
      else{ 
        order <-rmallows(listecomplete[[k]],beta2)$order     
        listecomplete[[k]] <- listecomplete[[k]][order]
        vecteurcomplet <- c(vecteurcomplet,listecomplete[[k]])
      }
    }   
    
    for(i in 1:length(vecteurcomplet)){
      a <- match(vecteurcomplet[i],ref)
      finalOrder <- c(finalOrder, a)
    }  
    return(list(order = finalOrder, new = ref[finalOrder]))
  }
  
  if(l == length(ref)){
    return(rmallows(ref,beta))
  }
  
  else{
    
    return(rmallows(ref,beta))}
} 

uniformBlocks <- function(ref, listblocks){
  l <- length(listblocks)
  finalOrder <- numeric()
  
  if(l>0){
    
    listecomplete <- list()
    for(i in 1:length(listblocks)){
      listecomplete[[i]] = listblocks[[i]]
    }
    
    for(i in 1:length(ref)){
      if(length(grep(ref[i],listblocks))==0){
        listecomplete=c(listecomplete,ref[i])
      }
    } 
    
      finalOrder <- numeric(0) 
      blocksOrder <- sample(1:length(listecomplete), length(listecomplete))
      listecomplete <- listecomplete[blocksOrder] 
      vecteurcomplet <- c()
      
      for(k in 1:length(listecomplete)){
        m <-length(listecomplete[[k]])
        if(m==1){
          vecteurcomplet <- c(vecteurcomplet,listecomplete[[k]]) 
        }
        else{ 
          order <- sample(1:m,m)     
          listecomplete[[k]] <- listecomplete[[k]][order]
          vecteurcomplet <- c(vecteurcomplet,listecomplete[[k]])
        }
      }
      
      for(j in 1:length(vecteurcomplet)){
        a <- match(vecteurcomplet[j],ref)
        finalOrder <- c(finalOrder, a)
         
    }
  }
  else{finalOrder <- sample(1:length(ref),length(ref), replace = FALSE)}
  
  return(list(order = finalOrder,new = ref[finalOrder]))
}

newMallowsProposal <- function(GBN, data, alpha, alpha2, lambda=0, 
                                listblocks = list(), 
                                str=matrix(1,length(GBN@resSigma),length(GBN@resSigma))) {
  ## alpha represents the Mallows proposal temperature
  W=GBN@WeightMatrix
  nodeOrder=colnames(W)
  newNodeOrder=rmallowsBlocks(nodeOrder, alpha,alpha2,listblocks)$order
  ## Permute nodes according to move above, graph is full <- 
  newGBN = GBNOrder(GBN, newNodeOrder) 
  newGBN@WeightMatrix = 0*newGBN@WeightMatrix; 
  newGBN@WeightMatrix[upper.tri(newGBN@WeightMatrix)] = 1;
  str = str[newNodeOrder,newNodeOrder]
  newGBN@WeightMatrix=newGBN@WeightMatrix* str;
  newGBN@resMean = 0*newGBN@resMean
  newGBN@resSigma = 0*newGBN@resSigma
  newdata = dataOrder(data, newNodeOrder)

  # Full network estimation
  newGBN = GBNmle(newGBN, newdata, lambda=lambda, sigmapre=GBN@resSigma)$GBN
  return(list(GBN=newGBN, data=newdata))
}


########################################
# Maximize parameters m, s, and W for a given structure with penality lambda
########################################
GBNmle <- function(GBN, data, lambda=0, sigmapre=rep(0,dim(data$x)[2]))
{
  ## Sanity check that all dimnames match, W is upper triangular:
  if(sum(colnames(data$x)!=colnames(data$int.nodes))>0 ||
       sum(colnames(data$x)!=colnames(GBN@WeightMatrix))>0 ||
       sum(colnames(data$x)!=names(GBN@resMean))>0 ||
       sum(colnames(data$x)!=names(GBN@resSigma))>0) stop("dimnames do not match")
  if(sum(diag(GBN@WeightMatrix))!=0 || sum(GBN@WeightMatrix[lower.tri(GBN@WeightMatrix)])!=0)
    stop("W must be upper triangular!") 
  struct = GBN@WeightMatrix!=0
  x = data$x
  p=dim(x)[2]; n=dim(x)[1];
  int = matrix(as.logical(data$int.nodes), nrow=n, ncol=p)
  if(length(which(data$int.means != 0)) > 0) stop("Only knockouts currently supported.")
  # informative sample size by variable
  N=apply(!int,2,sum); 
  # centered data
  y=array(NA,dim=c(dim(x),ncol(int)));
  for (j in 1:ncol(int)) y[,,j]=t(t(x)-apply(x[!int[,j],],2,mean));
  
  # solve W without using shortcut
  d=sum(struct);
  num=0*struct; num[struct]=1:d;
  b=rep(0,d);
  A=matrix(0,nrow=d,ncol=d);
  for (j in which(apply(struct,2,sum)>0)) {
    aux=t(y[!int[,j],,j])%*%y[!int[,j],,j]; #matrice p*p
    for (i in which(struct[,j])) {
      b[num[i,j]]=aux[i,j];
      for (ii in which(struct[,j])) {
        A[num[i,j],num[ii,j]]=aux[i,ii];
        if (ii==i) A[num[i,j],num[ii,j]]=A[num[i,j],num[ii,j]]+2*lambda*sigmapre[j]
      }
    }
  }
  # A<- as(A,"sparseMatrix")
  W=0*struct;
  W[struct]=solve(A,b);
  # then s
  I=diag(rep(1,p));
  s=rep(0,p);
  for (j in 1:p) s[j]=sqrt(apply((y[!int[,j],,j]%*%(I-W))^2,2,mean)[j])       
  # then m
  m=rep(0,p);
  for (j in 1:p) m[j]=apply(x[!int[,j],]%*%(I-W),2,mean)[j];       
  # fix dimnames
  rownames(W) = colnames(W) = names(m) = names(s) = colnames(GBN@WeightMatrix)
  GBNest = new("GBNetwork", WeightMatrix=W, resMean=m, resSigma=s)
  return(list(GBN=GBNest, A=A, b=b, y=y))
}


########################################
# Calculate log-likelihood given m, s, and W
# data = list(list(Data1,interv1,condit1),
#                list(Data2,interv2,condit2),...)
########################################
GBNlikelihood <- function(GBN, data)
{       
  ## Sanity check that all dimnames match:
  if(sum(colnames(data$x)!=colnames(data$int.nodes))>0 ||
       sum(colnames(data$x)!=colnames(GBN@WeightMatrix))>0 ||
       sum(colnames(data$x)!=names(GBN@resMean))>0 ||
       sum(colnames(data$x)!=names(GBN@resSigma))>0) stop("dimnames do not match")
  
  ## Other sanity check that W is upper triangular:
  if(sum(diag(GBN@WeightMatrix))!=0 || sum(GBN@WeightMatrix[lower.tri(GBN@WeightMatrix)])!=0)
    stop("W must be upper triangular!")
  
  struct = GBN@WeightMatrix!=0
  x = data$x
  p=dim(x)[2]; n=dim(x)[1];
  int = matrix(as.logical(data$int.nodes), nrow=n, ncol=p)
  if(length(which(data$int.means != 0)) > 0) stop("Only knockouts currently supported.")
  
  # informative sample size by variable
  N=apply(!int,2,sum);
  W=GBN@WeightMatrix
  m=GBN@resMean
  s=GBN@resSigma
  
  # calculate log likelihood
  first.term = -log(2*pi)/2 * sum(N)
  second.term = -sum(N*log(s))
  third.term = 0
  for(k in 1:n) {
    noint = int[k,]!=TRUE
    e = diag(1,p)
    third.term = third.term - (1/2)*sum(1/(s[noint])^2 * (x[k,noint] - x[k,]%*%W%*%e[,noint] - m[noint])^2)
  }
  LogVrais = first.term + second.term + third.term
  return(LogVrais)
}

########################################
# Simulate data from a GBN (intervention or conditional): Gregory's version
########################################
# L=inv(I-W) from parameter w
Lfct <- function(W) {
  d=dim(W)[1];
  I=diag(rep(1,d));
  L=I; tmp=I; i=0;
  while (i<d-1) {
    tmp=tmp%*%W; L=L+tmp; i=i+1;
  }
  return(L);
}
# covariance matrix from parameters s and w
Sfct = function(s,L) {
  return(t(L)%*%diag(s^2)%*%L);
}
# Simulate function
simulGBN = function(N,m,s,W,int=numeric(0),int_data=matrix(runif(N*length(int)),nrow=N),
                    seed) {
  set.seed(seed)
  s[int]=0; W[,int]=0;
  L=Lfct(W)
  Sigma=Sfct(s,L);
  nu=matrix(rep(m,N),nrow=N,byrow=TRUE);
  if (length(int)>0) nu[,int]=int_data;
  mu=nu%*%L;
  x=mvrnorm(N,mu=rep(0,length(m)),Sigma=Sigma)+mu;
  return(x);
}

########################################
## Calculate indirect and direct causal effects across all iterations  (alpha:Direct, beta:Indirect)
########################################
causalEffects = function(full.run)
{
  alphaRes = betaRes = matrix(NA, nrow=length(full.run), ncol=length(full.run[[1]]@WeightMatrix))
  for(i in 1:length(full.run))
  {
    alpha = full.run[[i]]@WeightMatrix
    beta = betaFromW(alpha)
    alphaRes[i,] = as.vector(alpha)
    betaRes[i,] = as.vector(beta)
  }
  return(list(alphaRes=alphaRes, betaRes=betaRes))
}

betaFromW <- function(Wgt)
{
  ## Permute nodes in former GBN using topOrder()
  ord = topOrder(abs(sign(Wgt)))
  Wo <- Wgt[ord,ord]
  beta <- Lfct(Wo)
  colnames(beta) <- colnames(Wo)
  rownames(beta) <- rownames(Wo)
  ## Put nodes back into correct order
  ord = order(as.numeric(substr(colnames(beta), 2, 10))) ## Assume names are "N#"
  beta <- beta[ord, ord]
  return(beta)
}


causalCI <- function(alphaRes, CIlb=0.05, CIub = 0.95) {
  qt = apply(alphaRes, 2, function(x) quantile(x, prob=c(CIlb,CIub)))
  zeroInCI = ifelse(qt[1,] < 0 & qt[2,] > 0 | qt[1,] == 0 & qt[2,] == 0, 0, 1)
  postMean = zeroInCI*apply(alphaRes, 2, mean)
  return(matrix(postMean, nrow=sqrt(dim(alphaRes)[2])))
}

########################################
## Examine posterior distribution of orders
########################################

posteriorOrder <- function(full.run) {
  p <- lapply(full.run, function(x) topOrder(abs(sign(x@WeightMatrix))))
  return(do.call("rbind", p))
}

########################################
## Plotting functions
########################################

plotCausalPaths = function(alphaRes, trueW=NULL) {
  par(mfcol = c(sqrt(ncol(alphaRes)),sqrt(ncol(alphaRes))), mar = c(0,0,0,0),
      oma = c(0,3,3,0))
  index = 1
  for(j in 1:sqrt(ncol(alphaRes))) {
    for(jj in 1:sqrt(ncol(alphaRes))) {
      plot(0,0, xlim=c(0,dim(alphaRes)[1]), ylim=c(min(alphaRes),max(alphaRes)),
           col="white", ylab="",xlab="", xaxt="n", yaxt="n")
      if(j == 1) mtext(outer=FALSE, side=2, line = 1, jj)
      if(jj == 1) mtext(outer=FALSE, side=3, line = 1, j)
      if(j != jj) {
        abline(h = 0, lty=2, col="grey")
        if(is.null(trueW[1])==FALSE) {
          lines(alphaRes[,index], col=ifelse(as.vector(trueW)[index]!=0,"red","black"))
          if(as.vector(trueW)[index]!=0) {
            legend("topright", bty="n", text.col="blue",
                   ifelse(sign(as.vector(trueW)[index])==1,"+","-"), cex=2)
          }
        }
        if(is.null(trueW[1])==TRUE) {
          lines(alphaRes[,index], col="black")
        }
      }
      index = index + 1
      if(j == jj) {
        rect(-100,-10,1000,10, col="black", bty="n")
      }
    }
  }
}

plotResPaths <- function(full.run, trueM=NULL, trueS=NULL) {
  par(mfrow = c(2, length(full.run[[1]]@resMean)), mar = c(0,0,0,0), oma = c(10,3,10,0))
  full.m <- do.call("rbind", lapply(full.run, function(x) x@resMean))
  full.s <- do.call("rbind", lapply(full.run, function(x) x@resSigma))
  for(j in 1:ncol(full.m)) {
    plot(full.m[,j], type = "l", ylim=c(min(full.m),max(full.m)), xaxt="n", yaxt="n")
    abline(h=0, lty = 2)
    if(is.null(trueM[1])==FALSE) {
      abline(h=trueM[j], lty=2, col="red")
    }
    if(j ==1) mtext(outer = FALSE, line = 1, side = 2, "M")
  }
  for(j in 1:ncol(full.s)) {
    plot(full.s[,j], type = "l", ylim=c(min(full.s),max(full.s)), xaxt="n", yaxt="n")
    abline(h=0, lty = 2)
    if(is.null(trueS[1])==FALSE) {
      abline(h=trueS[j], lty=2, col="red")
    }
    if(j ==1) mtext(outer = FALSE, line = 1, side = 2, "S")
  }
}

histCausalPaths <- function(alphaRes, trueW) {
  par(mfcol = c(sqrt(ncol(alphaRes)),sqrt(ncol(alphaRes))), mar = c(0,0,0,0),
      oma = c(0,3,3,0))
  index = 1
  for(j in 1:sqrt(ncol(alphaRes))) {
    for(jj in 1:sqrt(ncol(alphaRes))) {
      if(j!=jj) {
        hist(alphaRes[,index], xlim=c(min(alphaRes),max(alphaRes)),
             col=ifelse(as.vector(trueW)[index]!=0,"red","black"),
             main="", xlab="", ylab="", xaxt="n", yaxt="n")
        if(as.vector(trueW)[index]!=0) {
          legend("topright", bty="n", text.col="blue",
                 ifelse(sign(as.vector(trueW)[index])==1,"+","-"), cex=2)
        }
        
      }
      if(j==jj) {
        plot(0,0, xlim=c(0,dim(alphaRes)[1]), ylim=c(min(alphaRes),max(alphaRes)),
             col="white", ylab="",xlab="", xaxt="n", yaxt="n")
        rect(-100,-10,1000,10, col="black", bty="n")
      }
      if(j == 1) mtext(outer=FALSE, side=2, line = 1, jj)
      if(jj == 1) mtext(outer=FALSE, side=3, line = 1, j)
      index = index + 1
      if(j == jj) {
        rect(-100,-10,1000,10, col="black", bty="n")
      }
    }
  }
}

########################################
## Upper triangularization: CHECK THIS?
########################################

swap = function(mat, i, j)
  # Inverse les lignes et colonnes i et j d'une matrice, comme pour une permutation des indices i et j.
{
  n=length(mat[,1])
  v = c(1:(i-1),j,(i+1):(j-1),i,(j+1):n)
  if(i==1){v = v[-c(1,2)]}
  if(j==n){v = rev(rev(v)[-c(1,2)])}
  v = unique(v)
  mat = mat[v,v]
  return(mat)
}

detriangularisation = function(mat)
{       
  vect = as.numeric(substr(rownames(mat),2,10))
  n = length(vect)
  check = FALSE
  while(check==FALSE)
  {
    check = TRUE
    for(i in 1:n)
    {
      if(vect[i]<i)
      {
        mat = swap(mat,vect[i],i)
        vect = as.numeric(substr(rownames(mat),2,10))
        check = FALSE
      }
    }
  }
  return(mat)
}


#############################################################################
#############################################################################
# Fonction pour tracer les ordres
#############################################################################
#############################################################################

plotingPostOrder <- function(W,refNames,obs = "",GBN.all){  
  h  <- vector("list", length(GBN.all))
  for(i in 1:length(GBN.all)){
    direct <- W
    diag(direct) <- NA
    total <- betaFromW(W)
    diag(total) <- NA
    total <- total[refNames, refNames]
    nam <- factor(colnames(direct), levels=paste("N", 1:10, sep=""))
    ooo <- order(nam)
    direct <- direct[ooo,ooo]
    total <- total[ooo,ooo]  
    p <- nrow(direct)
    ## HEATMAP OF POSTERIOR ORDER DISTRIBUTION
    order.proportions <- matrix(0, nrow=10, ncol=10)
    colnames(order.proportions) <- 1:10
    rownames(order.proportions) <- refNames
    cat(i, "\n")
    choice <- c("total")
    tr <- choice
    truth <- get(tr)
    ## Read in results
    if(obs != "obsOnly") {
      GBN <- GBN.all[[i]]$full.run
    }
    if(obs == "obsOnly") {
      GBN <- GBN.all[[i]]
    }  
    ppp <- do.call("rbind", lapply(GBN, function(x) topOrder(abs(sign(x@WeightMatrix)))))
    ord <- t(apply(ppp, 1, order))
    colnames(ord) <- colnames(GBN[[1]]@WeightMatrix)
    ## Put nodes back in same order as reference
    ord <- ord[,refNames]
    for(p in 1:10) {
      tab <- table(ord[,p])/nrow(ord)
      order.proportions[p,names(tab)] <- order.proportions[p,names(tab)] + tab 
    }
    
    order.proportions <- order.proportions / 100
    ordprop <- melt(order.proportions)
    colnames(ordprop) <- c("Node", "Estimated", "value")
    ordprop$Node <- factor(ordprop$Node,levels=paste(refNames[10:1], sep=""))
    ordprop$Estimated <- factor(ordprop$Estimated,levels=1:10)
    ordprop <- arrange(ordprop,Node,Estimated)
    p <- ggplot(ordprop, aes(Estimated, Node)) + 
      geom_tile(aes(fill = value),
                colour = "white") + scale_fill_gradient(low = "white",
                                                        high = "steelblue")

    p <- p+theme(axis.title.x = element_text(size = 20),
                 axis.title.y = element_text(size = 20))
    h[[i]] <- p
    i <- i+1
  }
  
  return(h)
  
}

multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols)    # Number of rows needed
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
}
#############################################################################
#############################################################################
# Fonctions to create data
#############################################################################
#############################################################################


# nbData : number of observationnals data in X
# p : number of genes
# KO : list of wanted knock out vectors (allows multiple knock out)
# nbKO : vector of number of replicat wanted per knock out
# type : type of true graph (complet, BMC, Unique) 
# W : (if type not used) weighted matrix of true graph

dataCreate <- function(nbData,p,KO = list(), nbKO = c(),W = 1*upper.tri(matrix(0,p,p)),
                       m,sigma,seed){
  
  X<-simulGBN(nbData,m,sigma,W, seed = seed) 
  n <- sum(nbKO)
  intnode = matrix(0,nbData+n,p)
  
  if(length(KO)>0){
    for(i in 1:length(KO)){
      X <- rbind(X,simulGBN(nbKO[i],m,sigma,W,int = KO[[i]],
                            int_data = 0,seed = seed))
    }
    
    for(i in 1:length(KO)){
      if(i == 1){
        k <- nbData
      }
      else{
        k <- k+nbKO[i-1]
      } 
      intnode[(k+1):(k+nbKO[i]),KO[[i]]]=1
    }
  }
  
  colnames(X)=colnames(intnode)=colnames(W) = rownames(W)=paste("N",1:p,sep="")
  data<-dataFormat(X,intnode)
  
  return(list(data = data, X=X, W=W, intnode = intnode))
}
