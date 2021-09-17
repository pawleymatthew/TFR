prepPCA <- function(X){
  
  # description:
    # transform the raw data X to obtain Xtilde with Fréchet marginal distributions using rank-transform
  # inputs:
    # X : an nxd matrix of sample points (with arbitrary marginals)
  # outputs:
    # Xtilde : an nxd matrix of transformed sample points (with Fréchet marginals)
  
  # compute empirical CDFs 
  Fhat <- apply(X, 2, rank) / (nrow(X)+1)
  # transform to Fréchet marginals (mu=0, sigma=0, xi=2)
  frechet2 <- function(z) (-log(z))^(-1/2)
  Xtilde <- apply(Fhat, 2, frechet2)
  return(Xtilde)
}

estTPDM <- function(X, method = "v", u = 0.95) {
  
  # description:
    # estimate the TPDM for given data (with Fréchet marginals) using a specified threshold method and threshold quantile
  # inputs:
    # X : an nxd matrix of sample points (with Fréchet marginals)
    # method : which TPDM estimator to use, "v" for vector-based (default) or "p" for pairs-based
    # u : threshold quantile for choosing largest observations
  # outputs:
    # Sigma : a dxd matrix, the estimated TPDM
  
  # dimensions and preallocation
  n <- nrow(X)
  d <- ncol(X)
  Sigma <- matrix(NA, d, d)
  # compute TPDM estimate
  if (method == "v"){ # vector-based estimator
    r <- sqrt(rowSums(X^2))
    rstar <- quantile(r, u, na.rm = TRUE)
    bigr <- which(r > rstar)
    for (i in 1:d) {
      for (j in 1:d) {
        Sigma[i,j] <- d * sum(X[bigr,i]*X[bigr,j]/r[bigr]^2, na.rm = TRUE) / length(bigr)
      }
    }
  } else { # pairs-based estimator
    for (i in 1:d) {
      for (j in 1:d) {
        r <- sqrt(rowSums(X[,c(i,j)]^2))
        rstar <- quantile(r, u, na.rm = TRUE)
        bigr <- which(r > rstar)
        Sigma[i,j] <- 2 * sum(X[bigr,i]*X[bigr,j]/r[bigr]^2, na.rm = TRUE) / length(bigr) 
      }
    }
    Sigma <- nearPD(Sigma)[[1]] # find nearest pos. def. matrix
  }
  return(as.matrix(Sigma))
}

trans <- function(x){
  
  # description:
    # vectorised softplus function
  # inputs:
    # x : a vector of values in R
  # outputs:
    # y : a vector of values in (0,\infty), y_i = softplus(x_i)
  
  y <- log(1 + exp(x))
  y[is.infinite(y)] <- x[is.infinite(y)] # if exp caused Inf because x>>1, return original value, i.e. y_i=x_i
  return(y)
}

transinv <- function(x){
  
  # description:
    # vectorised inverse of softplus function
  # inputs:
  # x : a vector of values in (0,\infty)
  # outputs:
  # y : a vector of values in R, y_i = softplus^{-1}(x_i)
  
  y <- log(exp(x) - 1)
  y[is.infinite(y) && x > 1] <- x[is.infinite(y) && x > 1] # if exp caused Inf because x>>1, return original value, i.e. y_i=x_i
  y[is.infinite(y) && y < 0] <- -10^6 # if log caused Inf because x~=0, return y_i=-M for M>>0
  return(y)
}

extremePCA <- function(X, Sigma){
  
  # description:
    # compute the TPDM eigendecomposition and extremal principal components
  # inputs:
    # X : an nxd matrix of sample points (with Fréchet marginals)
    # Sigma : a dxd matrix, the estimated TPDM (must be positive definite)
  # outputs:
    # lambda : a d vector of the eigenvalues of Sigma
    # U : a dxd matrix whose columns are eigenvectors of Sigma (normalised to unit length and with U_{1,1}>0)
    # V : an nxd matrix of the extremal principal components 
  
  # dimensions
  n <- nrow(X)
  d <- ncol(X)
  # eigendecomposition of TPDM
  eigenSigma <- eigen(Sigma)
  U <- eigenSigma$vectors 
  if (U[1,1]<0) U <- -U
  lambda <- eigenSigma$values
  # extremal principal components
  V <- matrix(NA, nrow = n, ncol = d)
  for (t in 1:n){
    V[t,] <- t(U) %*% transinv(X[t,])
  }
  output <- list("lambda" = lambda, "U" = U, "V" = V)
  return(output)
}

reconstructX <- function(t, q, U, V){
  
  # description:
    # reconstruct event at time t based on the first q principal components
  # inputs:
    # t : time index of the event
    # q : number of principal components to be used
    # U : a dxd matrix whose columns are eigenvectors of Sigma
    # V : an nxd matrix of the extremal principal components 
  # outputs:
    # reconXt : a d vector of the reconstructed event at time t
  
  # compute weighted sum of vectors:
    # weights: vector of scalar weights
    # mat : matrix whose columns are the vectors
  weightedVecSum <- function(weights, mat) rowSums(t(t(mat) * weights))
  
  reconXt <- weightedVecSum(weights = V[t,1:q], mat = U[,1:q])
  reconXt <- trans(reconXt)
  return(reconXt)
}

simXtilde <- function(n, U, V, u, m){
  
  # description:
    # simulate event using Rohrbeck and Cooley (2021) algorithm
  # inputs:
    # n : number of events to simulate
    # U : a dxd matrix whose columns are eigenvectors of Sigma
    # V : an nxd matrix of the extremal principal components 
    # u : threshold quantile for choosing largest observations
    # m : number of principal components to be used in flexible model
  # outputs:
    # output : a list of outputs (kappa, W, Z, Wsample, Vsample, Zsample, Xtildesample)
  
  d <- ncol(V)
  # observed W vectors
  r <- sqrt(rowSums(V^2))
  rstar <- quantile(r, u, na.rm = TRUE)
  bigr <- which(r > rstar)
  W <- V[bigr,] / r[bigr] # W is #{extreme observations} x d matrix
  # observed Z vectors
  Z <- matrix(NA, nrow = nrow(W), ncol = m+1)
  Z[, 1:m] <- W[, 1:m]
  Z[, m+1] <- ifelse(W[, m+1] >= 0, sqrt(1-rowSums(W[,1:m]^2)), -sqrt(1-rowSums(W[,1:m]^2)))
  # fit Mises-Fisher mixture distribution to Z
  kappa <- Directional::vmfkde.tune(Z, low = 0.01, up = 1)[[1]]
  # generate n samples from the fitted Mises-Fisher model
  prob <- rep(1/nrow(Z), nrow(Z))
  k <- rep(1/kappa^2, nrow(Z))
  Zsample <- Directional::rmixvmf(n, prob, mu = Z, k)$x
  # extract indexes q=(q1,...,qn), nearest neighbours among the observations to each sample vector
  q <- max.col(Zsample %*% t(Z))
  # derive samples of W
  Wsample <- matrix(NA, nrow = n, ncol = d)
  Wsample[,1:m] <- Zsample[,1:m]
  Wsample[,-(1:m)] <- abs(Zsample[,m+1]/Z[q,m+1]) * W[q,-(1:m)]
  # generate the radii of the points
  rsample <- evd::rfrechet(n, loc = 0, scale = sqrt(d), shape = 2) # scale might be wrong, but shouldn't affect things
  # derive samples of V
  Vsample <- rsample * Wsample
  # derive sample of Xtilde
  Xtildesample <- matrix(NA, nrow = n, ncol = d)
  for (t in 1:n){
    Xtildesample[t,] <- reconstructX(t=t, q=d, U=U, V=Vsample)
  }
  # return outputs
  output <- list("W"=W, "Z"=Z, "kappa"=kappa, "Wsample"=Wsample, "Zsample"=Zsample, "Vsample"=Vsample, "Xtildesample"=Xtildesample)
  return(output)
}

simEventSet <- function(N, n, U, V, u, m, Sigmaobs){
  
  # description:
    # simulate multiple event sets using Rohrbeck and Cooley (2021)
  # inputs: 
    # N : number of events sets
    # n : number of events to simulate per set
    # U : a dxd matrix whose columns are eigenvectors of Sigma
    # V : an nxd matrix of the extremal principal components 
    # u : threshold quantile for choosing largest observations
    # m : number of principal components to be used in flexible model
    # Sigmaobs : the TPDM from the observed data
  # outputs:
    # diffSigma : a vector of Frobenius distances between the simulated and observed TPDMs for each set
  
  # vector to store the Frobenius distances
  diffSigma <- rep(NA, N)
  # simulate N sets and compute Frobenius distances
  for (i in 1:N){
    simXtildei <- simXtilde(n, U, V, u, m)
    Sigmai <- estTPDM(X = simXtildei$Xtildesample, method = "v", u = u)
    diffSigma[i] <- norm(Sigmai - Sigmaobs, type = "F")
  }
  return(diffSigma)
}

reconmerror <- function(U, V, u, m, Sigmaobs){
  
  # description:
    # compute Frobenius distance between reconstructed and observed TPDM
  # inputs: 
    # U : a dxd matrix whose columns are eigenvectors of Sigma
    # V : an nxd matrix of the extremal principal components 
    # u : threshold quantile for choosing largest observations
    # m : number of principal components to be used in flexible model
    # Sigmaobs : the TPDM from the observed data
  # outputs:
    # diffSigma : a vector of Frobenius distances between the simulated and reconstructed TPDMs for each set
  
  # vector to store the Frobenius distances
  reconset <- matrix(NA, nrow = nrow(V), ncol = ncol(V))
  # simulate N sets and compute Frobenius distances
  for (t in 1:nrow(V)){
    reconset[t,] <- reconstructX(t, m, U, V)
  }
  Sigmam <- estTPDM(X = reconset, method = "v", u = u)
  diffSigma <- norm(Sigmam - Sigmaobs, type = "F")
  return(diffSigma)
}