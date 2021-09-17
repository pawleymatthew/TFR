PAMcluster <- function(X, K){
  
  # DESCRIPTION:
  #   Computes K clusters from X using PAM with F-madogram distances
  # INPUTS:
  #   X : an nxd matrix, each column contains a time series for each location
  #   K : number of clusters
  # OUTPUTS:
  #   output : an object of class 'pam'
  
  d <- ncol(X) # number of sites
  n <- nrow(X) # number of measurements (denoted T in the paper)
  # F-madogram distance
  V = array(NaN, dim = c(n, d)) # initialise array to store empirical CDFs
  for(p in 1:d) {
    X.vec = as.vector(X[,p]) # time series for site p
    Femp = ecdf(X.vec)(X.vec) # compute empirical CDF
    V[,p] = Femp # store in V
  }
  D <- dist(t(V), method = "manhattan")
  # Clustering with PAM
  output = cluster::pam(x = D, k = K)
  return(output) 
}