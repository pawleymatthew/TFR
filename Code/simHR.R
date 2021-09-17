rHRclusters <- function(n, d, L = 10^3, beta=2){
  
  # description:
    # generate sample vectors from HR distribution with an angular measure supported on mutually perpendicular faces of specified dimensions
  # inputs:
    # n : number of sample vectors to be generated
    # d : a vector of the face dimensions, e.g. c(4,3,2) gives three faces and total dimension of sum(d)=9
    # L : a large constant to ensure approximate asymptotic independence between groups
  # outputs:
    # X : n x sum(d) matrix of the sample vectors
    # Gamma : the sum(d) x sum(d) variogram matrix
    # Chi : a sum(d) x sum(d) matrix containing the extremal dependence coefficients chi_ij
  
  total_d <- sum(d)
  cumsum_d <- cumsum(d)
  # generate the variogram
  H <- matrix(data = EnvStats::rpareto(total_d^2, location = 1, shape = 2.5), nrow = total_d, ncol = total_d)
  Gamma <- matrix(data = NA, ncol = total_d, nrow = total_d)
  for (i in 1:total_d){
    for (j in 1:total_d){
      # if i!=j are in same 'cluster', set Gamma_ij ~= 0 for asymptotic dependence
      # if i=j, set Gamma_ij = 0 for total asymptotic dependence
      if (sum(cumsum_d >= i) == sum(cumsum_d >= j)){ 
        Gamma[i,j] <- (beta/total_d) * sum((H[,i]-H[,j])^2)
        # otherwise, set Gamma_ij >> 0 for asymptotic independence
      } else
        Gamma[i,j] <- L
    }
  }
  # generate the sample vectors X1,...,Xn
  simdata <- graphicalExtremes::rmstable(n, model = "HR", d = total_d, par = Gamma)
  # compute matrix of pairwise extremal dependence coefficients
  Chi <- 2 * (1-pnorm(sqrt(Gamma)/2))
  # return list of outputs
  output <- list("X" = simdata, "Gamma" = Gamma, "Chi" = Chi)
  return(output)
}

meanChi <- function(Chi){
  
  # description:
    # compute mean of chi_{ij} for i!=j belonging to the same group (illustrates strength of asymptotic dependence within groups)
  # inputs:
    # Chi : matrix of pairwise extremal dependence coefficients (symmetric with ones of diagonal and zeroes off the block diagonals)
  # outputs:
    # mean of the relevant chi_ij values
  
  # make vector of lower triangular entries of Chi (excluding diagonal)
  Chi_vals <- as.numeric(Chi * lower.tri(Chi, diag = FALSE))
  # extract all non-zero elements
  Chi_vals <- Chi_vals[Chi_vals > 1e-8]
  # compute mean
  return(mean(Chi_vals))
}

