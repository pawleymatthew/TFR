plot3Dangles <- function(X, dims, u = 0.95){
  
  # description:
    # plot three dimensions of d-dimensional sample data in a 3D scatter plot. Each point is projected into the unit 2-sphere.
  # inputs:
    # X : an nxd matrix of sample points
    # dims : a vector of the three dimensions to be plotted, e.g. dims=c(1,31,51)
    # u : threshold for selecting extreme observations (based on the u quantile of the radii of X[dims])
  # outputs:
    # a 3D scatter plot
  
  # 3D scatter plot of extreme observations
  r <- sqrt(rowSums(X[,dims]^2))
  rstar <- quantile(r, u, na.rm = TRUE)
  bigr <- which(r > rstar)
  scatter3D(x = X[bigr,dims[1]] / r[bigr],
            y = X[bigr,dims[2]] / r[bigr],
            z = X[bigr,dims[3]] / r[bigr],
            col = "black", # change this to be indicator of radius?
            box=TRUE, bty="b2", axes=TRUE,
            xlim = c(0,1), ylim = c(0,1), zlim = c(0,1), 
            theta=135, phi=30, 
            pch=20, cex = 0.3,
            label=TRUE, xlab = paste0("X", dims[1]), ylab = paste0("X", dims[2]), zlab = paste0("X", dims[3]))
  # add unit sphere surface
  M <- mesh(seq(0, pi/2, length.out = 50), seq(0, pi/2, length.out = 50))
  phi <- M$x 
  theta <- M$y
  x <- cos(phi) * sin(theta)
  y <- sin(phi) * sin(theta)
  z <- cos(theta)
  surf3D(x, y, z, col = "black", alpha=0.1, add = TRUE)
}

levelplotmatrix <- function(M, xlab = "i", ylab = "j", main = NULL){
  
  # description:
    # make a levelplot for the entries of a matrix 
  # inputs:
    # M : matrix
    # xlab, ylab, main : axis labels and title
  # outputs:
    # a levelplot with colour key
  
  par(pty="s")
  levelplot(M, 
            cuts = 8, col.regions = brewer.pal(9, "Reds"),
            main = list(main, line=0.7),
            xlab = list(xlab, cex=0.8),
            ylab = list(ylab, cex=0.8),
            scales = list(x=list(cex=0.8), y=list(cex=0.8)),
            colorkey = list(labels=list(cex=0.8)))
}

plotfrevent <- function(t, Xtilde, qvals = c(2,5,10,20,45,92), PCA = PCAv, long = fr$long, lat = fr$lat){
  
  # description:
    # make exploratory extremal PCA plots
  # inputs:
    # t : time index of the event being explored
    # Xtilde : the transformed data
    # qvals : number of eigenvectors to use in the reconstructions
    # PCA : output of extremePCA function
    # long : longitude coordinates
    # lat : latitude coordinates
  # outputs:
    # a series of plots
  
  # setup grid layout of plots
  layout.matrix <- matrix(c(1,2,2,3,4,5,6,7,8,9,10,11), nrow = 4, ncol = 3, byrow = TRUE)
  layout(mat = layout.matrix,
         heights = c(1, 1, 1, 1), # Heights of the rows
         widths = c(1, 1, 1)) # Widths of the columns
  # plot map of the event
  col_loc <- brewer.pal(5, "Blues")
  temp <- seq(-1e-4, max(Xtilde[t,])+1e-4, length.out = 6)
  # prepare map
  map("worldHires", "France", col = "grey90", fill = T, xlim = c(-5.0,9.5), ylim = c(41.5,51), mar = c(2.5,2.5,2.5,2.5))
  title(bquote(tilde(bold(x))[.(t)]), line = 0.8)
  for (k in 1:92){
    temp_k <- sum(Xtilde[t,k] > temp)
    points(long[k], lat[k], pch=20, cex=0.8, col=col_loc[temp_k])
  }
  colkey(side = 1, col = col_loc, clim = c(min(temp), max(temp)), breaks = round(temp, 1), dist = -0.1, cex.axis = 0.8, mgp = c(3, 0.3, 0), add = T)
  # plot the values of the PCs
  # find min and max of V1(t),...,Vd(t) - used for y axis labels
  at_Vt <- c(min(PCA$V[t,]), max(0, min(PCA$V[t,])), max(PCA$V[t,]))
  # order evecs by strength of signal at time t
  evecs <- order(-abs(PCA$V[t,]))
  # set graphical parameters
  par(mar = c(2.1, 2, 2, 2), mgp = c(1.2, 0.3, 0))
  # plot ith PC and mark the largest observed events
  plot(PCA$V[t,], type = "h", yaxt = "n", xlab = "PC index, i", ylab="", cex.axis = 0.7, cex.lab = 0.7)
  title(ylab = bquote(V[.(t)][","][i]), cex.lab = 1)
  axis(side = 2, at = at_Vt, labels = round(at_Vt, digits = 0), cex.axis = 0.7)
  abline(h = 0, col = "black", lty = 1, lwd = 0.5)
  points(x = evecs[2:4], y = PCA$V[t,evecs[2:4]], pch = 20, cex = 0.5, col = "red")
  points(x = evecs[-(2:4)], y = PCA$V[t,evecs[-(2:4)]], pch = 20, cex = 0.5, col = "black")
  # plot eigenvectors
  col_loc <- c(brewer.pal(3,"Reds")[3:2], "white", brewer.pal(3,"Blues")[2:3])
  for (i in 2:4){
    v <- evecs[i]
    z <- max(abs(PCA$U[,v]))
    temp <- seq(-z-1e-4, z+1e-4, length.out = 6)
    # prepare map
    map("worldHires", "France", col = "grey85", fill = T, xlim = c(-5.0,9.5), ylim = c(41.5,51), mar = c(2.5,2.5,2.5,2.5))
    title(bquote(hat(bold(u))[.(v)] ~ "(" * V[.(t)][","][.(v)] ~ "=" ~ .(round(PCA$V[t,v],1)) * ")"), line = 0.8)
    for (k in 1:92){
      temp_k <- sum(PCAv$U[k,v] > temp)
      points(fr$long[k], fr$lat[k], pch=20, cex=0.8, col=col_loc[temp_k])
    }
    colkey(side = 1, col = col_loc, 
           clim = c(min(temp), max(temp)), breaks = round(temp, 2),
           dist = -0.1, cex.axis = 0.8, mgp = c(3, 0.3, 0), add = T)
  }
  # plot partial basis reconstructions
  col_loc <- brewer.pal(5, "Blues")
  for (q in qvals){
    reconXt <- reconstructX(t = t, q = q, U = PCA$U, V = PCA$V)
    temp <- seq(0, max(reconXt)+1e-4, length.out = 6)
    # prepare map
    map("worldHires", "France", col = "grey90", fill = T, xlim = c(-5.0,9.5), ylim = c(41.5,51), mar = c(2.5,2.5,2.5,2.5))
    title(paste(ifelse(q==92, "All", q), "eigenvectors"), line = 0.7)
    for (k in 1:92){
      temp_k <- sum(reconXt[k] > temp)
      points(long[k], lat[k], pch=20, cex=0.8, col=col_loc[temp_k])
    }
    colkey(side = 1, col = col_loc, clim = c(min(temp), max(temp)), breaks = round(temp, 1), dist = -0.1, cex.axis = 0.8, mgp = c(3, 0.3, 0), add = T)
  }
}