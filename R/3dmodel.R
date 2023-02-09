# Load the svgViewR package
library(svgViewR)

# Set number of points to draw
n <- 5

# Create a cloud of normally distributed 3D points
points3d <- cbind(rnorm(n, sd=1), rnorm(n, sd=1), rnorm(n, sd=1))

# Open a connection to .html file
svg.new(file='plot_static_points.html')

# Get distance of points from the center of point cloud
pdist <- sqrt(rowSums((points3d - matrix(colMeans(points3d), n, 3, byrow=TRUE))^2))

# Set color gradient between red and blue
colfunc <- colorRampPalette(c('red', 'blue'))

# Set desired number of colors along gradient
col_grad <- colfunc(50)

# Scale distances to indices and find corresponding color along gradient
col <- col_grad[(length(col_grad)-1)*(pdist - min(pdist)) / diff(range(pdist))+1]

# Add points to file
svg.points(points3d, col=col)

# Add coordinate axis planes around the points
svg_frame <- svg.frame(points3d)

# Close the file connection
svg.close()

#==================




# Set number of iterations
n_iter <- 100

# Create animated point array
points3da <- array(points3d, dim=c(dim(points3d), n_iter))

# Expand points from origin
for(iter in 0:(n_iter-1)){
  points3da[, , iter] <- points3da[, , iter] * 0.001 * iter^2
}






# move along direction
n <- 100
n_iter <- 2000

points3d <- matrix(runif(3*n, min=-1, max=1),nrow = n,ncol = 3)
points3da <- array(points3d, dim=c(dim(points3d), n_iter))
dir <- matrix(runif(3*n, min=-1, max=1),nrow = n,ncol = 3)
for(iter in 2:(n_iter)){
  if (iter%%10 == 0) {
    dir <- matrix(runif(3*n, min=-1, max=1),nrow = n,ncol = 3)
  }
  points3da[, , iter] <- points3da[, , iter-1] + dir
}

# Open a connection to .html file
svg.new(file='plot_animated_points.html')
svg.points(points3da, col=col)
svg_frame <- svg.frame(points3da)
svg.close()








#---------- border control
points3d <- matrix(runif(3*n, min=-1, max=1),nrow = n,ncol = 3)
points3da <- array(points3d, dim=c(dim(points3d), n_iter))
dir <- matrix(runif(3*n, min=-1, max=1),nrow = n,ncol = 3)
for(iter in 2:(n_iter)){
  if (iter%%20 == 0) {

    dir <- matrix(runif(3*n, min=-1, max=1),nrow = n,ncol = 3)
    #for each ind, get dir and border correct
    for(ind in 1:n){
      pos <- points3da[, , iter-1][ind,]
      dir[ind,] <- border.correct(pos, dir[ind,], 100)
    }
  }
  points3da[, , iter] <- points3da[, , iter-1] + dir
}

# Open a connection to .html file
svg.new(file='plot_animated_points.html')
svg.points(points3da, col=col)
svg_frame <- svg.frame(points3da)
svg.close()




#----- ATrrtaction
# this time each individual tumbles at different rates

#---------- border control
points3d <- matrix(runif(3*n, min=-1, max=1),nrow = n,ncol = 3)
points3da <- array(points3d, dim=c(dim(points3d), n_iter))
dir <- matrix(runif(3*n, min=-1, max=1),nrow = n,ncol = 3)
for(iter in 2:(n_iter)){
  if (iter%%20 == 0) {

    #dir <- matrix(runif(3*n, min=-1, max=1),nrow = n,ncol = 3)
    #for each ind, get dir and border correct
    for(ind in 1:n){
      pos <- points3da[, , iter-1][ind,]
      if(pos[1]>100){
        #tumble
        #dir[ind,] <- runif(3, min=-0.1, max=0.1)
      }else{
        dir[ind,] <- runif(3, min=-1, max=1)
      }

      dir[ind,] <- border.correct(pos, dir[ind,], 100)
    }
  }
  points3da[, , iter] <- points3da[, , iter-1] + dir
}

# Open a connection to .html file
svg.new(file='plot_animated_points.html')
svg.points(points3da, col=col)
svg_frame <- svg.frame(points3da)
svg.close()








#pos <- c(-100,0,0)
#dir <- c(1,1,1)
#lim<- 10
border.correct <- function(pos,dir,lim){
  x.min.d = -1
  x.max.d = 1
  y.min.d = -1
  y.max.d = 1
  z.min.d = -1
  z.max.d = 1
  x.l <- FALSE
  y.l <-FALSE
  z.l <- FALSE
  min.lim <- -lim
  max.lim <- lim
  if(pos[1]< min.lim){
    x.min.d <- 0
    x.l <- TRUE
  }

  if(pos[1]>max.lim){
    x.max.d <- 0
    x.l <- TRUE
  }
  if(pos[2]< min.lim){
    y.min.d <- 0
    y.l <- TRUE
  }

  if(pos[2]>max.lim){
    y.max.d <- 0
    y.l <- TRUE
  }

  if(pos[3]<min.lim){
    z.min.d <- 0
    z.l <- TRUE
  }

  if(pos[3]>max.lim){
    z.max.d <- 0
    z.l <- TRUE
  }


  out <- c(ifelse(x.l==TRUE, runif(1, min=x.min.d, max=x.max.d),dir[1]),
           ifelse(y.l==TRUE, runif(1, min=y.min.d, max=y.max.d),dir[2]),
           ifelse(z.l==TRUE, runif(1, min=z.min.d, max=z.max.d),dir[3]))
  return(out)
}




x.min.d






# boxing , limits border
lmts <- cbind(min=c(-50,-50,-50),max=c(50,50,50))





# assume attraction
# attractor at 50,50,50
# frequency of tumbling depends on concentation of an attractant







#single.dir <- runif(3, min=-1, max=1)









# Open a connection to .html file
svg.new(file='plot_animated_points.html')

# Add points to file
svg.points(points3da, col=col)

# Add coordinate axis planes around the points
svg_frame <- svg.frame(points3da)

# Close the file connection
svg.close()







