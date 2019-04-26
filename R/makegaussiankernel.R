makegaussiankernel <- function( sd, max.r=3*sd, cellsize=1, include.center.cell=TRUE) {
  #-------------------------------------------------------
  #  makegaussiankernel ~ function written by Ethan Plunkett Jan 2006
  #
  # This function creates create a gaussian kernel in which the sum of all cells is 1.
  #
  #	Parameters:
  #
  #	sd is the standard deviation of the distribution. If sd is set to zero the
  #		function will return a uniform kernel (all cells within the radius set to the same value).
  # 	max.r is the maximum radius of the kernel.  All cells beyond that distance
  # 		have a value of 0. This also determines the final dimensions of the kernel.
  # 	cellsize is the size of the cell in whatever units you care to use it is
  # 		assumed that cellsize, sd, and r are all in the same units.
  # 	include.center if FALSE the center cell has a value
  # 		of 0.  This is useful for my partiuclar application but probably to no one else.
  #  Modified Jan 2018 to eliminate inefficient looping over cells
  #-------------------------------------------------------

  max.r=max.r/cellsize  #convert radius to cell units
  sd=sd/cellsize
  size = ceiling(max.r)*2+1
  center = ceiling(max.r)+1
  kernel <- new('matrix', 0, size, size)
  rows <- row(kernel)
  cols <- col(kernel)
  dists <- sqrt( (rows - center)^2 + (cols - center)^2 )

  if (sd!=0){  #make a gaussian kernel
    kernel[ , ] <- dnorm(dists, mean=0, sd, log=FALSE)
    kernel[  dists > max.r] <- 0
   }
  else {		#sd == 0,  make uniform kernel
    kernel[ , ] <- as.numeric(dists <= max.r)
  }
  if (include.center.cell==FALSE) {
    kernel[center,center]<-0
  }
  #normalize so kernel sums to 1
  kernel <- kernel/sum(kernel)

  # Trim outer cells if they are unecessary
  if(all(kernel[1,]==0, kernel[,1]==0, kernel[nrow(kernel), ]==0, kernel[, ncol(kernel)]==0))
    kernel <- kernel[2:(nrow(kernel)-1), 2:(ncol(kernel)-1)]
  return(kernel)
}
