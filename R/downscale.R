downscale <- function(x, factor, match.extent=TRUE){
  # Arguments:
  #  x - a grid object
  #  factor - How many cells should the vertical and horizonal dimensions of 
  #    each original cell be sliced into.  A factor of 3 means each cell will
  #    be converted into 9 new cells.
  # Currently if the factor is odd then one row and one column is lost from the
  #   bottom and right of the final mantrix that could in theory be calculated.

  if( !isTRUE(all.equal(x$cellsize %% factor, 0))) 
    warning("Using a factor that doesn't result in an whole numbered cellsize is not recommended")
  
  even.factor <- isTRUE(all.equal(factor %/% 2,0))  # If the factor is even then
    # the large cell centers fall in the center of small cells;  otherwise they
    # fall on a vertex between small cells. Some calculations change.
  				
  # Identify the coordinates of the large cells centers in cell units of the small celled window
  if(even.factor) {
    corner.coords <- rbind(c(1,1), c( factor +1, 1), c(1, factor+1), c(factor+1, factor+1))
  } else {
    corner.coords <- rbind(c(0.5,0.5), c( factor +0.5, 0.5 ), c(0.5, factor+0.5), c(factor+0.5, factor+0.5))  ## Changed
  } 

  # Make weight matrix where each row is a cell in the downscaled window and
  # each column is the value of one of the four cells in the upscaled matrix. 
  # Matrix positions are converted to row number in the standard r way starting
  # in the top left and then progressing down the first column before jumping to
  # the top of the second column
  weights <- matrix(NA, factor * factor, 4)
  for(i in 1:nrow(weights)){
    r <- i %% factor
    if(r == 0) r <- factor
    c <- ((i - 1) %/% factor)+1     	
    distances <- abs(cbind(r-corner.coords[,1], c-corner.coords[,2]))
    w = 1 - distances/(factor)
    w<- w[,1]* w[,2]
    weights[i, ] <- w
  }  	 

  # Make an empty map to fill in with values
  if(match.extent){
     # Make empty map of the same extent as the original but higher pixel resolution
     # (This will result in a band of NA pixels around the edge
     nr <- x$nrow * factor
     nc <- x$ncol * factor
     xll <- x$xll
     yll <- x$yll 
     res <- list(m= matrix(NA,nr , nc,), nrow=nr, ncol=nc, xll=xll, yll=yll, cellsize=x$cellsize/factor)
     class(res) <- c("grid", class(res))
     offset <-  floor(factor/2) # number that needs to be added to 1 to get the first cell we will be writing too.
  } else {
    # Make empty map of just the cells that can be interpolated properly
    if(even.factor){
       nr <- round((x$nrow-1)*factor) ###
       nc <- round((x$ncol-1)*factor) ###
    }else{
      nr <- round((x$nrow-1)*factor+1)
      nc <- round((x$ncol-1)*factor+1)
    }
    xll <- x$xll + floor(.5*x$cellsize)
    yll <- x$yll + floor(.5*x$cellsize)
    res <- list(m= matrix(NA,nr , nc,), nrow=nr, ncol=nc, xll=xll, yll=yll, cellsize=x$cellsize/factor)
    class(res) <- c("grid", class(res))
    offset <- 0
  }

  # Do the interpolation as a matrix multiplication on each set of four corners in the lower res matrix
  for(r in 1:(x$nrow-1))for(c in 1:(x$ncol-1)){
  	corners <- as.numeric(x$m[r:(r+1), c:(c+1)]) # the values of the for corners of the window	
  	window <- matrix(weights %*% matrix(corners, length(corners), 1), factor, factor) # the downscaled window of pixels
	rows <- ( ((r-1)*factor+1):(r*factor) ) + offset
	cols <- ( ((c-1)*factor+1):(c*factor) ) + offset
  	res$m[rows, cols] <- window	
  }
  
  return(res)
}

