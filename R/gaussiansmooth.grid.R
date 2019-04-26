
gaussiansmooth.grid<- function(x, sd, max.r=3*sd, kernel.dim=21, max.p.na=1, sample.points, weights,
                                no.match=0, use.old=FALSE, na.value = NA, ...	){
 
  # This function performs a gaussiansmooth on a grid object. 
  #   For efficiency sake if the sd is much larger than the cellsize the grid is first upscaled then smoothed and then downscaled
  #	 
  #
  # Arguments 
  #	x : grid object to be smoothed
  #	sd: the standard deviation (in map units) of the gaussian kernel
  #     max.r : the maximum radius (in map units) of the guassian kernel
  #     kernel.dim : the target dimension of the kernel set to 0 to avoid upscaling 
  #	  otherwise x will be upscaled such that a kernel of kernel.dim is sufficient to reach max.r units out.
  #	  higher values result in more precise output at the expense of longer computation times.
  #	  max.p.na : is passed onto upscale (if upscaling is necessary)
  #  	  sample.points : a matrix or data.frame with x and y columns which specify points to sample
  #	  weights : a dataframe with columns "value" and "weight" the first containing grid values and the second containing the associated weight. If these columns are missing and the dataframe contains two columns it will be assumed that the first is "value" and second "weight"
  # Value : A grid of the same extent and cellsize as x which has been smoothed.
  
  # Setup for sampling points
  if(!missing(sample.points)){
    if(!all(c("x", "y") %in% names(sample.points)))
      stop("sample.points must have columns named \"x\" and \"y\"")
    sp <- data.frame(r=y2r(sample.points$y, x), c=x2c(sample.points$x, x))
    sv <- sp$r > 0 & sp$r <= x$nrow & sp$c > 0 & sp$c <= x$ncol  # selection vector for points in range
    sampled.values <- numeric(length(sv))
    sampled.values[] <- NA
    sample=TRUE
  } else {
    sample=FALSE
  }
  
  if(!missing(weights)){
    if(all(c("value", "weight") %in% colnames(weights))){
      x <- swap(x, weights$value, weights$weight, no.match=no.match, na.value = na.value)
    } else if(ncol(weights) == 2) {
      x <- swap(x, weights[1,], weights[,2], no.match=no.match, na.value = na.value)
    } else {
      stop("If weights are supplied and the dataframe has more than two columns it must have columns labled \"value\" and \"weight\".")
    }
  }
  
  
  current.window.dim <- max.r*2/x$cellsize # in cells
  o.grid.info <- x[names(x)!= "m"] # save the metadata from the orignal grid (drop the matrix)
  upscaled= FALSE # flag to keep track of upscaling
  if(kernel.dim != 0 & current.window.dim > kernel.dim){
    factor <- max(ceiling(current.window.dim/kernel.dim), 1)  # factor can't be less than 1
    print(paste("Upscaling to cells of ", factor*x$cellsize, "mapunits."))
    x <- upscale(x, factor=factor, max.p.na=max.p.na, use.old=use.old)
    # Expand the upscaled version of x by one cell
    x2 <- list(m=matrix(NA, x$nrow+2, x$ncol+2), nrow=x$nrow+2, ncol=x$ncol+2, xll=x$xll-x$cellsize, yll=x$yll-x$cellsize, cellsize=x$cellsize)
    x2$m[2:(x$nrow+1), 2:(x$ncol+1)] <- x$m
    class(x2) <- c("grid", class(x2))
    x <- x2
    rm(x2)
    gc()
    upscaled  <- TRUE 
  }
  k <- makegaussiankernel(sd=sd, max.r=max.r, cellsize=x$cellsize)
  print("Smoothing")
  x$m <- kernelsmooth(x$m, k, use.old = use.old)
  gc()
  if(upscaled){    
    print("Downscaling")
    print(paste("Factor =", factor))
    x<- downscale(x, factor=factor)
    gc()
    x <- matchextent(x, o.grid.info)
    gc()
  }

  if(sample){
    sampled.values[sv] <- x$m[as.matrix(sp[sv,])]
    return(list(grid=x, sample=sampled.values))
  }
  return(x)
}


