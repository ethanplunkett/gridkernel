kernelscalingdetails <- function( max.r, cellsize, tilesize = 2000, kernel.dim = 21, method = 3){
  # This function determines the upscaling factor, the realized tilesize, and the required buffer
  #  used by gaussiansmooth and kernelsmooth when coarsening grids to increase efficiency
  # Arguments:
  #  max.r this is the radius of the kernel in the units of the projection (mapunits)
  #    it should be set to the same value as in the call to the smoothing function.
  #    For gaussiansmooth the max.r defaults to 3 times the sd.
  #  cellsize : the resultion of the grid being smoothed
  #  tilesize : the target tilesize as set in the call to the smoothing functions.
  #    It defaults to 2000 both here and in the smoothing functions.
  #  kernel.dim : should be set the same as in the smoothing function, this
  #    is the target kerenl dimension (width and height) in cells in the upscaled
  #    grid.  Essentially this determines how coarse to make the gridl a larger kernel.dim
  #    means less upscaling (because more cells will be used to represent the kernel in the coarse grid)
  #    Here and in the smoothing function it defaults to 21 which means that the grid will be coarsened
  #    such that it takes 21 cells to represent a kernel with a radius of max.r.
  #  method: set the same as in the smoothing functions. The method affects the buffer as
  #   method 2 requires that the buffer be set big enough to encompass the max.r while
  #   method 1 reads the whole grid at once so buffer isn't relevant and method 3 requires
  #    no buffer.
  #
  #  This is used internally in kernelsmooth and gaussiansmoothto determine
  # how much upscaling (and thus approximation) to do.
  # It is provided to users as it is necessary to anticipate the realized
  #  tilesize  correctly to use the tiles arguments to those functions
  #
  #
  current.window.dim <- max.r*2/cellsize # in cells
  factor <- max(1, ceiling(current.window.dim/kernel.dim))
  # round tilesize up to the nearest (larger) multiple of factor
  tilesize <- ceiling(tilesize/factor)*factor
  if(method == 2){
    buffer <- ceiling(max.r/cellsize/factor)*factor+factor*2
  } else {
    buffer <-  0
  }
  return(list(factor = factor, tilesize = tilesize, buffer = buffer))
}




