calckernel.grid <- function(obj, kernel, x, y, row=y2r(y, list=obj), col=x2c(x, list=obj), ...){
#-------------------------------------------------------
#  calckernel.grid  -  Ethan Plunkett - Dec 2, 2009
#
#	Arguments:
#    obj - an object of class ("grid")
#	  kernel - the kernel to apply - should be a square, odd dimensioned matrix that sums to 1.
#    x, y - the x and y coordinates (in map units) of the focal cell
#    row, col - the row and column of the focal cell.
#  Details: 
#    This function calculates the sum of a kernel applied to a single cell of a grid object
#	  The focal cell can be specified in either 
#	Revised Feb 1 2010 to work with "grid" class objects and to accept either x.coord and y.coord, or row and col
#-------------------------------------------------------
if(row < 0 || row > nrow(obj$m)) stop("The specified row or y coordinate is not on the grid.  row = ", row)
if(col < 0 || col > ncol(obj$m)) stop("The specified col or x coordinate is not on the grid.  col = ", col)
calckernel(obj$m, kernel, row, col, ...) 
}

