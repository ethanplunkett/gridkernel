kernelsmooth.matrix <- function(x, kernel, use.old=FALSE, ...){
  #-------------------------------------------------------
  #	kernelsmooth.matrix  written by Ethan Plunkett Jan 2007
  #
  #	Applies a kernel to each cell in the matrix
  #	by repeatedly calling calckernel
  #	Returns a matrix in which each cell contains the sum of a kernel applied to the cells around that cell in the original matrix.
  #
  #	Parameters:
  #  	obj = a matrix with cells to which we want to apply the kernel
  #  	kernel = the kernel we wish to apply
  #-------------------------------------------------------
  if(!use.old){
    result <- kernelsmoothc(x=x, k=kernel)
  } else  {
    result <- matrix(NA, dim(x)[1], dim(x)[2])
    for( r in 1:dim(x)[1]) for( c in 1:dim(x)[2]) {
      result[r,c]<- calckernel(x,kernel, r, c, use.old=TRUE)
    } # end for each row and column
  }
  result[is.nan(result)] <- NA
  return(result)
} # end function definition
