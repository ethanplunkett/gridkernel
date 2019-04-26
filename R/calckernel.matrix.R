calckernel.matrix <- function(obj, kernel, row, col, use.old=FALSE, ...)  {
#-------------------------------------------------------
#	calckernel.matrix written by Ethan Plunkett Jan 2007
# 
# 	returns the value for that cell obtained by summing the product of the kernel and the matrix it overlaps.
#  
#  The function assumes the sum of the kernel is 1; and results will be incosistent near edges and NA's
#	if this assumption is violated. If this assumption is met than the result is a weighted average of the 
#	area around the target cell.  
#	Near edges and NA values the kernel is renormalized so the sum of the used cells is 1. 
#
#	The function could easily be rewritten to always divide by the sum
#  of the kernel but it's computationally more efficient to normalize the kernel once rather than
#  every time it's applied.  Use:  kernel <- kernel/sum(kernel) to normalize it.
#
#	Parameters:
#  	obj = a matrix with cells to which we want to apply the kernel
#  	kernel = the kernel we wish to apply
#  row = the row of the focal map cell
#  	col = the column of focal map cell
#
#
#  Modified slightly Feb 1, 2010 :
#		Made calckernel a generic function and converted this to a matrix method for that function
#-------------------------------------------------------
	kernel.r <- (dim(kernel)[1]-1)/2
	row.min <- row-kernel.r
	row.max <- row+kernel.r
	column.min <- col-kernel.r
	column.max <- col+kernel.r
	
	if( row < 1 || col < 1 || row > dim(obj)[1] || col > dim(obj)[2]){
		stop('Focal Cell isn\'t in matrix.')
	}

  if(!use.old)
    return(calckernelc(x=obj, k=kernel, row=row, col=col))

	# if the kernel hangs over the edge of the map
	# then take a cell by cell approach 
	if((row.min <1)||(column.min < 1) || (row.max > dim(obj)[1]) || column.max > dim(obj)[2]) {
		# warning(paste("Kernel centered at [", row, ",", col,"] extends beyond edge of matrix.", sep="" ))
		used.kernel.sum<-0
		weighted.sum<-0
		# for each row and column of the kernel:
		for(r in 1:dim(kernel)[1]) for( c in 1:dim(kernel)[2]){  
			# the correspoinding part of the matrix:
			matrix.r <- r-1+row.min
			matrix.c <- c-1+column.min
			# if that cell of the kernel does overlap the map
			if(matrix.r > 0 & matrix.c > 0 & matrix.r <= dim(obj)[1] & matrix.c <= dim(obj)[2]){
				if(is.na(obj[matrix.r, matrix.c])) next  # don't tally NA cells
				weighted.sum <- weighted.sum + kernel[r,c]*obj[matrix.r, matrix.c]
				used.kernel.sum <- used.kernel.sum + kernel[r,c]
			} #end if (the cell overlays)	
		} # end for each row and column
		
		mean <- weighted.sum/used.kernel.sum
		return(mean)
	} # end if the kernel overhangs
	
	# if the whole kernel can be used take a patch approach
	patch<- obj[(row-kernel.r):(row+kernel.r),(col-kernel.r):(col+kernel.r)] 
	mean<-sum(patch*kernel, na.rm=TRUE)
	if(any(is.na(patch)))
		mean <- mean/sum(kernel[!is.na(patch)])
	return(mean)
} # end function definition
