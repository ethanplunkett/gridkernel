'kernsource' <- function(x, bandwidth, cellsize = 30, search = 3, buffer = 0) {
# kernsource.r
# Source-based standard (non-resistant) Gaussian kernel estimator
# Arguments:
#	x			source grid, with >=0 values for source points.  Kernels are weighted
#				by values in x; use 1's for unweighted kernels.
#	bandwidth	kernel bandwidth (h), in meters
#	cellsize	cell size, in meters
#	search		search distance, in s.d.  search = 3 gives 99% of kernel; use values as low
#				as 2 for quick and sloppy results
#	buffer		edge buffer, in meters.  Source cells in buffer are ignored (presumably these
#				are dealt with in adjacent blocks).  Buffer must be either 0 or >= bandwidth 
#				* search / cellsize.
# Notes:
# 	- This version does not do edge corrections.  Results will be biased downward at
# 	hard data edges.  
# 	- Source data should be buffered by bandwidth * search / cellsize; write results in summation mode.
# B. Compton, 28 Apr 2011.  From kern.r, 16 Jan 2008; translated from APL function KERN, 25 Sep 2003.



	if(!(buffer == 0 | buffer >= (bandwidth * search / cellsize))) 
		stop('Error in call to kernsource: buffer too small.')
	
	v <- ceiling(bandwidth * search / cellsize)						# search distance, in cells
	e <- buffer
	x <- as.matrix(x)
	if(buffer == 0) {												# if buffer not supplied, 
		x <- rbind(matrix(0,v,dim(x)[2]),x,matrix(0,v,dim(x)[2]))	# 	expand to include edges around points
		x <- cbind(matrix(0,dim(x)[1],v),x,matrix(0,dim(x)[1],v))
		e <- v
	}
	z <- matrix(0, dim(x)[1], dim(x)[2])
		
	b <- 2 * v + 1
	b <- (ceiling(b / 2) - 1:(b)) ^ 2
	m <- outer(b, b, FUN='+') ^ .5
	m <- m * cellsize / bandwidth
	q <- dnorm(m)											# bivariate kernel
	q <- q * (q >= q[ceiling(dim(q)[1] / 2), 1])			# avoid square artifacts
	q <- q / (sum(q))										# standardize
	
	for(i in (e + 1):(dim(z)[1] - e)) {						# for each row,
		for(j in (e + 1):(dim(z)[2] - e)) {					#	for each column,
			if(x[i,j] == 0) next
			t <- (1:dim(q)[1]) - ceiling(dim(q)[1] / 2)		#	Make window
			k <- t + i
			l <- t + j
			z[k,l] <- z[k,l] + x[i,j] * q					#	Add kernel
		}
	}
	if(buffer == 0) {
		z <- z[(e + 1):(dim(z)[1] - e), (e + 1):(dim(z)[2] - e)]# drop edges
	}
	z
}
