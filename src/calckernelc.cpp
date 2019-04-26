#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double calckernelc(Rcpp::NumericMatrix x, Rcpp::NumericMatrix k, int row, int col){
  
   row = row - 1;  // switch to C++ indexing
   col = col - 1;
   int radius = (k.nrow() - 1)/2;
   int rmin = row - radius;
   int cmin = col - radius;
   double usedSum = 0;
   double weightSum =0;
   int krow = k.nrow();
   int kcol = k.ncol();
   int xrow = x.nrow();
   int xcol = x.ncol();
   for(int r = 0; r < krow; r++){
     int xr = r + rmin;
      if (xr < 0 || xr >= xrow)
          continue;
     for(int c = 0;  c < kcol; c++){
        int xc = c + cmin;
        if (xc < 0 || xc >= xcol )
          continue;
        double a = k(r, c);
        double b = a * x(xr, xc);      
        if(R_IsNA(b))
          continue;
        usedSum += a;
        weightSum += b;
        
     }
   }
   return weightSum/usedSum;
}


/*** R
k <- matrix(1, 3, 3)
k[2, 2] <- 3
k <- k/sum(k)

l <- matrix(2, 100, 100)


calckernelc(l, k, 50, 50)
calckernel(l, k, 50, 50)
rm(k, l)
*/

//
//  kernel.r <- (dim(kernel)[1]-1)/2
//	row.min <- row-kernel.r
//	row.max <- row+kernel.r
//	column.min <- col-kernel.r
//	column.max <- col+kernel.r
//	
//	if( row < 1 || col < 1 || row > dim(obj)[1] || col > dim(obj)[2]){
//		stop('Focal Cell isn\'t in matrix.')
//		}
//	
//	# if the kernel hangs over the edge of the map
//	# then take a cell by cell approach 
//	if((row.min <1)||(column.min < 1) || (row.max > dim(obj)[1]) || column.max > dim(obj)[2]) {
//		# warning(paste("Kernel centered at [", row, ",", col,"] extends beyond edge of matrix.", sep="" ))
//		used.kernel.sum<-0
//		weighted.sum<-0
//		# for each row and column of the kernel:
//		for(r in 1:dim(kernel)[1]) for( c in 1:dim(kernel)[2]){  
//			# the correspoinding part of the matrix:
//			matrix.r <- r-1+row.min
//			matrix.c <- c-1+column.min
//			# if that cell of the kernel does overlap the map
//			if(matrix.r > 0 & matrix.c > 0 & matrix.r <= dim(obj)[1] & matrix.c <= dim(obj)[2]){
//				if(is.na(obj[matrix.r, matrix.c])) next  # don't tally NA cells
//				weighted.sum <- weighted.sum + kernel[r,c]*obj[matrix.r, matrix.c]
//				used.kernel.sum <- used.kernel.sum + kernel[r,c]
//			} #end if (the cell overlays)	
//		} # end for each row and column
//		
//		mean <- weighted.sum/used.kernel.sum
//		return(mean)
//	} # end if the kernel overhangs
//	
//	# if the whole kernel can be used take a patch approach
//	patch<- obj[(row-kernel.r):(row+kernel.r),(col-kernel.r):(col+kernel.r)] 
//	mean<-sum(patch*kernel, na.rm=TRUE)
//	if(any(is.na(patch)))
//		mean <- mean/sum(kernel[!is.na(patch)])
//	return(mean)