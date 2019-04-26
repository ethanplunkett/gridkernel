#include <Rcpp.h>
using namespace Rcpp;

// upscalec implements matrix upscaling 
// it is called by the r function upscale


// [[Rcpp::export]]
Rcpp::NumericMatrix upscalec(Rcpp::NumericMatrix x, int factor, double maxPNA ){
  
  int xrow = x.nrow();
  int xcol = x.ncol();
  int urow = (int) ceil( (double) xrow / factor);
  int ucol = (int) ceil( (double) xcol / factor);
  
  
  // Rcout << "urow = " << urow << " " << ceil(xrow / factor) << std::endl;
  // Rcout << "ucol = " << ucol << std::endl;
   
  // NumericMatrix res = x;
  NumericMatrix res(urow, ucol);
  
//  upscale.cellsize <- x$cellsize*factor

  int maxNA = factor * factor * maxPNA;
  for(int r=0; r < urow; r++){  // loop through rows of the coarse grid
    for(int c=0; c < ucol; c++){ // loop through columns of the course grid
      
      // Define the start and end row in the fine grid that correponds to 
      //  the current cell of the coarse grid
      int startRow = r * factor;
      int endRow  = std::min( (r+1)*factor - 1, xrow-1);  
      int startCol = c * factor;
      int endCol = std::min( (c+1) *factor - 1, xcol-1);
      
//      Rcout << "r = " << r << std::endl;
//      Rcout << "c = " << c << std::endl;
//      Rcout << "startRow = " << startRow << std::endl;
//      Rcout << "endRow = " << endRow << std::endl;
//      Rcout << "startCol = " << startCol << std::endl;
//      Rcout << "endCol = " << endCol << std::endl;
      double usedSum = 0;
      int usedCount = 0;
     
      for(int rf = startRow; rf <= endRow; rf++){
        for(int cf = startCol; cf <= endCol; cf++){
          double a = x(rf, cf);
          if(!R_IsNA(a)){
            usedSum += a;
            usedCount ++;
          }
        } // end loop through fine columns (cf)
      }  // end lop through fine rows (rf)
      int naCount = (factor * factor) - usedCount;
      double cellValue;
      if(usedCount == 0 || naCount > maxNA){
        cellValue = NA_REAL;
      } else {
        cellValue = usedSum / (double) usedCount;
      }
      res(r, c) = cellValue;
        
    } // end col loop
  } // end row loop
  
  return res;
}


