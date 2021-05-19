#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
Rcpp::NumericMatrix kernelsmoothc(Rcpp::NumericMatrix x, Rcpp::NumericMatrix k){
   
   
   int krow = k.nrow();
   int kcol = k.ncol();
   int xrow = x.nrow();
   int xcol = x.ncol();
   int radius = (k.nrow() - 1)/2;
   
   // NumericMatrix res = x;
   NumericMatrix res(xrow, xcol);
   
   for(int row=0; row < xrow; row++){
     for(int col=0; col < xcol; col++){
       if (col % 100 == 0)
         Rcpp::checkUserInterrupt();
       int rmin = row - radius;
       int cmin = col - radius;
       double usedSum = 0;
       double weightSum =0;
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
        
        }  // end c loop (through kernel cols)
     }  // end r loop (through kernel rows)
    res(row, col) = weightSum/usedSum;
 } // end row loop  (through x rows)
 } // end col loop  (through x cols)
 return res;
}