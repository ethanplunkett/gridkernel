# gridkernel

This package performs kernel calculations on raster data.  It can perform
Gaussian and Pareto smoothing, as well as apply user defined kernels.

It uses edge and NA correction to eliminate edge effects, and  it can use 
upscaling to increase performence with very large bandwidths which would
otherwise be computationally expensive.

A typical workflow would be to read a raster file with the raster package, 
convert to a grid object in memory with as.grid(),  run either gaussiansmooth()
or kernelsmooth(),  and then convert back  to a raster with raster().  

## Installation

This package imports rasterprocess so inherits that package's dependence on 
windows. 

Use the code below to install gridprocess
 [gridprocess](https://github.com/ethanplunkett/gridprocess) and [gridkernel](https://github.com/ethanplunkett/gridkernel)
 
``` r
devtools::install_github("ethanplunkett/gridprocess")
devtools::install_github("ethanplunkett/gridkernel")

```
  
