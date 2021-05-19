# gridkernel

This package performs kernel calculations on raster data.  It can perform
Gaussian and Pareto smoothing, as well as apply user defined kernels.

It uses edge and NA correction to eliminate edge effects, and it can use 
upscaling to increase performence with very large bandwidths which would
otherwise be computationally expensive.

A typical workflow would be to read a raster file with the raster package, 
convert to a grid object in memory with as.grid(),  run either gaussiansmooth()
or kernelsmooth(),  and then convert back  to a raster with raster().  

## Installation


Use the code below to install gridprocess
 [gridprocess](https://github.com/ethanplunkett/gridprocess) and [gridkernel](https://github.com/ethanplunkett/gridkernel)
 
``` r
devtools::install_github("ethanplunkett/gridprocess")
devtools::install_github("ethanplunkett/gridkernel")

```

## Changelog

May 18, 2022 (0.1.1) Now "suggests" gridprocess rather than depends. As of today
I'm dropping the functions in gridkernel from gridio so that gridio now depends
on gridkernel.  Practically speaking it makes sesne to use this package either 
with gridio or with gridprocess and I don't want to import gridprocess when 
I'm using it with gridio.  Also fixed a bug in the kernel creation functions;
they used to work fine with R 3.4 but stopped with R 4.0.
  
