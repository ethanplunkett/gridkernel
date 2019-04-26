upscale <- function(x, factor, max.p.na = 0.2, use.old = FALSE){
  # Function to reduce the number of pixels in a grid object.
  # Arguments :
  #  x: a grid object
  #  factor: the number of rows and columns to use in x to make each new row in the result
  #  (esentialy factor is the amount the horizontal and vertical resolution will decrease)
  # max.p.na : As cells in x are averaged to make new cells NA's are dropped unless the proportion of NA's in a window exceeds this 
  # threshold in which case the new cell will have an NA value.
  
  
  
  upscale.cellsize <- x$cellsize*factor
  
  res <- matrix(NA, ceiling(nrow(x$m)/factor), ceiling(ncol(x$m)/factor))
  max.na <- max.p.na * factor * factor
  if(!use.old){
    res <- upscalec(x=x$m, factor=factor, maxPNA=max.p.na)
  } else {
    for(r in 1:nrow(res))for(c in 1:ncol(res)){
      row.index <- (1 + (r-1) * factor):min((r*factor),nrow(x$m))
      col.index <- (1 + (c-1) * factor):min(c*factor, ncol(x$m))
      if(sum(is.na(x$m[row.index, col.index])) > max.na){
        res[r,c] <- NA
      } else {
        res[r, c] <- mean(x$m[row.index,col.index], na.rm=TRUE)
      }
    }  
  } # end use old
  res[is.nan(res)] <- NA
  
  
  
  res <- list(m=res, nrow=nrow(res), ncol=ncol(res), xll=x$xll, yll=x$yll+(x$cellsize * x$nrow) - nrow(res)*upscale.cellsize, cellsize=upscale.cellsize)
  class(res) <- c("grid", class(res))
  return(res)
}

