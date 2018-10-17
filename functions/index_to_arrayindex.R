index_to_arrayindex <- function(index, arraysize, arrayindex_to_index = FALSE)
{
  # < CDD for FMRI >
  #   
  #   Copyright (C) 2018 Namgil Lee & Jong-Min Kim
  
  N = length(arraysize)  # number of dimensions of array
  
  cum_prod_sizes = cumprod(c(1,arraysize[-N]))  #(1, J{1}, J{1}J{2}, ..., J{1}***J{N-1})
  
  if (!arrayindex_to_index) {
    #***** Convert index to arrayindex ************#
    out = rep(0,N)
    rindex = index-1  #iterate with residuals
    for (n in N:1) {
      out[n] = rindex %/% cum_prod_sizes[n] + 1
      rindex = rindex %% cum_prod_sizes[n]
    }
  } else {
    #***** Convert arrayindex to index ************#
    # The  input  'index' is a vector of indices
    #out === SUM{ (j{n}-1)*J{1}J{2}*...*J{n-1} } + 1
    out = sum((index-1) *  cum_prod_sizes) + 1 
  }
  
  return(out)
  
}