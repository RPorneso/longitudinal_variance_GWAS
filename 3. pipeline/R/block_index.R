block_index = function(ind, nblocks, size, maxcol) {
  
  if (ind > nblocks) stop("Blockindex cannot be greater than total number of blocks")
  first = 1
  if (ind > 1) {
    first = first + (ind - 1) * size
  }
  last = first + size - 1
  if (ind == nblocks) {
    last = first + maxcol %% size - 1  ## changed last to first and added -1 
  }
  c(first, last)
}