At some point models were using considerably more memory than before, due to a memory leak (I suspect) 
involving calling R from SLiM to do some matrix transformations. I ended up cutting that and using a
different approach.