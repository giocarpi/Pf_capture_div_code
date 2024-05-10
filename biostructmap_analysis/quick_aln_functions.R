# ------------------------------ #
# quick functions for seeing if two sequences are the same
# and for aligning two sequences
# ------------------------------ #


library(tidyverse)
library(msa)


check_seq = function(a, b){
  seq1 = gsub('\n|\t|\r| ', '', a)
  seq2 = gsub('\n|\t|\r| ', '', b)
  r = seq1 == seq2
  print(r)
}

quick_aln = function(a, b){
  seq1 = gsub('\n|\t|\r| ', '', a)
  seq2 = gsub('\n|\t|\r| ', '', b)
  x = msaMuscle(c(seq1,seq2), type = 'protein', order = 'input') %>% suppressWarnings()
  y = msaConsensusSequence(x)
  y = gsub('\\?', ' ', y)
  names(y) = 'check'
  x = x %>% as.character()
  print(x)
  print(y)
  
  # split both and print where different
  s1 = strsplit(x[1],'')[[1]]
  s2 = strsplit(x[2],'')[[1]]
  
  miss = which(s1 != s2)
  
  print(miss)
}


