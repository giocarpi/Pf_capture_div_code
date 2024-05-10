# ------------------------------------------------------- #
#   Brad Broyles, Purdue University, He Lab - 7/7/23      #
#                                                         #
#   Nucleotide diversity is too low to display in pymol   #
#   I will scale up by taking log10(pi) for the sake of   #
#   vizualizing spectrum in pymol                         #
#                                                         #
#   Saving updated pdb involved extracting and rewriting  #
#     b-factor column                                     #
# ------------------------------------------------------- #

# load libraries
library(tidyverse)
library(bio3d)

# set wd
setwd('tajima_d_models/')
dir()

# update B-factor column for these pi pdbs (remove pfs230_seg1 and AMA1 for now)
p_names = c('csp','rh5','pfs25','pfs47',
            'pfs4845', 'pfs28', 'celtos', 'msp119')

p_names1 = paste(p_names, '_norsa_trim_pi.csv',sep='')
p_names2 = paste(p_names, '_norsa_trim_pi.pdb',sep='')

# loop through files -- load both csv and pdb
# use csv b values to fill in csv -- take log10 of b to get scaled b

for(num in 1:length(p_names1)){
  pdb = read.pdb(p_names2[num])
  csv = read_csv(p_names1[num])
  
  # trim off hetatm
  pdb = bio3d::trim.pdb(pdb, atom.select(pdb, type = 'ATOM'))
  
  # check lengths
  print(unique(pdb$atom$resno) %>% length() ==
          unique(csv$reference_residue) %>% length()
  )
  
  # scale up csv beta factors
  # if any scores go to -Inf set them as 0
  csv$score = as.numeric(csv$score)
  
  # scores of 0 wont be colored in plot
  # add a small psuedocount (1e-5) --- updated 2/14
  #csv$score = csv$score + 1e-5
  csv$score = log10(csv$score)
  csv$score = ifelse(csv$score == -Inf, 0, csv$score)
  
  # need residue map
  x = tibble(csv = unique(csv$reference_residue),
             pdb = unique(pdb$atom$resno))
  
  # replace b factor column in pdb
  for(i in 1:nrow(x)){
    ro = which(csv$reference_residue == x$csv[i])
    val = csv$score[ro]
    
    # if val is NA, I will set it to 1 since these are not displayed
    val = ifelse(is.na(val), 1, val)
    
    ro = which(pdb$atom$resno == x$pdb[i])
    pdb$atom$b[ro] = val
  }
  
  # save updated pdb
  fname = p_names2[num]
  fname = gsub('_pi', '_pi_adj', fname)
  
  # print these minimums and maximums
  print('minimum / maximum')
  print(min(pdb$atom$b))
  print(max(pdb$atom$b))
  
  write.pdb(pdb, file = fname)
}

# do ama1 and pfs230 ----
p_names = c('ama1_2q8a', 'ama1_full', 'ama1_no4g2_short', 'ama1_no4g2', 'ama1_short',
            'pfs230.1', 'pfs230.3', 'pfs230.5', 'pfs230.7', 'pfs230.9', 'pfs230.11', 'pfs230.13')

p_names1 = paste(p_names, '_norsa_trim_pi.csv',sep='')
p_names2 = paste(p_names, '_norsa_trim_pi.pdb',sep='')

for(num in 1:length(p_names1)){
  pdb = read.pdb(p_names2[num])
  csv = read_csv(p_names1[num])
  
  # trim off hetatm
  pdb = bio3d::trim.pdb(pdb, atom.select(pdb, type = 'ATOM'))
  
  # check lengths
  print(unique(pdb$atom$resno) %>% length() ==
          unique(csv$reference_residue) %>% length()
  )
  
  # scale up csv beta factors
  # if any scores go to -Inf set them as 0
  csv$score = as.numeric(csv$score)
  
  # scores of 0 wont be colored in plot
  # add a small psuedocount (1e-5) --- updated 2/14
  #csv$score = csv$score + 1e-5
  csv$score = log10(csv$score)
  csv$score = ifelse(csv$score == -Inf, 0, csv$score)
  
  # need residue map
  x = tibble(csv = unique(csv$reference_residue),
             pdb = unique(pdb$atom$resno))
  
  # replace b factor column in pdb
  for(i in 1:nrow(x)){
    ro = which(csv$reference_residue == x$csv[i])
    val = csv$score[ro]
    
    # if val is NA, I will set it to 1 since these are not displayed
    val = ifelse(is.na(val), 1, val)
    
    ro = which(pdb$atom$resno == x$pdb[i])
    pdb$atom$b[ro] = val
  }
  
  # save updated pdb
  fname = p_names2[num]
  fname = gsub('_pi', '_pi_adj', fname)
  
  # print these minimums and maximums
  print('minimum / maximum')
  print(min(pdb$atom$b))
  print(max(pdb$atom$b))
  
  write.pdb(pdb, file = fname)
}


