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

# update B-factor column for these pi pdbs
p_names = c('csp','ama1','rh5','pfs25',
            'pfs47_long','pfs47_short',
            'pfs4845','pfs230_seg1')

pfs230_names = c('pfs230_seg1','pfs230_seg2','pfs230_seg3','pfs230_seg4',
            'pfs230_seg5','pfs230_seg6','pfs230_seg7','pfs230_seg8',
            'pfs230_seg9','pfs230_seg10','pfs230_seg11','pfs230_seg12',
            'pfs230_seg13')

pfs230_names = paste('pfs230/',pfs230_names,sep='')

pnames = c(p_names, pfs230_names)

p_names1 = paste(p_names, '_trim_pi.csv',sep='')
p_names2 = paste(p_names, '_trim_pi.pdb',sep='')

## loop through p_names - load csv and pdb - scale up pi taking log10(pi) - save updated pdb ----
for(num in 1:length(p_names2)){
# read in pdb
pdb = read.pdb(p_names2[num])

# this will serve as a map between resno and line number
hi = pdb$atom %>% select(resno, b)

range1 = unique(hi$resno)

# read in csv
csv = read.csv(p_names1[num])

# if sequences have gaps in solved residues I needed this piece
csv$pdb_res = range1

for(i in 1:nrow(csv)){
  ro = which(hi$resno == csv$pdb_res[i])
  sc = csv$score[i]
  hi$b[ro] = sc
}

## round b column to max 7 positions
hi$b = as.numeric(hi$b)
hi$b = ifelse(hi$b==0,NA,hi$b)
hi$b = log10(hi$b)
hi$b =  hi$b %>% round(2)
hi$b = ifelse(is.na(hi$b),0,hi$b)

test = read_lines(p_names2[num])

# 63-77 captures b factor and spacing
for(i in 1:nrow(hi)){
capt1 = str_sub(test[i],1,62)
capt2 = str_sub(test[i],78,-1)
b = hi$b[i]
insert1 = format(b, scientific = F)
len = 15 - nchar(insert1)
insert1 = paste(insert1, paste(rep(' ',len),collapse=''),sep = '')
test[i] = paste(capt1, insert1, capt2,sep='')
}

fname = gsub('.pdb','_adj.pdb',p_names2[num])
cat(test,sep='\n',file = fname)
}


