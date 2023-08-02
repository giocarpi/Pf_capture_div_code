# ------------------------------------------------------- #
#   Brad Broyles, Purdue University, He Lab - 7/7/23      #
#                                                         #
#   Two jobs:  1. Trim fasta files to match PDB segments  #
#              2. Trim and reset PDB indices              #
#                                                         #
#   Notes on segmentation of Pfs230 --                    #
#     InterPro scan reveals 14 6C domains of Pfs230       #
#     For folding purposes we break pfs230 into 13        #
#       overlapping segments. Each segment captures two   #
#       6C domains                                        #
# ------------------------------------------------------- #

# load libraries
library(tidyverse)
library(bio3d)      # for parsing pdbs
library(phylotools) # for parsing fasta

## process csp -- 309 - 373 ----
hi = phylotools::read.fasta('field_sample_fastas/csp_w_ref.fasta')

pdb = read.pdb('solved_csp/7rxp.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))

mi = 309
ma = 373
range = mi:ma

mi_nuc = mi * 3 - 2
ma_nuc = ma * 3

hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

hi$seq.name[1]
hi$trim[1]
#cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
#    sep = '\n', file = 'Final fasta/csp_trimmed.fasta')

# reset pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

mod$atom$resno = mod$atom$resno - (mi - 1)

#write.pdb(mod, file = 'Final fasta/csp_7rxp_cleaned.pdb')

## process AMA1 -- 107 - 438 ----
hi = phylotools::read.fasta('field_sample_fastas/ama1_w_ref.fasta')

pdb = read.pdb('solved_ama1/2q8a.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))

mi = 107
ma = 438
range = mi:ma

mi_nuc = mi * 3 - 2
ma_nuc = ma * 3

hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

hi$seq.name[1]
hi$trim[1] %>% nchar()
#cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
#    sep = '\n', file = 'Final fasta/ama1_trimmed.fasta')

# reset pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

mod$atom$resno = mod$atom$resno - (mi - 1)

#write.pdb(mod, file = 'Final fasta/ama1_2q8a_cleaned.pdb')

## process rh5 -- 146 - 513 ----
hi = phylotools::read.fasta('field_sample_fastas/rh5_w_ref.fasta')

pdb = read.pdb('Final fasta/Rh5_4wat_cleaned.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))

mi = 146
ma = 513

range = mi:ma

mi_nuc = mi * 3 - 2
ma_nuc = ma * 3

hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

hi$seq.name[1]
hi$trim[1] %>% nchar()
#cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
#    sep = '\n', file = 'Final fasta/rh5_trimmed.fasta')

# reset pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

mod$atom$resno = mod$atom$resno - (mi - 1)

#which(!1:368 %in% unique(mod$atom$resno))

write.pdb(mod, file = 'Final fasta/Rh5_4wat_cleaned.pdb') #-- no reset needed

## process pfs25 -- 23 - 193 ----
hi = phylotools::read.fasta('field_sample_fastas/pfs25_w_ref.fasta')

pdb = read.pdb('af_models_six_recycles/Pfs25_aa_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))

mi = 23
ma = 193

mi_nuc = mi * 3 - 2
ma_nuc = ma * 3

hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

hi$seq.name[1]
hi$trim[1] %>% nchar()
#cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
#    sep = '\n', file = 'Final fasta/pfs25_trimmed.fasta')

# reset pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

mod$atom$resno = mod$atom$resno - (mi - 1)

#write.pdb(mod, file = 'Final fasta/pfs25_AF_cleaned.pdb')

## process pfs4845 -- 45 - 428 ----
hi = phylotools::read.fasta('field_sample_fastas/pfs4845_w_ref.fasta')

pdb = read.pdb('af_models_six_recycles/Pfs45-48_aa_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))

mi = 45
ma = 428

mi_nuc = mi * 3 - 2
ma_nuc = ma * 3

hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

hi$seq.name[1]
hi$trim[1] %>% nchar()
#cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
#    sep = '\n', file = 'Final fasta/pfs4845_trimmed.fasta')

# reset pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

mod$atom$resno = mod$atom$resno - (mi - 1)
#chainA$atom$resno = chainA$atom$resno - 106

#write.pdb(mod, file = 'Final fasta/pfs4845_AF_cleaned.pdb')


## process pfs47 -- 31 - 415 (short) ----
hi = phylotools::read.fasta('field_sample_fastas/pfs47_w_ref.fasta')

pdb = read.pdb('Final fasta/Pfs47_aa_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))

## two versions short and long
mi = 31
ma = 415

mi_nuc = mi * 3 - 2
ma_nuc = ma * 3

hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

hi$seq.name[1]
hi$trim[1] %>% nchar()
#cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
#    sep = '\n', file = 'Final fasta/pfs47_trimmed_short.fasta')

# reset pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

mod$atom$resno = mod$atom$resno - 30

#write.pdb(mod, file = 'Final fasta/pfs47_AF_cleaned_short.pdb')

## process pfs47 -- 22 - 419 (long) ----
hi = phylotools::read.fasta('field_sample_fastas/pfs47_w_ref.fasta')

pdb = read.pdb('Final fasta/Pfs47_aa_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))

## two versions short and long
mi = 22
ma = 419

mi_nuc = mi * 3 - 2
ma_nuc = ma * 3

hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

hi$seq.name[1]
hi$trim[1] %>% nchar()
#cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
#    sep = '\n', file = 'Final fasta/pfs47_trimmed_long.fasta')

# reset pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

mod$atom$resno = mod$atom$resno - 21

#write.pdb(mod, file = 'Final fasta/pfs47_AF_cleaned_long.pdb')

## process pfs230 -- 13 segments ----
trib = tribble(~seg,~mi,~ma,
  1, 552, 890,   
  2, 730, 1136,
  3, 915, 1278,
  4, 1133, 1435,
  5, 1282, 1563,
  6, 1432, 1910,
  7, 1691, 2038,
  8, 1907, 2202,
  9, 2049, 2377,
  10, 2201, 2666,
  11, 2445, 2830,
  12, 2663, 2982,
  13, 2828, 3116)

## load pfs230 fasta
hi = phylotools::read.fasta('field_sample_fastas/pfs230_w_ref.fasta')

## loop through each segmentation of pfs230 (domain doublets)
for(i in 1:nrow(trib)){
mi = trib$mi[i]
ma = trib$ma[i]

mi_nuc = mi * 3 - 2
ma_nuc = ma * 3

hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

fname = paste('Final fasta/pfs230_seg',i,'_trimmed.fasta',sep='')
cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = fname)

#update pdb
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

mod$atom$resno = mod$atom$resno - mi + 1
fname2 = paste('Final fasta/pfs230_segment',i,'_AF_cleaned.pdb',sep='')
write.pdb(mod, file = fname2)
}


