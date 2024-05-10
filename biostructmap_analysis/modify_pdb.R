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
#                                                         #
#   Make sure to filter out HETATM as well                #
# ------------------------------------------------------- #

# load libraries
library(tidyverse)
library(bio3d)      # for parsing pdbs
library(phylotools) # for parsing fasta
library(seqinr)
library(msa)


## process csp -- 309 - 373 ----
hi = phylotools::read.fasta('fasta_regions/PF3D7_0304600_CSP_n829_1194_C-terminal.fasta')

pdb = read.pdb('solved_csp/7rxp.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

fa_aa = s2c(hi$seq.text[1]) %>% seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[2], '[^-*]+')

mi = 33
ma = 97
range = mi:ma

mi_nuc = mi * 3 - 2
ma_nuc = ma * 3

hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

hi$seq.name[1]
pdb_aa
hi$trim[1] %>% s2c() %>% seqinr::translate() %>% paste(collapse = '')

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/csp_trimmed.fasta')

# check range
chainA$atom$resno %>% unique()

# extract structure
write.pdb(chainA, file = 'biostructmap_input/csp_7rxp_cleaned.pdb')

## process AMA1 -- 100 - 545 ----
hi = phylotools::read.fasta('fasta_regions/PF3D7_1133400_AMA1_n298_1635_D1+D2+D3.fasta')

pdb = read.pdb('solved_ama1/AF-C5HW96-F1-model_v4.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^-*]+')

mi = 87
ma = 532
range = mi:ma

#mi_nuc = mi * 3 - 2
#ma_nuc = ma * 3

#hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

#hi$seq.name[1]
#pdb_aa
#hi$trim[1] %>% s2c() %>%seqinr::translate() %>% paste(collapse = '')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/ama1_trimmed_full.fasta')

# check range
chainA$atom$resno %>% unique()

# trim pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

# check range
mod$atom$resno %>% unique()

# reset pdb indices
mod$atom$resno = mod$atom$resno + 13

write.pdb(mod, file = 'biostructmap_input/ama1_AF_cleaned_full.pdb')

# supply a couple variants of this structure 
var = mod

# trim off 4g2 loop (351 -- 390)
inds = atom.select(var, resno = c(100:350, 391:545))
var = trim.pdb(var, inds=inds)

# check range
var$atom$resno %>% unique()

write.pdb(var, file = 'biostructmap_input/ama1_AF_cleaned_no4g2.pdb')

# trim C term a bit (532 onwards)
var = mod
inds = atom.select(var, resno = c(100:350, 391:532))
var = trim.pdb(var, inds=inds)

# check range
var$atom$resno %>% unique()

write.pdb(var, file = 'biostructmap_input/ama1_AF_cleaned_no4g2_short.pdb')

# trim C term a bit (532 onwards)
var = mod
inds = atom.select(var, resno = c(100:532))
var = trim.pdb(var, inds=inds)

# check range
var$atom$resno %>% unique()

write.pdb(var, file = 'biostructmap_input/ama1_AF_cleaned_short.pdb')





## process AMA1 off 2q8a -- 100 - 532 ----
hi = phylotools::read.fasta('fasta_regions/PF3D7_1133400_AMA1_n298_1635_D1+D2+D3.fasta')

pdb = read.pdb('solved_ama1/2q8a.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[2], '[^-*]+')

mi = 8
ma = 339
range = mi:ma

mi_nuc = mi * 3 - 2
ma_nuc = ma * 3

hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

hi$seq.name[1]
pdb_aa
hi$trim[1] %>% s2c() %>%seqinr::translate() %>% paste(collapse = '')

#hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/ama1_trimmed_2q8a.fasta')

# check range
chainA$atom$resno %>% unique()

# trim pdb indices
#chainA$atom
#inds = atom.select(chainA, resno = mi:ma)
#mod = trim.pdb(chainA, inds=inds)

# check range
#mod$atom$resno %>% unique()

# reset pdb indices
#mod$atom$resno = mod$atom$resno + 13

write.pdb(chainA, file = 'biostructmap_input/ama1_2q8a_cleaned.pdb')







## process rh5 -- 146 - 513 ----
hi = phylotools::read.fasta('fasta_regions/PF3D7_0424100_RH5_n436_1539_structure.fasta')

pdb = read.pdb('solved_rh5/4wat.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[2], '[^-*]+')

#mi = 146
#ma = 513

#range = mi:ma

#mi_nuc = mi * 3 - 2
#ma_nuc = ma * 3

#hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

#hi$seq.name[1]
#pdb_aa
#hi$trim[1] %>% s2c() %>%seqinr::translate() %>% paste(collapse = '')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/rh5_trimmed.fasta')

# check range
chainA$atom$resno %>% unique()

# reset pdb indices

write.pdb(chainA, file = 'biostructmap_input/rh5_4wat_cleaned.pdb')

## process pfs25 -- 23 - 193 ----
hi = phylotools::read.fasta('fasta_regions/PF3D7_1031000_Pfs25_n67_579_structure.fasta')

pdb = read.pdb('af_models_six_recycles/Pfs25_aa_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^-*]+')

mi = 23
ma = 193

#mi_nuc = mi * 3 - 2
#ma_nuc = ma * 3

#hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

#hi$seq.name[1]
#pdb_aa
#hi$trim[1] %>% s2c() %>%seqinr::translate() %>% paste(collapse = '')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/pfs25_trimmed.fasta')

# trim pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

# check range
mod$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno - (mi - 1)

write.pdb(mod, file = 'biostructmap_input/pfs25_AF_cleaned.pdb')

## process pfs4845 -- 45 - 428 ----
hi = phylotools::read.fasta('fasta_regions/PF3D7_1346700_Pfs4548_n133_1284_structure.fasta')

pdb = read.pdb('af_models_six_recycles/Pfs45-48_aa_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^-*]+')

mi = 45
ma = 428

#mi_nuc = mi * 3 - 2
#ma_nuc = ma * 3

#hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

#hi$seq.name[1]
#pdb_aa
#hi$trim[1] %>% s2c() %>%seqinr::translate() %>% paste(collapse = '')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/pfs4845_trimmed.fasta')

# reset pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

# check range
mod$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno - (mi - 1)

write.pdb(mod, file = 'biostructmap_input/pfs4845_AF_cleaned.pdb')


## process pfs47 -- 31 - 415 (short) ----
hi = phylotools::read.fasta('fasta_regions/PF3D7_1346800_Pfs47_n64_1257_structure.fasta')

pdb = read.pdb('af_models_six_recycles/Pfs47_aa_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^*-]+')

## two versions short and long
#mi = 31
#ma = 415

#mi_nuc = mi * 3 - 2
#ma_nuc = ma * 3

#hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

#hi$seq.name[1]
#pdb_aa
#hi$trim[1] %>% s2c() %>%seqinr::translate() %>% paste(collapse = '')

#cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
#    sep = '\n', file = 'biostructmap_input/pfs47_trimmed_short.fasta')

# reset pdb indices
#chainA$atom
#inds = atom.select(chainA, resno = mi:ma)
#mod = trim.pdb(chainA, inds=inds)

# check range
#mod$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno - 30

#write.pdb(mod, file = 'biostructmap_input/pfs47_AF_cleaned_short.pdb')

## process pfs47 -- 22 - 419 (long) ----
hi = phylotools::read.fasta('fasta_regions/PF3D7_1346800_Pfs47_n64_1257_structure.fasta')

pdb = read.pdb('af_models_six_recycles/Pfs47_aa_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^*-]+')

## two versions short and long
mi = 22
ma = 419

#mi_nuc = mi * 3 - 2
#ma_nuc = ma * 3

#hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

#hi$seq.name[1]
#pdb_aa
#hi$trim[1] %>% s2c() %>%seqinr::translate() %>% paste(collapse = '')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/pfs47_trimmed_long.fasta')

# reset pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

# check range
mod$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno - 21

write.pdb(mod, file = 'biostructmap_input/pfs47_AF_cleaned_long.pdb')

## process pfs230 -- 13 segments ---- # now doing one by one ----
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

## load pfs230.1 -----
hi = phylotools::read.fasta('fasta_regions/PF3D7_0209000_Pfs230_n1654_2670_D1+D2.fasta')

pdb = read.pdb('af_models_six_recycles/pfs230 AF models/pfs230_segment1_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

# fasta needs cleaned up
# where is gaps
x = str_locate_all(hi$seq.text, '[^ATGC]+')

tib = tibble(i = 1:length(x), start = 0, stop = 0)

for(i in 1:length(x)){
  if(nrow(x[[i]]) < 1){
    print(i)
    print('no gaps')
    print('------')
    next
  }
  
  if(nrow(x[[i]]) > 1){
    print(i)
    print('two gaps')
    print('------')
    next
  }
  
  tib$start = x[[i]][1,1]
  tib$stop = x[[i]][1,2]
}

hi$seq.text[c(1, 18, 93)]

# drop 18 and 93
hi = hi[-c(18, 93),]
hi$seq.text = gsub('-', '', hi$seq.text)

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^*-]+')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/pfs230_seg1.fasta')

# reset pdb indices
chainA$atom$resno = chainA$atom$resno + 551
chainA$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno - 21
write.pdb(chainA, file = 'biostructmap_input/pfs230_seg1.pdb')



#files = list.files('af_models_six_recycles/pfs230 AF models/', full.names = T)

## loop through each segmentation of pfs230 (domain doublets)
#for(i in 1:nrow(trib)){
#  print(i)
#mi = trib$mi[i]
#ma = trib$ma[i]

#mi_nuc = mi * 3 - 2
#ma_nuc = ma * 3

#hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

#fname = paste('biostructmap_input/pfs230_seg',i,'_trimmed.fasta',sep='')
#cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
#    sep = '\n', file = fname)

# load pdb file
#a = grep(paste('segment',i, '_',sep=''), files)
#b = a[grep('.pdb', files[a])]
#c = b[grep('rank_001', files[b])]

#pdb = read.pdb(files[c])
#chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
#chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

#update pdb
#inds = atom.select(chainA)#, resno = mi:ma)
#mod = trim.pdb(chainA, inds=inds)

#mod$atom$resno = mod$atom$resno + (mi - 1)

#print(mod$atom$resno %>% unique())

#fname2 = paste('biostructmap_input/pfs230_segment',i,'_AF_cleaned.pdb',sep='')
#write.pdb(mod, file = fname2)
#}



## load pfs230.3 -----
hi = phylotools::read.fasta('fasta_regions/PF3D7_0209000_Pfs230_n2743_3834_D3+D4.fasta')

pdb = read.pdb('af_models_six_recycles/pfs230 AF models/pfs230_segment3_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

# fasta needs cleaned up
# where is gaps
x = str_locate_all(hi$seq.text, '[^ATGC]+')
x %>% unique()

tib = tibble(i = 1:length(x), start = 0, stop = 0)

for(i in 1:length(x)){
  if(nrow(x[[i]]) < 1){
    print(i)
    print('no gaps')
    print('------')
    next
  }
  
  if(nrow(x[[i]]) > 1){
    print(i)
    print('two gaps')
    print('------')
    next
  }
  
  tib$start[i] = x[[i]][1,1]
  tib$stop[i] = x[[i]][1,2]
}

tib %>% filter(start > 0)

# drop 18 and 93
hi = hi[-c(tib$i[which(tib$start>0)]),]
hi$seq.text = gsub('-', '', hi$seq.text)

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^*-]+')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/pfs230_seg3.fasta')

# reset pdb indices
chainA$atom$resno = chainA$atom$resno + 914
chainA$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno - 21
write.pdb(chainA, file = 'biostructmap_input/pfs230_seg3.pdb')

## load pfs230.5 -----
hi = phylotools::read.fasta('fasta_regions/PF3D7_0209000_Pfs230_n3844_4689_D5+D6.fasta')

pdb = read.pdb('af_models_six_recycles/pfs230 AF models/pfs230_segment5_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

# fasta needs cleaned up
# where is gaps
x = str_locate_all(hi$seq.text, '[^ATGC]+')
x %>% unique()

tib = tibble(i = 1:length(x), start = 0, stop = 0)

for(i in 1:length(x)){
  if(nrow(x[[i]]) < 1){
    print(i)
    print('no gaps')
    print('------')
    next
  }
  
  if(nrow(x[[i]]) > 1){
    print(i)
    print('two gaps')
    print('------')
    next
  }
  
  tib$start[i] = x[[i]][1,1]
  tib$stop[i] = x[[i]][1,2]
}

tib %>% filter(start > 0)

# drop 18 and 93
#hi = hi[-c(tib$i[which(tib$start>0)]),]
#hi$seq.text = gsub('-', '', hi$seq.text)

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^*-]+')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/pfs230_seg5.fasta')

# reset pdb indices
chainA$atom$resno = chainA$atom$resno + 1281
chainA$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno - 21
write.pdb(chainA, file = 'biostructmap_input/pfs230_seg5.pdb')

## load pfs230.7 -----
hi = phylotools::read.fasta('fasta_regions/PF3D7_0209000_Pfs230_n5071_6114_D7+D8.fasta')

pdb = read.pdb('af_models_six_recycles/pfs230 AF models/pfs230_segment7_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

# fasta needs cleaned up
# where is gaps
x = str_locate_all(hi$seq.text, '[^ATGC]+')
x %>% unique()

tib = tibble(i = 1:length(x), start = 0, stop = 0)

for(i in 1:length(x)){
  if(nrow(x[[i]]) < 1){
    print(i)
    print('no gaps')
    print('------')
    next
  }
  
  if(nrow(x[[i]]) > 1){
    print(i)
    print('two gaps')
    print('------')
    next
  }
  
  tib$start[i] = x[[i]][1,1]
  tib$stop[i] = x[[i]][1,2]
}

tib %>% filter(start > 0)

# drop 18 and 93
#hi = hi[-c(tib$i[which(tib$start>0)]),]
#hi$seq.text = gsub('-', '', hi$seq.text)

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^*-]+')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/pfs230_seg7.fasta')

# reset pdb indices
chainA$atom$resno = chainA$atom$resno + 1690
chainA$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno - 21
write.pdb(chainA, file = 'biostructmap_input/pfs230_seg7.pdb')

## load pfs230.9 -----
hi = phylotools::read.fasta('fasta_regions/PF3D7_0209000_Pfs230_n6145_7131_D9+D10.fasta')

pdb = read.pdb('af_models_six_recycles/pfs230 AF models/pfs230_segment9_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

# fasta needs cleaned up
# where is gaps
x = str_locate_all(hi$seq.text, '[^ATGC]+')
x %>% unique()

tib = tibble(i = 1:length(x), start = 0, stop = 0)

for(i in 1:length(x)){
  if(nrow(x[[i]]) < 1){
    print(i)
    print('no gaps')
    print('------')
    next
  }
  
  if(nrow(x[[i]]) > 1){
    print(i)
    print('two gaps')
    print('------')
    next
  }
  
  tib$start[i] = x[[i]][1,1]
  tib$stop[i] = x[[i]][1,2]
}

tib %>% filter(start == 0)

# drop 18 and 93
hi = hi[-c(tib$i[which(tib$start==0)]),]
hi$seq.text = gsub('-', '', hi$seq.text)

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^*-]+')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/pfs230_seg9.fasta')

# reset pdb indices
chainA$atom$resno = chainA$atom$resno + 2048
chainA$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno - 21
write.pdb(chainA, file = 'biostructmap_input/pfs230_seg9.pdb')

## load pfs230.11 -----
hi = phylotools::read.fasta('fasta_regions/PF3D7_0209000_Pfs230_n7333_8490_D11+D12.fasta')

pdb = read.pdb('af_models_six_recycles/pfs230 AF models/pfs230_segment11_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

# fasta needs cleaned up
# where is gaps
x = str_locate_all(hi$seq.text, '[^ATGC]+')
x %>% unique()

tib = tibble(i = 1:length(x), start = 0, stop = 0)

for(i in 1:length(x)){
  if(nrow(x[[i]]) < 1){
    print(i)
    print('no gaps')
    print('------')
    next
  }
  
  if(nrow(x[[i]]) > 1){
    print(i)
    print('two gaps')
    print('------')
    next
  }
  
  tib$start[i] = x[[i]][1,1]
  tib$stop[i] = x[[i]][1,2]
}

tib %>% filter(start == 0)

# drop 18 and 93
hi = hi[-c(tib$i[which(tib$start==0)]),]
hi$seq.text = gsub('-', '', hi$seq.text)

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^*-]+')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/pfs230_seg11.fasta')

# reset pdb indices
chainA$atom$resno = chainA$atom$resno + 2444
chainA$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno - 21
write.pdb(chainA, file = 'biostructmap_input/pfs230_seg11.pdb')

## load pfs230.13 -----
hi = phylotools::read.fasta('fasta_regions/PF3D7_0209000_Pfs230_n8482_9348_D13+D14.fasta')

pdb = read.pdb('af_models_six_recycles/pfs230 AF models/pfs230_segment13_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

# fasta needs cleaned up
# where is gaps
x = str_locate_all(hi$seq.text, '[^ATGC]+')
x %>% unique()

tib = tibble(i = 1:length(x), start = 0, stop = 0)

for(i in 1:length(x)){
  if(nrow(x[[i]]) < 1){
    print(i)
    print('no gaps')
    print('------')
    next
  }
  
  if(nrow(x[[i]]) > 1){
    print(i)
    print('two gaps')
    print('------')
    next
  }
  
  tib$start[i] = x[[i]][1,1]
  tib$stop[i] = x[[i]][1,2]
}

tib %>% filter(start == 0)

# drop 18 and 93
#hi = hi[-c(tib$i[which(tib$start==0)]),]
#hi$seq.text = gsub('-', '', hi$seq.text)

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^*-]+')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/pfs230_seg13.fasta')

# reset pdb indices
chainA$atom$resno = chainA$atom$resno + 2827
chainA$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno - 21
write.pdb(chainA, file = 'biostructmap_input/pfs230_seg13.pdb')

## process pfs28 -- 24 - 178 ----
hi = phylotools::read.fasta('fasta_regions/PF3D7_1030900_Pfs28_n70_534_structure.fasta')

pdb = read.pdb('solved_pfs28/8e1z.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^-*]+')

# need to trim off 2 aa residues (don't match)
mi = 3
ma = 157

#mi_nuc = mi * 3 - 2
#ma_nuc = ma * 3

#hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

#hi$seq.name[1]
#pdb_aa
#hi$trim[1] %>% s2c() %>%seqinr::translate() %>% paste(collapse = '')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/pfs28_trimmed.fasta')

# trim pdb indices
inds = atom.select(chainA, resno = 1:max(mod$atom$resno))
mod = trim.pdb(chainA, inds=inds)

mod$atom$resno = mod$atom$resno + 23

# check range
mod$atom$resno %>% unique()

write.pdb(mod, file = 'biostructmap_input/pfs28_8e1z_cleaned.pdb')

## process msp119 -- 1608 - 1699 ----
hi = phylotools::read.fasta('fasta_regions/PF3D7_0930300_MSP1_n4822_5097_D19.fasta')

pdb = read.pdb('new_AF_models/msp119/msp119_AF_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[2], '[^-*]+')

#mi = 1608
#ma = 1699

#mi_nuc = mi * 3 - 2
#ma_nuc = ma * 3

#hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

#hi$seq.name[1]
#pdb_aa
#hi$trim[1] %>% s2c() %>%seqinr::translate() %>% paste(collapse = '')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/msp119_trimmed.fasta')

# trim pdb indices
chainA$atom
inds = atom.select(chainA)
mod = trim.pdb(chainA, inds=inds)

mod$atom$resno = mod$atom$resno + 1607

# check range
mod$atom$resno %>% unique()

write.pdb(mod, file = 'biostructmap_input/msp119_AF_cleaned.pdb')

## process celtos 35 - 164 ---- 
hi = phylotools::read.fasta('fasta_regions/PF3D7_1216600_CelTOS_n103_492_structure.fasta')

pdb = read.pdb('final_pdbs/celTOS-AF-Q8I5P1-F1-model_v4.pdb')

chainA = trim.pdb(pdb,inds = atom.select(pdb, chain='A'))
chainA = trim.pdb(chainA, inds = atom.select(chainA, type = 'ATOM'))

fa_aa = s2c(hi$seq.text[1]) %>%seqinr::translate() %>% paste(collapse = '') %>% gsub('\\*','',.)
pdb_aa = pdbseq(chainA) %>% paste(collapse = '')

x = msaMuscle(c(fa_aa,pdb_aa), type = 'protein', order = 'input') %>% suppressWarnings()
y = x %>% as.character()
str_locate_all(y[1], '[^-*]+')

mi = 35
ma = 164

#mi_nuc = mi * 3 - 2
#ma_nuc = ma * 3

#hi$trim = str_sub(hi$seq.text, mi_nuc, ma_nuc)

#hi$seq.name[1]
#pdb_aa
#hi$trim[1] %>% s2c() %>%seqinr::translate() %>% paste(collapse = '')

hi$trim = hi$seq.text

cat(paste('>',hi$seq.name,'\n',hi$trim,sep=''),
    sep = '\n', file = 'biostructmap_input/celtos_trimmed.fasta')

# trim pdb indices
chainA$atom
inds = atom.select(chainA, resno = mi:ma)
mod = trim.pdb(chainA, inds=inds)

# check range
mod$atom$resno %>% unique()

#mod$atom$resno = mod$atom$resno + (mi - 1)
pdb_aa = pdbseq(mod) %>% paste(collapse = '')

write.pdb(mod, file = 'biostructmap_input/celtos_AF_cleaned.pdb')
