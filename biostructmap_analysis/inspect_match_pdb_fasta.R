# ---------------------------------------------------------------------------- #
# Pairing fasta sequences with PDB files                                       #
#   main goal is to trim fastas to fit the resolved residues of pdb file       #
#                                                                              #
#   this script does not output any files, only checks inputs files for        #
#   further analysis                                                           #
#                                                                              #
#   Brad Broyles -- Qixin He Lab -- Purdue Univerisity -- 3/4/24               #
# ---------------------------------------------------------------------------- #

# load libraries
library(tidyverse)
library(phylotools)
library(bio3d)
library(seqinr)

# load quick_aln function
source('quick_aln_functions.R')

# make show gaps function
show_gaps = function(x){
  x = gsub(' |\n','',x)
  print(str_locate_all(x, '-+'))
}

# set fasta and pdb directories
fa_dir = 'fasta_regions/'
pdb_dir = 'final_pdbs/'

# get names of fastas
fa_fi = list.files(fa_dir, full.names = T)

# get names of pdbs
pdb_fi = list.files(pdb_dir, full.names = T)

# CSP ----
i = grep('CSP|csp', fa_fi)
fa_fi[i]

j = grep('CSP|csp', pdb_fi)
pdb_fi[j]

fa = phylotools::read.fasta(fa_fi[i])
pdb = read.pdb(pdb_fi[j])

fa_aa = s2c(fa$seq.text[1]) %>% seqinr::translate() %>% paste(collapse = '') %>% gsub('*','',.)
pdb_aa = pdbseq(pdb) %>% paste(collapse = '')

quick_aln(fa_aa, pdb_aa)
show_gaps('--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------EEPSDKHIKEYLNKIQNSLSTEWSPCSVTCGNGIQVRIKPGSANKPKDELDYANDIEKKICKMEK-------------------------')

# AMA1 ----
i = grep('AMA|ama', fa_fi)
fa_fi[i]

j = grep('AMA|ama', pdb_fi)
pdb_fi[j]

fa = phylotools::read.fasta(fa_fi[i])
pdb = read.pdb(pdb_fi[j])

fa_aa = s2c(fa$seq.text[1]) %>% seqinr::translate() %>% paste(collapse = '') %>% gsub('*','',.)
pdb_aa = pdbseq(pdb) %>% paste(collapse = '')

quick_aln(fa_aa, pdb_aa)


# Rh5 ----
i = grep('RH5|Rh5', fa_fi)
fa_fi[i]

j = grep('Rh5|Rh5', pdb_fi)
pdb_fi[j]

fa = phylotools::read.fasta(fa_fi[i])
pdb = read.pdb(pdb_fi[j])

fa_aa = s2c(fa$seq.text[1]) %>% seqinr::translate() %>% paste(collapse = '') %>% gsub('*','',.)
pdb_aa = pdbseq(pdb) %>% paste(collapse = '')

quick_aln(fa_aa, pdb_aa)

# Pfs25 ----
i = grep('Pfs25|pfs25', fa_fi)
fa_fi[i]

j = grep('Pfs25|pfs25', pdb_fi)
pdb_fi[j]

fa = phylotools::read.fasta(fa_fi[i])
pdb = read.pdb(pdb_fi[j])

fa_aa = s2c(fa$seq.text[1]) %>% seqinr::translate() %>% paste(collapse = '') %>% gsub('*','',.)
pdb_aa = pdbseq(pdb) %>% paste(collapse = '')

quick_aln(fa_aa, pdb_aa)


# Pfs4845 ----
i = grep('Pfs45|pfs45', fa_fi)
fa_fi[i]

j = grep('Pfs48|pfs48', pdb_fi)
pdb_fi[j]

fa = phylotools::read.fasta(fa_fi[i])
pdb = read.pdb(pdb_fi[j])

fa_aa = s2c(fa$seq.text[1]) %>% seqinr::translate() %>% paste(collapse = '') %>% gsub('*','',.)
pdb_aa = pdbseq(pdb) %>% paste(collapse = '')

quick_aln(fa_aa, pdb_aa)

# Pfs47 ----
i = grep('Pfs47|pfs47', fa_fi)
fa_fi[i]

j = grep('Pfs47|pfs47', pdb_fi)
pdb_fi[j]

fa = phylotools::read.fasta(fa_fi[i])
pdb = read.pdb(pdb_fi[j[1]])

fa_aa = s2c(fa$seq.text[1]) %>% seqinr::translate() %>% paste(collapse = '') %>% gsub('*','',.)
pdb_aa = pdbseq(pdb) %>% paste(collapse = '')

quick_aln(fa_aa, pdb_aa)

# Pfs230 ----
i = grep('Pfs230|pfs230', fa_fi)
fa_fi[i]

j = grep('Pfs230|pfs230', pdb_fi)
pdb_fi[j]

fa = phylotools::read.fasta(fa_fi[i[13]])
pdb = read.pdb(pdb_fi[j[5]])

fa_aa = s2c(fa$seq.text[1]) %>% seqinr::translate() %>% paste(collapse = '') %>% gsub('*','',.)
pdb_aa = pdbseq(pdb) %>% paste(collapse = '')

quick_aln(fa_aa, pdb_aa)

# Pfs28 new ----
i = grep('Pfs28|pfs28', fa_fi)
fa_fi[i]

j = grep('Pfs28|pfs28', pdb_fi)
pdb_fi[j]

fa = phylotools::read.fasta(fa_fi[i])
pdb = read.pdb(pdb_fi[j[2]])

fa_aa = s2c(fa$seq.text[1]) %>% seqinr::translate() %>% paste(collapse = '') %>% gsub('*','',.)
pdb_aa = pdbseq(pdb) %>% paste(collapse = '')

quick_aln(fa_aa, pdb_aa)

# MSP1 new ----
i = grep('MSP1|msp1', fa_fi)
fa_fi[i]

j = grep('MSP1|msp1', pdb_fi)
pdb_fi[j]

fa = phylotools::read.fasta(fa_fi[i])
pdb = read.pdb(pdb_fi[j[8]])

fa_aa = s2c(fa$seq.text[1]) %>% seqinr::translate() %>% paste(collapse = '') %>% gsub('*','',.)
pdb_aa = pdbseq(pdb) %>% paste(collapse = '')

quick_aln(fa_aa, pdb_aa)

# CelTOS new ----
i = grep('CelTOS', fa_fi)
fa_fi[i]

j = grep('celtos|celTOS', pdb_fi)
pdb_fi[j]

fa = phylotools::read.fasta(fa_fi[i])
pdb = read.pdb(pdb_fi[j[1]])

fa_aa = s2c(fa$seq.text[1]) %>% seqinr::translate() %>% paste(collapse = '') %>% gsub('*','',.)
pdb_aa = pdbseq(pdb) %>% paste(collapse = '')

quick_aln(fa_aa, pdb_aa)
