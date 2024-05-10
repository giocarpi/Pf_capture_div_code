# ------------------------------------------------------- #
#   Brad Broyles, Purdue University, He Lab - 7/7/23      #
#                                                         #
#   Extract epitope interfaces from solved antibody-      #
#   antigen structures using a distance cutoff of 5 Ang   #
#   No structure for Pfs47, but there is linear epitope   #
#                                                         #
#   Save an updated pdb file with b column indicating     #
#   number of times an epitope was in a solved pdb        #
#   antibody-antigen interaction                          #
#                                                         #
#   Update 2/7/23 --- merge same antibody into one hit    #
#                                                         #
# ------------------------------------------------------- #

# load libraries
library(tidyverse)
library(bio3d)
library(seqinr)
library(phylotools)


# custom functions --- calc_dist, check_seq ----
calc_dist = function(pdb1, pdb2, cut=5){
  ## First filter out heteroatoms
  pdb1 = pdb1 %>% filter(type=='ATOM')
  pdb2 = pdb2 %>% filter(type=='ATOM')
  
  ## keep track of amino acid residues that fall within 'cut' angstroms
  dist = tibble(resn = pdb1$resno, within_dist = F)
  
  # compare each atom of pdb1 to every atom of pdb2, if there is a close distance 
  # update this 
  for (i in 1:nrow(pdb1)) {
    my_dist = sqrt((pdb1$x[i]-pdb2$x)^2 + (pdb1$y[i] - pdb2$y)^2 + (pdb1$z[i] - pdb2$z)^2)
    ro = which(my_dist <= cut)
    if(length(ro)>1){
      dist$within_dist[i] = T
    }
  }
  
  # return unique resn that have met the distance criteria
  interface = dist %>% filter(within_dist==T) %>% select(resn) %>% unlist() %>% unique()
  return(interface)
}

check_seq = function(pdb1, ref_seq, chainid='C'){
  ##############################################.
  # Function does 2 things
  # 1. Prints pdb_sequence and matching positions
  # 2. returns the range pdb sequence covered
  ##############################################.
  targ = pdb1$atom %>% filter(chain == chainid)
  hold = targ %>% filter(type == 'ATOM') %>% select(resno)
  mi = min(hold)
  ma = max(hold)
  my_range = paste(mi,'--',ma,sep='')
  
  ## take a pdb file grab the chain id of interest
  ## clean heteroatoms --- and grab sequence
  seq = trim.pdb(pdb1, inds = atom.select(pdb1, chain=chainid)) %>% 
    pdbseq()
  
  tib = tibble(resn = names(seq) %>% as.numeric(), aa = seq)
  tib$match = F
  
  # loop through tib and check aa vs reference seq
  for(i in 1:nrow(tib)){
    let = str_sub(ref_seq,tib$resn[i],tib$resn[i])
    if(let == tib$aa[i]){
      tib$match[i] = T
    }
  }
  
  #print(' ***** Reference Sequence ***** ')
  #print(ref_seq)
  #print('--------------------------------')
  print(' ******* PDB Sequence ******* ')
  print(paste(seq, collapse=''))
  print(' ***** Position Match ***** ')
  seq = ifelse(tib$match==T,'.','x') %>% paste(.,collapse = '')
  print(seq)
  
  
  
  return(c(my_range, seq))
}

wrap_fun = function(pdb1, targ.c, heav.c, ligh.c, dist_cut=5){
  # grab target sequence
  targ = pdb1$atom %>% filter(chain == targ.c)
  
  # grab antibody and paste together if necessary
  heav = pdb$atom %>% filter(chain == heav.c)
  ligh = pdb$atom %>% filter(chain == ligh.c)
  anti = rbind(heav, ligh)
  
  # run calc_dist function
  return(calc_dist(targ, anti, cut = dist_cut))
}

###############################################################
## CSP known antibody interactions ----

# load in 3d7 reference sequence -- copied from GenBank
csp = phylotools::read.fasta('field_sample_fastas/PF3D7_0304600_CSP_CDS.fasta')          # **** updated
csp = csp$seq.text[1]                                                        # **** updated
csp = gsub(' |\n','',csp)
csp = s2c(csp) %>% seqinr::translate() %>% paste(., collapse = '')          # **** updated

csp_tibble = tibble(targ_id = '', interface = '', 
                    range = '', missing_res = '', mis_match = NA,
                    h.c = '', l.c = '')          # **** updated

csp_tibble = csp_tibble[-1,]

  # 6B0S - Scally et al 2018: https://doi.org/10.1084/jem.20170869 ----
  pdb_id = '6b0s'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) # C ~ csp, (H,L)
  targ.c = 'C'
  heav.c = 'H'
  ligh.c = 'L'
  
  # grab range and check sequence
  my_range = check_seq(pdb,csp,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  csp_tibble[nrow(csp_tibble), 'h.c'] = h         # **** updated
  csp_tibble[nrow(csp_tibble), 'l.c'] = l         # **** updated
  
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    csp_tibble[nrow(csp_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 7rxi - Beutler et al 2022: https://doi.org/10.1371/journal.ppat.1010409 ----
  pdb_id = '7rxi'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = 'L'
  
  # grab range and check sequence
  my_range = check_seq(pdb,csp,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  csp_tibble[nrow(csp_tibble), 'h.c'] = h         # **** updated
  csp_tibble[nrow(csp_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    csp_tibble[nrow(csp_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7rxj - Beutler et al 2022 ----
  pdb_id = '7rxj'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)

  table(pdb$atom$chain)  # lots of chains -- inspect in pymol (A,B --> G), (H,L --> I)
  
  # first chain
  targ.c = 'G'
  heav.c = 'A'
  ligh.c = 'B'

  # grab range and check sequence
  my_range = check_seq(pdb,csp,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  csp_tibble[nrow(csp_tibble), 'h.c'] = h         # **** updated
  csp_tibble[nrow(csp_tibble), 'l.c'] = l         # **** updated

  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    csp_tibble[nrow(csp_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # ---- #
  # next chain 
  targ.c = 'I'
  heav.c = 'H'
  ligh.c = 'L'
  
  # grab range and check sequence
  my_range = check_seq(pdb,csp,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  csp_tibble[nrow(csp_tibble), 'h.c'] = h         # **** updated
  csp_tibble[nrow(csp_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    csp_tibble[nrow(csp_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7rxl - Beutler et al 2022 ----
  pdb_id = '7rxl'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) # lots of chains -- inspect in pymol (C,D --> F), (H,L --> E)
  
  # first chain
  targ.c = 'F'
  heav.c = 'C'
  ligh.c = 'D'

  # grab range and check sequence
  my_range = check_seq(pdb,csp,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  csp_tibble[nrow(csp_tibble), 'h.c'] = h         # **** updated
  csp_tibble[nrow(csp_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    csp_tibble[nrow(csp_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # next chain
  targ.c = 'E'
  heav.c = 'H'
  ligh.c = 'L'

  # grab range and check sequence
  my_range = check_seq(pdb,csp,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  csp_tibble[nrow(csp_tibble), 'h.c'] = h         # **** updated
  csp_tibble[nrow(csp_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    csp_tibble[nrow(csp_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7rxp - Beutler et al 2022 ----
  pdb_id = '7rxp'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) # A ~ csp, (H,L)
  
  # first chain
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = 'L'

  # grab range and check sequence
  my_range = check_seq(pdb,csp,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  csp_tibble[nrow(csp_tibble), 'h.c'] = h         # **** updated
  csp_tibble[nrow(csp_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    csp_tibble[nrow(csp_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 7s0x - Beutler et al 2022 ----
  pdb_id = '7s0x'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) # lots of chains -- inspect in pymol (C,D --> F), (H,L --> E)
        
  # first chain
  targ.c = 'F'
  heav.c = 'C'
  ligh.c = 'D'

  # grab range and check sequence
  my_range = check_seq(pdb,csp,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  csp_tibble[nrow(csp_tibble), 'h.c'] = h         # **** updated
  csp_tibble[nrow(csp_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    csp_tibble[nrow(csp_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ---- #
  
  # second chain
  targ.c = 'E'
  heav.c = 'H'
  ligh.c = 'L'

  # grab range and check sequence
  my_range = check_seq(pdb,csp,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  csp_tibble[nrow(csp_tibble), 'h.c'] = h         # **** updated
  csp_tibble[nrow(csp_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    csp_tibble[nrow(csp_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
## SAVE CSP antibody data ----
  saveRDS(csp_tibble, 'csp_antibody_tibble.rds')
###############################################################
## AMA1 known antibody interactions ----
  # load in 3d7 reference sequence -- copied from GenBank
  ama1 = phylotools::read.fasta('field_sample_fastas/PF3D7_1133400_AMA1_CDS.fasta')          # **** updated
  ama1 = ama1$seq.text[1]                                                        # **** updated
  ama1 = gsub(' |\n','',ama1)
  ama1 = s2c(ama1) %>% seqinr::translate() %>% paste(., collapse = '')          # **** updated
  
  ama1_tibble = tibble(targ_id = '', interface = '', 
                      range = '', missing_res = '', mis_match = NA,
                      h.c = '', l.c = '')          # **** updated

  ama1_tibble = ama1_tibble[-1,]
  
  # 2q8a - Coley et al 2007:  https://doi.org/10.1371/journal.ppat.0030172 ----
  pdb_id = '2q8a'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H,L)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = 'L'
  
  # grab range and check sequence
  my_range = check_seq(pdb,ama1,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into ama1_tibble
  ama1_tibble[nrow(ama1_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  ama1_tibble[nrow(ama1_tibble), 'interface'] = paste(interf,collapse = '+')
  ama1_tibble[nrow(ama1_tibble),'range'] = my_range
  ama1_tibble[nrow(ama1_tibble),'mis_match'] = seq
  ama1_tibble[nrow(ama1_tibble), 'h.c'] = h         # **** updated
  ama1_tibble[nrow(ama1_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    ama1_tibble[nrow(ama1_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 2q8b - Coley et al 2007 ----
  pdb_id = '2q8b'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H,L)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = 'L'
  
  my_range = check_seq(pdb,ama1,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into ama1_tibble
  ama1_tibble[nrow(ama1_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  ama1_tibble[nrow(ama1_tibble), 'interface'] = paste(interf,collapse = '+')
  ama1_tibble[nrow(ama1_tibble),'range'] = my_range
  ama1_tibble[nrow(ama1_tibble),'mis_match'] = seq
  ama1_tibble[nrow(ama1_tibble), 'h.c'] = h         # **** updated
  ama1_tibble[nrow(ama1_tibble), 'l.c'] = l         # **** updated
  
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    ama1_tibble[nrow(ama1_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
## SAVE AMA1 antibody data ----
  saveRDS(ama1_tibble, 'AMA1_antibody_tibble.rds')
###############################################################
## Rh5 known antibody interactions ----
  # load in 3d7 reference sequence -- copied from GenBank
  rh5 = phylotools::read.fasta('field_sample_fastas/PF3D7_0424100_RH5_CDS.fasta')          # **** updated
  rh5 = rh5$seq.text[1]                                                        # **** updated
  rh5 = gsub(' |\n','',rh5)
  rh5 = s2c(rh5) %>% seqinr::translate() %>% paste(., collapse = '')          # **** updated
  
  rh5_tibble = tibble(targ_id = '', interface = '', 
                       range = '', missing_res = '', mis_match = NA,
                       h.c = '', l.c = '')          # **** updated
  
  rh5_tibble = rh5_tibble[-1,]

  # 4u0r - Wright et al 2014: https://doi.org/10.1038/nature13715 ----
  pdb_id = '4u0r'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated
  
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 4u1g - Wright et al 2014 ----
  pdb_id = '4u1g'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (E,F --> D)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  targ.c = 'D'
  heav.c = 'E'
  ligh.c = 'F'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated
  
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
   
  # 5mi0 - Campeotto et al 2017: https://doi.org/10.1073/pnas.1616903114 ----
  pdb_id = '5mi0'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 6rcu - Alanine et al 2019: https://doi.org/10.1016/j.cell.2019.05.025 ----
  pdb_id = '6rcu'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (D,E --> A)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  
  targ.c = 'A'
  heav.c = 'D'
  ligh.c = 'E'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 6rcv - Alanine et al 2019 ----
  pdb_id = '6rcv'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (D,E --> A), (G,H --> F), (I,J --> F)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  targ.c = 'A'
  heav.c = 'D'
  ligh.c = 'E'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  
  targ.c = 'F'
  heav.c = 'G'
  ligh.c = 'H'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

# ----- #
  
  targ.c = 'F'
  heav.c = 'I'
  ligh.c = 'J'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7phu - Ragotte et al 2022: https://doi.org/10.1038/s41467-022-28601-4 ----
  pdb_id = '7phu'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (D,E --> A)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  
  targ.c = 'A'
  heav.c = 'D'
  ligh.c = 'E'
  
  # grab range and check sequence
  my_range = check_seq(pdb,rh5,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  rh5_tibble[nrow(rh5_tibble), 'h.c'] = h         # **** updated
  rh5_tibble[nrow(rh5_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    rh5_tibble[nrow(rh5_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
## SAVE Rh5 antibody data ----
  saveRDS(rh5_tibble, 'Rh5_antibody_tibble.rds')
  # load 4wat ----
  pdb_id = '4wat'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
###############################################################
## Pfs25 known antibody interactions ----
  # load in 3d7 reference sequence -- copied from GenBank
  pfs25 = phylotools::read.fasta('field_sample_fastas/PF3D7_1031000_Pfs25_CDS.fasta')          # **** updated
  pfs25 = pfs25$seq.text[1]                                                        # **** updated
  pfs25 = gsub(' |\n','',pfs25)
  pfs25 = s2c(pfs25) %>% seqinr::translate() %>% paste(., collapse = '')          # **** updated
  
  pfs25_tibble = tibble(targ_id = '', interface = '', 
                      range = '', missing_res = '', mis_match = NA,
                      h.c = '', l.c = '')          # **** updated
  
  pfs25_tibble = pfs25_tibble[-1,]
  # 6azz - McLeod/Scally et al 2017: https://doi.org/10.1038/s41467-017-01924-3  ----
  pdb_id = '6azz'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (E,F --> D)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  pdb$atom$resno = pdb$atom$resno + 21   # 2/7 why this is here
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  
  targ.c = 'D'
  heav.c = 'E'
  ligh.c = 'F'
  
  # already updated resn positions
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 6b0h - McLeod/Scally et al 2017  ----
  pdb_id = '6b0h'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(A,B --> J), (C,D --> I)
  targ.c = 'J'
  heav.c = 'A'
  ligh.c = 'B'
  
  pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

# ----- #
  
  targ.c = 'I'
  heav.c = 'C'
  ligh.c = 'D'
  
  #pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 6b0e - McLeod/Scally et al 2017  ----
  pdb_id = '6b0e'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(A,B --> E)
  targ.c = 'E'
  heav.c = 'A'
  ligh.c = 'B'
  
  pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 6b08 - McLeod/Scally et al 2017  ----
  pdb_id = '6b08'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 6b0a - McLeod/Scally et al 2017  ----
  pdb_id = '6b0a'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H,L --> A)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = 'L'
  
  pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 6b0g - McLeod/Scally et al 2017  ----
  pdb_id = '6b0g'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(C,D --> E)
  targ.c = 'E'
  heav.c = 'C'
  ligh.c = 'D'
  
  pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 6phb - McLeod et al 2019: https://doi.org/10.1038/s41467-019-11980-6  ----
  pdb_id = '6phb'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(A,B --> E), (C,D --> I)
  targ.c = 'E'
  heav.c = 'A'
  ligh.c = 'B'
  
  pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

# ----- #
  
  targ.c = 'I'
  heav.c = 'C'
  ligh.c = 'D'
  
  #pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 6phd - McLeod et al 2019  ----
  pdb_id = '6phd'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H,L --> C)
  targ.c = 'C'
  heav.c = 'H'
  ligh.c = 'L'
  
  pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 6phf - McLeod et al 2019  ----
  pdb_id = '6phf'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(A,B --> G), (C,D --> E)
  targ.c = 'G'
  heav.c = 'A'
  ligh.c = 'B'
  
  pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

# ----- #

  targ.c = 'E'
  heav.c = 'C'
  ligh.c = 'D'
  
  #pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 6phc - Mcleod et al 2019 ----
  pdb_id = '6phc'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(A,B --> I), (C,D --> E)
  targ.c = 'I'
  heav.c = 'A'
  ligh.c = 'B'
  
  pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # ----- #
  
  targ.c = 'E'
  heav.c = 'C'
  ligh.c = 'D'
  
  #pdb$atom$resno = pdb$atom$resno + 21
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  
  # 7txw - MacDonald et al 2023: https://doi.org/10.1038/s41541-023-00655-5  ----
  pdb_id = '7txw'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H,L --> A)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = 'L'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs25,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  pfs25_tibble[nrow(pfs25_tibble), 'h.c'] = h         # **** updated
  pfs25_tibble[nrow(pfs25_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs25_tibble[nrow(pfs25_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

## SAVE Pfs25 antibody data ----
  saveRDS(pfs25_tibble, 'Pfs25_antibody_tibble.rds')
################################################################
## pfs4845 known antibody interactions ----
  # load in 3d7 reference sequence -- copied from GenBank
  pfs4845 = phylotools::read.fasta('field_sample_fastas/PF3D7_1346700_Pfs4548_CDS.fasta')          # **** updated
  pfs4845 = pfs4845$seq.text[1]                                                        # **** updated
  pfs4845 = gsub(' |\n','',pfs4845)
  pfs4845 = s2c(pfs4845) %>% seqinr::translate() %>% paste(., collapse = '')          # **** updated
  
  pfs4845_tibble = tibble(targ_id = '', interface = '', 
                        range = '', missing_res = '', mis_match = NA,
                        h.c = '', l.c = '')          # **** updated
  
  pfs4845_tibble = pfs4845_tibble[-1,]
  
  # 7zwf - Ko et al 2022: doi: 10.1038/s41467-018-06340-9  ----
  pdb_id = '7zwf'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B --> A)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  
  # 7zwi - Ko et al 2022  ----
  pdb_id = '7zwi'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (E,F --> D)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

# ----- #  
  
  targ.c = 'D'
  heav.c = 'E'
  ligh.c = 'F'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
    
  # 7zxf - Ko et al 2022  ----
  pdb_id = '7zxf'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (D,E --> A)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # fix first 89 residues of chain A
  hold = trim.pdb(pdb, inds = atom.select(pdb, chain = 'A'))
  seq = pdbseq(hold) %>% paste(.,collapse = '')
  fix = str_sub(seq, 1, 89)
  
  #for (i in 1:85) {
  #  lil_seq = str_sub(fix, i, i+4)
  #  hi = str_locate_all(pfs4845,lil_seq)
  #  print('----------------------')
  #  print(paste(i, lil_seq, sep = ' ------- '))
  #  print(hi)
  #}
  
  # Pos 1-9 (need +2), pos 10-27 (need -2), pos 28-89 (need -2)
  
  ro = which(pdb$atom$chain=='A')
  res = unique(pdb$atom$resno[ro])
  ro2 = which(pdb$atom$resno %in% res[1:9])
  pdb$atom$resno[ro2] = pdb$atom$resno[ro2] + 2
  ro3 = which(pdb$atom$resno %in% res[10:89])
  pdb$atom$resno[ro3] = pdb$atom$resno[ro3] - 2
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # ----- #  
  
  targ.c = 'A'
  heav.c = 'D'
  ligh.c = 'E'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7zwm - Ko et al 2022  ----
  pdb_id = '7zwm'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (D,E --> A), (G,H --> F), (I,J --> F)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #  
  
  targ.c = 'A'
  heav.c = 'D'
  ligh.c = 'E'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  
  targ.c = 'F'
  heav.c = 'G'
  ligh.c = 'H'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  
  targ.c = 'F'
  heav.c = 'I'
  ligh.c = 'J'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
    
  # 7zxg - Ko et al 2022  ----
  pdb_id = '7zxg'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (E,F --> D)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ---- #
  
  targ.c = 'D'
  heav.c = 'E'
  ligh.c = 'F'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 6h5n - Lennartz et al 2018: https://doi.org/10.1038/s41467-018-06340-9  ----
  pdb_id = '6h5n'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (E,F --> D)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # ---- #
  
  targ.c = 'D'
  heav.c = 'E'
  ligh.c = 'F'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 7uxl - Fabra-Garcia et al 2023: https://doi.org/10.1016/j.immuni.2023.01.009  ----
  pdb_id = '7uxl'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(A,B --> R), (E,F --> R)
  targ.c = 'R'
  heav.c = 'A'
  ligh.c = 'B'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # ---- #
  
  targ.c = 'R'
  heav.c = 'E'
  ligh.c = 'F'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 7unb - McLeod et al 2022: https://doi.org/10.1016/j.immuni.2022.07.015 ----
  pdb_id = '7unb'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H,L --> R), (E,F --> R)
  targ.c = 'R'
  heav.c = 'H'
  ligh.c = 'L'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # ---- #
  
  targ.c = 'R'
  heav.c = 'E'
  ligh.c = 'F'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 6e63 - Kundu et al 2018: https://doi.org/10.1038/s41467-018-06742-9  ----
  pdb_id = '6e63'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (H,L --> P)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # ---- #
  
  targ.c = 'P'
  heav.c = 'H'
  ligh.c = 'L'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 6e62 - Kundu et al 2018  ----
  pdb_id = '6e62'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (H,L --> P)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # ---- #
  
  targ.c = 'P'
  heav.c = 'H'
  ligh.c = 'L'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs4845,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  pfs4845_tibble[nrow(pfs4845_tibble), 'h.c'] = h         # **** updated
  pfs4845_tibble[nrow(pfs4845_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs4845_tibble[nrow(pfs4845_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

## SAVE Pfs4845 antibody data ----
  saveRDS(pfs4845_tibble, 'Pfs4845_antibody_tibble.rds')
################################################################
## pfs230 known antibody interactions ----
  # load in 3d7 reference sequence -- copied from GenBank
  pfs230 = phylotools::read.fasta('field_sample_fastas/PF3D7_0209000_Pfs230_CDS.fasta')          # **** updated
  pfs230 = pfs230$seq.text[1]                                                        # **** updated
  pfs230 = gsub(' |\n','',pfs230)
  pfs230 = s2c(pfs230) %>% seqinr::translate() %>% paste(., collapse = '')          # **** updated
  
  pfs230_tibble = tibble(targ_id = '', interface = '', 
                          range = '', missing_res = '', mis_match = NA,
                          h.c = '', l.c = '')          # **** updated
  
  pfs230_tibble = pfs230_tibble[-1,]
  
  # 6ohg - Singh et al 2020:  https://doi.org/10.1038/s42003-020-01123-9  ----
  pdb_id = '6ohg'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7jum - Coehlo et al 2021: https://doi.org/10.1038/s41467-021-21955-1  ----
  pdb_id = '7jum'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H --> A, I --> B, J --> C)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ---- #
  
  targ.c = 'B'
  heav.c = 'I'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

# ---- #
  
  targ.c = 'C'
  heav.c = 'J'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7uvi - Ivanochko et al 2023: https://doi.org/10.1016/j.immuni.2023.01.013  ----
  pdb_id = '7uvi'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(A,B --> C), (D,E --> F)
  targ.c = 'C'
  heav.c = 'B'
  ligh.c = 'A'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

# ----- #
  
  targ.c = 'F'
  heav.c = 'D'
  ligh.c = 'E'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7uvq - Ivanochko et al 2023  ----
  pdb_id = '7uvq'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (D,E --> A)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = 'C'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  
  targ.c = 'A'
  heav.c = 'D'
  ligh.c = 'E'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7uvs - Ivanochko et al 2023  ----
  pdb_id = '7uvs'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(A,B --> C), (D,E --> F)
  targ.c = 'C'
  heav.c = 'B'
  ligh.c = 'A'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  
  targ.c = 'F'
  heav.c = 'D'
  ligh.c = 'E'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7uvo - Ivanochko et al 2023  ----
  pdb_id = '7uvo'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(A,B --> C)
  targ.c = 'C'
  heav.c = 'B'
  ligh.c = 'A'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7ufw - Tang et al 2023: https://doi.org/10.1016/j.immuni.2023.01.012  ----
  pdb_id = '7ufw'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H --> A), (I --> B)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  
  targ.c = 'B'
  heav.c = 'I'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7uvh - Ivanochko et al 2023  ----
  pdb_id = '7uvh'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(A,B --> C), (D,E --> F)
  targ.c = 'C'
  heav.c = 'B'
  ligh.c = 'A'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #
  
  targ.c = 'F'
  heav.c = 'D'
  ligh.c = 'E'
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
    
  # 7u9w - Tang et al 2023  ----
  pdb_id = '7u9w'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H --> A)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7ua2 - Tang et al 2023  ----
  pdb_id = '7ua2'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H --> A), (B -- > A)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

# ----- #  
  
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 7ua8 - Tang et al 2023  ----
  pdb_id = '7ua8'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H --> A), (I -- > B)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #  
  
  targ.c = 'B'
  heav.c = 'I'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

  # 7ubs - Tang et al 2023  ----
  pdb_id = '7ubs'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H --> A), (I -- > B), (J -- > C), (K --> D)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #  
  
  targ.c = 'B'
  heav.c = 'I'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

# ----- #  
  
  targ.c = 'C'
  heav.c = 'J'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

# ----- #  
  
  targ.c = 'D'
  heav.c = 'K'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7ucq - Tang et al 2023  ----
  pdb_id = '7ucq'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H --> A), (I -- > B), (J -- > C), (K --> D)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #  
  
  targ.c = 'B'
  heav.c = 'I'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #  
  
  targ.c = 'C'
  heav.c = 'J'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)

  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- #  
  
  targ.c = 'D'
  heav.c = 'K'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7uc8 - Tang et al 2023  ----
  pdb_id = '7uc8'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H --> A), (I -- > B), (J -- > C)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- # 
  
  targ.c = 'B'
  heav.c = 'I'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }

# ----- #
  
  targ.c = 'C'
  heav.c = 'J'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated

  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7u9e - Tang et al 2023  ----
  pdb_id = '7u9e'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B --> A), (H -- > A)
  targ.c = 'A'
  heav.c = 'B'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- # 
  
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
 
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 7ui1 - Tang et al 2023  ----
  pdb_id = '7ui1'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(H --> A), (I -- > B), (J -- > C), (K -- > D)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- # 
  
  targ.c = 'B'
  heav.c = 'I'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- # 
  
  targ.c = 'C'
  heav.c = 'J'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
# ----- # 
  
  targ.c = 'D'
  heav.c = 'K'
  ligh.c = ''
  
  # grab range and check sequence
  my_range = check_seq(pdb,pfs230,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  pfs230_tibble[nrow(pfs230_tibble), 'h.c'] = h         # **** updated
  pfs230_tibble[nrow(pfs230_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    pfs230_tibble[nrow(pfs230_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
    
  
## SAVE Pfs230 anitbody data ----
  saveRDS(pfs230_tibble, 'Pfs230_anitbody_tibble.rds')
########################################

# Update start positions to these
# cSP starts at 309
# AMA1 starts at 107
# Rh5 at 146
# Pfs25 at 23
# Pfs4845 at 45
# Pfs230 at 552

  # pfs47 ** Linear epitope from Canepa - 2018:  https://doi.org/10.1038/s41541-018-0065-5 ----

  # save antibody tibble
  hold = tibble(targ_id = 'linear', interface = paste(178:229, collapse='+'))
  saveRDS(hold, 'pfs47_antibody_tibble.rds')
  
###########################
  # msp119 known antibody interactions ----
  # load in 3d7 reference sequence -- copied from GenBank
  msp119 = phylotools::read.fasta('field_sample_fastas/PF3D7_0930300_MSP119_CDS.fasta')          # **** updated
  msp119 = msp119$seq.text[1]                                                        # **** updated
  msp119 = gsub(' |\n','',msp119)
  msp119 = s2c(msp119) %>% seqinr::translate() %>% paste(., collapse = '')          # **** updated
  
  msp119_tibble = tibble(targ_id = '', interface = '', 
                         range = '', missing_res = '', mis_match = NA,
                         h.c = '', l.c = '')          # **** updated
  
  msp119_tibble = msp119_tibble[-1,]
  
  # 8DFH Patel et al https://doi.org/10.1038/s41467-022-33336-3 ----
  pdb_id = '8dfh'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = 'L'
  
  # fix pdb resno column (says 2 but should say 1608)
  pdb$atom$resno = pdb$atom$resno + 1606
  
  # grab range and check sequence
  my_range = check_seq(pdb,msp119,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into msp119_tibble
  msp119_tibble[nrow(msp119_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  msp119_tibble[nrow(msp119_tibble), 'interface'] = paste(interf,collapse = '+')
  msp119_tibble[nrow(msp119_tibble),'range'] = my_range
  msp119_tibble[nrow(msp119_tibble),'mis_match'] = seq
  msp119_tibble[nrow(msp119_tibble), 'h.c'] = h         # **** updated
  msp119_tibble[nrow(msp119_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    msp119_tibble[nrow(msp119_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
    
  # 8DFI Patel et al ----
  pdb_id = '8dfi'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = 'L'
  
  # fix pdb resno column (says 2 but should say 1608)
  pdb$atom$resno = pdb$atom$resno + 1606
  
  # grab range and check sequence
  my_range = check_seq(pdb,msp119,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into msp119_tibble
  msp119_tibble[nrow(msp119_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  msp119_tibble[nrow(msp119_tibble), 'interface'] = paste(interf,collapse = '+')
  msp119_tibble[nrow(msp119_tibble),'range'] = my_range
  msp119_tibble[nrow(msp119_tibble),'mis_match'] = seq
  msp119_tibble[nrow(msp119_tibble), 'h.c'] = h         # **** updated
  msp119_tibble[nrow(msp119_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    msp119_tibble[nrow(msp119_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  # 8DFG Patel et al  ----
  pdb_id = '8dfg'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(HL > A, IM > B)
  targ.c = 'A'
  heav.c = 'H'
  ligh.c = 'L'
  
  # fix pdb resno column (says 2 but should say 1608)
  pdb$atom$resno = pdb$atom$resno + 1606
  
  # grab range and check sequence
  my_range = check_seq(pdb,msp119,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into msp119_tibble
  msp119_tibble[nrow(msp119_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  msp119_tibble[nrow(msp119_tibble), 'interface'] = paste(interf,collapse = '+')
  msp119_tibble[nrow(msp119_tibble),'range'] = my_range
  msp119_tibble[nrow(msp119_tibble),'mis_match'] = seq
  msp119_tibble[nrow(msp119_tibble), 'h.c'] = h         # **** updated
  msp119_tibble[nrow(msp119_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    msp119_tibble[nrow(msp119_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  ## and do second structure #(HL > A, IM > B)
  targ.c = 'B'
  heav.c = 'I'
  ligh.c = 'M'
  
  # fix pdb resno column (says 2 but should say 1608)
  #pdb$atom$resno = pdb$atom$resno + 1606
  
  # grab range and check sequence
  my_range = check_seq(pdb,msp119,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into msp119_tibble
  msp119_tibble[nrow(msp119_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  msp119_tibble[nrow(msp119_tibble), 'interface'] = paste(interf,collapse = '+')
  msp119_tibble[nrow(msp119_tibble),'range'] = my_range
  msp119_tibble[nrow(msp119_tibble),'mis_match'] = seq
  msp119_tibble[nrow(msp119_tibble), 'h.c'] = h         # **** updated
  msp119_tibble[nrow(msp119_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    msp119_tibble[nrow(msp119_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  
  
  # 1ob1 Pizzaro et al https://doi.org/10.1016/s0022-2836(03)00376-0 ----
  pdb_id = '1ob1'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(AB > C, DE > F)
  targ.c = 'C'
  heav.c = 'B'
  ligh.c = 'A'
  
  # fix pdb resno column (says 2 but should say 1608)
  pdb$atom$resno = pdb$atom$resno + 1606
  
  # grab range and check sequence
  my_range = check_seq(pdb,msp119,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into msp119_tibble
  msp119_tibble[nrow(msp119_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  msp119_tibble[nrow(msp119_tibble), 'interface'] = paste(interf,collapse = '+')
  msp119_tibble[nrow(msp119_tibble),'range'] = my_range
  msp119_tibble[nrow(msp119_tibble),'mis_match'] = seq
  msp119_tibble[nrow(msp119_tibble), 'h.c'] = h         # **** updated
  msp119_tibble[nrow(msp119_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    msp119_tibble[nrow(msp119_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  ## second one
  table(pdb$atom$chain) #(AB > C, DE > F)
  targ.c = 'F'
  heav.c = 'E'
  ligh.c = 'D'
  
  # fix pdb resno column (says 2 but should say 1608)
  #pdb$atom$resno = pdb$atom$resno + 1606
  
  # grab range and check sequence
  my_range = check_seq(pdb,msp119,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into msp119_tibble
  msp119_tibble[nrow(msp119_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  msp119_tibble[nrow(msp119_tibble), 'interface'] = paste(interf,collapse = '+')
  msp119_tibble[nrow(msp119_tibble),'range'] = my_range
  msp119_tibble[nrow(msp119_tibble),'mis_match'] = seq
  msp119_tibble[nrow(msp119_tibble), 'h.c'] = h         # **** updated
  msp119_tibble[nrow(msp119_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    msp119_tibble[nrow(msp119_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  
  
  
  # 6xqw Thouvenel et al https://doi.org/10.1084/jem.20200942 ----
  pdb_id = '6xqw'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A)
  targ.c = 'E'
  heav.c = 'H'
  ligh.c = 'L'
  
  # fix pdb resno column (says 2 but should say 1608)
  pdb$atom$resno = pdb$atom$resno + 1606
  
  # grab range and check sequence
  my_range = check_seq(pdb,msp119,targ.c)
  seq = my_range[2]
  my_range = my_range[1]
  seq = str_count(seq, 'x')
  
  # find interface
  interf = wrap_fun(pdb,targ.c,heav.c,ligh.c,5)
  
  # get heavy and light sequence
  h = trim.pdb(pdb, inds = atom.select(pdb, chain=heav.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  l = trim.pdb(pdb, inds = atom.select(pdb, chain=ligh.c)) %>% 
    pdbseq() %>% paste(collapse = '')             # **** updated
  
  # save into msp119_tibble
  msp119_tibble[nrow(msp119_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  msp119_tibble[nrow(msp119_tibble), 'interface'] = paste(interf,collapse = '+')
  msp119_tibble[nrow(msp119_tibble),'range'] = my_range
  msp119_tibble[nrow(msp119_tibble),'mis_match'] = seq
  msp119_tibble[nrow(msp119_tibble), 'h.c'] = h         # **** updated
  msp119_tibble[nrow(msp119_tibble), 'l.c'] = l         # **** updated
  
  # check missing residues
  hi = pdb$atom %>% filter(chain == targ.c, type == 'ATOM')
  bye = unique(hi$resno)
  my_r = min(bye):max(bye)
  ro = which(!my_r %in% bye)
  if(length(ro) > 0){
    msp119_tibble[nrow(msp119_tibble),"missing_res"] = paste(my_r[ro],collapse='+')
  }
  
  
## SVAE MSP119 antibody data ----
  saveRDS(msp119_tibble, 'msp119_antibody_tibble.rds')
  
############################################

## Summaries antibody data ----
  ## Need to collapse antibody data ----
  csp_tibble = read_rds('csp_antibody_tibble.rds')
  ama1_tibble = read_rds('AMA1_antibody_tibble.rds')
  rh5_tibble = read_rds('Rh5_antibody_tibble.rds')
  pfs25_tibble = read_rds('Pfs25_antibody_tibble.rds')
  pfs4845_tibble = read_rds('Pfs4845_antibody_tibble.rds')
  pfs230_tibble = read_rds('Pfs230_anitbody_tibble.rds')
  pfs47_tibble = read_rds('pfs47_antibody_tibble.rds')
  msp119_tibble = read_rds('msp119_antibody_tibble.rds')
  
  library(Biostrings)
  
  source('quick_aln_functions.R')
  
  # helper functions
  pairwise_compare = function(df){
    match = tibble(index = 1:nrow(df), lm = '', l2hm = '', hm = '')
    print('working on light to light and heavy to heavy')
    for(i in 1:(nrow(df) - 1)){
      for(j in (i+1):nrow(df)){
        # checking light chain similarity
        x = AAString(df$l.c[i])
        y = AAString(df$l.c[j])
        
        z = pairwiseAlignment(x, y)
        
        a = pid(z, type = 'PID1') %>% round(4)
        b = pid(z, type = 'PID2') %>% round(4)
        c = pid(z, type = 'PID3') %>% round(4)
        d = pid(z, type = 'PID4') %>% round(4)
        
        pos = which(c(a,b,c,d) > 98)
        
        if(T %in% (c(a,b,c,d)>90)){
        print('light chain compare')
        print(paste(i, j, a, b, c, d, sep = '|'))
        }
        
        if(length(pos) > 0){
          match$lm[i] = paste(match$lm[i], '+', j, sep = '')
        }
        
        # ------------------------------- #
        # checking heave chain similarity
        x = AAString(df$h.c[i])
        y = AAString(df$h.c[j])
        
        z = pairwiseAlignment(x, y)
        
        a = pid(z, type = 'PID1') %>% round(4)
        b = pid(z, type = 'PID2') %>% round(4)
        c = pid(z, type = 'PID3') %>% round(4)
        d = pid(z, type = 'PID4') %>% round(4)
        
        pos = which(c(a,b,c,d) > 98)

        if(T %in% (c(a,b,c,d)>90)){
          print('heavy chain compare')
          print(paste(i, j, a, b, c, d, sep = '|'))
        }
        
        if(length(pos) > 0){
          match$hm[i] = paste(match$hm[i], '+', j, sep = '')
        }
      }
    }
    
    print('working on light to heavy')
    for(i in 1:nrow(df)){
      for(j in 1:nrow(df)){
        # ------------------------------- #
        # checking light to heavy chain similarity
        x = AAString(df$l.c[i])
        y = AAString(df$h.c[j])
        
        z = pairwiseAlignment(x, y)
        
        a = pid(z, type = 'PID1') %>% round(4)
        b = pid(z, type = 'PID2') %>% round(4)
        c = pid(z, type = 'PID3') %>% round(4)
        d = pid(z, type = 'PID4') %>% round(4)
        
        pos = which(c(a,b,c,d) > 98)

        if(T %in% (c(a,b,c,d)>90)){
          print('light vs heavy chain compare')
          print(paste(i, j, a, b, c, d, sep = '|'))
        }
        
        if(length(pos) > 0){
          match$l2hm[i] = paste(match$l2hm[i], '+', j, sep = '')
        }
      }
    }
    
    return(match)
  }
  
  new_merged_ro = function(df, ros){
    interfs = df$interface[ros]
    ids = df$targ_id[ros]
    ranges = df$range[ros]
    missings = df$missing_res[ros]
    mismatches = df$mis_match[ros]
    
    interfs = paste(interfs, collapse = '+')
    interfs = str_split(interfs, '\\+')[[1]] %>% unique() %>% sort() %>% paste(collapse = '+')
    
    ids = paste(ids, collapse = ' | ')
    ranges = paste(ranges, collapse = ' | ')
    missings = paste(missings, collapse = ' | ')
    mismatches = paste(mismatches, collapse = ' | ')
    
    x = tibble(targ_id = ids, interface = interfs, range = ranges, missing_res = missings,
               mis_match = mismatches, h.c = df$h.c[ros[1]], l.c = df$l.c[ros[1]])
    
    return(x)
  }
  
# csp ----
  #match = pairwise_compare(csp_tibble)
  #match 
  
  # this time I will just manually merge based on literature review
  csp_tibble
  
  # dups 3-4, 5-6, 8-9
  ros1 = c(3,4)
  new_ro1 = new_merged_ro(csp_tibble, ros1)
  
  ros2 = c(5,6)
  new_ro2 = new_merged_ro(csp_tibble, ros2)
  
  ros3 = c(8,9)
  new_ro3 = new_merged_ro(csp_tibble, ros3)
  
  # add to csp
  csp_tibble = rbind(csp_tibble, new_ro1, new_ro2, new_ro3)
  
  # remove merged rows
  csp_tibble = csp_tibble[-c(ros1, ros2, ros3),]
  
  # add Ab name
  csp_tibble$ab_name = c('1710', 'Fab234', 'Fab1512', 'Fab236', 'Fab1488', 'Fab352')
  
  # drop the h.c and l.c columns
  csp_tibble = csp_tibble %>% select(-h.c, -l.c)
  
  # save updated table
  saveRDS(csp_tibble, 'collapsed_csp_antibody_tibble.rds')

  
# ama1 ----
  #match = pairwise_compare(ama1_tibble)
  #match 
  
  #quick_aln(ama1_tibble$h.c[1], ama1_tibble$h.c[2])
  
  # dups 1-2 (but heavy chain doesn't line up -- different at one AA pos, both labeled as 1F9)
  ros1 = c(1,2)
  new_ro1 = new_merged_ro(ama1_tibble, ros1)
  
  # add to ama1
  ama1_tibble = rbind(ama1_tibble, new_ro1)
  
  # remove merged rows
  ama1_tibble = ama1_tibble[-c(ros1),]
  
  # add Ab name
  ama1_tibble$ab_name = c('1F9')
  
  # drop the h.c and l.c columns
  ama1_tibble = ama1_tibble %>% select(-h.c, -l.c)
  
  # save updated table
  saveRDS(ama1_tibble, 'collapsed_ama1_antibody_tibble.rds')
  
# rh5 ----
  #match = pairwise_compare(rh5_tibble)
  #match 
  
  rh5_tibble
  
  # dups 1-4, 2-3, 6-8-10-12, 7-9
  ros1 = c(1,4)
  new_ro1 = new_merged_ro(rh5_tibble, ros1)
  
  ros2 = c(2,3)
  new_ro2 = new_merged_ro(rh5_tibble, ros2)
  
  ros3 = c(6,8,10,12)
  new_ro3 = new_merged_ro(rh5_tibble, ros3)
  
  ros4 = c(7,9)
  new_ro4 = new_merged_ro(rh5_tibble, ros4)
  
  # add to csp
  rh5_tibble = rbind(rh5_tibble, new_ro1, new_ro2, new_ro3, new_ro4)
  
  # remove merged rows
  rh5_tibble = rh5_tibble[-c(ros1, ros2, ros3, ros4),]
  
  # add Ab name
  rh5_tibble$ab_name = c('R5.004', 'R5.015', '9AD4', 'QA1', 'R5.016', 'R5.011')
  
  # drop the h.c and l.c columns
  rh5_tibble = rh5_tibble %>% select(-h.c, -l.c)
  
  # save updated table
  saveRDS(rh5_tibble, 'collapsed_rh5_antibody_tibble.rds')
  
# pfs25 ----
  #match = pairwise_compare(pfs25_tibble)
  #match 
  
  pfs25_tibble
  
  quick_aln(pfs25_tibble$l.c[2], pfs25_tibble$l.c[6])
  
  # dups 1-2, 3-4, 9-10, 12-13 (hc no match), 14-15
  ros1 = c(1,2)
  new_ro1 = new_merged_ro(pfs25_tibble, ros1)
  
  ros2 = c(3,4)
  new_ro2 = new_merged_ro(pfs25_tibble, ros2)
  
  ros3 = c(9,10)
  new_ro3 = new_merged_ro(pfs25_tibble, ros3)
  
  ros4 = c(12,13)
  new_ro4 = new_merged_ro(pfs25_tibble, ros4)
  
  ros5 = c(14,15)
  new_ro5 = new_merged_ro(pfs25_tibble, ros5)
  
  # add to csp
  pfs25_tibble = rbind(pfs25_tibble, 
                       new_ro1, new_ro2, new_ro3, new_ro4, new_ro5)
  
  # remove merged rows
  pfs25_tibble = pfs25_tibble[-c(ros1, ros2, ros3, ros4, ros5),]
  
  # add Ab name
  pfs25_tibble$ab_name = c('1260', '1276', '1269', '1245', '2586', '1G2',
                           '1190', '1262', '2530', '2587', '2544')
  
  # drop the h.c and l.c columns
  pfs25_tibble = pfs25_tibble %>% select(-h.c, -l.c)
  
  # save updated table
  saveRDS(pfs25_tibble, 'collapsed_pfs25_antibody_tibble.rds')
  
# pfs4845 ----
  #match = pairwise_compare(pfs4845_tibble)
  #match 
  
  pfs4845_tibble %>% head(n = 21)
  
  quick_aln(pfs4845_tibble$h.c[15], pfs4845_tibble$h.c[17])
  
  # dups (1,2,3,6,8) - (5,7,9,10,11) - (4,12,13,20,21) - (15,17) - 18,19 - 
  ros1 = c(1,2,3,6,8)
  new_ro1 = new_merged_ro(pfs4845_tibble, ros1)
  
  ros2 = c(5,7,9,10,11)
  new_ro2 = new_merged_ro(pfs4845_tibble, ros2)
  
  ros3 = c(4,12,13,20,21)
  new_ro3 = new_merged_ro(pfs4845_tibble, ros3)
  
  ros4 = c(18,19)
  new_ro4 = new_merged_ro(pfs4845_tibble, ros4)
  
  # add to csp
  pfs4845_tibble = rbind(pfs4845_tibble, 
                       new_ro1, new_ro2, new_ro3, new_ro4)
  
  # remove merged rows
  pfs4845_tibble = pfs4845_tibble[-c(ros1, ros2, ros3, ros4),]
  
  # add Ab name
  pfs4845_tibble$ab_name = c('RUPA-29','RUPA-44','RUPA-47','RUPA-117',
                             '32F3', '10D8', '85RF45.1', 'TB31F')
  
  # drop the h.c and l.c columns
  pfs4845_tibble = pfs4845_tibble %>% select(-h.c, -l.c)
  
  # save updated table
  saveRDS(pfs4845_tibble, 'collapsed_pfs4845_antibody_tibble.rds')
  
  
  
# pfs230 ----
  #match = pairwise_compare(pfs230_tibble)
  #match 
  
  pfs230_tibble
  
  quick_aln(pfs230_tibble$h.c[22], pfs230_tibble$h.c[33])
  
  # dups (2,3,4,12,13,18,32) - (5,6) - (9,10) - (14,15) - (19,20) 
  # - (21,22,23,24) - (25,26,27,28) - (29,30,31) - (34,35,36,37)
  ros1 = c(2,3,4,12,13,18,32)
  new_ro1 = new_merged_ro(pfs230_tibble, ros1)
  
  ros2 = c(5,6)
  new_ro2 = new_merged_ro(pfs230_tibble, ros2)
  
  ros3 = c(9,10)
  new_ro3 = new_merged_ro(pfs230_tibble, ros3)
  
  ros4 = c(14,15)
  new_ro4 = new_merged_ro(pfs230_tibble, ros4)
  
  ros5 = c(19,20)
  new_ro5 = new_merged_ro(pfs230_tibble, ros5)
  
  ros6 = c(21,22,23,24)
  new_ro6 = new_merged_ro(pfs230_tibble, ros6)
  
  ros7 = c(25,26,27,28)
  new_ro7 = new_merged_ro(pfs230_tibble, ros7)
  
  ros8 = c(29,30,31)
  new_ro8 = new_merged_ro(pfs230_tibble, ros8)
  
  ros9 = c(34,35,36,37)
  new_ro9 = new_merged_ro(pfs230_tibble, ros9)
  
  # add to csp
  pfs230_tibble = rbind(pfs230_tibble, 
                         new_ro1, new_ro2, new_ro3, new_ro4, 
                        new_ro5, new_ro6, new_ro7, new_ro8, new_ro9)
  
  # remove merged rows
  pfs230_tibble = pfs230_tibble[-c(ros1, ros2, ros3, ros4, 
                                   ros5, ros6, ros7, ros8, ros9),]
  
  # add Ab name
  pfs230_tibble$ab_name = c('4F12', 'RUPA-97', '15C5', 'RUPA-38', '230AS-88',
                            '230AL-18', '230AL-26', 'LMIV230-01', 'RUPA-55',
                            'LMIV230-02', 'RUPA-32', '230AL-20', '230AS-26',
                            '230AS-18', '230AS-73', '230AL-37')
  
  # drop the h.c and l.c columns
  pfs230_tibble = pfs230_tibble %>% select(-h.c, -l.c)
  
  # save updated table
  saveRDS(pfs230_tibble, 'collapsed_pfs230_antibody_tibble.rds')
  
## msp119 ----
  msp119_tibble
  
  ros1 = c(3,4)
  new_ro1 = new_merged_ro(msp119_tibble, ros1)
  
  ros2 = c(5,6)
  new_ro2 = new_merged_ro(msp119_tibble, ros2)
  
  # add to tibble
  msp119_tibble = rbind(msp119_tibble, 
                        new_ro1, new_ro2)
  
  # remove merged rows
  msp119_tibble = msp119_tibble[-c(ros1, ros2),]
  
  # add Ab name
  msp119_tibble$ab_name = c('42C3', '42C11', 'MaliM03',
                            '42D6', 'G17.12')
  
  # drop the h.c and l.c columns
  msp119_tibble = msp119_tibble %>% select(-h.c, -l.c)
  
  # save updated table
  saveRDS(msp119_tibble, 'collapsed_msp119_antibody_tibble.rds')
  
  
#############################

## write 4 level pdb epitope files ----

# csp----
df = read_rds('collapsed_csp_antibody_tibble.rds')
pdb = read.pdb('biostructmap_input/csp_7rxp_cleaned.pdb')
pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

# 0, 1, 2, 3+
y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

# I will just manually remove HETATM
write.pdb(pdb,'csp_interface_four.pdb')

# ama1 ----
df = read_rds('collapsed_ama1_antibody_tibble.rds')
pdb = read.pdb('biostructmap_input/ama1_AF_cleaned_no4g2_short.pdb')
pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

# remove HETATM manually
write.pdb(pdb,'ama1_interface_four.pdb')

# rh5 ----
df = read_rds('collapsed_rh5_antibody_tibble.rds')
pdb = read.pdb('biostructmap_input/rh5_4wat_cleaned.pdb')
pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

# remove HETATM manually
write.pdb(pdb,'rh5_interface_four.pdb')

# pfs25 ----
df = read_rds('collapsed_pfs25_antibody_tibble.rds')
pdb = read.pdb('biostructmap_input/pfs25_AF_cleaned.pdb')

pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

write.pdb(pdb,'pfs25_interface_four.pdb')

# pfs4845 ----
df = read_rds('collapsed_pfs4845_antibody_tibble.rds')
pdb = read.pdb('biostructmap_input/pfs4845_AF_cleaned.pdb')

pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

write.pdb(pdb,'pfs4845_interface_four.pdb')

# pfs230 ----
df = read_rds('collapsed_pfs230_antibody_tibble.rds')
pdb = read.pdb('biostructmap_input/pfs230_seg1.pdb')

pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

write.pdb(pdb,'pfs230_interface_four.pdb')

# pfs47 ----
pdb = read.pdb('biostructmap_input/pfs47_AF_cleaned_long.pdb')

pdb$atom$b = 0

# epitope 178 - 229
epi = 178:229

ro = which(pdb$atom$resno %in% epi)
pdb$atom$b[ro] = 1

write.pdb(pdb,'pfs47_interface_four.pdb')

# msp119 ----
df = read_rds('collapsed_msp119_antibody_tibble.rds')
pdb = read.pdb('biostructmap_input/msp119_AF_cleaned.pdb')

pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

write.pdb(pdb,'msp119_interface_four.pdb')
