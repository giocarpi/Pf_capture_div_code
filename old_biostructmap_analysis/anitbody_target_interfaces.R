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
# ------------------------------------------------------- #

# load libraries
library(tidyverse)
library(bio3d)
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
csp = 'MMRKLAILSVSSFLFVEALFQEYQCYGSSSNTRVLNELNYDNAG TNLYNELEMNYYGKQENWYSLKKNSRSLGENDDGNNEDNEKLRKPKHKKLKQPADGNP DPNANPNVDPNANPNVDPNANPNVDPNANPNANPNANPNANPNANPNANPNANPNANP NANPNANPNANPNANPNANPNANPNANPNANPNANPNVDPNANPNANPNANPNANPNA NPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNANPNKNN QGNGQGHNMPNDPNRNVDENANANSAVKNNNNEEPSDKHIKEYLNKIQNSLSTEWSPC SVTCGNGIQVRIKPGSANKPKDELDYANDIEKKICKMEKCSSVFNVVNSSIGLIMVLS FLFLN'
csp = gsub(' |\n','',csp)

csp_tibble = tibble(targ_id = '', interface = '', range = '', missing_res = '', mis_match = NA)
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
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  
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
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  
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
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  
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
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  
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
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  
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
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  
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
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  
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
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  
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
  
  # save into csp_tibble
  csp_tibble[nrow(csp_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  csp_tibble[nrow(csp_tibble), 'interface'] = paste(interf,collapse = '+')
  csp_tibble[nrow(csp_tibble),'range'] = my_range
  csp_tibble[nrow(csp_tibble),'mis_match'] = seq
  
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
  # load in 3d7 reference sequence
  ama1 = 'MRKLYCVLLLSAFEFTYMINFGRGQNYWEHPYQNSDVYRPINEH REHPKEYEYPLHQEHTYQQEDSGEDENTLQHAYPIDHEGAEPAPQEQNLFSSIEIVER SNYMGNPWTEYMAKYDIEEVHGSGIRVDLGEDAEVAGTQYRLPSGKCPVFGKGIIIEN SNTTFLTPVATGNQYLKDGGFAFPPTEPLMSPMTLDEMRHFYKDNKYVKNLDELTLCS RHAGNMIPDNDKNSNYKYPAVYDDKDKKCHILYIAAQENNGPRYCNKDESKRNSMFCF RPAKDISFQNYTYLSKNVVDNWEKVCPRKNLQNAKFGLWVDGNCEDIPHVNEFPAIDL FECNKLVFELSASDQPKQYEQHLTDYEKIKEGFKNKNASMIKSAFLPTGAFKADRYKS HGKGYNWGNYNTETQKCEIFNVKPTCLINNSSYIATTALSHPIEVENNFPCSLYKDEI MKEIERESKRIKLNDNDDEGNKKIIAPRIFISDDKDSLKCPCDPEMVSNSTCRFFVCK CVERRAEVTSNNEVVVKEEYKDEYADIPEHKPTYDKMKIIIASSAAVAVLATILMVYL YKRKGNAEKYDKMDEPQDYGKSNSRNDEMLDPEASFWGEEKRASHTTPVLMEKPYY '
  ama1 = gsub(' |\n','',ama1)
  
  ama1_tibble = tibble(targ_id = '', interface = '', range = '', missing_res = '', mis_match = NA)
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
  
  # save into ama1_tibble
  ama1_tibble[nrow(ama1_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  ama1_tibble[nrow(ama1_tibble), 'interface'] = paste(interf,collapse = '+')
  ama1_tibble[nrow(ama1_tibble),'range'] = my_range
  ama1_tibble[nrow(ama1_tibble),'mis_match'] = seq
  
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
  
  # save into ama1_tibble
  ama1_tibble[nrow(ama1_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  ama1_tibble[nrow(ama1_tibble), 'interface'] = paste(interf,collapse = '+')
  ama1_tibble[nrow(ama1_tibble),'range'] = my_range
  ama1_tibble[nrow(ama1_tibble),'mis_match'] = seq
  
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
  # load in 3d7 reference sequence
  rh5 = 'MIRIKKKLILTIIYIHLFILNRLSFENAIKKTKNQENNLTLLPI KSTEEEKDDIKNGKDIKKEIDNDKENIKTNNAKDHSTYIKSYLNTNVNDGLKYLFIPS HNSFIKKYSVFNQINDGMLLNEKNDVKNNEDYKNVDYKNVNFLQYHFKELSNYNIANS IDILQEKEGHLDFVIIPHYTFLDYYKHLSYNSIYHKSSTYGKCIAVDAFIKKINETYD KVKSKCNDIKNDLIATIKKLEHPYDINNKNDDSYRYDISEEIDDKSEETDDETEEVED SIQDTDSNHTPSNKKKNDLMNRTFKKMMDEYNTKKKKLIKCIKNHENDFNKICMDMKN YGTNLFEQLSCYNNNFCNTNGIRYHYDEYIHKLILSVKSKNLNKDLSDMTNILQQSEL LLTNLNKKMGSYIYIDTIKFIHKEMKHIFNRIEYHTKIINDKTKIIQDKIKLNIWRTF QKDELLKRILDMSNEYSLFITSDHLRQMLYNTFYSKEKHLNNIFHHLIYVLQMKFNDV PIKMEYFQTYKKNKPLTQ '
  rh5 = gsub(' |\n','',rh5)
  
  rh5_tibble = tibble(targ_id = '', interface = '', range = '', missing_res = '', mis_match = NA)
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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  
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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  
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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  
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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  
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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq

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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  
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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq

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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq

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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  
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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  
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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  
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
  
  # save into rh5_tibble
  rh5_tibble[nrow(rh5_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  rh5_tibble[nrow(rh5_tibble), 'interface'] = paste(interf,collapse = '+')
  rh5_tibble[nrow(rh5_tibble),'range'] = my_range
  rh5_tibble[nrow(rh5_tibble),'mis_match'] = seq
  
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
  # load in 3d7 reference sequence
  pfs25 = 'MNKLYSLFLFLFIQLSIKYNNAKVTVDTVCKRGFLIQMSGHLEC KCENDLVLVNEETCEEKVLKCDEKTVNKPCGDFSKCIKIDGNPVSYACKCNLGYDMVN NVCIPNECKNVTCGNGKCILDTSNPVKTGVCSCNIGKVPNVQDQNKCSKDGETKCSLK CLKENETCKAVDGIYKCDCKDGFIIDNESSICTAFSAYNILNLSIMFILFSVCFFIM '
  pfs25 = gsub(' |\n','',pfs25)
  
  pfs25_tibble = tibble(targ_id = '', interface = '', range = '', missing_res = '', mis_match = NA)
  pfs25_tibble = pfs25_tibble[-1,]
  
  # 6azz - McLeod/Scally et al 2017: https://doi.org/10.1038/s41467-017-01924-3  ----
  pdb_id = '6azz'
  url = paste("https://files.rcsb.org/download/",pdb_id,".pdb",sep='')
  pdb = read.pdb(url)
  
  table(pdb$atom$chain) #(B,C --> A), (E,F --> D)
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq

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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq 
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq 
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq 
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq
  
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
  
  # save into pfs25_tibble
  pfs25_tibble[nrow(pfs25_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs25_tibble[nrow(pfs25_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs25_tibble[nrow(pfs25_tibble),'range'] = my_range
  pfs25_tibble[nrow(pfs25_tibble),'mis_match'] = seq 
  
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
  # load in 3d7 reference sequence
  pfs4845 = 'MMLYISAKKAQVAFILYIVLVLRIISGNNDFCKPSSLNSEISGF IGYKCNFSNEGVHNLKPDMRERRSIFCTIHSYFIYDKIRLIIPKKSSSPEFKILPEKC FQKVYTDYENRVETDISELGLIEYEIEENDTNPNYNERTITISPFSPKDIEFFCFCDN TEKVISSIEGRSAMVHVRVLKYPHNILFTNLTNDLFTYLPKTYNESNFVSNVLEVELN DGELFVLACELINKKCFQEGKEKALYKSNKIIYHKNLTIFKAPFYVTSKDVNTECTCK FKNNNYKIVLKPKYEKKVIHGCNFSSNVSSKHTFTDSLDISLVDDSAHISCNVHLSEP KYNHLVGLNCPGDIIPDCFFQVYQPESEELEPSNIVYLDSQINIGDIEYYEDAEGDDK IKLFGIVGSIPKTTSFTCICKKDKKSAYMTVTIDSAYYGFLAKTFIFLIVAILLYI '
  pfs4845 = gsub(' |\n','',pfs4845)
  
  pfs4845_tibble = tibble(targ_id = '', interface = '', range = '', missing_res = '', mis_match = NA)
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq

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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq
  
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
  
  # save into pfs4845_tibble
  pfs4845_tibble[nrow(pfs4845_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs4845_tibble[nrow(pfs4845_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs4845_tibble[nrow(pfs4845_tibble),'range'] = my_range
  pfs4845_tibble[nrow(pfs4845_tibble),'mis_match'] = seq

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
  # load in 3d7 reference sequence
  pfs230 = 'MKKIITLKNLFLIILVYIFSEKKDLRCNVIKGNNIKDDEDKRFH LFYYSHNLFKTPETKEKKNKKECFYKNGGIYNLSKEIRMRKDTSVKIKQRTCPFHKEG SSFEMGSKNITCFYPIVGKKERKTLDTIIIKKNVTNDHVVSSDMHSNVQEKNMILIRN IDKENKNDIQNVEEKIQRDTYENKDYESDDTLIEWFDDNTNEENFLLTFLKRCLMKIF SSPKRKKTVVQKKHKSNFFINSSLKYIYMYLTPSDSFNLVRRNRNLDEEDMSPRDNFV IDDEEEEEEEEEEEEEEEEEEEEEEEEEYDDYVYEESGDETEEQLQEEHQEEVGAESS EESFNDEDEDSVEARDGDMIRVDEYYEDQDGDTYDSTIKNEDVDEEVGEEVGEEVGEE VGEEVGEEVGEEVGEEVGEEVGEEEGEEVGEGVGEEVGEEEGEEVGEEEGEYVDEKER QGEIYPFGDEEEKDEGGESFTYEKSEVDKTDLFKFIEGGEGDDVYKVDGSKVLLDDDT ISRVSKKHTARDGEYGEYGEAVEDGENVIKIIRSVLQSGALPSVGVDELDKIDLSYET TESGDTAVSEDSYDKYASNNTNKEYVCDFTDQLKPTESGPKVKKCEVKVNEPLIKVKI ICPLKGSVEKLYDNIEYVPKKSPYVVLTKEETKLKEKLLSKLIYGLLISPTVNEKENN FKEGVIEFTLPPVVHKATVFYFICDNSKTEDDNKKGNRGIVEVYVEPYGNKINGCAFL DEDEEEEKYGNQIEEDEHNEKIKMKTFFTQNIYKKNNIYPCYMKLYSGDIGGILFPKN IKSTTCFEEMIPYNKEIKWNKENKSLGNLVNNSVVYNKEMNAKYFNVQYVHIPTSYKD TLNLFCSIILKEEESNLISTSYLVYVSINEELNFSLFDFYESFVPIKKTIQVAQKNVN NKEHDYTCDFTDKLDKTVPSTANGKKLFICRKHLKEFDTFTLKCNVNKTQYPNIEIFP'
  pfs230 = gsub(' |\n','',pfs230)
  
  pfs230_tibble = tibble(targ_id = '', interface = '', range = '', missing_res = '', mis_match = NA)
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq

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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq

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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq

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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq

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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq

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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq

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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq

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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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
  
  # save into pfs230_tibble
  pfs230_tibble[nrow(pfs230_tibble)+1,'targ_id'] = paste(pdb_id, '.', heav.c, ',', ligh.c, '>', targ.c, sep='')
  pfs230_tibble[nrow(pfs230_tibble), 'interface'] = paste(interf,collapse = '+')
  pfs230_tibble[nrow(pfs230_tibble),'range'] = my_range
  pfs230_tibble[nrow(pfs230_tibble),'mis_match'] = seq
  
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

## Summaries antibody data ----
## write summary to b column

# csp----
df = read_rds('csp_antibody_tibble.rds')
pdb = read.pdb('Final fasta/csp_7rxp_cleaned.pdb')
pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

# 309 is adjustment from reseting pdb indices
y$res = y$res - 309 + 1
y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

write.pdb(pdb,'csp_interface_four.pdb')

# ama1 ----
df = read_rds('AMA1_antibody_tibble.rds')
pdb = read.pdb('Final fasta/ama1_2q8a_cleaned.pdb')
pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

y$res = y$res - 107 + 1
y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

write.pdb(pdb,'ama1_interface_four.pdb')

# rh5 ----
df = read_rds('Rh5_antibody_tibble.rds')
pdb = read.pdb('Final fasta/Rh5_4wat_cleaned.pdb')
pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

y$res = y$res - 146 + 1
y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

write.pdb(pdb,'rh5_interface_four.pdb')

# pfs25 ----
df = read_rds('Pfs25_antibody_tibble.rds')
pdb = read.pdb('Final fasta/pfs25_AF_cleaned.pdb')

pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

y$res = y$res - 23 + 1
y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

write.pdb(pdb,'pfs25_interface_four.pdb')

# pfs4845 ----
df = read_rds('Pfs4845_antibody_tibble.rds')
pdb = read.pdb('Final fasta/pfs4845_AF_cleaned.pdb')

pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

y$res = y$res - 45 + 1
y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

write.pdb(pdb,'pfs4845_interface_four.pdb')

# pfs230 ----
df = read_rds('Pfs230_anitbody_tibble.rds')
pdb = read.pdb('Final fasta/pfs230_segment1_AF_cleaned.pdb')

pdb$atom$b = 0

# unpack residue column
x = strsplit(df$interface,'\\+') %>% unlist() %>% table()
y = tibble(res = names(x) %>% as.numeric(), count = x %>% as.numeric())

y$res = y$res - 552 + 1
y$count = ifelse(y$count>2,3,y$count)

for(i in 1:nrow(y)){
  ro = which(pdb$atom$resno==y$res[i])
  pdb$atom$b[ro] = y$count[i]
}

write.pdb(pdb,'pfs230_interface_four.pdb')

# pfs47 ** Linear epitope from Canepa - 2018:  https://doi.org/10.1038/s41541-018-0065-5 ----
pdb = read.pdb('Final fasta/pfs47_AF_cleaned_long.pdb')

pdb$atom$b = 0

# epitope 178 - 229
epi = 178:229

# convert to adjusted start (22)
epi = epi - 21

ro = which(pdb$atom$resno %in% epi)
pdb$atom$b[ro] = 1

write.pdb(pdb,'pfs47_interface_four.pdb')

# also save antibody tibble
hold = tibble(targ_id = 'linear', interface = paste(178:229, collapse='+'))
saveRDS(hold, 'pfs47_antibody_tibble.rds')
