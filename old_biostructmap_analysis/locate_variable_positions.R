# ------------------------------------------------------- #
#   Brad Broyles, Purdue University, He Lab - 7/7/23      #
#                                                         #
#   Script to find variable positions in nucleotide seqs  #
#     and variable positions in aa sequences              #
#                                                         #
#   Save a dataframe of variable aa positions per protein #
#   Save a dataframe of aa frequencies per position       #
#   Update combined_taj_pi_epitope with polymorphic pos   #
#   Update antibody tibbles with polymorphic positions    #
# ------------------------------------------------------- #

# load libraries
library(seqinr)
library(tidyverse)
library(phylotools)

# load in FASTA datasets ----
ro = grep('.fasta',dir('field_sample_fastas/'))

# set up dataframe for counting snps
snp_loc_df = tibble(fname = dir('field_sample_fastas/')[ro],
                    pname = c('ama1','csp','pfs230','pfs25','pfs47','pfs4845','rh5'),
                    length_check = 0, nuc_count = 0, aa_count = 0, 
                    no_change = 0, aa_pos = '')

# will need aa vector for variable positions frequency 
aavector = strsplit('AVILMWYFSTNQCGPRHKDE*', '')[[1]] 

# Loop through file names ---- 
# Load fasta -- find nucleotide snps -- translate to aa -- find aa variants

for(p in 1:nrow(snp_loc_df)){
  
#get file name and load fasta
fname = paste('field_sample_fastas/', snp_loc_df$fname[p], sep='')
hi = phylotools::read.fasta(fname) %>% as_tibble()

# add an aa translation column
hi$aa_seq = lapply(strsplit(hi$seq.text,''), function(x){
  paste(seqinr::translate(x),collapse = '') %>% 
    return()
}) %>% unlist()

# loop through nucleotides and annotate positions that have variation
snp_loc = c()
for(i in 1:nchar(hi$seq.text[1])){
  ro = unique(str_sub(hi$seq.text,i,i)) %>% length()
  if(ro > 1){
    snp_loc = c(snp_loc,i)
  }
}

# loop through amino acids and annotate positions that have variation
aa_loc = c()
for(i in 1:nchar(hi$aa_seq[1])){
  ro = unique(str_sub(hi$aa_seq,i,i)) %>% length()
  if(ro > 1){
    aa_loc = c(aa_loc,i)
  }
}

# add nuc and aa variable position counts
snp_loc_df$nuc_count[p] = length(snp_loc)
snp_loc_df$aa_count[p] = length(aa_loc)

# convert nuc positions to aa position - check which nuc variation are missing in aa
nuc_snp_to_aa = ceiling(snp_loc/3)
ro = which(!nuc_snp_to_aa %in% aa_loc)
snp_loc_df$no_change[p] = length(ro)

# positions of polymorphic aa residues
snp_loc_df$aa_pos[p] = paste(aa_loc, collapse='+')

# are all nuc sequences the same length
snp_loc_df$length_check[p] = nchar(hi$seq.text) %>% unique() %>% length()

# --- BUILDING SNP_DF is done, now I want aa frequency at variable positions --- #
# set up a tibble
hold = tibble(pname = snp_loc_df$pname[p], aa_pos = aa_loc)

# add aa columns
for(i in aavector){
  col = tolower(i)
  hold[,col] = 0
}

# go to position, build table of aa counts, add to hold
for(i in 1:nrow(hold)){
  x = str_sub(hi$aa_seq, hold$aa_pos[i], hold$aa_pos[i])
  y = table(x)
  tib = tibble(aa = tolower(names(y)), count = y %>% as.numeric())
  for(j in 1:nrow(tib)){
    hold[i,tib$aa[j]] = tib$count[j]
  }
}

tib_name = paste( 'aa_freq_', snp_loc_df$pname[p], sep = '')
assign(tib_name, hold)

# ---- End of loop - print some statement ---- #
print('Finishing:')
print(snp_loc_df$pname[p])
print('-------')
}

# clean up env after loop 
rm(list = c('hi','hold','tib','aa_loc','col','fname','i','j',
            'nuc_snp_to_aa','p','ro','snp_loc','tib_name','x','y'))

# Explore SNP data with plots -----
snp_loc_df %>% pivot_longer(cols = c(nuc_count,aa_count,no_change)) %>% 
  ggplot(.,aes(pname,value,fill=name))+geom_col(position = 'fill')

snp_loc_df %>% pivot_longer(cols = c(nuc_count,aa_count,no_change)) %>% 
  ggplot(.,aes(pname,value,fill=name))+geom_col(position = 'dodge')

# Add mapping to pdb coordinates of trimmed fasta/pdb files ----
snp_loc_df$exp_start = c(107,309,552,23,22,45,146)
snp_loc_df$exp_end = c(438,373,890,193,419,428,513)
snp_loc_df$exp_aa_count = 0
snp_loc_df$exp_aa_pos = ''
snp_loc_df$reset = ''

unpack_pos = function(x){
  return(strsplit(x,'\\+')[[1]] %>% as.numeric())
}

# loop through proteins
# grab residues within structure range
# update aa_pos to exp_pos
for(i in 1:nrow(snp_loc_df)){
  resi = unpack_pos(snp_loc_df$aa_pos[i])
  mi = snp_loc_df$exp_start[i]
  ma = snp_loc_df$exp_end[i]
  # these residues are in solved structure
  keep = which(resi %in% mi:ma)
  snp_loc_df$exp_aa_count[i] = length(keep)
  keep = resi[keep]
  snp_loc_df$exp_aa_pos[i] = paste(keep,collapse = '+')
  # resetting residue index for easy pymol scripts
  keep = keep - mi + 1
  snp_loc_df$reset[i] = paste(keep,collapse = '+')
}

# pymol grab variable residues
paste('select vari, resi ',
      snp_loc_df$reset[7],
      sep = '')

# Unpack these aa_frequency datasets ----

# lets join all these together
aa_freq = rbind(aa_freq_ama1 %>% mutate(index = 1:nrow(aa_freq_ama1)),
                aa_freq_csp %>% mutate(index = 1:nrow(aa_freq_csp)),
                aa_freq_pfs230 %>% mutate(index = 1:nrow(aa_freq_pfs230)),
                aa_freq_pfs25 %>% mutate(index = 1:nrow(aa_freq_pfs25)),
                aa_freq_pfs47 %>% mutate(index = 1:nrow(aa_freq_pfs47)),
                aa_freq_pfs4845 %>% mutate(index = 1:nrow(aa_freq_pfs4845)),
                aa_freq_rh5 %>% mutate(index = 1:nrow(aa_freq_rh5)))

# turn counts into frequency
for(i in 1:nrow(aa_freq)){
  rs = sum(aa_freq[i,3:23])
  aa_freq[i,3:23] = aa_freq[i,3:23]/rs
}

# number of options for a column
aa_freq$options = 0
for(i in 1:nrow(aa_freq)){
 ro = which(aa_freq[i,3:23] > 0)
 aa_freq$options[i] = length(ro)
}

# add if position was in experimental structure
aa_freq$in_exp = 0
for(i in 1:nrow(snp_loc_df)){
  pos = unpack_pos(snp_loc_df$exp_aa_pos[i] %>% unlist())
  ro = which(aa_freq$pname == snp_loc_df$pname[i] & aa_freq$aa_pos %in% pos)
  aa_freq$in_exp[ro] = 1
}

plot_me = aa_freq %>% 
  filter(pname == 'pfs230', in_exp == 1) %>%
  pivot_longer(cols = tolower(aavector))

plot_me$name = factor(plot_me$name, levels = tolower(aavector))
plot_me$aa_pos = factor(plot_me$aa_pos)
plot_me$value = ifelse(plot_me$value == 0, NA, plot_me$value)

ggplot(plot_me, aes(name, aa_pos, fill=value))+
  geom_tile(color = 'black')+
  scale_fill_viridis_c(na.value = 'white')

# save snp_loc_df and aa_freq at polymorphic sites ----

#saveRDS(snp_loc_df,'snp_loc_df.rds')
#saveRDS(aa_freq,'aa_freq_at_snp.rds')

# adding polymorphic data to tajima and pi data ----
cb = read_rds('combined_pi_tajima_epitope2.rds')
cb$aa_variant=0

# edit pfs230 name
hold = snp_loc_df
hold$pname = ifelse(hold$pname=='pfs230','pfs230.1',hold$pname)

for(i in 1:nrow(hold)){
  pos = unpack_pos(hold$exp_aa_pos[i] %>% unlist())
  ro = which(cb$name == hold$pname[i] & cb$ref_pos %in% pos)
  cb$aa_variant[ro] = 1
}

# save an updated version of combined pi tajima ----
#saveRDS(cb, 'combined_pi_tajima_epitope2.rds')

# update antibody tables ----

# - AMA1 - # ----
ab = read_rds('AMA1_antibody_tibble.rds')
ab$snp_res = ''
ab$snp_count = 0

# grab snps for each epitope
for(i in 1:nrow(ab)){
x = ab$interface[i] %>% unpack_pos()
y = snp_loc_df$aa_pos[1] %>% unpack_pos()

snp_res = x[which(x %in% y)]
ab$snp_count[i] = length(snp_res)
ab$snp_res[i] = paste(snp_res, collapse = '+')
}

#saveRDS(ab, 'AMA1_antibody_tibble.rds')

# - CSP - # ----
ab = read_rds('csp_antibody_tibble.rds')
ab$snp_res = ''
ab$snp_count = 0

# grab snps for each epitope
for(i in 1:nrow(ab)){
  x = ab$interface[i] %>% unpack_pos()
  y = snp_loc_df$aa_pos[2] %>% unpack_pos()
  
  snp_res = x[which(x %in% y)]
  ab$snp_count[i] = length(snp_res)
  ab$snp_res[i] = paste(snp_res, collapse = '+')
}

#saveRDS(ab, 'csp_antibody_tibble.rds')

# - Rh5 - # ----
ab = read_rds('Rh5_antibody_tibble.rds')
ab$snp_res = ''
ab$snp_count = 0

# grab snps for each epitope
for(i in 1:nrow(ab)){
  x = ab$interface[i] %>% unpack_pos()
  y = snp_loc_df$aa_pos[7] %>% unpack_pos()
  
  snp_res = x[which(x %in% y)]
  ab$snp_count[i] = length(snp_res)
  ab$snp_res[i] = paste(snp_res, collapse = '+')
}

#saveRDS(ab, 'Rh5_antibody_tibble.rds')

# - Pfs25 - # ----
ab = read_rds('Pfs25_antibody_tibble.rds')
ab$snp_res = ''
ab$snp_count = 0

# grab snps for each epitope
for(i in 1:nrow(ab)){
  x = ab$interface[i] %>% unpack_pos()
  y = snp_loc_df$aa_pos[4] %>% unpack_pos()
  
  snp_res = x[which(x %in% y)]
  ab$snp_count[i] = length(snp_res)
  ab$snp_res[i] = paste(snp_res, collapse = '+')
}

#saveRDS(ab, 'Pfs25_antibody_tibble.rds')

# - Pfs4845 - # ----
ab = read_rds('Pfs4845_antibody_tibble.rds')
ab$snp_res = ''
ab$snp_count = 0

# grab snps for each epitope
for(i in 1:nrow(ab)){
  x = ab$interface[i] %>% unpack_pos()
  y = snp_loc_df$aa_pos[6] %>% unpack_pos()
  
  snp_res = x[which(x %in% y)]
  ab$snp_count[i] = length(snp_res)
  ab$snp_res[i] = paste(snp_res, collapse = '+')
}

#saveRDS(ab, 'Pfs4845_antibody_tibble.rds')

# - Pfs230 - # ----
ab = read_rds('Pfs230_anitbody_tibble.rds')
ab$snp_res = ''
ab$snp_count = 0

# grab snps for each epitope
for(i in 1:nrow(ab)){
  x = ab$interface[i] %>% unpack_pos()
  y = snp_loc_df$aa_pos[3] %>% unpack_pos()
  
  snp_res = x[which(x %in% y)]
  ab$snp_count[i] = length(snp_res)
  ab$snp_res[i] = paste(snp_res, collapse = '+')
}

#saveRDS(ab, 'Pfs230_anitbody_tibble.rds')

# - Pfs47 - # ----
ab = read_rds('pfs47_antibody_tibble.rds')
ab$snp_res = ''
ab$snp_count = 0

# grab snps for each epitope
for(i in 1:nrow(ab)){
  x = ab$interface[i] %>% unpack_pos()
  y = snp_loc_df$aa_pos[5] %>% unpack_pos()
  
  snp_res = x[which(x %in% y)]
  ab$snp_count[i] = length(snp_res)
  ab$snp_res[i] = paste(snp_res, collapse = '+')
}

#saveRDS(ab, 'pfs47_antibody_tibble.rds')
