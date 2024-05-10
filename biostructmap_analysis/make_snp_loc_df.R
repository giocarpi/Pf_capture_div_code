# SNP / AA polymorphic positions #
# ------------------------------ #

library(tidyverse)
library(phylotools)
library(seqinr)

# for now drop pfs230 and ama1

fa_f = list.files('biostructmap_input/', full.names = T)
fa_f = fa_f[grep('fasta', fa_f)]
fa_f = fa_f[-1]

ag = c('ama1', 'celtos', 'csp', 'msp119', 'pfs230.1', 'pfs230.11',
       'pfs230.13', 'pfs230.3', 'pfs230.5', 'pfs230.7', 'pfs230.9',
       'pfs25', 'pfs28', 'pfs47', 'pfs4845', 'rh5')


# need to add an adjustment (ex CSP nuc 1 is really nuc 919)
snp_loc_df = tibble(fname = fa_f, ag = ag, 
                    nuc_adj = c(),
                    aa_adj = c(100, 35, 309, 1608, 552, 2445,
                               2828, 915, 1282, 1691, 2049,
                               23, 24, 22, 45, 146),
                    length_check = 0, nuc_count = 0, aa_count = 0, 
                    no_change = 0, aa_pos = 0, nuc_pos = 0)

snp_loc_df$nuc_adj = (snp_loc_df$aa_adj * 3) - 2

# find variable positions ----

# will need aa vector for variable positions frequency 
aavector = strsplit('AVILMWYFSTNQCGPRHKDE*', '')[[1]] 

# Loop through file names ---- 
# Load fasta -- find nucleotide snps -- translate to aa -- find aa variants
for(p in 1:nrow(snp_loc_df)){
  hi = phylotools::read.fasta(snp_loc_df$fname[p])
  
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
  
  snp_loc2 = snp_loc + snp_loc_df$nuc_adj[p] - 1
  
  # loop through amino acids and annotate positions that have variation
  aa_loc = c()
  for(i in 1:nchar(hi$aa_seq[1])){
    ro = unique(str_sub(hi$aa_seq,i,i)) %>% length()
    if(ro > 1){
      aa_loc = c(aa_loc,i)
    }
  }
  
  aa_loc2 = aa_loc + snp_loc_df$aa_adj[p] - 1
  
  # add nuc and aa variable position counts
  snp_loc_df$nuc_count[p] = length(snp_loc)
  snp_loc_df$aa_count[p] = length(aa_loc)
  
  # convert nuc positions to aa position - check which nuc variation are missing in aa
  nuc_snp_to_aa = ceiling(snp_loc/3)
  ro = which(!nuc_snp_to_aa %in% aa_loc)
  snp_loc_df$no_change[p] = length(ro)
  
  # positions of polymorphic aa residues
  snp_loc_df$aa_pos[p] = paste(aa_loc2, collapse='+')
  snp_loc_df$nuc_pos[p] = paste(snp_loc2, collapse = '+')
  
  # are all nuc sequences the same length
  snp_loc_df$length_check[p] = nchar(hi$seq.text) %>% unique() %>% length()
  
  # --- BUILDING SNP_DF is done, now I want aa frequency at variable positions --- #
  # set up a tibble
  hold = tibble(ag = snp_loc_df$ag[p], aa_pos = aa_loc)
  
  # add aa columns
  for(i in aavector){
    col = tolower(i)
    hold[,col] = 0
  }
  
  # go to position, build table of aa counts, add to hold
  len = nrow(hi)
  for(i in 1:nrow(hold)){
    x = str_sub(hi$aa_seq, hold$aa_pos[i], hold$aa_pos[i])
    y = table(x)
    tib = tibble(aa = tolower(names(y)), count = y %>% as.numeric())
    for(j in 1:nrow(tib)){
      hold[i,tib$aa[j]] = tib$count[j] / len
    }
  }
  
  hold$aa_pos = hold$aa_pos + snp_loc_df$aa_adj[p] - 1
  
  tib_name = paste( 'aa_freq_', snp_loc_df$ag[p], sep = '')
  assign(tib_name, hold)
  
  # ---- End of loop - print some statement ---- #
  print('Finishing:')
  print(snp_loc_df$ag[p])
  print('-------')
}



# plot the data ----

plotlist = list()

for(i in 1:nrow(snp_loc_df)){
  x = paste('aa_freq_', snp_loc_df$ag[i], sep = '')
  x = get(x)
  
  ag = x$ag[1]
  
  plot_me = x %>% 
    pivot_longer(cols = tolower(aavector))
  
  plot_me$name = factor(plot_me$name, levels = tolower(aavector))
  plot_me$aa_pos = factor(plot_me$aa_pos)
  plot_me$value = ifelse(plot_me$value == 0, NA, plot_me$value)
  
  plotlist[[i]] = ggplot(plot_me, aes(name, aa_pos, fill=value))+
    geom_tile(color = 'black')+
    scale_fill_viridis_b(breaks = seq(0,1,0.1), na.value = 'white')+
    theme_bw()+
    ggtitle(ag)
}

pdf('snp_loc.pdf', height = 8, width = 5)

for(i in 1:length(plotlist)){
  print(plotlist[[i]])
}

dev.off()

# save the table -----
saveRDS(snp_loc_df, 'snp_loc_df.rds')


