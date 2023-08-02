# ------------------------------------------------------- #
#   Brad Broyles, Purdue University, He Lab - 7/7/23      #
#                                                         #
#   Unpack results from Biostructmap Tajima's D           #
#   Load in previously made dataset with pi and epitope   #
#   Add tajima scores to this dataset                     #
#                                                         #
#   Save a final dataset                                  #
# ------------------------------------------------------- #

# load libraries
library(tidyverse)
library(bio3d)
library(ggpubr)

## function for fixing pdb indices to original 3d7 location
fix_indices = function(csv, pdb, factor = 0){
bye = unique(pdb$atom$resno)

mi = min(bye)
ma = max(bye)

miss = which(!mi:ma %in% bye)

csv$updated = csv$reference_residue
for(i in miss){
  ro = which(csv$updated>=i)
  csv$updated[ro] = csv$updated[ro] + 1
}

csv$updated = csv$updated + factor

## I want to add the pdb sequence here as well
csv$seq = pdbseq(pdb)

return(csv)
}

# ama1 -- 107 start ----
csv = read.csv('tajima_d_models/ama1_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/ama1_trim_tajima.pdb')

ama1 = fix_indices(csv,pdb,factor = 106)

ama1$score = ama1$score %>% as.numeric()

# CSP -- 309 start ----
csv = read.csv('tajima_d_models/csp_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/csp_trim_tajima.pdb')

csp = fix_indices(csv,pdb,factor = 308)

csp$score = csp$score %>% as.numeric()

# Rh5 -- 146 start ----
csv = read.csv('tajima_d_models/rh5_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/rh5_trim_tajima.pdb')

rh5 = fix_indices(csv,pdb,factor = 145)

rh5$score = rh5$score %>% as.numeric()

# Pfs25 -- 23 start ----
csv = read.csv('tajima_d_models/Pfs25_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/Pfs25_trim_tajima.pdb')

pfs25 = fix_indices(csv,pdb,factor = 22)

pfs25$score = pfs25$score %>% as.numeric()

# Pfs4845 -- 45 start ----
csv = read.csv('tajima_d_models/Pfs4845_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/Pfs4845_trim_tajima.pdb')

pfs4845 = fix_indices(csv,pdb,factor = 44)

pfs4845$score = pfs4845$score %>% as.numeric()

# Pfs47 -- 22 start ----
csv = read.csv('tajima_d_models/Pfs47_long_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/Pfs47_long_trim_tajima.pdb')

pfs47 = fix_indices(csv,pdb,factor = 21)

pfs47$score = pfs47$score %>% as.numeric()

########################################################
# Pfs230.1 -- 552 start ----
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

csv = read.csv('tajima_d_models/pfs230/pfs230_seg1_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg1_trim_tajima.pdb')

pfs230.1 = fix_indices(csv,pdb,factor = 551)

pfs230.1$score = pfs230.1$score %>% as.numeric()

# Pfs230.3 -- 915 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg3_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg3_trim_tajima.pdb')

pfs230.3 = fix_indices(csv,pdb,factor = 914)

pfs230.3$score = pfs230.3$score %>% as.numeric()

# Pfs230.5 -- 1282 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg5_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg5_trim_tajima.pdb')

pfs230.5 = fix_indices(csv,pdb,factor = 1281)

pfs230.5$score = pfs230.5$score %>% as.numeric()

# Pfs230.7 -- 1691 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg7_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg7_trim_tajima.pdb')

pfs230.7 = fix_indices(csv,pdb,factor = 1690)

pfs230.7$score = pfs230.7$score %>% as.numeric()

# Pfs230.9 -- 2049 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg9_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg9_trim_tajima.pdb')

pfs230.9 = fix_indices(csv,pdb,factor = 2048)

pfs230.9$score = pfs230.9$score %>% as.numeric()

# Pfs230.11 -- 2445 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg11_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg11_trim_tajima.pdb')

pfs230.11 = fix_indices(csv,pdb,factor = 2444)

pfs230.11$score = pfs230.11$score %>% as.numeric()

# Pfs230.13 -- 2828 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg13_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg13_trim_tajima.pdb')

pfs230.13 = fix_indices(csv,pdb,factor = 2827)

pfs230.13$score = pfs230.13$score %>% as.numeric()

# Pfs230.2 -- 730 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg2_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg2_trim_tajima.pdb')

pfs230.2 = fix_indices(csv,pdb,factor = 729)

pfs230.2$score = pfs230.2$score %>% as.numeric()

# Pfs230.4 -- 1133 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg4_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg4_trim_tajima.pdb')

pfs230.4 = fix_indices(csv,pdb,factor = 1132)

pfs230.4$score = pfs230.4$score %>% as.numeric()

# Pfs230.6 -- 1432 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg6_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg6_trim_tajima.pdb')

pfs230.6 = fix_indices(csv,pdb,factor = 1431)

pfs230.6$score = pfs230.6$score %>% as.numeric()

# Pfs230.8 -- 1907 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg8_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg8_trim_tajima.pdb')

pfs230.8 = fix_indices(csv,pdb,factor = 1906)

pfs230.8$score = pfs230.8$score %>% as.numeric()

# Pfs230.10 -- 2201 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg10_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg10_trim_tajima.pdb')

pfs230.10 = fix_indices(csv,pdb,factor = 2200)

pfs230.10$score = pfs230.10$score %>% as.numeric()

# Pfs230.12 -- 2663 start ----
csv = read.csv('tajima_d_models/pfs230/pfs230_seg12_trim_tajima.csv')
pdb = read.pdb('tajima_d_models/pfs230/pfs230_seg12_trim_tajima.pdb')

pfs230.12 = fix_indices(csv,pdb,factor = 2662)

pfs230.12$score = pfs230.12$score %>% as.numeric()

#############################################################
# load in combined_pi and add tajima ----
df = read_rds('combined_pi.rds')

df$tajima = 0

df$tajima[which(df$name=='csp')] = csp$score
df$tajima[which(df$name=='ama1')] = ama1$score
df$tajima[which(df$name=='rh5')] = rh5$score
df$tajima[which(df$name=='pfs25')] = pfs25$score
df$tajima[which(df$name=='pfs4845')] = pfs4845$score
df$tajima[which(df$name=='pfs47')] = pfs47$score
df$tajima[which(df$name=='pfs230.1')] = pfs230.1$score
df$tajima[which(df$name=='pfs230.3')] = pfs230.3$score
df$tajima[which(df$name=='pfs230.5')] = pfs230.5$score
df$tajima[which(df$name=='pfs230.7')] = pfs230.7$score
df$tajima[which(df$name=='pfs230.9')] = pfs230.9$score
df$tajima[which(df$name=='pfs230.11')] = pfs230.11$score
df$tajima[which(df$name=='pfs230.13')] = pfs230.13$score
df$tajima[which(df$name=='pfs230.2')] = pfs230.2$score
df$tajima[which(df$name=='pfs230.4')] = pfs230.4$score
df$tajima[which(df$name=='pfs230.6')] = pfs230.6$score
df$tajima[which(df$name=='pfs230.8')] = pfs230.8$score
df$tajima[which(df$name=='pfs230.10')] = pfs230.10$score
df$tajima[which(df$name=='pfs230.12')] = pfs230.12$score

df = df %>% select(name, seq, updated, score, tajima, epitope)
colnames(df) = c('name','res','ref_pos','pi','tajima','epitope')

#saveRDS(df, 'combined_pi_tajima_epitope.rds')

## add one more column for known domain boundries ----
#df = read_rds('combined_pi_tajima_epitope.rds')

df$domain = ''
ro = grep('pfs230',df$name)
df$domain[ro[which(df$ref_pos[ro] %in% 589:730)]] = 'D1'
df$domain[ro[which(df$ref_pos[ro] %in% 733:887)]] = 'D2'
df$domain[ro[which(df$ref_pos[ro] %in% 918:1133)]] = 'D3'
df$domain[ro[which(df$ref_pos[ro] %in% 1136:1275)]] = 'D4'
df$domain[ro[which(df$ref_pos[ro] %in% 1285:1432)]] = 'D5'
df$domain[ro[which(df$ref_pos[ro] %in% 1435:1560)]] = 'D6'
df$domain[ro[which(df$ref_pos[ro] %in% 1694:1907)]] = 'D7'
df$domain[ro[which(df$ref_pos[ro] %in% 1910:2035)]] = 'D8'
df$domain[ro[which(df$ref_pos[ro] %in% 2052:2199)]] = 'D9'
df$domain[ro[which(df$ref_pos[ro] %in% 2204:2374)]] = 'D10'
df$domain[ro[which(df$ref_pos[ro] %in% 2448:2663)]] = 'D11'
df$domain[ro[which(df$ref_pos[ro] %in% 2666:2827)]] = 'D12'
df$domain[ro[which(df$ref_pos[ro] %in% 2831:2979)]] = 'D13'
df$domain[ro[which(df$ref_pos[ro] %in% 2982:3113)]] = 'D14'

# save updated RDS
hold = df %>% filter(name %in% c('ama1','csp','rh5',
'pfs25','pfs4845','pfs47','pfs230.1','pfs230.3','pfs230.5',
'pfs230.7','pfs230.9','pfs230.11','pfs230.13'))

hold$domain = ifelse(hold$domain=='',hold$name,hold$domain)

# save tajima combined with pi and epitope info ----
saveRDS(hold, 'combined_pi_tajima_epitope2.rds')

