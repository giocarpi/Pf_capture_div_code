# ------------------------------------------------------- #
#   Brad Broyles, Purdue University, He Lab - 7/7/23      #
#                                                         #
#   Gathering proportion of Tajima's D and nucleotide pi  #
#   scores in bins of low - moderate - high values.       #
#   Plotting full protein vs epitope residues             #
#            --- Supplemental figure ---                  #
# ------------------------------------------------------- #

# load libraries
library(tidyverse)
library(ggpubr)

# load pi and Tajima's combined data set
df = read_csv('antigen_selection_diversity_and_epitope.csv')

## To start gather binned data for 7 vaccine targets ----
# Get a list of proteins to summarize
pnames = unique(df$pname)[c(1:3, 7, 4, 5, 10, 6, 8, 9)]

# Make a dataframe to store bin values
tib = tibble(protein = pnames, `pi* >= 0.01` = 0, `0.001 <= pi* < 0.01` = 0, 
             `0.0001 <= pi* < 0.001` = 0, `pi* < 0.0001` = 0,
             `D* >= 3` = 0, `2 <= D* < 3` = 0, `1 <= D* < 2` = 0,
             `0 <= D* < 1` = 0, `D* < 0` = 0
             )

# filter out regions where pi is NA --- these were not included in the analysis
df = df %>% filter(!is.na(`pi*`))

# where tajima is NA -- means no diversity around so lets set to -1 so it is inluded in the purifying selection set
df$`D*` = ifelse(is.na(df$`D*`), -1, df$`D*`)

# loop through protein table, count residues in each of the predifined pi and tajima bins
for(i in 1:nrow(tib)){
  #subset df for this protein only
  hold = df %>% filter(pname == tib$protein[i])
  
  # count different nucleotide diversity levels
  tib$`pi* >= 0.01`[i] = which(hold$`pi*`>=0.01) %>% length()
  tib$`0.001 <= pi* < 0.01`[i] = which(hold$`pi*`<0.01 & hold$`pi*`>=0.001) %>% length()
  tib$`0.0001 <= pi* < 0.001`[i] = which(hold$`pi*`<0.001 & hold$`pi*`>=0.0001) %>% length()
  tib$`pi* < 0.0001`[i] = which(hold$`pi*`<0.0001) %>% length()   # what if I turn of > 0 req
  
  # count diferent tajima's bins
  tib$`D* >= 3`[i] = which(hold$`D*`>=3) %>% length()
  tib$`2 <= D* < 3`[i] = which(hold$`D*`<3 & hold$`D*`>=2) %>% length()
  tib$`1 <= D* < 2`[i] = which(hold$`D*`<2 & hold$`D*`>=1) %>% length()
  tib$`0 <= D* < 1`[i] = which(hold$`D*`<1 & hold$`D*`>=0) %>% length()
  tib$`D* < 0`[i] = which(hold$`D*`<0) %>% length()  # what if I add 0 (D <= 0)
}

## adjust protein names
tib$protein[which(tib$protein=='csp')] = 'CSP'
tib$protein[which(tib$protein=='ama1')] = 'AMA1'
tib$protein[which(tib$protein=='rh5')] = 'Rh5'
tib$protein[which(tib$protein=='msp119')] = 'MSP119'
tib$protein[which(tib$protein=='pfs25')] = 'Pfs25'
tib$protein[which(tib$protein=='pfs4845')] = 'Pfs4845'
tib$protein[which(tib$protein=='pfs47')] = 'Pfs47'
tib$protein[which(tib$protein=='pfs230.1')] = 'Pfs230.D1,D2'
tib$protein[which(tib$protein=='pfs28')] = 'Pfs28'
tib$protein[which(tib$protein=='celtos')] = 'CelTOS'


## Now gather data for just epitope residues ----
csp = read_rds('collapsed_csp_antibody_tibble.rds')
ama1 = read_rds('collapsed_AMA1_antibody_tibble.rds')
rh5 = read_rds('collapsed_Rh5_antibody_tibble.rds')
pfs25 = read_rds('collapsed_Pfs25_antibody_tibble.rds')
pfs4845 = read_rds('collapsed_Pfs4845_antibody_tibble.rds')
pfs230.seg1 = read_rds('collapsed_pfs230_antibody_tibble.rds')
msp119 = read_rds('collapsed_msp119_antibody_tibble.rds')

# Use linear epitope of Pfs47 explained in Cruz 2022 (aa 178 - 229)
pfs47 = tibble(interface = paste(178:229, collapse='+'))

## summarize these datasets into any residue found in epitope
tib2 = tibble(protein = c('csp','ama1','rh5','pfs25','pfs4845','pfs230.seg1','pfs47', 'msp119'),
                 interface = '', `pi* >= 0.01` = 0, `0.001 <= pi* < 0.01` = 0, 
              `0.0001 <= pi* < 0.001` = 0, `pi* < 0.0001` = 0,
              `D* >= 3` = 0, `2 <= D* < 3` = 0, `1 <= D* < 2` = 0,
              `0 <= D* < 1` = 0, `D* < 0` = 0)

# grab residues that fall in any of the solved epitope residues
for(i in 1:nrow(tib2)){
  hi = get(tib2$protein[i])
  interf = strsplit(hi$interface,'\\+') %>% unlist() %>% unique()
  tib2$interface[i] = paste(interf,collapse = '+')
}

# loop through protein table, count residues in each of the predifined pi and tajima bins
for(i in 1:nrow(tib2)){
  interf = strsplit(tib2$interface[i], '\\+') %>% unlist() %>% as.numeric()
  
  #subset df for this protein only
  hold = df %>% filter(pname == tib2$protein[i], reference_residue %in% interf)
  
  # count different nucleotide diversity levels
  tib2$`pi* >= 0.01`[i] = which(hold$`pi*`>=0.01) %>% length()
  tib2$`0.001 <= pi* < 0.01`[i] = which(hold$`pi*`<0.01 & hold$`pi*`>=0.001) %>% length()
  tib2$`0.0001 <= pi* < 0.001`[i] = which(hold$`pi*`<0.001 & hold$`pi*`>=0.0001) %>% length()
  tib2$`pi* < 0.0001`[i] = which(hold$`pi*`<0.0001) %>% length()   # could turn off >0
  
  # count diferent tajima's bins
  tib2$`D* >= 3`[i] = which(hold$`D*`>=3) %>% length()
  tib2$`2 <= D* < 3`[i] = which(hold$`D*`<3 & hold$`D*`>=2) %>% length()
  tib2$`1 <= D* < 2`[i] = which(hold$`D*`<2 & hold$`D*`>=1) %>% length()
  tib2$`0 <= D* < 1`[i] = which(hold$`D*`<1 & hold$`D*`>=0) %>% length()
  tib2$`D* < 0`[i] = which(hold$`D*`<0) %>% length()    # could add 0
}

# edit protein names
## adjust protein names
tib2$protein[which(tib2$protein=='csp')] = 'CSP'
tib2$protein[which(tib2$protein=='ama1')] = 'AMA1'
tib2$protein[which(tib2$protein=='rh5')] = 'Rh5'
tib2$protein[which(tib2$protein=='msp119')] = 'MSP119'
tib2$protein[which(tib2$protein=='pfs25')] = 'Pfs25'
tib2$protein[which(tib2$protein=='pfs4845')] = 'Pfs4845'
tib2$protein[which(tib2$protein=='pfs47')] = 'Pfs47'
tib2$protein[which(tib2$protein=='pfs230.seg1')] = 'Pfs230.D1,D2'

# should we create empty columns for pfs28 and celtos
tib2[nrow(tib2)+1,] = list('Pfs28', NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tib2[nrow(tib2)+1,] = list('CelTOS', NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

# drop interface 
tib2 = tib2 %>% select(-interface)

# adjust pnames
tib2$protein = paste('epi: ', tib2$protein, sep = '')

## Combine whole protein and epitope data ----
combined = rbind(tib,tib2)

## create null rows for ggplot spacing
combined[nrow(combined)+1,'protein'] = ''
combined[nrow(combined)+1,'protein'] = ' '
combined[nrow(combined)+1,'protein'] = '  '
combined[nrow(combined)+1,'protein'] = '   '
combined[nrow(combined)+1,'protein'] = '    '
combined[nrow(combined)+1,'protein'] = '     '
combined[nrow(combined)+1,'protein'] = '      '
combined[nrow(combined)+1,'protein'] = '       '
combined[nrow(combined)+1,'protein'] = '        '

# fix pfs230 name
combined$protein = ifelse(combined$protein == 'pfs230.seg1', 'Pfs230.D1,D2', combined$protein)

## plot pi ----
plot_me = combined %>% pivot_longer(cols = contains('pi'))

# better >= sign
plot_me$name = gsub('<=', '\u2264', plot_me$name)
plot_me$name = gsub('>=', '\u2265', plot_me$name)

# swap in pi character
plot_me$name = gsub('pi', '\u03c0', plot_me$name)

# set order of color scale
plot_me$name = factor(plot_me$name,
                      levels = unique(plot_me$name))

# set order of proteins
plot_me$protein = factor(plot_me$protein,
                         levels = unique(plot_me$protein)[c(1, 11, 21,
                                                            2, 12, 22,
                                                            3, 13, 23,
                                                            4, 18, 24,
                                                            5, 14, 25,
                                                            6, 15, 26,
                                                            7, 16, 27,
                                                            8, 17, 28,
                                                            9, 19, 29,
                                                            10, 20)])

p1 = ggplot(plot_me, aes(protein, value, fill = name))+
    geom_col(position = 'fill', color = 'black')+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title= element_blank())+
    ylab('Proportion')+xlab('Protein')+
    scale_fill_brewer(palette = 17, direction = -1)

## plot tajima ----
plot_me = combined %>% pivot_longer(cols = contains('D*'))

# better >= sign
plot_me$name = gsub('<=', '\u2264', plot_me$name)
plot_me$name = gsub('>=', '\u2265', plot_me$name)


# set order of color scale
plot_me$name = factor(plot_me$name,
                      levels = unique(plot_me$name))

# set order of proteins
plot_me$protein = factor(plot_me$protein,
                         levels = unique(plot_me$protein)[c(1, 11, 21,
                                                            2, 12, 22,
                                                            3, 13, 23,
                                                            4, 18, 24,
                                                            5, 14, 25,
                                                            6, 15, 26,
                                                            7, 16, 27,
                                                            8, 17, 28,
                                                            9, 19, 29,
                                                            10, 20)])

p2 = ggplot(plot_me, aes(protein, value, fill = name))+
  geom_col(position = 'fill', color = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title= element_blank())+
  ylab('Proportion')+xlab('Protein')+
  scale_fill_brewer(palette = 17, direction = -1)

## arrange plots -----
p1 = ggpubr::ggarrange(p1+ggtitle('Nucleotide Diversity')+xlab(''),
          p2+ggtitle("Tajima's D"),ncol=1, align = 'v',
          labels = c('A','B'))

  ## save this plot as figure
ggsave('supp_fig_5.png', plot = p1, height = 5.21, width = 6.5,
       dpi = 800)

ggsave('supp_fig_5.pdf', plot = p1, height = 5.21, width = 6.5,
       dpi = 800)

## save csv ----
# lets add total residues, %  > 0.01 pi, % > 0.001 pi, % > 2 D, % > 1 D
combined = combined %>% mutate(
  total = `pi* >= 0.01` + `pi* < 0.0001` + `0.001 <= pi* < 0.01` + `0.0001 <= pi* < 0.001`, 
  p0.01 = `pi* >= 0.01` / total * 100,
  p0.001 = (`pi* >= 0.01` + `0.001 <= pi* < 0.01`) / total * 100,
  p2 = (`D* >= 3` + `2 <= D* < 3`) / total * 100,
  p1 = (`D* >= 3` + `2 <= D* < 3` + `1 <= D* < 2`) / total * 100,
  p_pi_lowest = `pi* < 0.0001` / total * 100,
  p_d_lowest = `D* < 0` / total * 100
)

write_csv(combined, '3d_stats_bins.csv')

#######################################
## same strategy for pfs230 segments ----

# load pi and Tajima's combined data set
df = read_csv('antigen_selection_diversity_and_epitope.csv')

## To start gather binned data for 7 vaccine targets ----
# Get a list of proteins to summarize
pnames = unique(df$pname)[10:16]

# Make a dataframe to store bin values
tib = tibble(protein = pnames, `pi* >= 0.01` = 0, `0.001 <= pi* < 0.01` = 0, 
             `0.0001 <= pi* < 0.001` = 0, `pi* < 0.0001` = 0,
             `D* >= 3` = 0, `2 <= D* < 3` = 0, `1 <= D* < 2` = 0,
             `0 <= D* < 1` = 0, `D* < 0` = 0
)

# filter out regions where pi is NA --- these were not included in the analysis
df = df %>% filter(!is.na(`pi*`))

# where tajima is NA -- means no diversity around so lets set to -1 so it is inluded in the purifying selection set
df$`D*` = ifelse(is.na(df$`D*`), -1, df$`D*`)


# loop through protein table, count residues in each of the predifined pi and tajima bins
for(i in 1:nrow(tib)){
  #subset df for this protein only
  hold = df %>% filter(pname == tib$protein[i])
  
  # count different nucleotide diversity levels
  tib$`pi* >= 0.01`[i] = which(hold$`pi*`>=0.01) %>% length()
  tib$`0.001 <= pi* < 0.01`[i] = which(hold$`pi*`<0.01 & hold$`pi*`>=0.001) %>% length()
  tib$`0.0001 <= pi* < 0.001`[i] = which(hold$`pi*`<0.001 & hold$`pi*`>=0.0001) %>% length()
  tib$`pi* < 0.0001`[i] = which(hold$`pi*`<0.0001) %>% length()   # what if I turn of > 0 req
  
  # count diferent tajima's bins
  tib$`D* >= 3`[i] = which(hold$`D*`>=3) %>% length()
  tib$`2 <= D* < 3`[i] = which(hold$`D*`<3 & hold$`D*`>=2) %>% length()
  tib$`1 <= D* < 2`[i] = which(hold$`D*`<2 & hold$`D*`>=1) %>% length()
  tib$`0 <= D* < 1`[i] = which(hold$`D*`<1 & hold$`D*`>=0) %>% length()
  tib$`D* < 0`[i] = which(hold$`D*`<0) %>% length()  # what if I add 0 (D <= 0)
}

## adjust protein names
tib$protein[which(tib$protein=='pfs230.seg1')] = 'Pfs230.D1,D2'
tib$protein[which(tib$protein=='pfs230.seg3')] = 'Pfs230.D3,D4'
tib$protein[which(tib$protein=='pfs230.seg5')] = 'Pfs230.D5,D6'
tib$protein[which(tib$protein=='pfs230.seg7')] = 'Pfs230.D7,D8'
tib$protein[which(tib$protein=='pfs230.seg9')] = 'Pfs230.D9,D10'
tib$protein[which(tib$protein=='pfs230.seg11')] = 'Pfs230.D11,D12'
tib$protein[which(tib$protein=='pfs230.seg13')] = 'Pfs230.D13,D14'


## create null rows for ggplot spacing
tib[nrow(tib)+1,'protein'] = ''
tib[nrow(tib)+1,'protein'] = ' '
tib[nrow(tib)+1,'protein'] = '  '
tib[nrow(tib)+1,'protein'] = '   '
tib[nrow(tib)+1,'protein'] = '    '
tib[nrow(tib)+1,'protein'] = '     '

## plot pi ----
plot_me = tib %>% pivot_longer(cols = contains('pi'))

# better >= sign
plot_me$name = gsub('<=', '\u2264', plot_me$name)
plot_me$name = gsub('>=', '\u2265', plot_me$name)

# swap in pi character
plot_me$name = gsub('pi', '\u03c0', plot_me$name)

# set order of color scale
plot_me$name = factor(plot_me$name,
                      levels = unique(plot_me$name))

# set order of proteins
plot_me$protein = factor(plot_me$protein,
                         levels = unique(plot_me$protein)[c(1, 8,
                                                            2, 9,
                                                            3, 10,
                                                            4, 11,
                                                            5, 12,
                                                            6, 13, 7)])

p1 = ggplot(plot_me, aes(protein, value, fill = name))+
  geom_col(position = 'fill', color = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title= element_blank())+
  ylab('Proportion')+xlab('Protein')+
  scale_fill_brewer(palette = 17, direction = -1)

## plot tajima ----
plot_me = tib %>% pivot_longer(cols = contains('D*'))

# better >= sign
plot_me$name = gsub('<=', '\u2264', plot_me$name)
plot_me$name = gsub('>=', '\u2265', plot_me$name)


# set order of color scale
plot_me$name = factor(plot_me$name,
                      levels = unique(plot_me$name))

# set order of proteins
plot_me$protein = factor(plot_me$protein,
                         levels = unique(plot_me$protein)[c(1, 8,
                                                            2, 9,
                                                            3, 10,
                                                            4, 11,
                                                            5, 12,
                                                            6, 13, 7)])

p2 = ggplot(plot_me, aes(protein, value, fill = name))+
  geom_col(position = 'fill', color = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title= element_blank())+
  ylab('Proportion')+xlab('Protein')+
  scale_fill_brewer(palette = 17, direction = -1)

## arrange plots -----
ggpubr::ggarrange(p1+ggtitle('Nucleotide Diversity')+xlab(''),
                  p2+ggtitle("Tajima's D"),ncol=1, align = 'v',
                  labels = c('A','B'))

## save this plot as figure
ggsave('supp_fig_pfs230.png', height = 5.21, width = 6.5,
       dpi = 800)

#df[grep('pfs230', df$pname),] %>% view()

## save csv ----
# lets add total residues, %  > 0.01 pi, % > 0.001 pi, % > 2 D, % > 1 D
combined = tib

combined = combined %>% mutate(
  total = `pi* >= 0.01` + `pi* < 0.0001` + `0.001 <= pi* < 0.01` + `0.0001 <= pi* < 0.001`, 
  p0.01 = `pi* >= 0.01` / total * 100,
  p0.001 = (`pi* >= 0.01` + `0.001 <= pi* < 0.01`) / total * 100,
  p2 = (`D* >= 3` + `2 <= D* < 3`) / total * 100,
  p1 = (`D* >= 3` + `2 <= D* < 3` + `1 <= D* < 2`) / total * 100,
  p_pi_lowest = `pi* < 0.0001` / total * 100,
  p_d_lowest = `D* < 0` / total * 100
)

write_csv(combined, 'pfs230_3d_stats_bins.csv')

#######################################
## same idea but limit the plot to surface exposed residues ----
# load pi and Tajima's combined data set
df = read_csv('antigen_selection_diversity_and_epitope.csv')
df = df %>% filter(exposed_yes_no == T)

## To start gather binned data for 7 vaccine targets ----
# Get a list of proteins to summarize
pnames = unique(df$pname)[c(1:3, 7, 4, 5, 10, 6, 8, 9)]

# Make a dataframe to store bin values
tib = tibble(protein = pnames, `pi* >= 0.01` = 0, `0.001 <= pi* < 0.01` = 0, 
             `0.0001 <= pi* < 0.001` = 0, `pi* < 0.0001` = 0,
             `D* >= 3` = 0, `2 <= D* < 3` = 0, `1 <= D* < 2` = 0,
             `0 <= D* < 1` = 0, `D* < 0` = 0
)

# filter out regions where pi is NA --- these were not included in the analysis
df = df %>% filter(!is.na(`pi*`))

# where tajima is NA -- means no diversity around so lets set to -1 so it is inluded in the purifying selection set
df$`D*` = ifelse(is.na(df$`D*`), -1, df$`D*`)

# loop through protein table, count residues in each of the predifined pi and tajima bins
for(i in 1:nrow(tib)){
  #subset df for this protein only
  hold = df %>% filter(pname == tib$protein[i])
  
  # count different nucleotide diversity levels
  tib$`pi* >= 0.01`[i] = which(hold$`pi*`>=0.01) %>% length()
  tib$`0.001 <= pi* < 0.01`[i] = which(hold$`pi*`<0.01 & hold$`pi*`>=0.001) %>% length()
  tib$`0.0001 <= pi* < 0.001`[i] = which(hold$`pi*`<0.001 & hold$`pi*`>=0.0001) %>% length()
  tib$`pi* < 0.0001`[i] = which(hold$`pi*`<0.0001) %>% length()   # what if I turn of > 0 req
  
  # count diferent tajima's bins
  tib$`D* >= 3`[i] = which(hold$`D*`>=3) %>% length()
  tib$`2 <= D* < 3`[i] = which(hold$`D*`<3 & hold$`D*`>=2) %>% length()
  tib$`1 <= D* < 2`[i] = which(hold$`D*`<2 & hold$`D*`>=1) %>% length()
  tib$`0 <= D* < 1`[i] = which(hold$`D*`<1 & hold$`D*`>=0) %>% length()
  tib$`D* < 0`[i] = which(hold$`D*`<0) %>% length()  # what if I add 0 (D <= 0)
}

## adjust protein names
tib$protein[which(tib$protein=='csp')] = 'CSP'
tib$protein[which(tib$protein=='ama1')] = 'AMA1'
tib$protein[which(tib$protein=='rh5')] = 'Rh5'
tib$protein[which(tib$protein=='msp119')] = 'MSP119'
tib$protein[which(tib$protein=='pfs25')] = 'Pfs25'
tib$protein[which(tib$protein=='pfs4845')] = 'Pfs4845'
tib$protein[which(tib$protein=='pfs47')] = 'Pfs47'
tib$protein[which(tib$protein=='pfs230.1')] = 'Pfs230.D1,D2'
tib$protein[which(tib$protein=='pfs28')] = 'Pfs28'
tib$protein[which(tib$protein=='celtos')] = 'CelTOS'


## Now gather data for just epitope residues ----
csp = read_rds('collapsed_csp_antibody_tibble.rds')
ama1 = read_rds('collapsed_AMA1_antibody_tibble.rds')
rh5 = read_rds('collapsed_Rh5_antibody_tibble.rds')
pfs25 = read_rds('collapsed_Pfs25_antibody_tibble.rds')
pfs4845 = read_rds('collapsed_Pfs4845_antibody_tibble.rds')
pfs230.seg1 = read_rds('collapsed_pfs230_antibody_tibble.rds')
msp119 = read_rds('collapsed_msp119_antibody_tibble.rds')

# Use linear epitope of Pfs47 explained in Cruz 2022 (aa 178 - 229)
pfs47 = tibble(interface = paste(178:229, collapse='+'))

## summarize these datasets into any residue found in epitope
tib2 = tibble(protein = c('csp','ama1','rh5','pfs25','pfs4845','pfs230.seg1','pfs47', 'msp119'),
              interface = '', `pi* >= 0.01` = 0, `0.001 <= pi* < 0.01` = 0, 
              `0.0001 <= pi* < 0.001` = 0, `pi* < 0.0001` = 0,
              `D* >= 3` = 0, `2 <= D* < 3` = 0, `1 <= D* < 2` = 0,
              `0 <= D* < 1` = 0, `D* < 0` = 0)

# grab residues that fall in any of the solved epitope residues
for(i in 1:nrow(tib2)){
  hi = get(tib2$protein[i])
  interf = strsplit(hi$interface,'\\+') %>% unlist() %>% unique()
  tib2$interface[i] = paste(interf,collapse = '+')
}

# loop through protein table, count residues in each of the predifined pi and tajima bins
for(i in 1:nrow(tib2)){
  interf = strsplit(tib2$interface[i], '\\+') %>% unlist() %>% as.numeric()
  
  #subset df for this protein only
  hold = df %>% filter(pname == tib2$protein[i], reference_residue %in% interf)
  
  # count different nucleotide diversity levels
  tib2$`pi* >= 0.01`[i] = which(hold$`pi*`>=0.01) %>% length()
  tib2$`0.001 <= pi* < 0.01`[i] = which(hold$`pi*`<0.01 & hold$`pi*`>=0.001) %>% length()
  tib2$`0.0001 <= pi* < 0.001`[i] = which(hold$`pi*`<0.001 & hold$`pi*`>=0.0001) %>% length()
  tib2$`pi* < 0.0001`[i] = which(hold$`pi*`<0.0001) %>% length()   # could turn off >0
  
  # count diferent tajima's bins
  tib2$`D* >= 3`[i] = which(hold$`D*`>=3) %>% length()
  tib2$`2 <= D* < 3`[i] = which(hold$`D*`<3 & hold$`D*`>=2) %>% length()
  tib2$`1 <= D* < 2`[i] = which(hold$`D*`<2 & hold$`D*`>=1) %>% length()
  tib2$`0 <= D* < 1`[i] = which(hold$`D*`<1 & hold$`D*`>=0) %>% length()
  tib2$`D* < 0`[i] = which(hold$`D*`<0) %>% length()    # could add 0
}

# edit protein names
## adjust protein names
tib2$protein[which(tib2$protein=='csp')] = 'CSP'
tib2$protein[which(tib2$protein=='ama1')] = 'AMA1'
tib2$protein[which(tib2$protein=='rh5')] = 'Rh5'
tib2$protein[which(tib2$protein=='msp119')] = 'MSP119'
tib2$protein[which(tib2$protein=='pfs25')] = 'Pfs25'
tib2$protein[which(tib2$protein=='pfs4845')] = 'Pfs4845'
tib2$protein[which(tib2$protein=='pfs47')] = 'Pfs47'
tib2$protein[which(tib2$protein=='pfs230.seg1')] = 'Pfs230.D1,D2'

# should we create empty columns for pfs28 and celtos
tib2[nrow(tib2)+1,] = list('Pfs28', NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
tib2[nrow(tib2)+1,] = list('CelTOS', NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

# drop interface 
tib2 = tib2 %>% select(-interface)

# adjust pnames
tib2$protein = paste('epi: ', tib2$protein, sep = '')

## Combine whole protein and epitope data ----
combined = rbind(tib,tib2)

## create null rows for ggplot spacing
combined[nrow(combined)+1,'protein'] = ''
combined[nrow(combined)+1,'protein'] = ' '
combined[nrow(combined)+1,'protein'] = '  '
combined[nrow(combined)+1,'protein'] = '   '
combined[nrow(combined)+1,'protein'] = '    '
combined[nrow(combined)+1,'protein'] = '     '
combined[nrow(combined)+1,'protein'] = '      '
combined[nrow(combined)+1,'protein'] = '       '
combined[nrow(combined)+1,'protein'] = '        '

# fix pfs230 name
combined$protein = ifelse(combined$protein == 'pfs230.seg1', 'Pfs230.D1,D2', combined$protein)

## plot pi ----
plot_me = combined %>% pivot_longer(cols = contains('pi'))

# better >= sign
plot_me$name = gsub('<=', '\u2264', plot_me$name)
plot_me$name = gsub('>=', '\u2265', plot_me$name)

# swap in pi character
plot_me$name = gsub('pi', '\u03c0', plot_me$name)

# set order of color scale
plot_me$name = factor(plot_me$name,
                      levels = unique(plot_me$name))

# set order of proteins
plot_me$protein = factor(plot_me$protein,
                         levels = unique(plot_me$protein)[c(1, 11, 21,
                                                            2, 12, 22,
                                                            3, 13, 23,
                                                            4, 18, 24,
                                                            5, 14, 25,
                                                            6, 15, 26,
                                                            7, 16, 27,
                                                            8, 17, 28,
                                                            9, 19, 29,
                                                            10, 20)])

p1 = ggplot(plot_me, aes(protein, value, fill = name))+
  geom_col(position = 'fill', color = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title= element_blank())+
  ylab('Proportion')+xlab('Protein')+
  scale_fill_brewer(palette = 17, direction = -1)

## plot tajima ----
plot_me = combined %>% pivot_longer(cols = contains('D*'))

# better >= sign
plot_me$name = gsub('<=', '\u2264', plot_me$name)
plot_me$name = gsub('>=', '\u2265', plot_me$name)


# set order of color scale
plot_me$name = factor(plot_me$name,
                      levels = unique(plot_me$name))

# set order of proteins
plot_me$protein = factor(plot_me$protein,
                         levels = unique(plot_me$protein)[c(1, 11, 21,
                                                            2, 12, 22,
                                                            3, 13, 23,
                                                            4, 18, 24,
                                                            5, 14, 25,
                                                            6, 15, 26,
                                                            7, 16, 27,
                                                            8, 17, 28,
                                                            9, 19, 29,
                                                            10, 20)])

p2 = ggplot(plot_me, aes(protein, value, fill = name))+
  geom_col(position = 'fill', color = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title= element_blank())+
  ylab('Proportion')+xlab('Protein')+
  scale_fill_brewer(palette = 17, direction = -1)

## arrange plots -----
ggpubr::ggarrange(p1+ggtitle('Nucleotide Diversity')+xlab(''),
                  p2+ggtitle("Tajima's D"),ncol=1, align = 'v',
                  labels = c('A','B'))

## save this plot as figure
ggsave('supp_fig_5_surface_only.png', height = 5.21, width = 6.5,
       dpi = 800)


## save csv ----
# lets add total residues, %  > 0.01 pi, % > 0.001 pi, % > 2 D, % > 1 D
combined = combined %>% mutate(
  total = `pi* >= 0.01` + `pi* < 0.0001` + `0.001 <= pi* < 0.01` + `0.0001 <= pi* < 0.001`, 
  p0.01 = `pi* >= 0.01` / total * 100,
  p0.001 = (`pi* >= 0.01` + `0.001 <= pi* < 0.01`) / total * 100,
  p2 = (`D* >= 3` + `2 <= D* < 3`) / total * 100,
  p1 = (`D* >= 3` + `2 <= D* < 3` + `1 <= D* < 2`) / total * 100,
  p_pi_lowest = `pi* < 0.0001` / total * 100,
  p_d_lowest = `D* < 0` / total * 100
)

write_csv(combined, '3d_surface_stats_bins_.csv')

#######################################
## same strategy for pfs230 segments ----

# load pi and Tajima's combined data set
df = read_csv('antigen_selection_diversity_and_epitope.csv')
df = df %>% filter(exposed_yes_no == T)

## To start gather binned data for 7 vaccine targets ----
# Get a list of proteins to summarize
pnames = unique(df$pname)[10:16]

# Make a dataframe to store bin values
tib = tibble(protein = pnames, `pi* >= 0.01` = 0, `0.001 <= pi* < 0.01` = 0, 
             `0.0001 <= pi* < 0.001` = 0, `pi* < 0.0001` = 0,
             `D* >= 3` = 0, `2 <= D* < 3` = 0, `1 <= D* < 2` = 0,
             `0 <= D* < 1` = 0, `D* < 0` = 0
)

# filter out regions where pi is NA --- these were not included in the analysis
df = df %>% filter(!is.na(`pi*`))

# where tajima is NA -- means no diversity around so lets set to -1 so it is inluded in the purifying selection set
df$`D*` = ifelse(is.na(df$`D*`), -1, df$`D*`)


# loop through protein table, count residues in each of the predifined pi and tajima bins
for(i in 1:nrow(tib)){
  #subset df for this protein only
  hold = df %>% filter(pname == tib$protein[i])
  
  # count different nucleotide diversity levels
  tib$`pi* >= 0.01`[i] = which(hold$`pi*`>=0.01) %>% length()
  tib$`0.001 <= pi* < 0.01`[i] = which(hold$`pi*`<0.01 & hold$`pi*`>=0.001) %>% length()
  tib$`0.0001 <= pi* < 0.001`[i] = which(hold$`pi*`<0.001 & hold$`pi*`>=0.0001) %>% length()
  tib$`pi* < 0.0001`[i] = which(hold$`pi*`<0.0001) %>% length()   # what if I turn of > 0 req
  
  # count diferent tajima's bins
  tib$`D* >= 3`[i] = which(hold$`D*`>=3) %>% length()
  tib$`2 <= D* < 3`[i] = which(hold$`D*`<3 & hold$`D*`>=2) %>% length()
  tib$`1 <= D* < 2`[i] = which(hold$`D*`<2 & hold$`D*`>=1) %>% length()
  tib$`0 <= D* < 1`[i] = which(hold$`D*`<1 & hold$`D*`>=0) %>% length()
  tib$`D* < 0`[i] = which(hold$`D*`<0) %>% length()  # what if I add 0 (D <= 0)
}

## adjust protein names
tib$protein[which(tib$protein=='pfs230.seg1')] = 'Pfs230.D1,D2'
tib$protein[which(tib$protein=='pfs230.seg3')] = 'Pfs230.D3,D4'
tib$protein[which(tib$protein=='pfs230.seg5')] = 'Pfs230.D5,D6'
tib$protein[which(tib$protein=='pfs230.seg7')] = 'Pfs230.D7,D8'
tib$protein[which(tib$protein=='pfs230.seg9')] = 'Pfs230.D9,D10'
tib$protein[which(tib$protein=='pfs230.seg11')] = 'Pfs230.D11,D12'
tib$protein[which(tib$protein=='pfs230.seg13')] = 'Pfs230.D13,D14'


## create null rows for ggplot spacing
tib[nrow(tib)+1,'protein'] = ''
tib[nrow(tib)+1,'protein'] = ' '
tib[nrow(tib)+1,'protein'] = '  '
tib[nrow(tib)+1,'protein'] = '   '
tib[nrow(tib)+1,'protein'] = '    '
tib[nrow(tib)+1,'protein'] = '     '

## plot pi ----
plot_me = tib %>% pivot_longer(cols = contains('pi'))

# better >= sign
plot_me$name = gsub('<=', '\u2264', plot_me$name)
plot_me$name = gsub('>=', '\u2265', plot_me$name)

# swap in pi character
plot_me$name = gsub('pi', '\u03c0', plot_me$name)

# set order of color scale
plot_me$name = factor(plot_me$name,
                      levels = unique(plot_me$name))

# set order of proteins
plot_me$protein = factor(plot_me$protein,
                         levels = unique(plot_me$protein)[c(1, 8,
                                                            2, 9,
                                                            3, 10,
                                                            4, 11,
                                                            5, 12,
                                                            6, 13, 7)])

p1 = ggplot(plot_me, aes(protein, value, fill = name))+
  geom_col(position = 'fill', color = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title= element_blank())+
  ylab('Proportion')+xlab('Protein')+
  scale_fill_brewer(palette = 17, direction = -1)

## plot tajima ----
plot_me = tib %>% pivot_longer(cols = contains('D*'))

# better >= sign
plot_me$name = gsub('<=', '\u2264', plot_me$name)
plot_me$name = gsub('>=', '\u2265', plot_me$name)


# set order of color scale
plot_me$name = factor(plot_me$name,
                      levels = unique(plot_me$name))

# set order of proteins
plot_me$protein = factor(plot_me$protein,
                         levels = unique(plot_me$protein)[c(1, 8,
                                                            2, 9,
                                                            3, 10,
                                                            4, 11,
                                                            5, 12,
                                                            6, 13, 7)])

p2 = ggplot(plot_me, aes(protein, value, fill = name))+
  geom_col(position = 'fill', color = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title= element_blank())+
  ylab('Proportion')+xlab('Protein')+
  scale_fill_brewer(palette = 17, direction = -1)

## arrange plots -----
ggpubr::ggarrange(p1+ggtitle('Nucleotide Diversity')+xlab(''),
                  p2+ggtitle("Tajima's D"),ncol=1, align = 'v',
                  labels = c('A','B'))

## save this plot as figure
ggsave('supp_fig_pfs230_surface_only.png', height = 5.21, width = 6.5,
       dpi = 800)

#df[grep('pfs230', df$pname),] %>% view()

## save csv ----
# lets add total residues, %  > 0.01 pi, % > 0.001 pi, % > 2 D, % > 1 D
combined = tib

combined = combined %>% mutate(
  total = `pi* >= 0.01` + `pi* < 0.0001` + `0.001 <= pi* < 0.01` + `0.0001 <= pi* < 0.001`, 
  p0.01 = `pi* >= 0.01` / total * 100,
  p0.001 = (`pi* >= 0.01` + `0.001 <= pi* < 0.01`) / total * 100,
  p2 = (`D* >= 3` + `2 <= D* < 3`) / total * 100,
  p1 = (`D* >= 3` + `2 <= D* < 3` + `1 <= D* < 2`) / total * 100,
  p_pi_lowest = `pi* < 0.0001` / total * 100,
  p_d_lowest = `D* < 0` / total * 100
)

write_csv(combined, 'pfs230_3d_surface_stats_bins.csv')


#######################################
## same strategy for individual antibody epitopes ----
# load pi and Tajima's combined data set
df = read_csv('antigen_selection_diversity_and_epitope.csv')

# filter out regions where pi is NA --- these were not included in the analysis
df = df %>% filter(!is.na(`pi*`))

# where tajima is NA -- means no diversity around so lets set to -1 so it is included in the purifying selection set
df$`D*` = ifelse(is.na(df$`D*`), -1, df$`D*`)

## Now gather data for just epitope residues ----
csp = read_rds('collapsed_csp_antibody_tibble.rds')
ama1 = read_rds('collapsed_AMA1_antibody_tibble.rds')
rh5 = read_rds('collapsed_Rh5_antibody_tibble.rds')
pfs25 = read_rds('collapsed_Pfs25_antibody_tibble.rds')
pfs4845 = read_rds('collapsed_Pfs4845_antibody_tibble.rds')
pfs230.seg1 = read_rds('collapsed_pfs230_antibody_tibble.rds')
msp119 = read_rds('collapsed_msp119_antibody_tibble.rds')

# Use linear epitope of Pfs47 explained in Cruz 2022 (aa 178 - 229)
pfs47 = tibble(interface = paste(178:229, collapse='+'))
pfs47 = tibble(targ_id = 'linear', interface = pfs47$interface, range = '178--229', 
               missing_res = NA, mis_match = NA, ab_name = 'linear', ag = 'pfs47')

# add ag to antibody tables
csp$ag = 'csp'
ama1$ag = 'ama1'
rh5$ag = 'rh5'
pfs25$ag = 'pfs25'
pfs4845$ag = 'pfs4845'
pfs230.seg1$ag = 'pfs230.seg1'
msp119$ag = 'msp119'


ab_df = rbind(csp, ama1, rh5, pfs25, pfs4845, pfs230.seg1, msp119, pfs47)

# drop columns not needed
ab_df = ab_df %>% select(ag, ab_name, interface)

# add bin columns
ab_df = ab_df %>% mutate(polymorphic = 0, `pi* >= 0.01` = 0, `0.001 <= pi* < 0.01` = 0, 
                         `0.0001 <= pi* < 0.001` = 0, `pi* < 0.0001` = 0,
                         `D* >= 3` = 0, `2 <= D* < 3` = 0, `1 <= D* < 2` = 0,
                         `0 <= D* < 1` = 0, `D* < 0` = 0)

# loop through protein table, count residues in each of the predifined pi and tajima bins
for(i in 1:nrow(ab_df)){
  interf = strsplit(ab_df$interface[i], '\\+') %>% unlist() %>% as.numeric()
  
  #subset df for this protein only
  hold = df %>% filter(pname == ab_df$ag[i], reference_residue %in% interf)
  
  # count different nucleotide diversity levels
  ab_df$`pi* >= 0.01`[i] = which(hold$`pi*`>=0.01) %>% length()
  ab_df$`0.001 <= pi* < 0.01`[i] = which(hold$`pi*`<0.01 & hold$`pi*`>=0.001) %>% length()
  ab_df$`0.0001 <= pi* < 0.001`[i] = which(hold$`pi*`<0.001 & hold$`pi*`>=0.0001) %>% length()
  ab_df$`pi* < 0.0001`[i] = which(hold$`pi*`<0.0001) %>% length()   # could turn off >0
  
  # count diferent tajima's bins
  ab_df$`D* >= 3`[i] = which(hold$`D*`>=3) %>% length()
  ab_df$`2 <= D* < 3`[i] = which(hold$`D*`<3 & hold$`D*`>=2) %>% length()
  ab_df$`1 <= D* < 2`[i] = which(hold$`D*`<2 & hold$`D*`>=1) %>% length()
  ab_df$`0 <= D* < 1`[i] = which(hold$`D*`<1 & hold$`D*`>=0) %>% length()
  ab_df$`D* < 0`[i] = which(hold$`D*`<0) %>% length()    # could add 0
  
  ab_df$polymorphic[i] = which(hold$polymorphic_yes_no == 1) %>% length()
}

# edit protein names
## adjust protein names
ab_df$ag[which(ab_df$ag=='csp')] = 'CSP'
ab_df$ag[which(ab_df$ag=='ama1')] = 'AMA1'
ab_df$ag[which(ab_df$ag=='rh5')] = 'Rh5'
ab_df$ag[which(ab_df$ag=='msp119')] = 'MSP119'
ab_df$ag[which(ab_df$ag=='pfs25')] = 'Pfs25'
ab_df$ag[which(ab_df$ag=='pfs4845')] = 'Pfs4845'
ab_df$ag[which(ab_df$ag=='pfs47')] = 'Pfs47'
ab_df$ag[which(ab_df$ag=='pfs230.seg1')] = 'Pfs230.D1,D2'

# adjust pnames
ab_df$ab_name = paste(ab_df$ag, ": ",ab_df$ab_name, sep = '')


## create null rows for ggplot spacing
ab_df[nrow(ab_df)+1,'ab_name'] = ''
ab_df[nrow(ab_df)+1,'ab_name'] = ' '
ab_df[nrow(ab_df)+1,'ab_name'] = '  '
ab_df[nrow(ab_df)+1,'ab_name'] = '   '
ab_df[nrow(ab_df)+1,'ab_name'] = '    '
ab_df[nrow(ab_df)+1,'ab_name'] = '     '
ab_df[nrow(ab_df)+1,'ab_name'] = '      '


## plot pi ----
plot_me = ab_df %>% pivot_longer(cols = contains('pi'))

# better >= sign
plot_me$name = gsub('<=', '\u2264', plot_me$name)
plot_me$name = gsub('>=', '\u2265', plot_me$name)

# swap in pi character
plot_me$name = gsub('pi', '\u03c0', plot_me$name)

# set order of color scale
plot_me$name = factor(plot_me$name,
                      levels = unique(plot_me$name))

# set order of proteins
plot_me$ab_name = factor(plot_me$ab_name,
                         levels = unique(plot_me$ab_name)[c(1:6, 55,
                                                            7, 56,
                                                            8:13, 57,
                                                            14:24, 58,
                                                            25:32, 59,
                                                            33:48, 60,
                                                            49:53, 61,
                                                            54)])

p1 = ggplot(plot_me, aes(ab_name, value, fill = name))+
  geom_col(position = 'fill', color = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title= element_blank())+
  ylab('Proportion')+xlab('Protein')+
  scale_fill_brewer(palette = 17, direction = -1)

## plot tajima ----
plot_me = ab_df %>% pivot_longer(cols = contains('D*'))

# better >= sign
plot_me$name = gsub('<=', '\u2264', plot_me$name)
plot_me$name = gsub('>=', '\u2265', plot_me$name)


# set order of color scale
plot_me$name = factor(plot_me$name,
                      levels = unique(plot_me$name))

# set order of proteins
plot_me$ab_name = factor(plot_me$ab_name,
                         levels = unique(plot_me$ab_name)[c(1:6, 55,
                                                            7, 56,
                                                            8:13, 57,
                                                            14:24, 58,
                                                            25:32, 59,
                                                            33:48, 60,
                                                            49:53, 61,
                                                            54)])

p2 = ggplot(plot_me, aes(ab_name, value, fill = name))+
  geom_col(position = 'fill', color = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title= element_blank())+
  ylab('Proportion')+xlab('Protein')+
  scale_fill_brewer(palette = 17, direction = -1)

## arrange plots -----
ggpubr::ggarrange(p1+ggtitle('Nucleotide Diversity')+xlab(''),
                  p2+ggtitle("Tajima's D"),ncol=1, align = 'v',
                  labels = c('A','B'))

## save this plot as figure
ggsave('ab_diversity_selection.png', height = 8, width = 12,
       dpi = 800)

## save csv ----
# lets add total residues, %  > 0.01 pi, % > 0.001 pi, % > 2 D, % > 1 D
ab_df = ab_df %>% mutate(
  total = `pi* >= 0.01` + `pi* < 0.0001` + `0.001 <= pi* < 0.01` + `0.0001 <= pi* < 0.001`, 
  p0.01 = `pi* >= 0.01` / total * 100,
  p0.001 = (`pi* >= 0.01` + `0.001 <= pi* < 0.01`) / total * 100,
  p2 = (`D* >= 3` + `2 <= D* < 3`) / total * 100,
  p1 = (`D* >= 3` + `2 <= D* < 3` + `1 <= D* < 2`) / total * 100,
  p_pi_lowest = `pi* < 0.0001` / total * 100,
  p_d_lowest = `D* < 0` / total * 100
)

write_csv(ab_df, '3d_epitopes_stats_bins.csv')
