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
df = read_rds('combined_pi_tajima_epitope2.rds')

## To start gather binned data for 7 vaccine targets ----
# Get a list of proteins to summarize
pnames = unique(df$name)[1:7]

# Make a dataframe to store bin values
tib = tibble(protein = pnames, `pi* >= 0.01` = 0, `0.001 <= pi* < 0.01` = 0, 
             `0.0001 <= pi* < 0.001` = 0, `pi* < 0.0001` = 0,
             `D* >= 3` = 0, `2 <= D* < 3` = 0, `1 <= D* < 2` = 0,
             `0 <= D* < 1` = 0, `D* < 0` = 0
             )

# loop through protein table, count residues in each of the predifined pi and tajima bins
for(i in 1:nrow(tib)){
  #subset df for this protein only
  hold = df %>% filter(name == tib$protein[i])
  
  # count different nucleotide diversity levels
  tib$`pi* >= 0.01`[i] = which(hold$pi>=0.01) %>% length()
  tib$`0.001 <= pi* < 0.01`[i] = which(hold$pi<0.01 & hold$pi>=0.001) %>% length()
  tib$`0.0001 <= pi* < 0.001`[i] = which(hold$pi<0.001 & hold$pi>=0.0001) %>% length()
  tib$`pi* < 0.0001`[i] = which(hold$pi<0.0001 & hold$pi>0) %>% length()
  
  # count diferent tajima's bins
  tib$`D* >= 3`[i] = which(hold$tajima>=3) %>% length()
  tib$`2 <= D* < 3`[i] = which(hold$tajima<3 & hold$tajima>=2) %>% length()
  tib$`1 <= D* < 2`[i] = which(hold$tajima<2 & hold$tajima>=1) %>% length()
  tib$`0 <= D* < 1`[i] = which(hold$tajima<1 & hold$tajima>=0) %>% length()
  tib$`D* < 0`[i] = which(hold$tajima<0) %>% length()
}

## adjust protein names
ro = grep('csp|ama1', tib$protein)
tib$protein[ro] = toupper(tib$protein[ro])
tib$protein[which(tib$protein=='rh5')] = 'Rh5'
tib$protein[which(tib$protein=='pfs25')] = 'Pfs25'
tib$protein[which(tib$protein=='pfs4845')] = 'Pfs4845'
tib$protein[which(tib$protein=='pfs47')] = 'Pfs47'
tib$protein[which(tib$protein=='pfs230.1')] = 'Pfs230.D1,D2'

## Now gather data for just epitope residues ----
csp = read_rds('csp_antibody_tibble.rds')
ama1 = read_rds('AMA1_antibody_tibble.rds')
rh5 = read_rds('Rh5_antibody_tibble.rds')
pfs25 = read_rds('Pfs25_antibody_tibble.rds')
pfs4845 = read_rds('Pfs4845_antibody_tibble.rds')
pfs230.1 = read_rds('Pfs230_anitbody_tibble.rds')

# Use linear epitope of Pfs47 explained in Cruz 2022 (aa 178 - 229)
pfs47 = tibble(interface = paste(178:229, collapse='+'))

## summarize these datasets into any residue found in epitope
tib2 = tibble(protein = c('csp','ama1','rh5','pfs25','pfs4845','pfs230.1','pfs47'),
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
  hold = df %>% filter(name == tib2$protein[i], ref_pos %in% interf)
  
  # count different nucleotide diversity levels
  tib2$`pi* >= 0.01`[i] = which(hold$pi>=0.01) %>% length()
  tib2$`0.001 <= pi* < 0.01`[i] = which(hold$pi<0.01 & hold$pi>=0.001) %>% length()
  tib2$`0.0001 <= pi* < 0.001`[i] = which(hold$pi<0.001 & hold$pi>=0.0001) %>% length()
  tib2$`pi* < 0.0001`[i] = which(hold$pi<0.0001 & hold$pi>0) %>% length()
  
  # count diferent tajima's bins
  tib2$`D* >= 3`[i] = which(hold$tajima>=3) %>% length()
  tib2$`2 <= D* < 3`[i] = which(hold$tajima<3 & hold$tajima>=2) %>% length()
  tib2$`1 <= D* < 2`[i] = which(hold$tajima<2 & hold$tajima>=1) %>% length()
  tib2$`0 <= D* < 1`[i] = which(hold$tajima<1 & hold$tajima>=0) %>% length()
  tib2$`D* < 0`[i] = which(hold$tajima<0) %>% length()
}

# edit protein names
ro = grep('csp|ama1', tib2$protein)
tib2$protein[ro] = toupper(tib2$protein[ro])
tib2$protein[which(tib2$protein=='rh5')] = 'Rh5'
tib2$protein[which(tib2$protein=='pfs25')] = 'Pfs25'
tib2$protein[which(tib2$protein=='pfs4845')] = 'Pfs4845'
tib2$protein[which(tib2$protein=='pfs47')] = 'Pfs47'
tib2$protein[which(tib2$protein=='pfs230.1')] = 'Pfs230.D1,D2'

# append 'epi' to these protein names and drop interface column
tib2$protein = paste('epi: ', tib2$protein, sep = '')

tib2 = tib2 %>% select(-interface)

## Combine whole protein and epitope data ----
combined = rbind(tib,tib2)

## create 6 null rows for ggplot spacing
combined[nrow(combined)+1,'protein'] = ''
combined[nrow(combined)+1,'protein'] = ' '
combined[nrow(combined)+1,'protein'] = '  '
combined[nrow(combined)+1,'protein'] = '   '
combined[nrow(combined)+1,'protein'] = '    '
combined[nrow(combined)+1,'protein'] = '     '

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
                         levels = unique(plot_me$protein)[c(2,8,15,1,9,16,
                                                            3,10,17,4,11,18,
                                                            5,12,19,7,13,20,6,14)])

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
                         levels = unique(plot_me$protein)[c(2,8,15,1,9,16,
                                                            3,10,17,4,11,18,
                                                            5,12,19,7,13,20,6,14)])

p2 = ggplot(plot_me, aes(protein, value, fill = name))+
  geom_col(position = 'fill', color = 'black')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title= element_blank())+
  ylab('Proportion')+xlab('Protein')+
  scale_fill_brewer(palette = 17, direction = -1)

## arrange plots -----
ggarrange(p1+ggtitle('Nucleotide Diversity')+xlab(''),
          p2+ggtitle("Tajima's D"),ncol=1, align = 'v',
          labels = c('A','B'))

  ## save this plot as figure