## Just saving rds files as excel files for others to utilize
library(tidyverse)

rds = grep('.rds$', dir())

fnames = dir()[rds]

for(i in fnames){
  save_as = gsub('.rds','.csv',i)
  save_as = paste('Generated Datasets/', save_as, sep='')
  obj = read_rds(i)
  write.csv(obj, file = save_as, row.names = F)
}
