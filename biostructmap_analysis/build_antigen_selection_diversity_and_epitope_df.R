# --------------------------------------------- #
# combine all antigens into one table
# pname | resno | exposed_yn | pi* | D* | #ab structures | polymorphic
# --------------------------------------------- #

pnames = c('csp', 'ama1', 'rh5', 'celtos', 
           'pfs25', 'pfs28', 'pfs4845', 'pfs47',
           'pfs230', 'msp119')

library(tidyverse)
library(seqinr)

# load csp ----
st = 309

pi = read_csv('tajima_d_models/csp_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/csp_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/csp_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
csp = inner_join(taj, pi, by = 'reference_residue')

# reset indices
csp$reference_residue = csp$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_0304600_CSP_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
csp = right_join(seq, csp, by = 'reference_residue')

# load ama1 ----
st = 100

pi = read_csv('tajima_d_models/ama1_no4g2_short_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/ama1_no4g2_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/ama1_no4g2_short_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
ama1 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
ama1$reference_residue = ama1$reference_residue + (st - 1)
unique(ama1$reference_residue)

# indices need adjustments
ama1$reference_residue[which(ama1$reference_residue > 350)] = ama1$reference_residue[which(ama1$reference_residue > 350)] + 40

# fill in these missing positions now
tib = tibble(reference_residue = 351:390, `D*` = NA, `pi*` = NA, exposed_yes_no = NA)
ama1 = rbind(ama1, tib) %>% arrange(reference_residue)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_1133400_AMA1_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
ama1 = right_join(seq, ama1, by = 'reference_residue')


# load rh5 ----
st = 146

pi = read_csv('tajima_d_models/rh5_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/rh5_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/rh5_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
rh5 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
rh5$reference_residue = rh5$reference_residue + (st - 1)

# indices need adjustments
rh5$reference_residue[which(rh5$reference_residue > 257)] = rh5$reference_residue[which(rh5$reference_residue > 257)] + 36

# add back missing residues
tib = tibble(reference_residue = 258:293, `D*` = NA, `pi*` = NA, exposed_yes_no = NA)
rh5 = rbind(rh5, tib) %>% arrange(reference_residue)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_0424100_RH5_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
rh5 = right_join(seq, rh5, by = 'reference_residue')


# load celtos ----
st = 35

pi = read_csv('tajima_d_models/celtos_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/celtos_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/celtos_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
celtos = inner_join(taj, pi, by = 'reference_residue')

# reset indices
celtos$reference_residue = celtos$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_1216600_CelTOS_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
celtos = right_join(seq, celtos, by = 'reference_residue')


# load pfs25 ----
st = 22

pi = read_csv('tajima_d_models/pfs25_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/pfs25_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/pfs25_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
pfs25 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
pfs25$reference_residue = pfs25$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_1031000_Pfs25_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
pfs25 = right_join(seq, pfs25, by = 'reference_residue')


# load pfs28 ----
st = 24

pi = read_csv('tajima_d_models/pfs28_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/pfs28_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/pfs28_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
pfs28 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
pfs28$reference_residue = pfs28$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_1030900_Pfs28_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
pfs28 = right_join(seq, pfs28, by = 'reference_residue')

# load pfs4845 ----
st = 45

pi = read_csv('tajima_d_models/pfs4845_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/pfs4845_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/pfs4845_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
pfs4845 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
pfs4845$reference_residue = pfs4845$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_1346700_Pfs4548_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
pfs4845 = right_join(seq, pfs4845, by = 'reference_residue')


# load pfs47_long ----
st = 22

pi = read_csv('tajima_d_models/pfs47_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/pfs47_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/pfs47_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
pfs47_long = inner_join(taj, pi, by = 'reference_residue')

# reset indices
pfs47_long$reference_residue = pfs47_long$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_1346800_Pfs47_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
pfs47_long = right_join(seq, pfs47_long, by = 'reference_residue')

# load pfs47_short ----
#st = 31

#pi = read_csv('tajima_d_models/pfs47_short_trim_pi.csv')
#taj = read_csv('tajima_d_models/pfs47_short_trim_tajima.csv')

# exposed or not from pi
#pi$score = as.numeric(pi$score)
#taj$score = as.numeric(taj$score)

#pi$exposed_yes_no = !is.na(pi$score)

# update colnames
#colnames(pi)[3] = 'pi*'
#colnames(taj)[3] = 'D*'

# drop chain column
#pi = pi %>% select(-chain)
#taj = taj %>% select(-chain)

# join together
#pfs47_short = inner_join(taj, pi, by = 'reference_residue')

# reset indices
#pfs47_short$reference_residue = pfs47_short$reference_residue + (st - 1)

# add aa
#fa = phylotools::read.fasta('field_sample_fastas/PF3D7_1346800_Pfs47_CDS.fasta')
#seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
#seq$reference_residue = 1:nrow(seq)
#colnames(seq)[1] = 'ref_3d7_AA'
#pfs47_short = right_join(seq, pfs47_short, by = 'reference_residue')


# load pfs230.seg1 ----
st = 552

pi = read_csv('tajima_d_models/pfs230.1_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/pfs230.1_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/pfs230.1_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
pfs230.seg1 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
pfs230.seg1$reference_residue = pfs230.seg1$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_0209000_Pfs230_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
pfs230.seg1 = right_join(seq, pfs230.seg1, by = 'reference_residue')


# load pfs230.seg3 ----
st = 915

pi = read_csv('tajima_d_models/pfs230.3_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/pfs230.3_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/pfs230.3_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
pfs230.seg3 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
pfs230.seg3$reference_residue = pfs230.seg3$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_0209000_Pfs230_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
pfs230.seg3 = right_join(seq, pfs230.seg3, by = 'reference_residue')

# load pfs230.seg5 ----
st = 1282

pi = read_csv('tajima_d_models/pfs230.5_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/pfs230.5_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/pfs230.5_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
pfs230.seg5 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
pfs230.seg5$reference_residue = pfs230.seg5$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_0209000_Pfs230_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
pfs230.seg5 = right_join(seq, pfs230.seg5, by = 'reference_residue')


# load pfs230.seg7 ----
st = 1691

pi = read_csv('tajima_d_models/pfs230.7_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/pfs230.7_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/pfs230.7_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
pfs230.seg7 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
pfs230.seg7$reference_residue = pfs230.seg7$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_0209000_Pfs230_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
pfs230.seg7 = right_join(seq, pfs230.seg7, by = 'reference_residue')


# load pfs230.seg9 ----
st = 2049

pi = read_csv('tajima_d_models/pfs230.9_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/pfs230.9_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/pfs230.9_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
pfs230.seg9 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
pfs230.seg9$reference_residue = pfs230.seg9$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_0209000_Pfs230_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
pfs230.seg9 = right_join(seq, pfs230.seg9, by = 'reference_residue')

# load pfs230.seg11 ----
st = 2445

pi = read_csv('tajima_d_models/pfs230.11_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/pfs230.11_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/pfs230.11_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
pfs230.seg11 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
pfs230.seg11$reference_residue = pfs230.seg11$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_0209000_Pfs230_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
pfs230.seg11 = right_join(seq, pfs230.seg11, by = 'reference_residue')

# load pfs230.seg13 ----
st = 2828

pi = read_csv('tajima_d_models/pfs230.13_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/pfs230.13_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/pfs230.13_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
pfs230.seg13 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
pfs230.seg13$reference_residue = pfs230.seg13$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_0209000_Pfs230_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
pfs230.seg13 = right_join(seq, pfs230.seg13, by = 'reference_residue')


# load msp119 ----
st = 1608

pi = read_csv('tajima_d_models/msp119_norsa_trim_pi.csv')
taj = read_csv('tajima_d_models/msp119_norsa_trim_tajima.csv')

nof = read_csv('tajima_d_models/rsa_filtered/msp119_trim_pi.csv')
nof$score = as.numeric(nof$score)
nof$exposed_yes_no = !is.na(nof$score)
nof = nof %>% select(-score, -chain)

pi = left_join(pi, nof, by = 'reference_residue')

# exposed or not from pi
pi$score = as.numeric(pi$score)
taj$score = as.numeric(taj$score)

# update colnames
colnames(pi)[3] = 'pi*'
colnames(taj)[3] = 'D*'

# drop chain column
pi = pi %>% select(-chain)
taj = taj %>% select(-chain)

# join together
msp119 = inner_join(taj, pi, by = 'reference_residue')

# reset indices
msp119$reference_residue = msp119$reference_residue + (st - 1)

# add aa
fa = phylotools::read.fasta('field_sample_fastas/PF3D7_0930300_MSP119_CDS.fasta')
seq = fa$seq.text[1] %>% s2c() %>% seqinr::translate() %>% as_tibble()
seq$reference_residue = 1:nrow(seq)
colnames(seq)[1] = 'ref_3d7_AA'
msp119 = right_join(seq, msp119, by = 'reference_residue')

# load in antibody data ----
csp_anti = read_rds('collapsed_csp_antibody_tibble.rds')
ama1_anti = read_rds('collapsed_ama1_antibody_tibble.rds')
rh5_anti = read_rds('collapsed_rh5_antibody_tibble.rds')
pfs25_anti = read_rds('collapsed_pfs25_antibody_tibble.rds')
pfs4845_anti = read_rds('collapsed_pfs4845_antibody_tibble.rds')
pfs230_anti = read_rds('collapsed_pfs230_antibody_tibble.rds')
pfs47_anti = read_rds('pfs47_antibody_tibble.rds')
msp119_anti = read_rds('collapsed_msp119_antibody_tibble.rds')

# now get all residues that were every targetted -- csp
x = paste(csp_anti$interface, collapse = '+')
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
csp$epitope_yes_no = FALSE
csp$epitope_yes_no[which(csp$reference_residue %in% x)] = T

# now get all residues that were every targetted -- ama1
x = paste(ama1_anti$interface, collapse = '+')
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
ama1$epitope_yes_no = FALSE
ama1$epitope_yes_no[which(ama1$reference_residue %in% x)] = T

# now get all residues that were every targetted -- rh5
x = paste(rh5_anti$interface, collapse = '+')
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
rh5$epitope_yes_no = FALSE
rh5$epitope_yes_no[which(rh5$reference_residue %in% x)] = T

# now get all residues that were every targetted -- pfs25
x = paste(pfs25_anti$interface, collapse = '+')
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs25$epitope_yes_no = FALSE
pfs25$epitope_yes_no[which(pfs25$reference_residue %in% x)] = T

# now get all residues that were every targetted -- pfs4845
x = paste(pfs4845_anti$interface, collapse = '+')
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs4845$epitope_yes_no = FALSE
pfs4845$epitope_yes_no[which(pfs4845$reference_residue %in% x)] = T

# now get all residues that were every targetted -- pfs230
x = paste(pfs230_anti$interface, collapse = '+')
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs230.seg1$epitope_yes_no = FALSE
pfs230.seg1$epitope_yes_no[which(pfs230.seg1$reference_residue %in% x)] = T

# now get all residues that were every targetted -- pfs47_long
x = paste(pfs47_anti$interface, collapse = '+')
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs47_long$epitope_yes_no = FALSE
pfs47_long$epitope_yes_no[which(pfs47_long$reference_residue %in% x)] = T

# now get all residues that were every targetted -- pfs47_short
#x = paste(pfs47_anti$interface, collapse = '+')
#x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
#pfs47_short$epitope_yes_no = FALSE
#pfs47_short$epitope_yes_no[which(pfs47_short$reference_residue %in% x)] = T

# now get all residues that were every targetted -- msp119
x = paste(msp119_anti$interface, collapse = '+')
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
msp119$epitope_yes_no = FALSE
msp119$epitope_yes_no[which(msp119$reference_residue %in% x)] = T

# just put false for pfs28 and celtos and rest of pfs230
celtos$epitope_yes_no = F
pfs28$epitope_yes_no = F
pfs230.seg3$epitope_yes_no = F
pfs230.seg5$epitope_yes_no = F
pfs230.seg7$epitope_yes_no = F
pfs230.seg9$epitope_yes_no = F
pfs230.seg11$epitope_yes_no = F
pfs230.seg13$epitope_yes_no = F

# load in polymorphic data ----
snps = read_rds('snp_loc_df.rds')

# now get all residues that were polymorphic -- csp
x = snps %>% filter(ag == 'csp') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
csp$polymorphic_yes_no = FALSE
csp$polymorphic_yes_no[which(csp$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- ama1
x = snps %>% filter(ag == 'ama1') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
ama1$polymorphic_yes_no = FALSE
ama1$polymorphic_yes_no[which(ama1$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- rh5
x = snps %>% filter(ag == 'rh5') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
rh5$polymorphic_yes_no = FALSE
rh5$polymorphic_yes_no[which(rh5$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- pfs25
x = snps %>% filter(ag == 'pfs25') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs25$polymorphic_yes_no = FALSE
pfs25$polymorphic_yes_no[which(pfs25$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- pfs4845
x = snps %>% filter(ag == 'pfs4845') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs4845$polymorphic_yes_no = FALSE
pfs4845$polymorphic_yes_no[which(pfs4845$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- pfs47_long
x = snps %>% filter(ag == 'pfs47') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs47_long$polymorphic_yes_no = FALSE
pfs47_long$polymorphic_yes_no[which(pfs47_long$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- msp119
x = snps %>% filter(ag == 'msp119') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
msp119$polymorphic_yes_no = FALSE
msp119$polymorphic_yes_no[which(msp119$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- celtos
x = snps %>% filter(ag == 'celtos') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
celtos$polymorphic_yes_no = FALSE
celtos$polymorphic_yes_no[which(celtos$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- pfs28
x = snps %>% filter(ag == 'pfs28') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs28$polymorphic_yes_no = FALSE
pfs28$polymorphic_yes_no[which(pfs28$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- pfs47_short
#x = snps %>% filter(ag == 'pfs47') %>% select(aa_pos) %>% unlist()
#x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
#pfs47_short$polymorphic_yes_no = FALSE
#pfs47_short$polymorphic_yes_no[which(pfs47_short$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- pfs230
# first lets paste all pfs230 segments together

# now get all residues that were polymorphic -- seg1
x = snps %>% filter(ag == 'pfs230.1') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs230.seg1$polymorphic_yes_no = FALSE
pfs230.seg1$polymorphic_yes_no[which(pfs230.seg1$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- seg3
x = snps %>% filter(ag == 'pfs230.3') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs230.seg3$polymorphic_yes_no = FALSE
pfs230.seg3$polymorphic_yes_no[which(pfs230.seg3$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- seg5
x = snps %>% filter(ag == 'pfs230.5') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs230.seg5$polymorphic_yes_no = FALSE
pfs230.seg5$polymorphic_yes_no[which(pfs230.seg5$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- seg7
x = snps %>% filter(ag == 'pfs230.7') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs230.seg7$polymorphic_yes_no = FALSE
pfs230.seg7$polymorphic_yes_no[which(pfs230.seg7$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- seg9
x = snps %>% filter(ag == 'pfs230.9') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs230.seg9$polymorphic_yes_no = FALSE
pfs230.seg9$polymorphic_yes_no[which(pfs230.seg9$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- seg11
x = snps %>% filter(ag == 'pfs230.11') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs230.seg11$polymorphic_yes_no = FALSE
pfs230.seg11$polymorphic_yes_no[which(pfs230.seg11$reference_residue %in% x)] = T

# now get all residues that were polymorphic -- seg13
x = snps %>% filter(ag == 'pfs230.13') %>% select(aa_pos) %>% unlist()
x = str_split(x, '\\+')[[1]] %>% unique() %>% as.numeric() %>% sort()
pfs230.seg13$polymorphic_yes_no = FALSE
pfs230.seg13$polymorphic_yes_no[which(pfs230.seg13$reference_residue %in% x)] = T

# lets add segment also
pfs230.seg1$pname = 'pfs230.seg1'
pfs230.seg3$pname = 'pfs230.seg3'
pfs230.seg5$pname = 'pfs230.seg5'
pfs230.seg7$pname = 'pfs230.seg7'
pfs230.seg9$pname = 'pfs230.seg9'
pfs230.seg11$pname = 'pfs230.seg11'
pfs230.seg13$pname = 'pfs230.seg13'

pfs230_all = rbind(pfs230.seg1, pfs230.seg3, pfs230.seg5,
                   pfs230.seg7, pfs230.seg9, pfs230.seg11,  pfs230.seg13)

# merge into one final table ----
csp$pname = 'csp'
ama1$pname = 'ama1'
rh5$pname = 'rh5'
pfs25$pname = 'pfs25'
pfs4845$pname = 'pfs4845'
pfs47_long$pname = 'pfs47'
msp119$pname = 'msp119'
pfs28$pname = 'pfs28'
celtos$pname = 'celtos'

# pfs230 already done
antigen_df = rbind(csp, ama1, rh5, pfs25, pfs4845, 
                   pfs47_long, msp119, pfs28, celtos, pfs230_all)

antigen_df = antigen_df %>% select(pname, reference_residue, ref_3d7_AA, everything())

write_csv(antigen_df, 'antigen_selection_diversity_and_epitope.csv')

table(antigen_df$epitope_yes_no, antigen_df$polymorphic_yes_no, antigen_df$pname)

