library(tidyverse)
metadata <- read_tsv("inptu/meta_info4fasta.tsv")
moi <- read_tsv("input/Pf_capture1092_final.singleton.moimix_fws.tsv")

newMeta <-metadata %>% left_join(moi)

newMeta <- newMeta %>% mutate(moi = cut(fws,c(0,0.8,0.95,1)))

newMeta$moi_count <- factor(newMeta$moi, levels=c("(0.95,1]","(0.8,0.95]","(0,0.8]"), 
                            labels=c("MOI = 1", "MOI = 2", "MOI >2"))

newMetaSum <- newMeta %>% count(country,moi_count)%>%ungroup()%>%group_by(country)%>%
  mutate(moi_perc=n/sum(n))


p<-ggplot(newMetaSum, aes(country,moi_perc))+
  geom_col(aes(fill=moi_count),position="stack")+
  scale_fill_viridis_d(name="Multiplicity of infection")+
  ylab("Percentage")+xlab("Country")+
  theme_bw(18)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave("figS2.png", p, width=8, height=6)
