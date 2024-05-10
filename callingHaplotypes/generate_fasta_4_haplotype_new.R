library(readr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(tibble)
library(SeqArray)
library(dplyr)
library(tidyr)
library(reshape)
library(ape)
library(tidyverse)

## Output folder 
# Dir <- "path2output/pf3k_fasta_multi"  
# Create output directoy
#if(!dir.exists(Dir)) dir.create(Dir, recursive=T)
# setwd(Dir)

# set output file path
outputDir_1 = "vcf"
if(!dir.exists(outputDir_1)) dir.create(outputDir_1, recursive=T)
outputDir_2 = "fasta"
if(!dir.exists(outputDir_2)) dir.create(outputDir_2, recursive=T)
outputDir_3 = "pedmap"
if(!dir.exists(outputDir_3)) dir.create(outputDir_3, recursive=T)


#source("Naung_et_al.Haplotype_phasing.R")
source("pf3k_antigen.R")
#vcf_file = "Pf_capture1092_final.snp_filter_singleton.antigen.annotation.v5.1.vcf"
#vcf_file <- "Pf_capture1092_final.snp_filter_singleton_allAntigens_feb12024_chrR.ann.vcf"
vcf_file <- "input/Pf_capture1092_HapoCaller_qual_antigen_ann.vcf"
gds_file = "input/Pf_capture1092_HapoCaller_qual_antigen_ann.gds"
#SeqArray::seqVCF2GDS(vcf_file, gds_file)

pf3kv5 <- SeqArray::seqOpen(gds_file, readonly = TRUE)
seqResetFilter(pf3kv5)
metadata <- read.delim("input/meta_info4fasta.tsv", header=TRUE)

fws_info = "input/Pf_capture1092_final.singleton.moimix_fws.tsv"
#' ## Read moimix results
moi <- read.table(fws_info, header = T)
colnames(moi) = c("sample.id", "fws")
# modify the standard for moi=1 is 0.95s
moi1_samples <- moi$sample.id[moi$fws>0.95]
moi2_samples <- moi$sample.id[moi$fws<=0.95 & moi$fws>0.80] 
keep_samples <- c(moi1_samples, moi2_samples)
clonal_props <- data.frame(samples=as.character(moi2_samples))


#' # Parameters


BIN_SIZE <- 30 # for SNP density plots
MOI1_MIN_DEPTH <- 5 # minimum depth of coverage per site for sample with MOI 1
MOI2_MIN_DEPTH <- 10 # minimum depth of coverage per site for sample with MOI 2
HET_THRESH <- 0.2   # minimum minor allele proportion for a heterozygous call to be made
MINOR_PROP_THRESH <- 0.4 # we differentiate the major and minor haplotypes by
# read depth; the mean minor allele frequency must be below
# this threshold for the major and minor haplotypes to be
# differentiated properly

# target antigens
antigens <- read.delim("input/vaccine_candidate_antigen.tsv", header=TRUE, stringsAsFactors=FALSE)

#antigens <- read.delim("test_antigen.tsv", header=TRUE, stringsAsFactors=FALSE)
# get CDS
cds_seqs_tt<-read.dna("input/PlasmoDB-66_Pfalciparum3D7_AnnotatedCDSs.fasta",format="fasta",as.character = T)

selectedAntigen <- data.frame(RowNo = 1:length(cds_seqs_tt),label=names(cds_seqs_tt),Gene.ID=NA)%>%mutate(Gene.ID = str_split_i(label,"\\.",1))
antigens <- antigens %>% left_join(selectedAntigen)
rownames(antigens) <- antigens$Gene.ID

cds_seqs <- cds_seqs_tt[antigens$RowNo]


# generate variant annotation table
ANNinfo <- seqGetData(pf3kv5, 
                      c("chromosome","position","variant.id","allele", "annotation/info/ANN"),.padNA=T)

annfullInfo <- data.frame(chr=character(),pos=integer(),variant_id=character(),
                          ref=character(),allele=integer(),alt=character(),
                          annotation=character(),impact=character(),gene_name=character(),
                          gene_id=character(),transcript_biotype=character(),
                          nu2=character(),nu3=character(),nu4=character(),
                          HGVS.c=character(),HGVS.p=character(),cDNA=character(),CDS=character(),AA=character())
annBasicInfo <- data.frame(head(ANNinfo,n=4))
annList<- ANNinfo$`annotation/info/ANN`
annList$count <-cumsum(annList$length)

annfullInfo1<-data.frame(chr=character(),pos=integer(),variant_id=character(),
                         ref=character(),allele=integer(),alt=character())

annfullInfo2<-data.frame(chr=character(),pos=integer(),alt=character(),
                         annotation=character(),impact=character(),gene_name=character(),
                         gene_id=character(),transcript_biotype=character(),
                         nu2=character(),nu3=character(),nu4=character(),
                         HGVS.c=character(),HGVS.p=character(),cDNA=character(),CDS=character(),AA=character())

for (i in 1:length(annList$length)){
  k<-2
  basicRow <-annBasicInfo[i,1:3]
  tl <- str_split(annBasicInfo[i,4],",")[[1]]
  for(k in 2:length(tl)){
    tempRow <- c(basicRow, tl[1],k,tl[k])
    names(tempRow)<-colnames(annfullInfo1)
    annfullInfo1 <- rbind(annfullInfo1,tempRow)
  }
}


j<-1
for (i in 1:length(annList$length)){
  tl <- str_split(annBasicInfo[i,4],",")[[1]]
  while(j <= annList$count[i]){
    tempRow <-annBasicInfo[i,1:2]
    ann <- str_split(annList$data[j],"[(|)]")[[1]][1:14]
    tempRow <- c(tempRow,ann)
    names(tempRow)<-colnames(annfullInfo2)
    annfullInfo2<-rbind(annfullInfo2,tempRow)
    j<-j+1
  }
  if(("*" %in% tl)&(annList$length[i]>0)){
    ann[1]<-"*"
    if(str_detect(ann[2],"intron")){
      ann[2]<-"intron_spandeletion"
    }else{
      ann[2]<-"spanning_deletion"      
    }

    ann[3]<-"other"
    tempRow<-c(annBasicInfo[i,1:2],ann)
    names(tempRow)<-colnames(annfullInfo2)
    annfullInfo2<-rbind(annfullInfo2,tempRow)
  }
}

annfullInfo<-annfullInfo1%>%left_join(annfullInfo2,by=c("chr","pos","alt"))

#ANNinfo<-as.data.frame(ANNinfo)
#ANNinfo<-ANNinfo %>% tidyr::separate(allele,into=c("ref","alt"),sep=",",remove=T)
#colnames(ANNinfo)<-c("chr","pos","variant_id","ref","alt","ANN")
#ANNinfo<-ANNinfo %>% 
#  tidyr::separate(ANN,into=c("nu1","annotation","impact","gene_id",
#                  "feature_type","transcript_biotype","nu2","nu3","nu4",
#                  "HGVS.c","HGVS.p","cDNA","CDS","AA"), sep="[(|)]",remove=T,convert=F)
#ANNinfo<-ANNinfo %>% dplyr::select(-starts_with("nu"))
annfullInfo<-annfullInfo %>% dplyr::select(-starts_with("nu"))
#write_tsv(annfullInfo,"pf3k_functional_annotation.antigen.tsv")

annfullInfo$start<-NA
annfullInfo$end<-NA
annfullInfo$label<-NA
annfullInfo$seq1<-NA
annfullInfo$seq2<-NA

for (i in 1:nrow(annfullInfo)){
  annot<-annfullInfo$annotation[i]
  hgvs<-annfullInfo$HGVS.c[i]
  #test1<-splitAnnot(annot,hgvs)
  annfullInfo[i,17:21]<-splitAnnot(annot,hgvs)
}

#debug(splitAnnot)
# process gff file

gff<-read_tsv("input/PlasmoDB-66_Pfalciparum3D7.gff",col_names=F,comment="#")
colnames(gff)<-c("seqname","source","feature","start","end","score","strand","frame","attribute")
gff<-gff%>%tidyr::extract(attribute,regex="ID=([a-zA-Z0-9._-]+);\\S+=",
                          into=c("ID"),convert=F,remove=F)


#debug(getCDSpos)
#getCDSpos(370442,"Pf3D7_02_v3",gff)

head(annfullInfo)

#tt<- func_annotations[1:3,]
ttt<-mapply(getCDSpos,snp_pos=annfullInfo$pos,snp_chr=annfullInfo$chr,
            MoreArgs = list(gff_file = gff))
ttt<-t(ttt)
newAnnot<-as_tibble(ttt)
newAnnot<-type_convert(newAnnot)
colnames(newAnnot) <- c("chr","pos","gene_id","start","end","cds_pos","cds_length")
newAnnot<-newAnnot%>%distinct(chr,pos,cds_pos,cds_length)

annfullInfo<-annfullInfo%>%left_join(newAnnot,by=c("chr","pos"))


write_tsv(annfullInfo,"output/pf3k_functional_annotation.antigen.tsv")

func_annotations <- read.table("output/pf3k_functional_annotation.antigen.tsv", header=T, sep = "\t")

# get all regions of genes
gene_regions <- read_csv("input/vaccine_segments.csv")
#gene_regions <- gene_regions %>% left_join(antigens)
# Functional annotation
#func_annotations <- func_annotations %>% tidyr::extract(CDS, into=c("pos_CDS", "len_CDS"), regex="([0-9]+)/([0-9]+)", convert=TRUE)
gene<-antigens$Gene.ID[1]
filterAnnot <- c("frameshift_variant","frameshift_variant&stop_gained","stop_gained")
for (gene in antigens$Gene.ID) {
  cat("\n******  \n")
  cat(paste0("\n###", gene, "  \n"))
  # extract coding SNPs in gene of interest
  snps <- func_annotations %>% 
    subset(gene_id==gene & 
             ! is.na(label))%>%group_by(variant_id)%>%
    mutate(remove = ifelse(sum(annotation %in% filterAnnot)>0,T,F))%>%
    ungroup()%>%filter(remove==F)
  if (nrow(snps) > 0) {  
    # visualise the distribution of SNPs across the gene
    #gene_var <- snps%>%filter(variant_id==72)%>%distinct(chr, pos,ref,variant_id,cds_pos)
    gene_var <- snps%>%distinct(chr, pos,variant_id,cds_pos)
    
    # take a subset of data -> samples with MOI 1 or 2 + coding variants in gene 
    # of interest only
    seqSetFilter(pf3kv5, variant.id = gene_var$variant_id, sample.id = as.character(keep_samples))
    cat("  \n")
    # write data to a vcf file
    seqGDS2VCF(pf3kv5, paste0("vcf/", gene, "_coding.vcf"),
               fmt.var=c("AD", "DP", "GQ", "PL", "RGQ", "SB"))
    
    # extract haplotypes for samples with MOI 1
    #debug(moi1_sequences)
    #undebug(moi1_sequences)
    moi1_haps <- moi1_sequences(gene_var$variant_id, as.character(moi1_samples), pf3kv5, MOI1_MIN_DEPTH)
    
    # extract haplotypes for samples with MOI 2
    cat("MOI2")
    moi2_haps <- moi2_haplotypes(gene_var$variant_id, as.character(moi2_samples), pf3kv5, HET_THRESH, MOI2_MIN_DEPTH)
    moi2_seqs <- moi2_sequences(gene_var$variant_id, as.character(moi2_samples), moi2_haps, MINOR_PROP_THRESH)
    
    # combine the retained seqs
    retained_genotypes <- prep_retained_genotypes(moi1_haps,moi2_seqs,moi2_haps)
    
    # extract reference coding sequence for gene from fasta file
    #cds_seq <- read.dna(paste0("ref", gene, ".1.fasta"), format="fasta") %>% as.character %>% as.vector %>% toupper
    cds_seq <- cds_seqs[antigens[antigens$Gene.ID==gene,"label"]][[1]]%>% toupper
    
    # generate country files
    countryList <- unique(metadata$country)
    #countryFastaList <- paste0("fasta/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"],"_",countryList, "_CDS.fasta")
    #countryFastaList <- as.list(countryFastaList)
    #names(countryFastaList)<-countryList
    # generate a fasta file containing coding sequences for each clone 
    
    #debug(generate_fasta)
    #generate_fasta(cds_seq, gene_var$cds_pos, moi1_haps, moi2_seqs, 
    #               paste0("fasta/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"], "_CDS.fasta"), 
    #               antigens[gene, ]$dir,countryFastaList = countryFastaList,countryList = countryList)
    #generate_fasta(cds_seq, gene_var$cds_pos, moi1_haps, moi2_seqs, gene_var$ref,
    #               paste0("fasta/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"], "_CDS.fasta"), 
    #               antigens[gene, ]$dir)
    #generate_fasta_new(cds_seq, gene_var,snps, retained_genotypes,
    #                   paste0("fasta/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"], "_CDS.fasta"),
    #                   filterN = F)
    
    fullSeqs <- generate_final_seq_block(cds_seq, gene_var,snps, retained_genotypes)
    write_csv(as.data.frame(fullSeqs), paste0("pedmap/", gene, "_perSiteSeq.csv"))
    
    
    # write the full length cds first
    fastName <- paste0("fasta_fullCDS/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"], "_fullCDS.fasta")
    #debug(write_fasta)
    #undebug(write_fasta)
    write_fasta(fullSeqs,incSamples=NULL,1,length(cds_seq),fastName)
    for (country in countryList) {
      tempFas <- paste0("fasta_fullCDS/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"],"_",country, "_fullCDS.fasta")
      samples<- metadata$sample[metadata$country==country]
      write_fasta(fullSeqs,samples,1,length(cds_seq),tempFas)
    }
    cat("\n******  \n write regions \n")
    
    # get regions
    out_regions <- gene_regions%>%filter(Gene.ID==gene)
    for (i in 1:nrow(out_regions)){
    #for (i in 1:1){
      nstart<-out_regions$nuc_start[i]
      nend<-out_regions$nuc_end[i]
      if(nend<=nstart){
        stop("region end must be larger than region start")
      }
      domain <- out_regions$Domains[i]
      cat(paste(nstart,nend,domain,"\n"))
      #all samples
      fastName <- paste0("fasta_regions/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"], 
                        "_n",nstart,"_",nend,"_",domain,".fasta")
      write_fasta(fullSeqs,incSamples=NULL,nstart,nend,fastName)
      for (country in countryList) {
        cat(paste(country,"\n"))
        tempFas <- paste0("fasta_regions/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"],
                          "_n",nstart,"_",nend,"_",domain,"_",country,".fasta")
        samples<- metadata$sample[metadata$country==country]
        write_fasta(fullSeqs,samples,nstart,nend,tempFas)
      }
    }
    
    # keep track of clonal proportions for samples with MOI 2
    clonal_props <- cbind(moi2_seqs[["clones"]]$min_prop, clonal_props)
    colnames(clonal_props)[1] <- gene
    
    # graph a boxplot of minor clone proportions for the gene
    #boxplot(moi2_seqs[["clones"]]$min_prop, horizontal=TRUE, ylim=c(0, 0.5), xlab="Minor Clone Proportion (mean BAF across het loci)", main=paste0("Minor Clone Proportions for ", gene)) 
    
    # generate a pedmap file if there are more than 10 variants across the gene
    if (nrow(gene_var)>10) {
      generate_pedmap(gene_var, moi1_haps, moi2_haps, metadata, paste0("pedmap/", gene, "_coding.ped"), paste0("pedmap/", gene, "_coding.map"))
    }
    
    seqResetFilter(pf3kv5)
  }
}


calSites <- func_annotations%>%filter(!is.na(label), alt!="*")%>%group_by(variant_id)%>%
  summarize(keep_snpOnly=ifelse(sum(label != "SNP")>0,F,T))
func_annotations <- func_annotations%>%left_join(calSites)
# only SNPs
for (gene in antigens$Gene.ID) {
  cat("\n******  \n")
  cat(paste0("\n###", gene, "  \n"))
  # extract coding SNPs in gene of interest
  
  snps <- func_annotations %>% 
    subset(gene_id==gene & keep_snpOnly==T)
  if (nrow(snps) > 0) {  
    # visualise the distribution of SNPs across the gene
    #gene_var <- snps%>%filter(variant_id==72)%>%distinct(chr, pos,ref,variant_id,cds_pos)
    gene_var <- snps%>%distinct(chr, pos,variant_id,cds_pos)
    
    
    # take a subset of data -> samples with MOI 1 or 2 + coding variants in gene 
    # of interest only
    seqSetFilter(pf3kv5, variant.id = gene_var$variant_id, sample.id = as.character(keep_samples))
    cat("  \n")
    # write data to a vcf file
    seqGDS2VCF(pf3kv5, paste0("vcf/", gene, "_coding.vcf"),
               fmt.var=c("AD", "DP", "GQ", "PL", "RGQ", "SB"))
    
    # extract haplotypes for samples with MOI 1
    undebug(moi1_sequences)
    moi1_haps <- moi1_sequences(gene_var$variant_id, as.character(moi1_samples), pf3kv5, MOI1_MIN_DEPTH)
    
    # extract haplotypes for samples with MOI 2
    cat("MOI2")
    moi2_haps <- moi2_haplotypes(gene_var$variant_id, as.character(moi2_samples), pf3kv5, HET_THRESH, MOI2_MIN_DEPTH)
    moi2_seqs <- moi2_sequences(gene_var$variant_id, as.character(moi2_samples), moi2_haps, MINOR_PROP_THRESH)
    
    # combine the retained seqs
    retained_genotypes <- prep_retained_genotypes(moi1_haps,moi2_seqs,moi2_haps)
    
    # extract reference coding sequence for gene from fasta file
    #cds_seq <- read.dna(paste0("ref", gene, ".1.fasta"), format="fasta") %>% as.character %>% as.vector %>% toupper
    cds_seq <- cds_seqs[antigens[antigens$Gene.ID==gene,"label"]][[1]]%>% toupper
    
    # generate country files
    #countryList <- unique(metadata$country)
    #countryFastaList <- paste0("fasta/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"],"_",countryList, "_CDS.fasta")
    #countryFastaList <- as.list(countryFastaList)
    #names(countryFastaList)<-countryList
    # generate a fasta file containing coding sequences for each clone 
    
    #debug(generate_fasta)
    #generate_fasta(cds_seq, gene_var$cds_pos, moi1_haps, moi2_seqs, 
    #               paste0("fasta/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"], "_CDS.fasta"), 
    #               antigens[gene, ]$dir,countryFastaList = countryFastaList,countryList = countryList)
    #generate_fasta(cds_seq, gene_var$cds_pos, moi1_haps, moi2_seqs, gene_var$ref,
    #               paste0("fasta/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"], "_CDS.fasta"), 
    #               antigens[gene, ]$dir)
    generate_fasta_new(cds_seq, gene_var,snps, retained_genotypes,
                       paste0("fasta/", gene,"_",antigens[antigens$Gene.ID==gene,"Symbol"], "_CDS.fasta"),
                       filterN=F)
    
    
    
    # keep track of clonal proportions for samples with MOI 2
    clonal_props <- cbind(moi2_seqs[["clones"]]$min_prop, clonal_props)
    colnames(clonal_props)[1] <- gene
    
    # graph a boxplot of minor clone proportions for the gene
    boxplot(moi2_seqs[["clones"]]$min_prop, horizontal=TRUE, ylim=c(0, 0.5), xlab="Minor Clone Proportion (mean BAF across het loci)", main=paste0("Minor Clone Proportions for ", gene)) 
    
    # generate a pedmap file if there are more than 10 variants across the gene
    if (nrow(gene_var)>10) {
      generate_pedmap(gene_var, moi1_haps, moi2_haps, metadata, paste0("pedmap/", gene, "_coding.ped"), paste0("pedmap/", gene, "_coding.map"))
    }
    
    seqResetFilter(pf3kv5)
  }
}


#' Session info
sessionInfo()

