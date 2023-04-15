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

## Output folder 
Dir <- "path2output/pf3k_fasta_multi"  
# Create output directoy
if(!dir.exists(Dir)) dir.create(Dir, recursive=T)
setwd(Dir)

# set output file path
outputDir_1 = "vcf"
if(!dir.exists(outputDir_1)) dir.create(outputDir_1, recursive=T)
outputDir_2 = "fasta"
if(!dir.exists(outputDir_2)) dir.create(outputDir_2, recursive=T)
outputDir_3 = "pedmap"
if(!dir.exists(outputDir_3)) dir.create(outputDir_3, recursive=T)


source("Naung_et_al.Haplotype_phasing.R")

vcf_file = "Pf_capture1092_final.snp_filter_singleton.antigen.annotation.v5.1.vcf"
gds_file = "Pf_capture1092_final.snp_filter_singleton.antigen.annotation.gds"
SeqArray::seqVCF2GDS(vcf_file, gds_file)

pf3kv5 <- SeqArray::seqOpen(gds_file, readonly = TRUE)
metadata <- read.delim("meta_info4fasta.tsv", header=TRUE)

fws_info = "Pf_capture1092_final.singleton.moimix_fws.tsv"
#' ## Read moimix results
moi <- read.table(fws_info, header = T)
colnames(moi) = c("sample.id", "fws")
moi1_samples <- moi$sample.id[moi$fws>0.90]
moi2_samples <- moi$sample.id[moi$fws<=0.90 & moi$fws>0.80] 
keep_samples <- c(moi1_samples, moi2_samples)


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
antigens <- read.delim("../vaccine_candidate_antigen.tsv", header=TRUE, stringsAsFactors=FALSE)
#antigens <- read.delim("test_antigen.tsv", header=TRUE, stringsAsFactors=FALSE)


rownames(antigens) <- antigens$Gene.ID


clonal_props <- data.frame(samples=as.character(moi2_samples))
func_annotations <- read.table("pf3k_functional_annotation.antigen.tsv", header=T, sep = "\t")
# Functional annotation

for (gene in antigens$Gene.ID) {
  cat("\n******  \n")
  cat(paste0("\n###", gene, "  \n"))
  
  # extract coding SNPs in gene of interest
  snps <- func_annotations %>% 
    subset(gene_id==gene & 
             (annotation=="missense_variant" | annotation=="synonymous_variant" |
                annotation=="missense_variant&splice_region_variant" |
                annotation=="splice_region_variant&synonymous_variant" |
                annotation=="stop_gained" | annotation=="stop_lost" |
                annotation=="stop_lost&splice_region_variant" |
                annotation=="start_lost" | annotation=="stop_gained&splice_region_variant" |
                annotation=="splice_region_variant&stop_retained_variant" |
                annotation=="stop_retained_variant"))
  if (nrow(snps) > 0) {  
    # visualise the distribution of SNPs across the gene
    gene_var <- snps %>% select(chr, variant_id, pos, alt, CDS) %>% distinct %>%
      group_by(chr, variant_id, pos, CDS) %>% summarise(alleles=(n()+1)) %>% 
      extract(CDS, into=c("pos_CDS", "len_CDS"), regex="([0-9]+)/([0-9]+)", convert=TRUE) %>%
      as.data.frame
    snp_bins <- snps %>% mutate(bins=ceiling((pos-antigens[gene, "start"])/BIN_SIZE)) %>%
      group_by(bins) %>% summarise(count=n())
    snp_pos <- ggplot() + 
      geom_segment(data=gene_var, aes(x=pos, xend=pos, y=0, yend=1, 
                                      colour=as.character(alleles)), alpha=0.8) +
      scale_colour_manual(values=c("red", "green", "blue"), 
                          breaks=c("2", "3", "4"), name="Number of alleles") +
      xlab("Position (bp)") + xlim(c(antigens[gene, "start"], antigens[gene, "end"])) +
      ggtitle(paste0("SNP Loci for ", gene)) + theme_classic() + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title.y = element_blank(), axis.line.y = element_blank(),
            plot.title=element_text(hjust=0.5, face="bold"))
    # summary of coding SNPs across gene
    cat(paste0("Number of SNPs: ", nrow(gene_var)))
    cat("\n")
    print(kable(snps, caption = paste0("Coding SNPs for ", gene)) %>% 
            kable_styling() %>% scroll_box(width="100%", height="300px"))
    cat("\n  ")
    print(snp_pos)
    
    # take a subset of data -> samples with MOI 1 or 2 + coding variants in gene 
    # of interest only
    seqSetFilter(pf3kv5, variant.id = gene_var$variant_id, sample.id = as.character(keep_samples))
    cat("  \n")
    # write data to a vcf file
    seqGDS2VCF(pf3kv5, paste0("vcf/", gene, "_coding.vcf"),
               fmt.var=c("AD", "DP", "GQ", "PL", "RGQ", "SB"))
    
    # extract haplotypes for samples with MOI 1
    moi1_haps <- moi1_sequences(gene_var$variant_id, as.character(moi1_samples), pf3kv5, MOI1_MIN_DEPTH)
    
    # extract haplotypes for samples with MOI 2
    cat("MOI2")
    moi2_haps <- moi2_haplotypes(gene_var$variant_id, as.character(moi2_samples), pf3kv5, HET_THRESH, MOI2_MIN_DEPTH)
    moi2_seqs <- moi2_sequences(gene_var$variant_id, as.character(moi2_samples), moi2_haps, MINOR_PROP_THRESH)
    
    # extract reference coding sequence for gene from fasta file
    cds_seq <- read.dna(paste0("ref", gene, ".1.fasta"), format="fasta") %>% as.character %>% as.vector %>% toupper
    # generate a fasta file containing coding sequences for each clone 
    generate_fasta(cds_seq, gene_var$pos_CDS, moi1_haps, moi2_seqs, paste0("fasta/", gene, "_CDS.fasta"), antigens[gene, ]$dir)
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

