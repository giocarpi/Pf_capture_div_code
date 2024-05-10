
## ----hap----------------------------------------------------------------------------------------------------------
# EXTRACT MAJOR AND MINOR HAPLOTYPES FOR SAMPLES WITH MOI==2
# returns a list of dataframes, one for each sample, encoding haplotypes
# over a given set of loci; maj_clone and min_clone give the allele calls
# for each clone, while maf and baf denote the major and b-allele freqs
# het_thresh is the minimum b allele frequency for a true polyclonal call
# loci with read depth < min_depth are filtered out
moi2_haplotypes <- function(var_subset, sample_subset, mgen, 
                            het_thresh, min_depth) {
  seqResetFilter(mgen)
  seqSetFilter(mgen, variant.id=var_subset, sample.id=sample_subset)
  haplotypes <- list()
  sample_hap <- data.frame(variant.id=seqGetData(mgen, "variant.id"),
                           maj_clone=NA, min_clone=NA,
                           maj_clone_bin=NA, min_clone_bin=NA, maf=NA, baf=NA,
                           het=FALSE)
  
  # initialise list of haplotypes, one df for each sample
  for (sample in sample_subset) {
    haplotypes[[sample]] <- sample_hap
  }
  
  # pick out allele depths + allele calls
  all_depths <- seqGetData(mgen, "annotation/format/AD")
  all_calls <- seqGetData(mgen, "allele") %>% strsplit(",")
  num_per_locus <- all_depths$length
  depth_per_locus <- all_depths$data
  
  variant_num <- 1
  curr <- 0
  
  # iterate over each variant
  while (variant_num <= length(num_per_locus)) {
    num <- num_per_locus[variant_num]

    # extract major/minor clones and maf/baf for each sample
    for (i in 1:length(sample_subset)) {
      depths <- depth_per_locus[i,(curr+1):(curr+num)]
      calls <- order(-depths)
      if (sum(depths)>min_depth) {
        haplotypes[[i]][variant_num, "maj_clone"] <- all_calls[[variant_num]][calls[1]]
        haplotypes[[i]][variant_num, "min_clone"] <- all_calls[[variant_num]][calls[2]]
        haplotypes[[i]][variant_num, "maj_clone_bin"] <- calls[1]
        haplotypes[[i]][variant_num, "min_clone_bin"] <- calls[2]
        haplotypes[[i]][variant_num, "maf"] <- depths[calls[1]]/sum(depths)
        haplotypes[[i]][variant_num, "baf"] <- depths[calls[2]]/sum(depths)
        if (haplotypes[[i]][variant_num, "baf"]<het_thresh) {
          haplotypes[[i]][variant_num, "min_clone"] <- haplotypes[[i]][variant_num, "maj_clone"]
          haplotypes[[i]][variant_num, "min_clone_bin"] <- haplotypes[[i]][variant_num, "maj_clone_bin"]
        } else {
          haplotypes[[i]][variant_num, "het"] <- TRUE
        }
      } else {
        haplotypes[[i]][variant_num, "min_clone_bin"] <- haplotypes[[i]][variant_num, "maj_clone_bin"] <- 0
      }
    }
    
    variant_num <- variant_num + 1
    curr <- curr+num  
  }
  
  return(haplotypes)
}

# EXTRACT SEQUENCES AND CLONAL PROPORTIONS FOR SAMPLES WITH MOI==2
# returns a list of dataframes; clones contains the major and minor
# clone proportions for each sample, while seqs contains allele sequences
# for the major and (if applicable) minor clones in each sample
# we determine the major and minor haplotypes using read depths; the
# minor clone proportion must be below min_prop_thresh for haplotypes
# to be pieced together appropriately
moi2_sequences <- function(variants, samples, haplotypes, min_prop_thresh) {
  n_samp <- length(samples)
  clones <- data.frame(row.names=samples, maj_prop=rep(1, n_samp), 
                       min_prop=rep(0, n_samp), retained=FALSE)
  seqs <- data.frame(row.names=variants)
  
  # for true het calls, we calculate mean(maf) = maj_prop and mean(baf) = minor_prop
  # and output the allele sequences for both the major and minor clones
  # if a sample has no true het calls, we output the major clone only
  for (sample in samples) {
    if (any(haplotypes[[sample]]$het)) {
      het_calls <- haplotypes[[sample]] %>% subset(het==TRUE)
      clones[sample, "maj_prop"] <- mean(het_calls$maf)
      clones[sample, "min_prop"] <- mean(het_calls$baf)
      if (clones[sample, "min_prop"]<min_prop_thresh){
        min_name <- paste0(sample, "_min")
        seqs <- cbind(haplotypes[[sample]]$maj_clone,
                      haplotypes[[sample]]$min_clone, seqs)
        colnames(seqs)[1] <- paste0(sample, "_maj")
        colnames(seqs)[2] <- paste0(sample, "_min")
        clones[sample, "retained"] <- TRUE
      }
    } else {
      seqs <- cbind(haplotypes[[sample]]$maj_clone, seqs)
      colnames(seqs)[1] <- paste0(sample, "_maj")
    }
  }
  return(list(clones=clones, seqs=seqs))
}

#==========================================================================

# EXTRACT SEQUENCES FOR SAMPLES WITH MOI==1
# returns a dataframe of allele calls for a given set of variants and
# samples, given a per-site constraint on the depth of coverage
moi1_sequences <- function(variants, samples, mgen, min_depth) {
  n_sample = length(samples)
  n_var = length(variants)
  seqResetFilter(mgen)
  seqSetFilter(mgen, variant.id=variants, sample.id=samples)
  
  seqs <- gen_bin <- setNames(data.frame(matrix(ncol=n_sample, nrow=n_var), 
                                         row.names=seqGetData(mgen, "variant.id")), 
                              seqGetData(mgen, "sample.id"))
  

  seqSetFilter(mgen, variant.id=variants, sample.id=samples)

  # extract genotypes (0=ref, 1,2... are alt alleles)
  gen1 <- seqGetData(mgen, "genotype") 
  gen <- gen1[1,,] %>% matrix(ncol=n_var, nrow=n_sample)
  
  # pick out allele calls
  all_calls <- seqGetData(mgen, "allele") %>% strsplit(",")
  
  # pick out read depth at each site
  cov_depth <- seqGetData(mgen, "annotation/format/DP")
  cov_depth[is.na(cov_depth)] <- 0
  
  # determine the allele called for each sample, given the depth
  # of coverage at any given site is at least min_depth
  for (var_num in 1:n_var) {
    for (samp_num in 1:n_sample) {
      if (cov_depth[samp_num, var_num]>=min_depth) {
         seqs[var_num, samp_num] <- all_calls[[var_num]][gen[samp_num,var_num]+1]
         gen_bin[var_num, samp_num] <- gen[samp_num,var_num]+1
      } else {
         gen_bin[var_num, samp_num] <- 0
         seqs[var_num, samp_num] <- "N"
      }
    }
  }

  return(list(seqs=seqs, genotypes=gen_bin))
}

#==========================================================================
# FORMAT CONVERSION

# generate a ped file ped_out (with family ID given by country of origin, and
# maternal ID, paternal ID and phenotype all set to 0) and a map_file for
# samples with MOI==1 or MOI==2 across a subset of variants
generate_pedmap <- function(variants, moi1_haps, moi2_haps, metadata, ped_out, map_out) {
  cat("", file=ped_out, append=FALSE)
  for (sample in names(moi2_haps)) {
    loc <- moi2_haps[[sample]] %>% transmute(loc=paste0("\t", maj_clone_bin, "\t", min_clone_bin))
    cat(paste0(metadata[sample,]$country, "\t", sample, "\t", "0\t0\t2\t0", 
               paste(loc$loc, collapse=""), "\n"), file=ped_out, append=TRUE)
  }
  for (sample in colnames(moi1_haps[["genotypes"]])) {
    cat(paste0(metadata[sample,]$country, "\t", sample, "\t", "0\t0\t1\t0\t", 
               paste(rep(moi1_haps[["genotypes"]][, sample], each=2), collapse="\t"), 
               "\n"), file=ped_out, append=TRUE)
  }
  variants %>% mutate(cM=pos/15000) %>% select(chr, variant_id, cM, pos) %>% 
    write.table(file=map_out, row.names = FALSE, col.names = FALSE, sep="\t", quote=FALSE)
}

# we take the complementary base to the allele called in the VCE when annotating the CDS fasta file
comp_base <- function(base) {
  if (toupper(base)=="A") return("T")
  else if (toupper(base)=="T") return("A")
  else if (toupper(base)=="C") return("G")
  else if (toupper(base)=="G") return("C")
  else if (base=="-") return("-")
  else return("N")
}

getCDSpos<-function(snp_pos,snp_chr,gff_file){
  #print(paste(snp_chr,snp_pos))
  geneRegion<-gff_file%>%filter(seqname==snp_chr,feature=="protein_coding_gene")
  # get the whole gene range
  for (i in 1:nrow(geneRegion)){
    if (snp_pos>=geneRegion[i,"start"] & snp_pos<=geneRegion[i,"end"]){
      gene_start <- as.integer(geneRegion[i,"start"])
      gene_end <- as.integer(geneRegion[i,"end"])
      gene_ID <- as.character(geneRegion[i,"ID"])
      break
    }
  }

  subRegion<- gff_file%>%filter(seqname==snp_chr,feature=="CDS",
                                start>=gene_start,
                                end<=gene_end)
  
  subRegion<-subRegion%>%mutate(dur = end-start+1)
  tl<-sum(subRegion$dur)
  cum_len<-0
  cds_pos<-0
  for (i in 1:nrow(subRegion)){
    if (snp_pos>=as.integer(subRegion[i,"start"]) & snp_pos<=as.integer(subRegion[i,"end"])){
      cds_pos<-cum_len+(snp_pos-as.integer(subRegion[i,"start"])+1)
      if (subRegion[i,"strand"]=="-") {
        cds_pos<-tl-cds_pos+1
      }
      break 
    }else{
      cum_len <- cum_len + as.integer(subRegion[i,"dur"])
    }
    
    
  }
  return(c(snp_chr,snp_pos,gene_ID,gene_start,gene_end,cds_pos,tl))
}


# generate a fasta file fasta_out, containing both major and minor haplotypes
generate_fasta <- function(cds_seq, cds_pos, moi1_haps, moi2_seqs, ref,fasta_out, direction,
                           include3D7=T,countryFastaList=NULL,countryList=NULL) {
  cat("", file=fasta_out, append=FALSE)
  #if(include3D7){
  #  cat(paste0(">pf3d7\n"), file=fasta_out, append=TRUE)
  #  cat(paste0(paste(cds_seq, collapse=""), "\n"), file=fasta_out, append=TRUE)
  #}
  if(length(countryFastaList)>0){
    for (cntry in  countryList){
      cat("", file=countryFastaList[[cntry]], append=FALSE)
    }
  }
  seqs <- cbind(moi1_haps[["seqs"]], moi2_seqs[["seqs"]]) %>% apply(1:2, as.character)
  #seqs<-moi1_haps[["seqs"]]%>% apply(1:2, as.character)
  seqs[is.na(seqs)] <- "N"
  seqs<-seqs[,colSums(seqs=="N")==0]
  seqs[seqs=="*"] <-"-"
  seqs<-cbind(rep("N",nrow(seqs)),seqs)
  colnames(seqs)[1]<-"pf3d7"
  for (i in 1:length(cds_pos)) {
    seqs[i,1] <- cds_seq[cds_pos[i]]
    if(as.character(direction)=="NEGATIVE"){
      seqs[i,1]<-comp_base(seqs[i,1])
    }
  }
  
  for (i in 1:nrow(seqs)){
    print(i)
    if(len(ref[i])==1){
      padLen<-max(nchar(seqs[i,]))
      for(j in 1:ncol(seqs)){
        seqs[i,j]<-paste0(seqs[i,j],paste(rep("-",padLen-nchar(seqs[i,j])),collapse=""))
        if(as.character(direction)=="NEGATIVE"){
          seqs[i,j]<-comp_base(seqs[i,j])
        }
      }
    }else{
      
    }
  }

  for (sample in colnames(seqs)) {
    sample_seq <- cds_seq
    country <- metadata[metadata$sample==str_split_i(sample,"_",1),"country"]
    #print(c(country,sample))

    for (i in 1:length(cds_pos)) {
      sample_seq[cds_pos[i]] <- seqs[i, sample]
    }
    # output everything instead of filtering samples with N
    #if (sum(sample_seq=="N")==0){
      cat(paste0(">", sample, "\n"), file=fasta_out, append=TRUE)
      cat(paste0(paste(sample_seq, collapse=""), "\n"), file=fasta_out, append=TRUE)
      # add to country fasta if required
      if(length(countryFastaList)>0){
        
        cat(paste0(">", sample, "\n"), file=countryFastaList[[country]], append=TRUE)
        cat(paste0(paste(sample_seq, collapse=""), "\n"), file=countryFastaList[[country]], append=TRUE)
      }
    #}
  }
}


splitAnnot <- function(annot,hgvs){
  patdel<-"c\\.(\\d+)(_(\\d+))?(del|ins|dup)([ATCG]+)?"
  patvar <-"c\\.(\\d+)([ATCG])>([ATCG])"
  if(is.na(annot)|annot=="spanning_deletion"){
    start<-NA
    end<-NA
    label<-"spanning_deletion"
    seq1<-NA
    seq2<-NA
    
  }else if(str_detect(annot,"intron")){
    start<-NA
    end<-NA
    label<-NA
    seq1<-NA
    seq2<-NA    
  }else if(str_detect(annot,"deletion|insertion|frame")){
    if(str_detect(hgvs,patdel)){
      start<-as.integer(str_extract(hgvs,patdel,group=1))
      end<-as.integer(str_extract(hgvs,patdel,group=3))
      label<-str_extract(hgvs,patdel,group=4)
      seq1<-str_extract(hgvs,patdel,group=5)
      seq2<-seq1
      
    }else{
      print(paste(hgvs,annot))
    }
  }else if(str_detect(annot,"missense|synonymous|stop")){
    if(str_detect(hgvs,patvar)){
      start<-as.integer(str_extract(hgvs,patvar,group=1))
      end<-start
      label<-"SNP"
      seq1<-str_extract(hgvs,patvar,group=2)
      seq2<-str_extract(hgvs,patvar,group=3)
    }else{
      print(paste(hgvs,annot))
    }
  }else{
    start<-NA
    end<-NA
    label<-NA
    seq1<-NA
    seq2<-NA
  }
  return(list(start,end,label,seq1,seq2))
  
}


prep_retained_genotypes <- function(moi1_haps, moi2_seqs,moi2_haps){
  genotypes_moi2<-as.data.frame(lapply(moi2_haps,function(x) x$maj_clone_bin))
  colnames(genotypes_moi2)<-paste0(colnames(genotypes_moi2),"_maj")
  genotypes_moi2min<-as.data.frame(lapply(moi2_haps,function(x) x$min_clone_bin))
  genotypes_moi2min<-genotypes_moi2min[,moi2_seqs$clones$retained]
  colnames(genotypes_moi2min)<-paste0(colnames(genotypes_moi2min),"_min")
  return(cbind(moi1_haps$genotypes,genotypes_moi2,genotypes_moi2min))
  
}


generate_final_seq_block <- function(cds_seq, gene_var, snps, genotypes,
                           include3D7=T,filterN = F) {
  
  
  outMat<-matrix(cds_seq,nrow=length(cds_seq),ncol=ncol(genotypes))
  colnames(outMat)<-colnames(genotypes)
  
  for (i in 1:nrow(genotypes)){
    variant.id<-gene_var$variant_id[i]
    cds_pos<-gene_var$cds_pos[i]
    for(j in 1:ncol(genotypes)) {
      al<-genotypes[i,j]
      if(is.na(al)|al==0){
        outMat[cds_pos,j]<-"N"
      }else{
        if (al>1){
          #print(paste(i,j))
          siteInfo<-snps%>%filter(variant_id==variant.id,allele==al)%>%
            dplyr::select(start, end, label, seq1,seq2)%>%as.vector()
          if(siteInfo$label=="SNP"){
            outMat[siteInfo$start,j]<-siteInfo$seq2
          }else if(siteInfo$label=="ins"){
            outMat[siteInfo$start,j]<-paste0(outMat[siteInfo$start,j],siteInfo$seq1)
          }else if(siteInfo$label=="dup"){
            outMat[siteInfo$end,j]<-paste0(outMat[siteInfo$end,j],siteInfo$seq1)
          }else if(siteInfo$label=="del"){
            outMat[siteInfo$start:siteInfo$end,j]<-"-"
          }else{
            #print(paste(i,j,siteInfo$label))
          }
        }
      }
      }
  }
  
  # add pf3D7
  if(include3D7){
    outMat<-cbind(cds_seq,outMat)
    colnames(outMat)[1]<-"Pf3D7"
  }
  
  if(filterN){
    rsN<- (colSums(outMat=="N")==0)
    outMat <- outMat[,rsN]
  }
  
  # add -- to inserted regions
  padLenSites<-apply(nchar(outMat),1,max,na.rm=T)
  for(i in 1:nrow(outMat)){
    if(padLenSites[i]>1){
      for(j in 1:ncol(outMat)){
        outMat[i,j]<-paste0(outMat[i,j],paste(rep("-",padLenSites[i]-nchar(outMat[i,j])),collapse=""))
      }
      
      
    }
   }
   
  return(outMat)
}



write_fasta <- function(seqblock, incSamples=NULL, incRegionStart, incRegionEnd, fasta_out, filterN=T){
  if (is.null(colnames(seqblock))){
    stop("columns need to have sample names")
  }
  
  outMat <- seqblock[incRegionStart:incRegionEnd,]
  
  if(! is.null(incSamples)){
    colind<- (str_split_i(colnames(seqblock),"_",1) %in% incSamples)
    outMat <- outMat[,colind]
    
  }

  if(filterN){
    rsN<- (colSums(outMat=="N")==0)
    outMat <- outMat[,rsN]
  }
  
  
  cat("", file=fasta_out, append=FALSE)
  for(j in 1:ncol(outMat)) {
    
    cat(paste0(">", colnames(outMat)[j], "\n"), file=fasta_out, append=TRUE)
    cat(paste0(paste(outMat[,j], collapse=""), "\n"), file=fasta_out, append=TRUE)
    
  }
  
}



