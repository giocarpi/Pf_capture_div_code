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
  gen <- seqGetData(mgen, "genotype") %>% matrix(ncol=n_var, nrow=n_sample)
  
  # pick out allele calls
  all_calls <- seqGetData(mgen, "allele") %>% strsplit(",")
  
  # pick out read depth at each site
  cov_depth <- seqGetData(mgen, "annotation/format/DP") #$data
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
  else if (toupper(base)=="N") return("N")
  else return(NA)
}
# generate a fasta file fasta_out, containing both major and minor haplotypes
generate_fasta <- function(cds_seq, cds_pos, moi1_haps, moi2_seqs, fasta_out, direction) {
  cat("", file=fasta_out, append=FALSE)
  seqs <- cbind(moi1_haps[["seqs"]], moi2_seqs[["seqs"]]) %>% apply(1:2, as.character)
  seqs[is.na(seqs)] <- "N"
  for (sample in colnames(seqs)) {
    cat(paste0(">", sample, "\n"), file=fasta_out, append=TRUE)
    sample_seq <- cds_seq
    #cat(as.character(direction), "\n")
    #for (i in 1:length(cds_pos)) sample_seq[cds_pos[i]] <- ifelse(as.character(direction)=="-", comp_base(seqs[i, sample]), seqs[i, sample])
    for (i in 1:length(cds_pos)){
        sample_seq[cds_pos[i]] <- ifelse(as.character(direction)=="-", comp_base(seqs[i, sample]), seqs[i, sample])
        if (is.na(comp_base(seqs[i, sample]))){
            cat("CDS_position NA: ", cds_pos[i], " ", sample_seq[cds_pos[i]], " ", seqs[i, sample], " ", sample, "\n")
        }
    }
    cat(paste0(paste(sample_seq, collapse=""), "\n"), file=fasta_out, append=TRUE)
  }
}

