#two different huge functions have written, and results are extracted by number which acts like mini functions
#vaxpact_input mainly does the job of accepting input and core calculation.
#vaxpack_output does extract the results, and carry out graphical output.
vaxpack_output_params <- function (cutoff_maf_pct = 5, win_size = 50, step_size = 5) {

  if (exists("vp.data") != TRUE)
    stop("Please run the input function 'vaxpack_input()' first")
#assign everything needed into vp.data
  vp.data <- vp.data
  vp.AA.VARIANT.TABLE <- vp.data[[1]]
  vp.AAs.EACH.VARIABLE.CODON <- vp.data[[2]]
  vp.REF.EACH.VARIABLE.CODON <- vp.data[[3]]
  vp.BASIC.RESULTS <- vp.data[[4]] #table 1
  vp.MAF.NUC.TABLE <- vp.data[[5]]
  vp.MAF.AA.TABLE <- vp.data[[6]]
  vp.GLOBAL.SNPS <- vp.data[[7]]
  vp.SEQ.NUM <- vp.data[[8]]
  vp.SEQ.LENGTH <- vp.data[[9]]
  vp.GENE.NAME <- vp.data[[10]]
  vp.SEG.CODONS <- vp.data[[11]]
  vp.POP.NUM <- vp.data[[12]]
  vp.MAX.AA.VARIANTS <- vp.data[[13]]
  vp.A1 <- vp.data[[14]]
  vp.E1 <- vp.data[[15]]
  vp.E2 <- vp.data[[16]]
  vp.TD.SIM <- vp.data[[17]]
  vp.ALL.SEQUENCES.MATRIX <- vp.data[[18]]
  vp.REF.MATRIX <- vp.data[[19]]
  vp.SEG.SITES <- vp.data[[20]]
  vp.total.row.num <- nrow(vp.BASIC.RESULTS)

  cat("\014")
  cat ("Choose your output:\n")
  cat ("TABLES\n")
  cat ("1  - Results Table \n")
  cat ("2  - Haplotype Table\n")
  cat ("3  - Minor Allele Frequency Table (Nucleotides)\n")
  cat ("4  - Minor Allele Frequency Table (Amino Acids)\n")
  cat ("\n")
  cat ("GRAPHS\n")
  cat ("5  - Haplotype Population Pie Chart\n")
  cat ("6  - AA Variant Percentage Column Graph\n")
  cat ("7  - Phylogenetic Tree\n")
  cat ("8  - Haplotype Accumulation Plot\n")
  cat ("\n")
  cat ("SLIDING SCALE\n")
  cat ("9  - Polymorphism\n")
  cat ("10 - Tajima's D\n")
  cat ("11 - Nucleotide Diversity\n")
  cat ("12 - All Sliding Scale Graphs Overlapped \n")
  cat ("\n")
  cat ("13 - Sliding Scale data table \n")
  cat ("\n")

  #xie186#choice <- as.numeric( readline ("Enter a number from the selection above -  ")) #readline does read the user defined number
  #if (choice %in% c(1:13) == FALSE) stop("There is no output option for that number")
  cat ("\n")

  #RESULTS TABLE-------------------------------------------------------------------------
  #xie186#if (choice == 1) {
    cat("Minimal haplotypes will be calculated using only the segregation sites where polymorphism\n")
    cat("is found in at least 'x' percent of the population at that site")
    threshold <- cutoff_maf_pct #xie186#as.numeric(readline ("'x' <-  ")) 

    variant.codons.over.threshold <- vp.AA.VARIANT.TABLE[,!vp.AA.VARIANT.TABLE[1,]<=threshold/100]
    if (threshold == 0) variant.codons.over.threshold <- vp.AA.VARIANT.TABLE
    colnames(vp.AAs.EACH.VARIABLE.CODON) <- vp.SEG.CODONS #segregation sites AA
    AAs.each.variable.codon.over.threshold <- vp.AAs.EACH.VARIABLE.CODON
    AAs.each.variable.codon.over.threshold <- AAs.each.variable.codon.over.threshold[ ,
      colnames(AAs.each.variable.codon.over.threshold) %in% colnames(variant.codons.over.threshold)] #this return amino acid sites that have MAF over user specified values
    min.hap.strings <- matrix(nrow = vp.SEQ.NUM)
    #AA_based Hap calculation
    for (i in 1:vp.SEQ.NUM){
      min.hap.strings[i,1] <- paste(AAs.each.variable.codon.over.threshold[i,], collapse = "")}
    min.hap.count.col <- matrix(nrow = vp.POP.NUM + 1)
    min.hap.div.col <- matrix(nrow = vp.POP.NUM + 1)
    counter <- 0
    for (i in 1:vp.POP.NUM){
      min.hap.count.col[i,1] <- length(
        unique(min.hap.strings[(1+counter):((vp.BASIC.RESULTS[i,1]) + counter), 1]))
      min.hap.AAs.col <- matrix(nrow = min.hap.count.col[i,1])
      #the following code will give you unique hap above specified threshold
      min.hap.AAs.col[ ,1] <- unique(
        min.hap.strings[(1+counter):((vp.BASIC.RESULTS[i,1]) + counter), ])
      min.hap.table <- match(min.hap.strings[(1+counter):((vp.BASIC.RESULTS[i,1]) + counter), ],
                             min.hap.AAs.col)
      min.hap.freq <- matrix(ncol = 2, nrow = min.hap.count.col[i,1]) #Hap based on AA and its frequency
      min.hap.freq[ ,1] <- c(1:min.hap.count.col[i,1])
      for (j in 1:min.hap.count.col[i,1]){
        AA.hap.counter <- 0
        for (k in 1:vp.BASIC.RESULTS[i,1]){
          if (min.hap.table[k] == j) {AA.hap.counter <- AA.hap.counter + 1}}
        min.hap.freq[j,2] <- AA.hap.counter}
      min.hap.div.col[i,1] <-
        (vp.BASIC.RESULTS[i,1]/
           (vp.BASIC.RESULTS[i,1]-1)) * (1-(sum(((min.hap.freq[,2])/vp.BASIC.RESULTS[i,1])^2))) #hap diversity statistic
      counter <- counter + (vp.BASIC.RESULTS[i,1])}
#for total as above for result table option 1
      min.hap.count.col[vp.total.row.num,1] <- length(unique(min.hap.strings[,1]))
      min.hap.AAs.col <- matrix(nrow = min.hap.count.col[vp.total.row.num,1])
      min.hap.AAs.col[ ,1] <- unique(min.hap.strings[ ,1])
      min.hap.table <- match(min.hap.strings[ ,1], min.hap.AAs.col)
      min.hap.freq <- matrix(ncol = 2, nrow = min.hap.count.col[vp.total.row.num,1])
      min.hap.freq[ ,1] <- c(1:min.hap.count.col[vp.total.row.num,1])
      for (j in 1:min.hap.count.col[vp.total.row.num,1]){
        AA.hap.counter <- 0
        for (k in 1:vp.BASIC.RESULTS[vp.total.row.num,1]){
          if (min.hap.table[k] == j) {AA.hap.counter <- AA.hap.counter + 1}}
        min.hap.freq[j,2] <- AA.hap.counter}
      min.hap.div.col[vp.total.row.num,1] <-
        (vp.BASIC.RESULTS[vp.total.row.num,1]/
           (vp.BASIC.RESULTS[vp.total.row.num,1]-1)) * (1-(sum(((min.hap.freq[,2])/
                vp.BASIC.RESULTS[vp.total.row.num,1])^2)))

    vp.RESULTS.TABLE <<- cbind(vp.BASIC.RESULTS,
                               round(min.hap.count.col, 3), round(min.hap.div.col, 3))
    colnames(vp.RESULTS.TABLE)[12:13] <- c("Minimal AA h", "Minimal AA Hd")
    cat("\n")
    cat("Saved as \"vp.RESULTS.TABLE\", use write.csv() to save to excel \n")
    cat("e.g. write.csv(vp.RESULTS.TABLE, file = \"my.results.table.in.excel\")")
    #View(vp.RESULTS.TABLE)
    #xie186#}


  #POLYMORHPISM GRAPH (S)-----------------------------------------------------------------------------
  #xie186#if (choice == 9) {
    window.size <- win_size  #xie186# as.numeric(readline ("Sliding window size <-  "))
    step.size <- step_size  #xie186# as.numeric(readline ("Step size <-  "))
    number.of.windows <- (ceiling((vp.SEQ.LENGTH / step.size) - (window.size / step.size)))
    sliding.window.S <- matrix(ncol = number.of.windows, nrow = 1)
    windowstart <- 1
    windowend <- window.size

    label.jump <- ceiling((vp.SEQ.LENGTH / 100) / 7) * 100
    label.scaler <- vp.SEQ.LENGTH / number.of.windows
    xlabels.values <- matrix(ncol = 6)
    for (i in 1:6){
      xlabels.values[1,i] <- label.jump + (i-1)*label.jump
    }
    xlabels <- paste(xlabels.values)
    xlabels.values <- xlabels.values / label.scaler
     xlabels.values <- seq(from = xlabels.values[1], to = xlabels.values[6], by = xlabels.values[1])

    for (i in 1:number.of.windows){
      snps <- vp.GLOBAL.SNPS[windowstart:windowend]
      seg.sites <- matrix(0)
      for (j in 1:window.size){
        if (snps[j] != 0) {seg.sites <- cbind(seg.sites, as.matrix(paste(j)))}}
      seg.sites.num <- length(seg.sites) - 1
      sliding.window.S[1,i] <- seg.sites.num
      windowstart <- windowstart + step.size
      windowend <- windowend + step.size
      cat("Calculating window number", i, "of", number.of.windows, "\r")}

    sliding.window.S <- as.data.frame(t(sliding.window.S))


    vp.S.Graph <<- ggplot(data = sliding.window.S,
                          aes(y = V1, x = 1:number.of.windows), na.rm=TRUE)+
      geom_line(colour = "Gold")+
      geom_area(fill = "Yellow")+
      theme_classic()+
      scale_x_continuous(breaks = xlabels.values,
                         limits = c(0, (label.jump+50)*6/label.scaler),
                         labels = xlabels)+
      labs(subtitle = vp.GENE.NAME, y = "S", x = "Nucleotide", title = "Polymorphic Sites")+
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    cat("\n")
    cat("\n")
    cat("Saved as \"vp.S.Graph\"\n")
    #return(suppressMessages(vp.S.Graph))
  #xie186#}

  #TAJIMA'S D GRAPH---------------------------------------------------------------------------
  #if (choice == 10) {
    window.size <- win_size  #xie186# as.numeric(readline ("Sliding window size <-  "))
    step.size <- step_size  #xie186# as.numeric(readline ("Step size <-  "))
    number.of.windows <- (ceiling((vp.SEQ.LENGTH / step.size) - (window.size / step.size)))
    sliding.window.TD <- matrix(ncol = number.of.windows, nrow = 1)
    windowstart <- 1
    windowend <- window.size

    label.jump <- ceiling((vp.SEQ.LENGTH / 100) / 7) * 100
    label.scaler <- vp.SEQ.LENGTH / number.of.windows
    xlabels.values <- matrix(ncol = 6)
    for (i in 1:6){
      xlabels.values[1,i] <- label.jump + (i-1)*label.jump
    }
    xlabels <- paste(xlabels.values)
    xlabels.values <- xlabels.values / label.scaler
    xlabels.values <- seq(from = xlabels.values[1], to = xlabels.values[6], by = xlabels.values[1])

    for (i in 1:number.of.windows){
      snps <- vp.GLOBAL.SNPS[windowstart:windowend]
      seg.sites <- matrix(0)
      for (j in 1:window.size){
        if (snps[j] != 0) {seg.sites <- cbind(seg.sites, as.matrix(paste(j)))}}
      seg.sites.num <- length(seg.sites) - 1
      pi.per.nuc <- ((snps /  (vp.SEQ.NUM)^2) / window.size) * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      tpi <- sum(snps) /  (vp.SEQ.NUM) / vp.SEQ.NUM  * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      pi <- sum(pi.per.nuc)
      td <- ((tpi) -  (seg.sites.num / vp.A1)) /
        (sqrt(  (vp.E1*seg.sites.num)  + ( (vp.E2*seg.sites.num) * (seg.sites.num-1) ) ))
      sliding.window.TD[1,i] <- td + 0
      windowstart <- windowstart + step.size
      windowend <- windowend + step.size
      cat("Calculating window number", i, "of", number.of.windows, "\r")}

    sliding.window.TD <- as.data.frame(t(sliding.window.TD))
    sliding.window.TD[is.na(sliding.window.TD)] <- 0

    td.sim.for.count <- matrix(ncol = 9)
    for (i in 1:72){
      if (vp.TD.SIM[i,1] <= vp.SEQ.NUM)
        if (vp.TD.SIM[(i+1),1] > vp.SEQ.NUM)
          td.sim.for.count <- vp.TD.SIM[i, ]
    }
    if (vp.SEQ.NUM >= 1000) td.sim.for.count <- vp.TD.SIM[73, ]

    vp.TD.Graph <<- suppressMessages(ggplot(data = sliding.window.TD,
                           aes(y = V1, x = 1:number.of.windows), na.rm=TRUE)+
      geom_line(colour = "DarkRed")+
      geom_area(fill = "Red")+
      theme_classic()+
      scale_x_continuous(breaks = xlabels.values,
                         limits = c(0, (label.jump+50)*6/label.scaler),
                         labels = xlabels)+
      scale_y_continuous(sec.axis = sec_axis(~.+0, breaks = (td.sim.for.count[2:9]),
                                             labels = c("0.1", "0.1",
                                                        "0.05", "0.05",
                                                        "0.01", "0.01",
                                                        "0.001", "0.001"),
                                             name = "p value for Tajima's D"))+
      geom_hline(yintercept = (td.sim.for.count[2:9]),
                 colour = "DarkRed", linetype = "dashed", size = 0.4)+
      labs(subtitle = vp.GENE.NAME, y = "Tajima's D", x = "Nucleotide", title = "Tajima's D")+
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)))
    cat("\n")
    cat("Saved as \"vp.TD.Graph\"\n")
    #return(suppressMessages(vp.TD.Graph))
  #xie186#}

  #NUCLEOTIDE DIVERSITY GRAPH-----------------------------------------------------------------------
  #xie186#if (choice == 11) {
    window.size <- win_size  #xie186# as.numeric(readline ("Sliding window size <-  "))
    step.size <- step_size  #xie186# as.numeric(readline ("Step size <-  "))
    number.of.windows <- (ceiling((vp.SEQ.LENGTH / step.size) - (window.size / step.size)))
    sliding.window.pi <- matrix(ncol = number.of.windows, nrow = 1)
    windowstart <- 1
    windowend <- window.size

    label.jump <- ceiling((vp.SEQ.LENGTH / 100) / 7) * 100
    label.scaler <- vp.SEQ.LENGTH / number.of.windows
    xlabels.values <- matrix(ncol = 6)
    for (i in 1:6){
      xlabels.values[1,i] <- label.jump + (i-1)*label.jump
    }
    xlabels <- paste(xlabels.values)
    xlabels.values <- xlabels.values / label.scaler
    xlabels.values <- seq(from = xlabels.values[1], to = xlabels.values[6], by = xlabels.values[1])

    for (i in 1:number.of.windows){
      snps <- vp.GLOBAL.SNPS[windowstart:windowend]
      seg.sites <- matrix(0)
      for (j in 1:window.size){
        if (snps[j] != 0) {seg.sites <- cbind(seg.sites, as.matrix(paste(j)))}}
      seg.sites.num <- length(seg.sites) - 1
      pi.per.nuc <- ((snps /  (vp.SEQ.NUM)^2) / window.size) * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      pi <- sum(pi.per.nuc)
      sliding.window.pi[1,i] <- pi + 0
      windowstart <- windowstart + step.size
      windowend <- windowend + step.size
      cat("Calculating window number", i, "of", number.of.windows, "\r")}

    sliding.window.pi <- as.data.frame(t(sliding.window.pi))
    sliding.window.pi[is.na(sliding.window.pi)] <- 0

    vp.pi.Graph <<- ggplot(data = sliding.window.pi,
                           aes(y = V1, x = 1:number.of.windows), na.rm=TRUE)+
      geom_line(colour = "DarkBlue")+
      geom_area(fill = "Blue")+
      theme_classic()+
      scale_x_continuous(breaks = xlabels.values,
                         limits = c(0, (label.jump +50) *6/label.scaler),
                         labels = xlabels)+
      labs(subtitle = vp.GENE.NAME, y = "\u03a0", x = "Nucleotide", title = "Nucleotide Diversity")+
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    cat("\n")
    cat("Saved as \"vp.pi.Graph\"\n")
    #return(suppressMessages(vp.pi.Graph))
  #xie186#}


  #SLIDING SCALE DATA TABLE---------------------------------------------------------------------
  #xie186#if(choice == 13) {
    window.size <- win_size  #xie186# as.numeric(readline ("Sliding window size <-  "))
    step.size <- step_size  #xie186# as.numeric(readline ("Step size <-  "))
    number.of.windows <- (ceiling((vp.SEQ.LENGTH / step.size) - (window.size / step.size)))
    sliding.window.stats <- matrix(ncol = number.of.windows, nrow = 3)
    windowstart <- 1
    windowend <- window.size

    label.jump <- ceiling((vp.SEQ.LENGTH / 100) / 7) * 100
    label.scaler <- vp.SEQ.LENGTH / number.of.windows
    xlabels.values <- matrix(ncol = 6)
    for (i in 1:6){
      xlabels.values[1,i] <- label.jump + (i-1)*label.jump
    }
    xlabels <- paste(xlabels.values)
    xlabels.values <- xlabels.values / label.scaler

    for (i in 1:number.of.windows){
      snps <- vp.GLOBAL.SNPS[windowstart:windowend]
      seg.sites <- matrix(0)
      for (j in 1:window.size){
        if (snps[j] != 0) {seg.sites <- cbind(seg.sites, as.matrix(paste(j)))}}
      seg.sites.num <- length(seg.sites) - 1
      pi.per.nuc <- ((snps /  (vp.SEQ.NUM)^2) / window.size) * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      tpi <- sum(snps) /  (vp.SEQ.NUM) / vp.SEQ.NUM  * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      pi <- sum(pi.per.nuc)
      td <- ((tpi) -  (seg.sites.num / vp.A1)) /
        (sqrt(  (vp.E1*seg.sites.num)  + ( (vp.E2*seg.sites.num) * (seg.sites.num-1) ) ))
      sliding.window.stats[1,i] <- seg.sites.num
      sliding.window.stats[2,i] <- pi + 0
      sliding.window.stats[3,i] <- td + 0
      windowstart <- windowstart + step.size
      windowend <- windowend + step.size
      cat("Calculating window number", i, "of", number.of.windows, "\r")}
    windowstart <- 1
    windowsites <- paste(windowstart, "-", (windowstart + window.size - 1), sep = "")
    for (i in 1:number.of.windows){
      windowsites <- c(windowsites, paste(windowstart, "-", (windowstart + window.size - 1), sep = ""))
      windowstart <- windowstart + step.size}
    windowsites <- windowsites[-1]
    colnames(sliding.window.stats) <- windowsites
    rownames(sliding.window.stats) <- c("S", "Pi", "TD")

    sliding.window.stats <- as.data.frame(t(sliding.window.stats))
    sliding.window.stats[is.na(sliding.window.stats)] <- 0

    td.sim.for.count <- matrix(ncol = 9)
    for (i in 1:72){
      if (vp.TD.SIM[i,1] <= vp.SEQ.NUM)
        if (vp.TD.SIM[(i+1),1] > vp.SEQ.NUM)
          td.sim.for.count <- vp.TD.SIM[i, ]
    }
    if (vp.SEQ.NUM >= 1000) td.sim.for.count <- vp.TD.SIM[73, ]

    TD_p_values <- matrix(nrow = number.of.windows)
    for (i in 1:number.of.windows){
      if (sliding.window.stats[i,3] < td.sim.for.count[2])
        TD_p_values[i,1] <- 0.1
      if (sliding.window.stats[i,3] > td.sim.for.count[3])
        TD_p_values[i,1] <- 0.1
      if (sliding.window.stats[i,3] < td.sim.for.count[4])
        TD_p_values[i,1] <- 0.05
      if (sliding.window.stats[i,3] > td.sim.for.count[5])
        TD_p_values[i,1] <- 0.05
      if (sliding.window.stats[i,3] < td.sim.for.count[6])
        TD_p_values[i,1] <- 0.01
      if (sliding.window.stats[i,3] > td.sim.for.count[7])
        TD_p_values[i,1] <- 0.01
      if (sliding.window.stats[i,3] < td.sim.for.count[8])
        TD_p_values[i,1] <- 0.001
      if (sliding.window.stats[i,3] > td.sim.for.count[9])
        TD_p_values[i,1] <- 0.001
    }
    
    vp.Sliding.Window.Stats <<- cbind(sliding.window.stats, TD_p_values)
    colnames(vp.Sliding.Window.Stats)[4] <- ("TD p value")
    vp.Sliding.Window.Stats[is.na(vp.Sliding.Window.Stats)] <- ">0.1"
    vp.Sliding.Window.Stats <<- vp.Sliding.Window.Stats

  #xie186#View(vp.Sliding.Window.Stats)
  cat("\n")
  cat("\n")
  cat("Confidence limits of Tajima's D determined by comparison to original \n")
  cat("simulation of beta distribution by F. Tajima (1989)\n")
  cat("\n")
  cat("Saved as \"vp.Sliding.Window.Stats\", use write.csv() to save to excel\n")
  cat("e.g. write.csv(vp.Sliding.Window.Stats, file = \"my.sliding.window.stats.in.excel\")")
  #xie186#}
}
