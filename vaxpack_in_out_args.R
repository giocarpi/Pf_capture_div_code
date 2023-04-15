argv = commandArgs(trailingOnly = T)

path <- argv[1]
ref  <- argv[2]
genename = argv[3]
prefix  <- argv[4]


library(ggplot2)

ploidycount = 1
source("vaxpack_input.R")
source("vaxpack_output.R")

vaxpack_input_params(path = path, ref = ref, genename = genename, ploidycount = ploidycount)
vaxpack_output_params(cutoff_maf_pct = 5, win_size = 50, step_size = 5)

csv_out = paste0(prefix, ".div.csv")
tsv_out = paste0(prefix, ".div.tsv")
pdf_out = paste0(prefix, ".div.pdf")

write.csv(vp.RESULTS.TABLE, csv_out, row.names = T)
write.table(vp.RESULTS.TABLE, tsv_out, row.names = T, sep="\t", quote=F)


pdf(pdf_out, width = 5, height = 5)
print(vp.S.Graph)
print(vp.TD.Graph)
print(vp.pi.Graph)
dev.off()

csv_out_stat = paste0(prefix, ".win50_step5.Sliding.Window.Stats.csv")
write.csv(vp.Sliding.Window.Stats, file = csv_out_stat)

