args <- commandArgs(trailingOnly = TRUE)

setwd<-args[2]
dset<-file.path(args[2],"/LDWeaver_Results")
aln_path<-file.path(args[2],"/core.full.aln")
gbk_path<-args[1]

LDWeaver::LDWeaver(dset = dset, 
aln_path = aln_path, 
gbk_path = gbk_path, 
save_additional_outputs = T)
