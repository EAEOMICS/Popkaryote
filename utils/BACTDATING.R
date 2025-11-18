args <- commandArgs(trailingOnly = TRUE)

library(ape)
library(BactDating)
library(coda)
clonalframe<- file.path(args[2],"cf_out")
t=loadCFML(clonalframe)
DIR=file.path(args[2],"bactdating")
dir.create(DIR, recursive = TRUE)

dic <- read.table(args[1], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(dic) <- c("sample", "date")
dates_dict <- setNames(dic$date, dic$sample)
d <- t$tip.label
d <- dates_dict[d]

if (args[3]== 'True'){ 
	samples <- t$tip.label
	years <- dates_dict[samples]
	d_range <- cbind(years, years + 1)
  d_mean <- rowMeans(d_range, na.rm = TRUE)
  rooted=initRoot(t,d_mean)
  ROOT=file.path(args[2],"bactdating/roottotip_plot.png")
  png(ROOT, width = 800, height = 600)
  r=roottotip(rooted,d_mean)
  dev.off()
  res=bactdate(rooted,d_mean,nbIts=1e6)
  trace=file.path(args[2],"bactdating/trace.png")
  png(trace, width = 800, height = 600)
  plot(res, 'trace')
  dev.off()
  
  mcmc = as.mcmc.resBactDating(res)
  
  TREE=file.path(args[2],"bactdating/tree.png")
  png(TREE,width=800,height=600)
  plot(res,'treeCI')
  dev.off()
  
  ess <- effectiveSize(mcmc)
  ES=file.path(args[2],"bactdating/effective_size.txt")
  write.table(ess, file = ES, quote = FALSE, sep = "\t", col.names = NA)
} else {
  rooted=initRoot(t,d)
  ROOT=file.path(args[2],"bactdating/roottotip_plot.png")
  png(ROOT, width = 800, height = 600)
  r=roottotip(rooted,d)
  dev.off()
  res=bactdate(rooted,d,nbIts=1e6)

  trace=file.path(args[2],"bactdating/trace.png")
  png(trace, width = 800, height = 600)
  plot(res, 'trace')
  dev.off()

  mcmc = as.mcmc.resBactDating(res)

  TREE=file.path(args[2],"bactdating/tree.png")
  png(TREE,width=800,height=600)
  plot(res,'treeCI')
  dev.off()

  ess <- effectiveSize(mcmc)

  ES=file.path(args[2],"bactdating/effective_size.txt")
  write.table(ess, file = ES, quote = FALSE, sep = "\t", col.names = NA)
}
