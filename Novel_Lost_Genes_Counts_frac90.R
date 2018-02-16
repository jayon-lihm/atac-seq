## 2018.02.15
## Count lost genes and novel genes at 90% sampling
## For Figure 3(C)

s_info_file <- "./SRA_sample_info.color.txt"
sampleTab <- read.table(s_info_file, header=T, as.is=T, sep="\t")

geneFile <- "./Merged_mm9_refSeq_TSS_1kb.txt"
geneTab <- read.table(geneFile, header=T, as.is=T, sep="\t") ## Selected TSS based on strand, merged set (within same name2 field)

#############
res <- data.frame(paper=rep(sampleTab$paper, each=10), sample_alias=rep(sampleTab$sample_alias, each=10), ninetyPct_rep=rep(1:10, dim(sampleTab)[1]), orig_numGenes=-9, ninetyPct_numGenes=-9, novelGenes=-9, lostGenes=-9)
frac <- 90
for (i in 1:dim(sampleTab)[1]){
  cat(i, ",")
  paper <- sampleTab$paper[i]
  sample_alias <- sampleTab$sample_alias[i]

  file <- paste("./", sample_alias, ".TSS_calls_frac1090.byGene.txt", sep="")
  sample.df <- read.table(file, header=T, as.is=T, sep="\t")

  for (r in 1:10){
    rep.name <- paste("f", frac, "_r", r, sep="")
    a <- sample.df[,c(sample_alias, rep.name)]
    ## a. number of original genes
    res$orig_numGenes[which(res$sample_alias==sample_alias & res$ninetyPct_rep==r)] <- length(which(a[,sample_alias]==1))
    ##a. number of genes
    res$ninetyPct_numGenes[which(res$sample_alias==sample_alias & res$ninetyPct_rep==r)] <- length(which(a[,rep.name]==1))
    ##b. number of novel genes
    res$novelGenes[which(res$sample_alias==sample_alias & res$ninetyPct_rep==r)] <- length(which(a[,sample_alias]==0 & a[,rep.name]==1))
    ##c. number of lost genes
    res$lostGenes[which(res$sample_alias==sample_alias & res$ninetyPct_rep==r)] <- length(which(a[,sample_alias]==1 & a[,rep.name]==0))
  }
}

res$pctNovel <- res$novelGenes / res$ninetyPct_numGenes
res$pctLost <- res$lostGenes / res$orig_numGenes

write.table(res, "./Number_novel_lost_genes_90pct.txt", col.names=T, row.names=F, quote=F, sep="\t")
