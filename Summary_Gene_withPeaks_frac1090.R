## 2018.02.15
## Figure 3(B) & 3(C)
## 10% & 90% sampling
## output: gene table + original + 10pct_r1 + ... + 10pct_r10 + 90pct_r1 + ... + 90pct_r10
## 1 if the TSS or GENE is called

args <- commandArgs(TRUE)
i <- as.numeric(args[1]) ## # of row to be directed in sampleTab
tmpdir <- args[2]

if (tmpdir=="") {tmpdir <- "."}

infofile <- "./SRA_sample_info.color_20170913.txt"
sampleTab <- read.table(infofile, header=T, as.is=T, sep="\t")

geneFile <- paste(tmpdir, "/Merged_mm9_refSeq_TSS_1kb.txt", sep="")
geneTab <- read.table(geneFile, header=T, as.is=T, sep="\t") ## Selected TSS based on strand, merged set (within same name2 field)

geneNameFile <- paste(tmpdir, "/mm9_genename.txt", sep="") ## 24,361 gene names
genename <- read.table(geneNameFile, header=T, as.is=T, sep="\t")

orig_geneHit <- read.table(paste(tmpdir, "/Called_TSS_original_peaks.txt.gz", sep=""), header=T, as.is=T, sep="\t")

paper <- as.character(sampleTab$paper[i])
sample_alias <- as.character(sampleTab$sample_alias[i])

######## OUTPUT FILES###############
outfile_tss <- paste("./", sample_alias, ".TSS_calls_frac1090.byTSS.txt", sep="")
outfile_gene <- paste("./", sample_alias, ".TSS_calls_frac1090.byGene.txt", sep="")
####################################

geneRes <- orig_geneHit[,c("chrom", "merged_promoter.start", "merged_promoter.end", "name2", sample_alias)]
gene_ID <- paste(geneRes$chrom, geneRes$merged_promoter.start, geneRes$merged_promoter.end, geneRes$name2, sep="_")

## 1. 10%
frac <- 10
for (rep in 1:10){
  cat(rep, ",")
  rep.name <- paste("f", frac, "_r", rep, sep="")
  geneRes[,rep.name] <- 0

  #### PEAK FILE for 10% sampling that needs to be specified #####
  file10 <- paste("./10pct_", sample_alias, ".random_", rep, ".broad.macs2_peaks.broadPeak.gz", sep="")
  ################################################################
  cmd <- paste("cat ", geneFile, " | sed '1d' | bedtools intersect -a stdin -b ", file10, " -wa -u | wc -l", sep="")
  numL <- as.numeric(system(cmd, intern=TRUE))
  if (numL>0){
    cmd <- paste("cat ", geneFile, " | sed '1d' | bedtools intersect -a stdin -b ", file10, " -wa -u", sep="")
    a <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
    colnames(a) <- colnames(geneTab)
    a_ID <- paste(a$chrom, a$merged_promoter.start, a$merged_promoter.end, a$name2, sep="_")
    geneRes[which(gene_ID %in% a_ID), rep.name] <- 1
  }
}

## 2. 90%
frac <- 90
for (rep in 1:10){
  cat(rep, ",")
  rep.name <- paste("f", frac, "_r", rep, sep="")
  geneRes[,rep.name] <- 0

  ###### PEAK FILE for 90% sampling that needs to be specified ####
  file90 <- paste("./", sample_alias, ".f", frac, "_r", rep, ".broad.macs2_peaks.broadPeak.gz", sep="")
  ##################################################################
  cmd <- paste("cat ", geneFile, " | sed '1d' | bedtools intersect -a stdin -b ", file90, " -wa -u | wc -l", sep="")
  numL <- as.numeric(system(cmd, intern=TRUE))
  if (numL>0){
    cmd <- paste("cat ", geneFile, " | sed '1d' | bedtools intersect -a stdin -b ", file90, " -wa -u", sep="")
    a <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
    colnames(a) <- colnames(geneTab)
    a_ID <- paste(a$chrom, a$merged_promoter.start, a$merged_promoter.end, a$name2, sep="_")
    geneRes[which(gene_ID %in% a_ID), rep.name] <- 1
  }
}
                           

############# sum up to gene level ###############
## genename:: geneTab's name2 & chrom
namevector <- colnames(geneRes)[5:dim(geneRes)[2]]
tmp <- genename
for (i in namevector){
  tmp[,i] <- 0
  i_genename <- unique(geneRes$name2[which(geneRes[,i]==1)])
  a <- which(tmp$name2 %in% i_genename)
  if (length(a)>0){
    tmp[which(tmp$name2 %in% i_genename), i] <- 1
  }
}
                           
write.table(geneRes, outfile_tss, col.names=T, row.names=F, quote=F, sep="\t")
write.table(tmp, outfile_gene, col.names=T, row.names=F, quote=F, sep="\t")

