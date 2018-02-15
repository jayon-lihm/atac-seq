## 2018.02.15
## Count TSS-accessible genes
## Generate summary tables
sampleFile <- "./SRA_sample_info.color_20170913.txt"
sampleTab <- read.table(sampleFile, header=T, as.is=T, sep="\t")

geneFile <- "./Merged_mm9_refSeq_TSS_1kb.txt"
geneTab <- read.table(geneFile, header=T, as.is=T, sep="\t") ## Merged TSS regions (+/- 1kb) (within the same name2 field)

## Output: paper, sample_alias, number of TSS,  number of genes
o_res <- sampleTab[,1:2]

for (i in 1:dim(sampleTab)[1]){
  paper <- sampleTab$paper[i]
  sample_alias <- sampleTab$sample_alias[i]
  ########## PEAK FILE NAME needs to be set ##############
  peakFile <- paste("./", sample_alias, ".merged_genome.macs2_peaks.broadPeak", sep="")
  ########################################################
  cmd <- paste("cat ", geneFile, " | sed '1d' | bedtools intersect -a stdin -b ", peakFile, " -wa -u", sep="")
  a <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
  colnames(a) <- colnames(geneTab)
  o_res$Gene_TSS_covered[i] <- length(unique(a$name2)) ## Number of Genes
  o_res$TSS_covered[i] <- dim(a)[1] ## Number of TSS
}

write.table(o_res, "./Num_TSS_accessible_genes.txt", col.names=T, row.names=F, quote=F, sep="\t")

## Output: gene table + sample1 + ... + sample193 + Total number of samples with peaks
## 1 if that TSS is hit
geneRes <- geneTab[,1:4]
namevector <- sampleTab$sample_alias
geneRes[,namevector] <- 0
gene_ID <- paste(geneRes$chrom, geneRes$merged_promoter.start, geneRes$merged_promoter.end, geneRes$name2, sep="_")

for (i in 1:dim(sampleTab)[1]){
  paper <- sampleTab$paper[i]
  sample_alias <- sampleTab$sample_alias[i]
  ########## PEAK FILE NAME needs to be set ##############
  peakFile <- paste("./", sample_alias, ".merged_genome.macs2_peaks.broadPeak", sep="")
  ########################################################
  cmd <- paste("cat ", geneFile, " | sed '1d' | bedtools intersect -a stdin -b ", peakFile, " -wa -u", sep="")
  a <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
  colnames(a) <- colnames(geneTab)
  a_ID <- paste(a$chrom, a$merged_promoter.start, a$merged_promoter.end, a$name2, sep="_")
  geneRes[which(gene_ID %in% a_ID), sample_alias] <- 1
}

geneRes$num_samples <- apply(as.matrix(geneRes[,5:dim(geneRes)[2]]), 1, "sum")
write.table(geneRes, "./TSS_numSamples.txt", col.names=T, row.names=F, quote=F, sep="\t")

## Aggreagate at gene level
## Output: gene name + chromosome + sample 1 + ... + sample 193
## 1 if any of TSS in a gene is hit by peaks
tmp <- data.frame(name2=unique(geneTab$name2), stringsAsFactors=FALSE)
tmp$chrom <- sapply(tmp$name2, function(x){ return( paste(unique(geneTab$chrom[which(geneTab$name2==x)]), collapse=",") ) })
tmp[,namevector] <- 0

for (i in namevector){
  sampleName <- i
  cat(i, "\n")
  tmp[,sampleName] <- sapply(tmp$name2, function(x){ a <- sum(geneRes[which(geneRes$name2==x), sampleName]); if (a>0){ return(1) } else if (a==0){ return(0) } })
}

write.table(tmp, "./Gene_numSamples.txt", col.names=T, row.names=F, quote=F, sep="\t")

