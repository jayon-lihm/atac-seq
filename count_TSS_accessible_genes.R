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
  ########## PEAK FILE NAME needs to be set ######################################                                                                                                
  peakFile1 <- paste("./", sample_alias, ".original_default.broadPeak", sep="")
  peakFile2 <- paste("./", sample_alias, ".original_default.narrowPeak", sep="")
  ################################################################################                                                                                                

  cmd <- paste("cat ", geneFile, " | sed '1d' | bedtools intersect -a stdin -b ", peakFile1, " -wa -u", sep="")
  a <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
  colnames(a) <- colnames(geneTab)
  o_res$Broad_GeneTSS_covered[i] <- length(unique(a$name2)) ## number of "genes"                                                                                                  
  o_res$Broad_TSS_covered[i] <- dim(a)[1] ## number of "TSS regions"                                                                                                              

  cmd <- paste("cat ", geneFile, " | sed '1d' | bedtools intersect -a stdin -b ", peakFile2, " -wa -u", sep="")
  a <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
  colnames(a) <- colnames(geneTab)
  o_res$Narrow_GeneTSS_covered[i] <- length(unique(a$name2))
  o_res$Narrow_TSS_covered[i] <- dim(a)[1]
}

write.table(o_res, "./Num_TSS_accessible_genes_BrNa.txt", col.names=T, row.names=F, quote=F, sep="\t")

## Output: gene table + sample1 + ... + sample193 + Total number of samples with peaks
## 1 if that TSS is hit
geneRes_br <- geneTab[,1:4]
geneRes_na <- geneTab[,1:4]

namevector <- sampleTab$sample_alias
geneRes_br[,namevector] <- 0
geneRes_na[,namevector] <- 0
gene_ID <- paste(geneRes_br$chrom, geneRes_br$merged_promoter.start, geneRes_br$merged_promoter.end, geneRes_br$name2, sep="_")

for (i in 1:dim(sampleTab)[1]){
  cat(i, ",")
  paper <- sampleTab$paper[i]
  sample_alias <- sampleTab$sample_alias[i]

  ########## PEAK FILE NAME needs to be set ######################################                                                                                                
  peakFile1 <- paste("./", sample_alias, ".original_default.broadPeak", sep="")
  peakFile2 <- paste("./", sample_alias, ".original_default.narrowPeak", sep="")
  ################################################################################                                                                                                

  peakFile <- peakFile1
  cmd <- paste("cat ", geneFile, " | sed '1d' | bedtools intersect -a stdin -b ", peakFile, " -wa -u", sep="")
  a <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
  colnames(a) <- colnames(geneTab)
  a_ID <- paste(a$chrom, a$merged_promoter.start, a$merged_promoter.end, a$name2, sep="_")
  geneRes_br[which(gene_ID %in% a_ID), sample_alias] <- 1

  peakFile <- peakFile2
  cmd <- paste("cat ", geneFile, " | sed '1d' | bedtools intersect -a stdin -b ", peakFile, " -wa -u", sep="")
  a <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
  colnames(a) <- colnames(geneTab)
  a_ID <- paste(a$chrom, a$merged_promoter.start, a$merged_promoter.end, a$name2, sep="_")
  geneRes_na[which(gene_ID %in% a_ID), sample_alias] <- 1

}
cat("\n")

geneRes_br$num_samples <- apply(as.matrix(geneRes_br[, sampleTab$sample_alias]), 1, "sum")
geneRes_na$num_samples <- apply(as.matrix(geneRes_na[, sampleTab$sample_alias]), 1, "sum")

write.table(geneRes_br, "./TSS_AllSamples_Br.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(geneRes_na, "./TSS_AllSAmples_Na.txt", col.names=T, row.names=F, quote=F, sep="\t")
### uploaded files on GITHUB are gzipped                                                                                                                                          

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

