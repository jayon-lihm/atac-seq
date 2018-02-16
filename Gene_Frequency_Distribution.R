## 2018.02.15
## Probability of a gene exhibiting a TSS-accessible peak across 193 samples.
## Black line is for peaks overlapping with TSS, green line is for peaks overlapping in TSS regions, and blue line is for peaks overlapping with genic regions (+/- 1kb of gene body).
## 50% or more samples have peaks in roughly 55% of genic regions in RefSeq.

## Overlapping TSS
## Overlapping TSS +/- 1kb
## Overlapping Genebody +/- 1kb

sampleFile <- "./SRA_sample_info.color_20170913.txt"
tss_file <- "./Merged_mm9_refSeq_TSS.txt"
tss1k_file <- "./Merged_mm9_refSeq_TSS_1kb.txt"
gene1k_file <- "./Merged_mm9_refSeq_Gene_1kb.txt"
geneNameFile <- "./mm9_genename.txt" ## 24,361 gene names

#### output: chrom, start, end, name2, sample1, ..., sample193
gene_peak_overlaps <- function(sample.tab, gene.file){
  ## sample.tab=sampleTab
  ## gene.file=geneFile
  ## gene.res=tss.tab[,1:4] or tss1k.tab[,1:4] or gene1k.tab[,1:4]
  namevector <- sample.tab$sample_alias
  header <- unlist(strsplit(system(paste("cat ", gene.file, " | head -1", sep=""), intern=TRUE), "\t"))
  if (header[1]=="chrom"){
    gene.res <- read.table(gene.file, header=T, as.is=T, sep="\t")
  } else {
    gene.res <- read.table(gene.file, header=F, as.is=T, sep="\t")
  }
  gene.res[,namevector] <- 0
  gene_ID <- paste(gene.res[,1], gene.res[,2], gene.res[,3], gene.res[,4], sep="_") ## chrom, start, end, name2
  
  for (i in 1:dim(sample.tab)[1]){
    cat(i, ":")
    a <- NULL
    a_ID <- NULL
    paper <- sample.tab$paper[i]
    sample_alias <- sample.tab$sample_alias[i]
    ##### PEAK FILE needs to be specified #######
    peakFile <- paste("./", sample_alias, ".merged_genome.broad.macs2_peaks.broadPeak", sep="")
    #############################################
    
    if (header[1]=="chrom"){
      cmd <- paste("cat ", gene.file, " | sed '1d' | bedtools intersect -a stdin -b ", peakFile, " -wa -u", sep="")
    } else {
      cmd <- paste("cat ", gene.file, " | bedtools intersect -a stdin -b ", peakFile, " -wa -u", sep="")
    }
    lines <- as.numeric(system(paste(cmd, " | wc -l", sep=""), intern=TRUE))
    cat(lines, ",")
    if (lines>0){
      a <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
      a_ID <- paste(a[,1], a[,2], a[,3], a[,4], sep="_")
      gene.res[which(gene_ID %in% a_ID), sample_alias] <- 1
    }
  }
  return(gene.res)
}

gene_level_count <- function(genename, res){
  tmp <- genename
  for (i in namevector){
    tmp[,i] <- 0
    i_genename <- unique(res$name2[which(res[,i]==1)])
    a <- which(tmp$name2 %in% i_genename)
    if (length(a)>0){
      tmp[which(tmp$name2 %in% i_genename), i] <- 1
    }
  }
  return(tmp)
}


sampleTab <- read.table(sampleFile, header=T, as.is=T, sep="\t")
namevector <- sampleTab$sample_alias
genename <- read.table(geneNameFile, header=T, as.is=T, sep="\t")
## 1. TSS
tss_res <- gene_peak_overlaps(sampleTab, tss_file2)
tss_res.gene <- gene_level_count(genename, res=tss_res)

## 2. TSS 1kb
tss1k_res <- gene_peak_overlaps(sampleTab, tss1k_file2)
tss1k_res.gene <- gene_level_count(genename, res=tss1k_res)

## 3. Gene 1kb
gene1k_res <- gene_peak_overlaps(sampleTab, gene1k_file2)
gene1k_res.gene <- gene_level_count(genename, res=gene1k_res)

### cumulative distribution
## tmp ==> row: gene, col: tss, tss1k, gene1k
tmp <- genename
tmp$tss_samples <- 0
tmp$tss1k_samples <- 0
tmp$gene1k_samples <- 0

tmp$tss_samples <- rowSums(tss_res.gene[,namevector])
tmp$tss1k_samples <- rowSums(tss1k_res.gene[,namevector])
tmp$gene1k_samples <- rowSums(gene1k_res.gene[,namevector])

write.table(tmp, "./peaks_overlap_with_Genes_perGeneSummary.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(tmp[which(tmp$tss1k_samples==dim(sampleTab)[1]),], "./CommonGenes_451.txt", col.names=T, row.names=F, quote=F, sep="\t")

## Reverse cumulative ( P(X>=x) )
pdf("./RevCumulative_Histogram_gene_TSS_peaks.pdf")
x <- tmp$tss_samples/dim(sampleTab)[1]
cum_y1 <- ecdf(x)
plot(sort(x), 1-cum_y1(sort(x)), type="l", main="% RefSeq Genes with Peaks ", xlab=">= x% Samples with peaks", ylab="% Genes", xlim=c(0,1), ylim=c(0,1), lwd=1.5, xaxt="n", yaxt="n")

x <- tmp$tss1k_samples/dim(sampleTab)[1]
cum_y2 <- ecdf(x)
lines(sort(x), 1-cum_y2(sort(x)), col="green", lwd=1.5)

x <- tmp$gene1k_samples/dim(sampleTab)[1]
cum_y3 <- ecdf(x)
lines(sort(x), 1-cum_y3(sort(x)), col="blue", lwd=1.5)

legend("topright", legend=c("TSS", "TSS +/- 1Kbp", "Gene +/- 1Kbp"), col=c("black", "green", "blue"), lty=c(1,1,1), lwd=c(1.5, 1.5, 1.5))

lines(c(-0.1, 0.5), c(1-cum_y1(0.5), 1-cum_y1(0.5)), col="red", lty="dashed")
lines(c(0.5, 0.5), c(1-cum_y1(0.5), -0.03), col="red", lty="dashed", xpd=TRUE)

lines(c(-0.1, 0.5), c(1-cum_y2(0.5), 1-cum_y2(0.5)), col="red", lty="dashed")
lines(c(0.5, 0.5), c(1-cum_y2(0.5), -0.03), col="red", lty="dashed", xpd=TRUE)

lines(c(-0.1, 0.5), c(1-cum_y3(0.5), 1-cum_y3(0.5)), col="red", lty="dashed")
lines(c(0.5, 0.5), c(1-cum_y3(0.5), -0.03), col="red", lty="dashed", xpd=TRUE)

axis(1, at=seq(0, 1, by=0.2), labels=seq(0, 1, by=0.2)*100)
axis(2, at=seq(0, 1, by=0.2), labels=seq(0, 1, by=0.2)*100)
axis(2, at=c(1-cum_y1(0.5), 1-cum_y2(0.5), 1-cum_y3(0.5)), labels=round(c(1-cum_y1(0.5), 1-cum_y2(0.5), 1-cum_y3(0.5)),3)*100, col.axis="red", cex.axis=0.5, las=2, mgp=c(3, 0.1, 0), tcl=0)

dev.off()

####


