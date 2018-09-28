## Figure 4(A)

### INPUT FILES
## files for coordinates of TSS locations, TSS +/- 1kb regions, Gene body +/- 1kb regions
tss_file <- "./Merged_mm9_refSeq_TSS.txt"
tss1k_file <- "./Merged_mm9_refSeq_TSS_1kb.txt"
gene1k_file <- "./Merged_mm9_refSeq_Gene_1kb.txt"

sample_file <- "./SRA_sample_info.color_20170913.txt"
sampleTab <- read.table(sample_file, header=T, as.is=T, sep="\t")

geneNameFile <- "./mm9_genename.txt" ## 24,361 gene names
genename <- read.table(geneNameFile, header=T, as.is=T, sep="\t")

####
## FUNCTION
gene_peak_overlaps <- function(sample.tab, gene.file, peak_surfix){
  ## sample.tab=sampleTab
  ## gene.file=geneFile
  ## gene.res=tss.tab[,1:4] or tss1k.tab[,1:4] or gene1k.tab[,1:4]
  ## peak_surfix= if common surfix is used for peak file name across your samples (OR specify below)
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

    ######## Orignal peak file name ########
    peakFile <- paste("./", sample_alias, ".", peak_surfix, sep="")
    ########################################
    
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

#####
## CHANGE "peak_surfix" depending on your peak file names
## 1. TSS
tss_res <- gene_peak_overlaps(sampleTab, tss_file, peak_surfix="original_default.narrow_peaks.narrowPeak")
tss_res.gene <- gene_level_count(genename, res=tss_res)

## 2. TSS 1kb
tss1k_res <- gene_peak_overlaps(sampleTab, tss1k_file, peak_surfix="original_default.narrow_peaks.narrowPeak")
tss1k_res.gene <- gene_level_count(genename, res=tss1k_res)

## 3. Gene 1kb
gene1k_res <- gene_peak_overlaps(sampleTab, gene1k_file, "original_default.narrow_peaks.narrowPeak")
gene1k_res.gene <- gene_level_count(genename, res=gene1k_res)

#######################
### cumulative distribution
## tmp ==> row: gene, col: tss, tss1k, gene1k
tmp <- genename
tmp$tss_samples <- 0
tmp$tss1k_samples <- 0
tmp$gene1k_samples <- 0

tmp$tss_samples <- rowSums(tss_res.gene[,namevector])
tmp$tss1k_samples <- rowSums(tss1k_res.gene[,namevector])
tmp$gene1k_samples <- rowSums(gene1k_res.gene[,namevector])

write.table(tmp, "./Default_Narrow_peaks_overlap_with_Genes_perGeneSummary.txt", col.names=T, row.names=F, quote=F, sep="\t")

############################
### Calculate CDF: Reverse cumulative ( P(X>=x) 
### With narrow peaks
tmp_na <- read.table("./Default_Narrow_peaks_overlap_with_Genes_perGeneSummary.txt", header=T, as.is=T, sep="\t")
x_tss <- tmp_na$tss_samples/dim(sampleTab)[1]
y_tss <- sapply(1:193, function(x, vec_in){ frac <- x/193; return(length(which(vec_in>=frac))) }, vec_in=x_tss) / dim(tmp_na)[1]

x_tss1k <- tmp_na$tss1k_samples/dim(sampleTab)[1]
y_tss1k <- sapply(1:193, function(x, vec_in){ frac <- x/193; return(length(which(vec_in>=frac))) }, vec_in=x_tss1k) / dim(tmp_na)[1]

x_gene1k <- tmp_na$gene1k_samples/dim(sampleTab)[1]
y_gene1k <- sapply(1:193, function(x, vec_in){ frac <- x/193; return(length(which(vec_in>=frac))) }, vec_in=x_gene1k) / dim(tmp_na)[1]


midpoint <- 193/2
y_tss_mid <- y_tss[97]
y_tss1k_mid <- y_tss1k[97]
y_gene1k_mid <- y_gene1k[97]


pdf("./RevCumulative_Histogram_gene_TSS_narrow_peaks.pdf")
plot( (1:193)/193, y_tss, type="l", main="% RefSeq TSS-accessible Genes (narrow peaks)", xlab=">= x% Samples with peaks", ylab="% Genes", xlim=c(0,1), ylim=c(0,1), lwd=1.5, xaxt="n", yaxt="n")
lines( (1:193)/193, y_tss1k, col="green", lwd=1.5)
lines( (1:193)/193, y_gene1k, col="blue", lwd=1.5)

legend("topright", legend=c("TSS", "TSS +/- 1Kbp", "Gene +/- 1Kbp"), col=c("black", "green", "blue"), lty=c(1,1,1), lwd=c(1.5, 1.5, 1.5))

lines(c(-0.1, 0.5), c(y_tss_mid, y_tss_mid), col="red", lty="dashed")
lines(c(-0.1, 0.5), c(y_tss1k_mid, y_tss1k_mid), col="red", lty="dashed")
lines(c(-0.1, 0.5), c(y_gene1k_mid, y_gene1k_mid), col="red", lty="dashed")
lines(c(0.5, 0.5), c(y_gene1k_mid, -0.03), col="red", lty="dashed", xpd=TRUE)

axis(1, at=seq(0, 1, by=0.2), labels=seq(0, 1, by=0.2)*100)
axis(2, at=seq(0, 1, by=0.2), labels=seq(0, 1, by=0.2)*100)
axis(2, at=c(y_tss_mid, y_tss1k_mid, y_gene1k_mid), labels=round(c(y_tss_mid, y_tss1k_mid, y_gene1k_mid),3)*100, col.axis="red", cex.axis=0.5, las=2, mgp=c(3, 0.1, 0), tcl=0)

dev.off()

