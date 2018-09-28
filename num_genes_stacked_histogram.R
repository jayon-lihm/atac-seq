## 2018.02.15
## Figure 2 (C)
## Plot stacked histogram of TSS-accessible genes
## Color assigned by study group

sampleFile <- "./SRA_sample_info.color_20170913.txt"
TSS_File <- "./Num_TSS_accessible_genes_BrNa.txt"

sampleTab <- read.table(sampleFile, header=T, as.is=T, sep="\t")
tssTab <- read.table(TSS_File, header=T, as.is=T, sep="\t")

## Broad
res <- NULL
brk_bins <- seq(2000, 15500, by=500)
sample_counts <- NULL
for (i in 1:length(unique(sampleTab$paper))){
  paper <- unique(sampleTab$paper)[i]
  tmp <- tssTab[which(tssTab$paper==paper),]
  sample_counts <- c(sample_counts, dim(tmp)[1])
  i_count <- as.matrix(table(cut(tmp$Broad_GeneTSS_covered, breaks=brk_bins)))
  res <- cbind(res, i_count)
  cat(paper, mean(tmp$Broad_GeneTSS_covered), "\n")
}

colnames(res) <- c(unique(sampleTab$paper))
br.res <- res

## Narrow
res <- NULL
brk_bins <- seq(2000, 15500, by=500)
sample_counts <- NULL
for (i in 1:length(unique(sampleTab$paper))){
  paper <- unique(sampleTab$paper)[i]
  tmp <- tssTab[which(tssTab$paper==paper),]
  sample_counts <- c(sample_counts, dim(tmp)[1])
  i_count <- as.matrix(table(cut(tmp$Narrow_GeneTSS_covered, breaks=brk_bins)))
  res <- cbind(res, i_count)
  cat(paper, mean(tmp$Narrow_GeneTSS_covered), "\n")
}

colnames(res) <- c(unique(sampleTab$paper))
na.res <- res


## assign color per study
color.vec <- NULL
for (i in colnames(res)){
  color.vec <- c(color.vec, unique(sampleTab$color[which(sampleTab$paper==i)]))
}

## PLOT
## Initiate Bar plot
pdf("Number_of_Genes_Distribution_histogram_BrNa.pdf")
a <- barplot(t(br.res), main="Number of Genes with Broad Peaks at TSS", xlab="# Genes", ylab="# Samples", col=color.vec, legend=colnames(br.res),
             args.legend = list(x="topleft"), xaxt="n")
xlabel <- brk_bins
xlabel2 <- format(xlabel, big.mark=",",scientific=FALSE)
bar_interval <- unique(round(diff(a), 1))
tck_loc <- c(0, a+bar_interval/2)
text(tck_loc, rep(-0.5, length(xlabel2)), srt=45, adj=1, xpd=TRUE, labels=xlabel2, cex=0.7)

a <- barplot(t(na.res), main="Number of Genes with Narrow Peaks at TSS", xlab="# Genes", ylab="# Samples", col=color.vec, legend=colnames(na.res),
             args.legend = list(x="topleft"), xaxt="n")
xlabel <- brk_bins
xlabel2 <- format(xlabel, big.mark=",",scientific=FALSE)
bar_interval <- unique(round(diff(a), 1))
tck_loc <- c(0, a+bar_interval/2)
text(tck_loc, rep(-0.5, length(xlabel2)), srt=45, adj=1, xpd=TRUE, labels=xlabel2, cex=0.7)

dev.off()
