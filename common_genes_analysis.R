### Figure 4B
### Mean expression of commonly accessible genes
library(ggplot2)

sampleFile <- "./SRA_sample_info.color_20170913.txt"
sampleTab <- read.table(sampleFile, header=T, as.is=T, sep="\t")

geneFile <- "./Merged_mm9_refSeq_TSS_1kb.txt"
geneTab <- read.table(geneFile, header=T, as.is=T, sep="\t") ## Selected TSS based on strand, merged set (within same name2 field)

geneNameFile <- "./mm9_genename.txt" ## 24,361 gene names
genename <- read.table(geneNameFile, header=T, as.is=T, sep="\t")

## 2. Mean Expression data
expTable <- read.table("mouse.ranked_mean_sd_expression.txt.gz", header=T, as.is=T, sep=" ") ##38,223
expTable2 <- expTable[which(expTable$symbol %in% geneTab$name2),] ## 22,402
## 0: low expression, 1: high expression
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

expTable2$r.means_2 <- range01(rank(expTable2$r.means))
expTable2$r.sd_2 <- range01(rank(expTable2$r.sd))

## all genes
allGenes <- genename ## 24,361
tmp <- expTable2[which(expTable2$symbol %in% genename$name2),]
tmp2 <- merge(allGenes, tmp[,c("symbol", "r.means_2", "r.sd_2")], by.x="name2", by.y="symbol", all.x=TRUE)

## some gene symbols have multiple expression values
a <- as.data.frame(table(tmp2$name2))
tmp2[which(tmp2$name2 %in% a$Var1[which(a$Freq>1)]),]
allGenes <- tmp2[-which(tmp2$name2 %in% a$Var1[which(a$Freq>1)]),]

filt.allGenes <- allGenes[which(!is.na(allGenes$r.means_2) & !is.na(allGenes$r.sd_2)),]

#########
na_common_file <- "./CommonAccessibleGenes_164_narrow.txt"
na_commonGenes <- read.table(na_common_file, header=T, as.is=T, sep="\t")
na_filt.common <- na_commonGenes[which( !is.na(na_commonGenes$r.means_2) & !is.na(na_commonGenes$r.sd_2) ), ]

## insert marginal histogram
## PLOT FUNCTION from "http://www.statmethods.net/advgraphs/layout.html"
filt.common <- na_filt.common
png("./All_and_CommonGenes_mean_expression_vs_sd.narrow.png")
par(fig=c(0,0.8,0,0.8), new=TRUE)
plot(filt.allGenes$r.sd_2, filt.allGenes$r.means_2, xlab="ranked SD of ranked expression", ylab="average ranked expression")
points(filt.common$r.sd_2, filt.common$r.means_2, col="cadetblue2", pch=19)

par(fig=c(0,0.8,0.55,1), new=TRUE)
xhist <- hist(filt.common$r.sd_2, plot=FALSE, breaks=20)
top <- max(xhist$counts)
x_plot <- c(xhist$counts, 0, 0)
barplot(x_plot, axes=FALSE, ylim=c(0, top), space=0, col="cadetblue2")

par(fig=c(0.65,1,0,0.8),new=TRUE)
yhist <- hist(filt.common$r.means_2, plot=FALSE, breaks=20)
top <- max(yhist$counts)
y_plot <- c(0, 0, 0, yhist$counts)
barplot(y_plot, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, col="cadetblue2")
mtext("Ranked mean expression vs Ranked sd (Narrow)", side=3, outer=TRUE, line=-3)
dev.off()
