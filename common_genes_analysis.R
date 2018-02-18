## 2018.02.16
## 1. Common genes' expression profiles
## Plot mean expression vs mean sd for all genes (in black dots)
## Plot mean expression vs mean sd for common genes (in red dots)

## 2. Control analysis
## Select matching control of 451 genes and compare the distribution of SD.

## 3. GO enrichment analysis
## GO enrichment: Mann Whitney test between mean expression of those genes with GO overlaps and those genes that don't overlap ==> per GO term

sampleFile <- "./SRA_sample_info.color_20170913.txt"
geneFile <- "./Merged_mm9_refSeq_TSS_1kb.txt"
in_file <-"./Gene_numSamples.txt.gz"

sampleTab <- read.table(sampleFile, header=T, as.is=T, sep="\t")
geneTab <- read.table(geneFile, header=T, as.is=T, sep="\t") ## Selected TSS based on strand, merged set (within same name2 field)

####### 1. Mean Expression data #########
expTable <- read.table("./mouse.ranked_mean_sd_expression.txt.gz", header=T, as.is=T, sep=" ") ##38,223

## Select genes that intersect between expression tables and our tables
expTable2 <- expTable[which(expTable$symbol %in% geneTab$name2),] ## 22,402

## Rescale from 0 to 1 (0: low expression, 1: high expression)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
expTable2$r.means_2 <- range01(rank(expTable2$r.means))
expTable2$r.sd_2 <- range01(rank(expTable2$r.sd))

nsample <- 193
df <- read.table(in_file, header=T, as.is=T, sep="\t")
df$num_samples <- apply(as.matrix(df[,3:(3+nsample-1)]), 1, "sum")

## Summary number of samples
df.sum <- data.frame(num_samples=0:nsample, numGene=rep(0, (nsample+1)), cum.numGene=rep(0, (nsample+1)))
for (i in 0:nsample){
  df.sum$numGene[i+1] <- length(which(df$num_samples==i))
  if (i==0){
    df.sum$cum.numGene[i+1] <- df.sum$numGene[i+1]
  } else {
    df.sum$cum.numGene[i+1] <- df.sum$cum.numGene[i] + df.sum$numGene[i+1]
  }
}

## Merge expression table and our gene frequency table
allGenes <- df[,c("name2", "chrom")] ## 24,361
tmp <- expTable2[which(expTable2$symbol %in% df$name2),]
tmp2 <- merge(allGenes, tmp[,c("symbol", "r.means_2", "r.sd_2")], by.x="name2", by.y="symbol", all.x=TRUE)

## Remove gene symbols with multiple expression values
a <- as.data.frame(table(tmp2$name2))
##tmp2[which(tmp2$name2 %in% a$Var1[which(a$Freq>1)]),]
allGenes <- tmp2[-which(tmp2$name2 %in% a$Var1[which(a$Freq>1)]),]

filt.allGenes <- allGenes[which(!is.na(allGenes$r.means_2) & !is.na(allGenes$r.sd_2)),]

## Extract commonly accessible genes hit by all 193 samples
commonGenes <- df[which(df$num_samples==dim(sampleTab)[1]), c("name2", "chrom") ]
commonGenes$r.means_2 <- -9
commonGenes$r.sd_2 <- -9

for (i in 1:dim(commonGenes)[1]){
  gene <- commonGenes$name2[i]
  if ( gene %in% expTable2$symbol ){
    commonGenes$r.means_2[i] <- expTable2$r.means_2[which(expTable2$symbol==gene)]
    commonGenes$r.sd_2[i] <- expTable2$r.sd_2[which(expTable2$symbol==gene)]
  }
}

filt.common <- commonGenes[which(commonGenes$r.means_2!=-9 & commonGenes$r.sd_2!=-9),]

commonGenes$r.means_2[which(commonGenes$r.means_2==-9)] <- NA
commonGenes$r.sd_2[which(commonGenes$r.sd_2==-9)] <- NA

write.table(commonGenes, "./CommonAccessibleGenes_451_RankedGeneExpression.txt", col.names=T, row.names=F, quote=F, sep="\t")

## Insert marginal histogram
## PLOT FUNCTION from "http://www.statmethods.net/advgraphs/layout.html"
png("./All_and_CommonGenes_mean_expression_vs_sd.png")
par(fig=c(0,0.8,0,0.8), new=TRUE)
plot(filt.allGenes$r.sd_2, filt.allGenes$r.means_2, xlab="ranked SD of ranked expression", ylab="average ranked expression")
points(filt.common$r.sd_2, filt.common$r.means_2, col="red", pch=19)

par(fig=c(0,0.8,0.55,1), new=TRUE)
xhist <- hist(filt.common$r.sd_2, plot=FALSE, breaks=30)
top <- max(xhist$counts)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0, col="darkred")

par(fig=c(0.65,1,0,0.8),new=TRUE)
yhist <- hist(filt.common$r.means_2, plot=FALSE, breaks=30)
top <- max(yhist$counts)
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, col="darkred")
mtext("Ranked mean expression vs Ranked sd", side=3, outer=TRUE, line=-3)
dev.off()


#### 2. Control Test ####
## SELECT matching controls of 451 genes
std_expTab <- filt.allGenes[order(filt.allGenes$r.means_2),]
std_expTab$sd_rank <- rank(std_expTab$r.sd_2)
std_expTab$common <- 0
std_expTab$common[which(std_expTab$name2 %in% filt.common$name2)] <- 1

idx <- which(std_expTab$common==1)
idx_minus <- rep(NA, length(idx))
idx_plus <- rep(NA, length(idx))

for (i in 1:length(idx)){
  flag <- 0
  ctrl_idx <- idx[i] - 1
  while (flag==0){
    if ( ctrl_idx %in% unique(c(idx, idx_minus)) ){
      ctrl_idx <- ctrl_idx - 1
    } else {
      flag <- 1
      cat(i, ",", "index:", idx[i], ",", "control index:", ctrl_idx, "\n")
    }
  }
  idx_minus[i] <- ctrl_idx
}


for (i in 1:length(idx)){
  flag <- 0
  ctrl_idx <- idx[i] + 1
  while (flag==0){
    if ( ctrl_idx %in% unique(c(idx, idx_plus, idx_minus)) ){
      ctrl_idx <- ctrl_idx + 1
    } else {
      flag <- 1
      cat(i, ",", "index:", idx[i], ",", "control index:", ctrl_idx, "\n")
    }
  }
  idx_plus[i] <- ctrl_idx
}

all_ctrl_idx <- c(idx_plus, idx_minus)
std_expTab$control <- 0
std_expTab$control[all_ctrl_idx] <- 1

random_ass <- runif(length(idx), 0, 1)
matching_ctrl <- sort(c(idx_minus[which(random_ass<=0.5)], idx_plus[which(random_ass>0.5)]))

which(matching_ctrl %in% idx)

mean(std_expTab$r.means_2[idx])
mean(std_expTab$r.means_2[matching_ctrl])

mean(std_expTab$r.sd_2[idx])
mean(std_expTab$r.sd_2[matching_ctrl])

wilcox.test(std_expTab$r.sd_2[matching_ctrl], std_expTab$r.sd_2[idx])
wilcox.test(std_expTab$r.sd_2[idx], std_expTab$r.sd_2[matching_ctrl])


###### 3. Enrichment Analysis
### EGAD (https://www.bioconductor.org/packages/release/bioc/html/EGAD.html) data
load("EGADlite.RData")
load("biogrid.RData")
load("GO.mouse.RData")
load("GO.voc.RData")
source("gene_set_enrichment.r")

geneTab <- read.table(geneFile, header=T, as.is=T, sep="\t") ## Selected TSS based on strand, merged set (within same name2 field)
gene_list <- unique(geneTab$name2)

## Remove "IEA"
filt <- GO.mouse[,4] != "IEA"

## Select intersection of gene_list & GO.mouse
temp <- GO.mouse[filt,]
filt2 <- temp$name %in% gene_list
temp2 <- temp[filt2,]

## Make table
go.terms <- unique(temp2[,3])
gene.symb <- unique(temp2[,1])
annotations <- make_annotations(temp2[ ,c(1,3)], gene.symb, go.terms ) ## row: gene, column: GO

## Generate annotations: 17,045 genes & 17,724 GO terms
labels.counts <- colSums(annotations)
filt.GO_names <- colnames(annotations)[which(labels.counts>=10 & labels.counts<=300)]## 5,579 GO terms 
annotations2 <- annotations[,filt.GO_names]

## 357 common genes overlap with "temp2" GO involved genes
common.genes <- commonGenes$name2[which(commonGenes$name2 %in% temp2$name)]

###### Mann Whiteney U Test
df.rank <- df[which(df$name2 %in% rownames(annotations)),c("name2", "num_samples")]
df.rank$rank <- rank(df.rank$num_samples)
go.pval <- data.frame(GO=colnames(annotations2), pval=-9, stringsAsFactors=FALSE)

df.rank.s <- df.rank[order(df.rank$name2),] ## 17,045 genes
annotations2.s <- annotations2[order(rownames(annotations2)),] ## 17,045 genes

Y <- df.rank.s$rank
tempres <- apply(as.matrix(annotations2), 2, function(x, y){
   res <- wilcox.test(y~x)
   return(res$p.value)
}, y=Y)

go.pval2 <- go.pval[order(go.pval$GO),]
tempres2 <- tempres[order(names(tempres))]

go.pval2$pval <- as.vector(tempres2)
go.pval2$pval.adj <- p.adjust(go.pval2$pval, method="BH")
go.pval2$numGenes <- colSums(annotations2)[order(colnames(annotations2))]
go.pval2 <- go.pval2[order(go.pval2$pval.adj),]
sig.go.pval <- go.pval2[which(go.pval2$pval.adj<0.05),]
sig.go.pval2 <- merge(sig.go.pval, GO.voc, by.y="GOID", by.x="GO", all.x=TRUE)
sig.go.pval2 <- sig.go.pval2[order(sig.go.pval2$pval.adj),]

write.table(sig.go.pval2, "./MannWhitney_rankTest_CommonGenes_GO_enrichment.txt", col.names=T, row.names=F, quote=F, sep="\t")

## HYPERGEOMETRIC TEST => No significant GO term
## with your genes of interest, run gene set enrichment. Make sure you are consistent between gene identifiers with your gene list of interest.
enrichment_results <- gene_set_enrichment(common.genes, annotations2, GO.voc)
res <- enrichment_results[order(enrichment_results[,"pvals.adj[v.f]"]),]
sigGO <- res[,"pvals.adj[v.f]"]<0.05
sig.res <- res[sigGO,]

filt.GO_names <- colnames(annotations)[which(labels.counts>=5 & labels.counts<=100)]## 7,388 GO terms
annotations2 <- annotations[,filt.GO_names]
enrichment_results <- gene_set_enrichment(common.genes, annotations2, GO.voc)
res <- enrichment_results[order(enrichment_results[,"pvals.adj[v.f]"]),]
sigGO <- res[,"pvals.adj[v.f]"]<0.05
sig.res <- res[sigGO,]

filt.GO_names <- colnames(annotations)[which(labels.counts<=500)]## 17,291 GO terms
annotations2 <- annotations[,filt.GO_names]
enrichment_results <- gene_set_enrichment(common.genes, annotations2, GO.voc)
res <- enrichment_results[order(enrichment_results[,"pvals.adj[v.f]"]),]
sigGO <- res[,"pvals.adj[v.f]"]<0.05
sig.res <- res[sigGO,]


