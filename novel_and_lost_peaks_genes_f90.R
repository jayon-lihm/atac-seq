### Figure 3 B & C
## 1. Count novel and lost "Genes"
## 2. Generate table with row: genes, column: 1 if overlaps with original peak, 90% 1st replicate, ..., 10th replicate
## 3. Average the novel and lost peaks over 10 replicates
## 4. And plot (Figure 3b&c)

library("openintro")
library("gdata")
sampleFile <- "./SRA_sample_info.color_20170913.txt"
sampleTab <- read.table(sampleFile, header=T, as.is=T, sep="\t")

geneFile <- "./Merged_mm9_refSeq_TSS_1kb.txt"
geneTab <- read.table(geneFile, header=T, as.is=T, sep="\t") ## Selected TSS based on strand, merged set (within same name2 field)
gene_ID <- paste(geneTab$chrom, geneTab$merged_promoter.start, geneTab$merged_promoter.end, geneTab$name2, sep="_")

geneNameFile <- "./mm9_genename.txt" ## 24,361 gene names
genename <- read.table(geneNameFile, header=T, as.is=T, sep="\t")

fd <- toupper(as.hexmode(150)) # 64 in hexa == 100 in decimal

path <- "./"

f <- 90
frac <- 90
sum_file <- "./Num_TSS_accessible_genes_BrNa.txt"
geneSumTab <- read.table(sum_file, header=T, as.is=T, sep="\t")

BrFile <- "./TSS_AllSamples_Br.txt.gz"
NaFile <- "./TSS_AllSamples_Na.txt.gz"

Br_orig_genes <- read.table(BrFile, header=T, as.is=T, sep="\t")
Na_orig_genes <- read.table(NaFile, header=T, as.is=T, sep="\t")

#### FUNCTION
peaks_ovl_TSS <- function(peakFile, out_data, new_col_name, gene_file=geneFile){
  gene_tab <- read.table(gene_file, header=T, as.is=T, sep="\t")
  gene_ID <- paste(geneTab$chrom, geneTab$merged_promoter.start, geneTab$merged_promoter.end, geneTab$name2, sep="_")
  out_data[,new_col_name] <- 0
  cmd <- paste("cat ", gene_file, " | sed '1d' | bedtools intersect -a stdin -b ", peakFile, " -wa -u | wc -l", sep="")
  numL <- as.numeric(system(cmd, intern=TRUE))
  if (numL>0){
    cmd <- paste("cat ", gene_file, " | sed '1d' | bedtools intersect -a stdin -b ", peakFile, " -wa -u", sep="")
    a <- read.table(pipe(cmd), header=F, as.is=T, sep="\t")
    colnames(a) <- colnames(gene_tab)
    a_ID <- paste(a$chrom, a$merged_promoter.start, a$merged_promoter.end, a$name2, sep="_")
    out_data[which(gene_ID %in% a_ID), new_col_name] <- 1
  }
  return(out_data)
}

TSS_to_GENE <- function(data_in, name_vector, gene_tab=genename){
  ## data_in: table of chrom, start, end, gene name (name2), samples
  ## name_vector: colnames to aggregate TSS to gene level
  data_out <- genename
  data_out[,name_vector] <- 0
  for (i in name_vector){
    data_out[,i] <- 0
    i_genename <- unique(data_in$name2[which(data_in[,i]==1)])
    a <- which(data_out$name2 %in% i_genename)
    if (length(a)>0){
      data_out[which(data_out$name2 %in% i_genename), i] <- 1
    }
  }
  return(data_out)
}

### 1. Count novel/lost genes per sample, per replicate
for (j in 1:dim(sampleTab)[1]){
  paper <- sampleTab$paper[j]
  sample_alias <- sampleTab$sample_alias[j]

  cat(paper, sample_alias)
  Br_tssRes <- Br_orig_genes[,c("chrom", "merged_promoter.start", "merged_promoter.end", "name2")]
  Na_tssRes <- Na_orig_genes[,c("chrom", "merged_promoter.start", "merged_promoter.end", "name2")]

  ###### Original peak files need to be specified here ###########################################
  br_originalFile <- paste("./", sample_alias, ".original_default.broadPeak", sep="")
  na_originalFile <- paste("./", sample_alias, ".original_default.narrowPeak", sep="")
  #############################################################################################

  ## Broad & Narrow
  Br_tssRes <- peaks_ovl_TSS(br_originalFile, out_data=Br_tssRes, new_col_name=sample_alias)
  Na_tssRes <- peaks_ovl_TSS(na_originalFile, out_data=Na_tssRes, new_col_name=sample_alias)
  
  cat(":: ")
  for (rep in 1:10){
    cat(rep, ",")
    rep.name <- paste("f", frac, "_r", rep, sep="")

    ###### 90% sampling peak files need to be specified here ###########################################
    Br_f90_file <- paste("./", sample_alias, ".f", frac, "_r", rep, ".DefaultBroad.macs2_broadPeak.gz", sep="")
    Na_f90_file <- paste("./", sample_alias, ".f", frac, "_r", rep, ".DefaultNarrow.macs2_narrowPeak.gz", sep="")
    #############################################################################################
    
    Br_tssRes <- peaks_ovl_TSS(Br_f90_file, out_data=Br_tssRes, new_col_name=rep.name)
    Na_tssRes <- peaks_ovl_TSS(Na_f90_file, out_data=Na_tssRes, new_col_name=rep.name)
  }

  ## Sum up to gene level
  Br_geneRes <- TSS_to_GENE(Br_tssRes, name_vector=colnames(Br_tssRes)[5:dim(Br_tssRes)[2]])
  Na_geneRes <- TSS_to_GENE(Na_tssRes, name_vector=colnames(Na_tssRes)[5:dim(Na_tssRes)[2]])

  Br_outfile_tss <- paste("./", sample_alias, ".TSS_overlaps_Broad.byTSS.txt", sep="")
  Br_outfile_gene <- paste("./", sample_alias, ".TSS_overlaps_Broad.byGene.txt", sep="")
  
  Na_outfile_tss <- paste("./", sample_alias, ".TSS_overlaps_Narrow.byTSS.txt", sep="")
  Na_outfile_gene <- paste("./", sample_alias, ".TSS_overlaps_Narrow.byGene.txt", sep="")

  ############ OUTPUT ####################
  ## TSS-level output: chrom + merged_promoter.start + merged_promoter.end + gene name + 1/0 from original reads + f90_r1 + ... + f90_r10
  write.table(Br_tssRes, Br_outfile_tss, col.names=T, row.names=F, quote=F, sep="\t")
  write.table(Na_tssRes, Na_outfile_tss, col.names=T, row.names=F, quote=F, sep="\t")

  ## Gene-level output: chrom + gene name + 1/0 from original reads + f90_r1 + ... + f90_r10     
  write.table(Br_geneRes, Br_outfile_gene, col.names=T, row.names=F, quote=F, sep="\t")
  write.table(Na_geneRes, Na_outfile_gene, col.names=T, row.names=F, quote=F, sep="\t")
  cat("\n")
  
}

##########################################################
###2. Average number of Novel/Lost genes
lostSumFile <- "./NumPeaks_NovelLost_AllReps.txt.gz"
lostSumTab <- read.table(lostSumFile, header=T, as.is=T, sep="\t")

res <- lostSumTab[ which(lostSumTab$frac==90), ]

res$numGenes_def_br <- -9
res$numGenes_def_na <- -9

res$numNovelGenes_def_br <- -9
res$numNovelGenes_def_na <- -9

res$numLostGenes_def_br <- -9
res$numLostGenes_def_na <- -9

res$orig.numPeaks_def_br <- -9
res$orig.numPeaks_def_na <- -9

res$orig.numGenes_def_br <- -9
res$orig.numGenes_def_na <- -9

for (i in 1:dim(sampleTab)[1]){
  cat(i, "::")
  paper <- sampleTab$paper[i]
  sample_alias <- sampleTab$sample_alias[i]

  ## original numbers
  ## peaks

  ###### Original peak files need to be specified here ###########################################
  br_originalFile <- paste("./", sample_alias, ".original_default.broadPeak", sep="")
  na_originalFile <- paste("./", sample_alias, ".original_default.narrowPeak", sep="")
  #############################################################################################
  
  cmd <- paste("cat ", br_originalFile, " | wc -l", sep="")
  res$orig.numPeaks_def_br[which(res$paper==paper & res$sample_alias==sample_alias)] <-    as.numeric(system(cmd, intern=TRUE))

  cmd <- paste("cat ", na_originalFile, " | wc -l", sep="")
  res$orig.numPeaks_def_na[which(res$paper==paper & res$sample_alias==sample_alias)] <-    as.numeric(system(cmd, intern=TRUE))

  ## genes

  ###### Output files from part #1 needs to be pasted here ##############################
  Br_outfile_gene <- paste("./", sample_alias, ".TSS_overlaps_Broad.byGene.txt", sep="")
  Na_outfile_gene <- paste("./", sample_alias, ".TSS_overlaps_Narrow.byGene.txt", sep="")
  #############################################################################################
  
  br_df <- read.table(Br_outfile_gene, header=T, as.is=T, sep="\t")
  na_df <- read.table(Na_outfile_gene, header=T, as.is=T, sep="\t")
  
  res$orig.numGenes_def_br[which(res$paper==paper & res$sample_alias==sample_alias)]  <- sum(br_df[,sample_alias])
  res$orig.numGenes_def_na[which(res$paper==paper & res$sample_alias==sample_alias)]  <- sum(na_df[,sample_alias])

  for (r in 1:10){
    loc <- which(res$sample_alias==sample_alias & res$repl==r)
    if ( res$orig.numGenes_def_br[loc] != sum(br_df[,sample_alias]) | res$orig.numGenes_def_na[loc] != sum(na_df[,sample_alias]) ){ cat("ERROR in original number of genes. BREAK. \n"); break }
    rep_name <- paste("f90_r", r, sep="")
    res$numGenes_def_br[loc] <- sum(br_df[,rep_name])
    res$numGenes_def_na[loc] <- sum(na_df[,rep_name])

    res$numNovelGenes_def_br[loc] <- length(which(br_df[,rep_name]==1 & br_df[,sample_alias]==0))
    res$numNovelGenes_def_na[loc] <- length(which(na_df[,rep_name]==1 & na_df[,sample_alias]==0))
    
    res$numLostGenes_def_br[loc] <- length(which(br_df[,rep_name]==0 & br_df[,sample_alias]==1))
    res$numLostGenes_def_na[loc] <- length(which(na_df[,rep_name]==0 & na_df[,sample_alias]==1))
  }
}

write.table(res,"./Novel_lost_genes_f90_AllReps.txt", col.names=T, row.names=F, quote=F, sep="\t")


## 3. Summary over 10 replicates
res <- read.table("./Novel_lost_genes_f90_AllReps.txt", header=T, as.is=T, sep="\t") ## output file from part #2 above
res$br.frac_novel_peaks <- -9
res$br.frac_novel_genes <- -9
res$br.frac_lost_peaks <- -9
res$br.frac_lost_genes <- -9

res$na.frac_novel_peaks <- -9
res$na.frac_novel_genes <- -9
res$na.frac_lost_peaks <- -9
res$na.frac_lost_genes <- -9

res$br.frac_novel_peaks <- res$numNovelPeaks_def_br / res$numPeaks_def_br
res$na.frac_novel_peaks <- res$numNovelPeaks_def_na / res$numPeaks_def_na

res$br.frac_novel_genes <- res$numNovelGenes_def_br / res$numGenes_def_br
res$na.frac_novel_genes <- res$numNovelGenes_def_na / res$numGenes_def_na

res$br.frac_lost_peaks <- res$numLostPeaks_def_br / res$orig.numPeaks_def_br
res$na.frac_lost_peaks <- res$numLostPeaks_def_na / res$orig.numPeaks_def_na

res$br.frac_lost_genes <- res$numLostGenes_def_br / res$orig.numGenes_def_br
res$na.frac_lost_genes <- res$numLostGenes_def_na / res$orig.numGenes_def_na

tmp <- aggregate( . ~ paper + sample_alias + frac, data=res, FUN="mean")
tmp <- tmp[,which(colnames(tmp) != "repl")]

write.table(tmp, "./Average_novel_lost_genes_f90.txt", col.names=T, row.names=F, quote=F, sep="\t")

## 4. Plot
avg_res <- read.table("./Average_novel_lost_genes_f90.txt", header=T, as.is=T, sep="\t") ## output file from part #3 above

pdf("./histogram_frac_lost_peaks_genes_f90.pdf")
par(mfrow=c(4,1))
hist(log(avg_res$br.frac_novel_peaks, base=10), main="Fraction of novel peaks at 90% sampling (broad)", xlab="log10(fraction)", ylab="number of samples", breaks=20, xlim=c(-3, 1))
hist(log(avg_res$br.frac_novel_genes, base=10), main="Fraction of novel genes at 90% sampling (broad)", xlab="log10(fraction)", ylab="number of samples", breaks=20, xlim=c(-3, 1))
hist(log(avg_res$br.frac_lost_peaks, base=10), main="Fraction of lost peaks at 90% sampling (broad)", xlab="log10(fraction)", ylab="number of samples", breaks=20, xlim=c(-3, 1))
hist(log(avg_res$br.frac_lost_genes, base=10), main="Fraction of lost genes at 90% sampling (broad)", xlab="log10(fraction)", ylab="number of samples", breaks=20, xlim=c(-3, 1))

par(mfrow=c(4,1))
hist(log(avg_res$na.frac_novel_peaks, base=10), main="Fraction of novel peaks at 90% sampling (narrow)", xlab="log10(fraction)", ylab="number of samples", breaks=20, xlim=c(-3, 1))
hist(log(avg_res$na.frac_novel_genes, base=10), main="Fraction of novel genes at 90% sampling (narrow)", xlab="log10(fraction)", ylab="number of samples", breaks=20, xlim=c(-3, 1))
hist(log(avg_res$na.frac_lost_peaks, base=10), main="Fraction of lost peaks at 90% sampling (narrow)", xlab="log10(fraction)", ylab="number of samples", breaks=20, xlim=c(-3, 1))
hist(log(avg_res$na.frac_lost_genes, base=10), main="Fraction of lost genes at 90% sampling (narrow)", xlab="log10(fraction)", ylab="number of samples", breaks=20, xlim=c(-3, 1))

dev.off()

