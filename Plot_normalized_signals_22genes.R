## 2018.02.18
## 1. Normalized Read Counts
## 2. Relative position to TSS +/- 1kb, normalized signal plot (with random 22 genes as control)
s_info_file <- "MouseBrain_SampleInfo.txt"
sampleTab <- read.table(s_info_file, header=T, as.is=T, sep="\t")

filter_samples <- c("AMYGDALA_FEAR_2", "AMYGDALA_ERRB4_WT_2")
valid_samples <- setdiff(sampleTab$merged_ID, filter_samples)

amyg <- grep("AMYGDALA", valid_samples)
cortex <- grep("CORTEX", valid_samples)

all_HG_test <- read.table("hyp_geo_test.TISSUE.txt", header=T, as.is=T, sep="\t")
gene22 <- all_HG_test$name2[which(all_HG_test$yx_p.value_adj<0.05)]

###### 1. Normalized read counts of amygdala and cortex
org_reads <- read.table("MouseBrain_Read_Counts_TSS1kb.txt.gz", header=T, as.is=T, sep="\t")

Sample_List <- colnames(org_reads)[7:dim(org_reads)[2]]
tot_reads <- matrix(colSums(org_reads[,Sample_List]), ncol=1)
rownames(tot_reads) <- Sample_List

org_reads.merged <- aggregate( . ~ name2, data=org_reads[,c("name2", Sample_List)], FUN="sum")
tot_reads.mat <- tot_reads %*% matrix(1, nrow=1, ncol=dim(org_reads.merged)[1])
tot_reads.mat <- t(tot_reads.mat)

norm_reads.merged <- data.frame(name2=org_reads.merged[,"name2"], (org_reads.merged[,Sample_List] / tot_reads.mat), stringsAsFactors=FALSE)
colnames(norm_reads.merged) <- c("name2", Sample_List)

df <- norm_reads.merged[which(norm_reads.merged$name2 %in% gene22),]
data_sum <- apply(as.matrix(df[,Sample_List]), 2, "sum")
data_sum2 <- data_sum[valid_samples]

plot_data <- data.frame(normalized_reads=data_sum2)
plot_data$group <- "."
plot_data$group[grep("AMYGDALA", rownames(plot_data))] <- "Amygdala"
plot_data$group[grep("CORTEX", rownames(plot_data))] <- "Cortex"
plot_data$color <- "red" # cortex
plot_data$color[which(plot_data$group=="Amygdala")] <- "blue" # amygdala

plot_data$num.group <- 1.2 # Cortex
plot_data$num.group[which(plot_data$group=="Amygdala")] <- 1.1 #Amygdala

pdf("Normalized_read_counts.22genes.pdf")
plot(jitter(plot_data$num.group), plot_data$normalized_reads, col=plot_data$color, pch=16, main="Normalized number of reads in TSS-regions of 22 genes", xlab="Brain Regions", ylab="Sum of Normalized number of reads" , xlim=c(1.05, 1.25), xaxt="n")
axis(side=1, at=c(1.1,1.2), labels=c("Amygdala", "Cortex"))
dev.off()


t.test(plot_data$normalized_reads[which(plot_data$group=="Amygdala")], plot_data$normalized_reads[which(plot_data$group=="Cortex")], ) ## p-value = 0.0005578
t.test(plot_data$normalized_reads[which(plot_data$group=="Amygdala")], plot_data$normalized_reads[which(plot_data$group=="Cortex")], alternative="less") ##  p-value = 0.0002789



####### 2. Relative position to TSS
## Random 22 genes
for (rep in 1:10){
  random22.idx <- sample(which(all_HG_test$yx_p.value_adj>=0.05), 22) ## Select 22 genes with adjusted p-value greater than 0.05
  random22 <- all_HG_test$name2[random22.idx]
  randRes <- data.frame(index=1:2000)
  for (i in 1:length(valid_samples)){
    cat("rep", rep, ",", "sample", i, "\n")
    sampleID <- valid_samples[i]
    #### DEPTH FILE within 1kb of TSS ######
    dp_file <- paste("./", sampleID, ".merged_genome.sorted.bam.TSS1k_depth.gz", sep="")
    ########################################
    
    cmd <- paste("zcat ", dp_file, " | awk '", paste("$8==\"", random22, "\"", sep="", collapse=" || "), "'", sep="") ## Select random 22 genes from the depth file
    con <- pipe(cmd, open="r")
    tmp22 <- read.csv(con, header=F, as.is=T, sep="\t")
    colnames(tmp22) <- c("chrom.depth", "start.depth", "end.depth", "depth", "chrom.gene", "start.gene", "end.gene", "gene", "X", "Y")

    cmd2 <- paste("zcat ", dp_file, " | awk 'BEGIN{sum=0;count=0}{sum+=$4; count+=1}END{print sum, count}'", sep="") ## Total depth within 1kb of TSS of all genes
    a <- unlist(strsplit(system(cmd2, intern=TRUE), " "))
    tot_depth_TSS <- as.numeric(a[1])
    tot_covered_pos_TSS <- as.numeric(a[2])
    
    tmp22$norm.depth <- tmp22$depth/tot_depth_TSS

    tmp.mat <- data.frame(index=1:2000)
    for (j in 1:length(random22)){
      oneGene <- tmp22[which(tmp22$gene==random22[j]),]
      tmp.mat[,random22[j]] <- 0
      for (k in unique(oneGene$start.gene)){
        oneTSS <- oneGene[which(oneGene$start.gene==k),]
        idx <- oneTSS$end.depth - oneTSS$start.gene
        oneTSS_2 <- oneTSS[which(idx>=1 & idx<=2000),]
        idx <- oneTSS_2$end.depth - oneTSS_2$start.gene
        tmp.mat[idx, random22[j]] <- tmp.mat[idx,random22[j]] + oneTSS_2$norm.depth
      }
    }
    
    tmp.mat$sum <- apply(as.matrix(tmp.mat[,-1]), 1, "sum")
    randRes[,sampleID] <- tmp.mat$sum
  }

  #### Normalized signals per replicate
  write.table(randRes, paste("./Normalized_signals_random22genes.rep", rep, ".txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

  pdf(paste("./Normalized_signals_random22genes.TSScentered.rep", rep, ".pdf", sep=""))
  x <- randRes$index
  ya <- apply(as.matrix(randRes[,valid_samples[amyg]]), 1, "mean")
  yc <- apply(as.matrix(randRes[,valid_samples[cortex]]), 1, "mean")
  
  y_range <- c(min(ya,yc), max(ya, yc))
  
  plot(x,yc, col="red", type="l", main="Normalized signals of random 22 genes centered at TSS", xlab="Relative position to TSS", ylab="Normalized signal", xaxt="n", ylim=y_range)
  points(x,ya, col="blue", type="l")
  axis(side=1, at=c(1, 500, 1000, 1500, 2000), labels=c("-1000", "-500", "0", "+500", "+1000"))
  legend("topright", legend=c("Cortex", "Amygdala"), lty=1, col=c("red", "blue"))
  dev.off()
  
}


## Plot the difference of normlaized signals between amygdala and cortex
gene22_norm <- read.table("Normalized_signals_TSScenetered_sig22genes.txt", header=T, as.is=T, sep="\t")
gene_x <- gene22_norm$index
gene_ya <- apply(as.matrix(gene22_norm[,valid_samples[amyg]]), 1, "mean")
gene_yc <- apply(as.matrix(gene22_norm[,valid_samples[cortex]]), 1, "mean")
gene_diff <- gene_yc - gene_ya

rand_diff.mat <- NULL
for (rep in 1:10){
  ###### Normalized signals TSS centered per replicate (random 22 genes) ##########
  rand22_norm <- read.table(paste("./Normalized_signals_random22genes.rep", rep, ".txt", sep=""), header=T, as.is=T, sep="\t")
  #################################################################################
  
  rand_x <- rand22_norm$index
  rand_ya <- apply(as.matrix(rand22_norm[,valid_samples[amyg]]), 1, "mean")
  rand_yc <- apply(as.matrix(rand22_norm[,valid_samples[cortex]]), 1, "mean")
  rand_diff.mat <- cbind(rand_diff.mat, rand_yc - rand_ya)
}

rand_diff <- apply(as.matrix(rand_diff.mat), 1, "mean")
rand_diff.sd <- apply(as.matrix(rand_diff.mat), 1, "sd")

library("ggplot2")
df <- as.data.frame(cbind(rand_diff, rand_diff.sd))
df$x <- 1:2000
df$under <- df$rand_diff - df$rand_diff.sd
df$over <- df$rand_diff + df$rand_diff.sd

rand_ss <- smooth.spline(rand_x, rand_diff)
rand_under_ss <- smooth.spline(df$x, df$under)
rand_over_ss <- smooth.spline(df$x, df$over)
gene_ss <- smooth.spline(gene_x, gene_diff)

df$rand_diff_smoothed_x <- rand_ss$x
df$rand_diff_smoothed_y <- rand_ss$y
df$rand_diff_smoothed_y_over <- rand_over_ss$y
df$rand_diff_smoothed_y_under <- rand_under_ss$y

df$gene22_smoothed_x <- gene_ss$x
df$gene22_smoothed_y <- gene_ss$y

pdf("Difference_normalized_signals_22genes_wRand10reps_wSDbar.pdf")
ggplot(df, aes(x=x)) +
    geom_line(aes(y=rand_diff_smoothed_y), colour="black") +
    geom_line(aes(y=gene22_smoothed_y), colour="red") +
    geom_ribbon(aes(x=x, ymin= rand_diff_smoothed_y_under, ymax=rand_diff_smoothed_y_over), alpha=0.2, inherit.aes=FALSE, fill="grey") +
    theme_bw() +
    scale_x_continuous(name = "Relative position to TSS\n", breaks=c(1, 500, 1000, 1500, 2000), labels=c("-1000", "-500", "0", "+500", "+1000")) +
    scale_y_continuous(name = "Difference in normalized signal\n") +
    ggtitle("Difference between cortex and amygdala\n") +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    scale_colour_discrete(labels=c("22 genes (adj.p<0.05)", "random 22 genes (10 reps)"))
dev.off()
