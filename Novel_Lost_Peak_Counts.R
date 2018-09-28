## Calculate percentage of novel and lost peaks
## Percentage of novel peaks: (# new peaks) / (# sampled peaks)
## Percentage of lost peaks: (# peaks disappeared) / (# of original peaks)

library("openintro")
library("gdata")
sampleFile <- "./SRA_sample_info.color_20170913.txt"
sampleTab <- read.table(sampleFile, header=T, as.is=T, sep="\t")

fd <- toupper(as.hexmode(150)) # 64 in hexa == 100 in decimal
path <- "./"

################################
## Original number of peaks
outFile <- paste(path, "Original_number_peaks.txt", sep="")
orig_peaks <- sampleTab[, c("paper", "sample_alias")]
orig_peaks$default_br <- 0
orig_peaks$default_na <- 0

for (i in 1:dim(sampleTab)[1]){
  cat(i, ",")
  paper <- sampleTab$paper[i]
  sample_alias <- sampleTab$sample_alias[i]

  ###### BROAD PEAK FILE PATH needs to be specified here #######
  peakFile <- paste("./", sample_alias, ".default_broad.macs2_peaks.broadPeak", sep="")
  ###############################################################################
  cmd <- paste("cat ", peakFile, " | wc -l", sep="")
  orig_peaks$default_br[i] <- as.numeric(system(cmd, intern=TRUE))

  ###### NARROW PEAK FILE PATH needs to be specified here #######
  peakFile <- paste("./", sample_alias, ".default_narrow.macs2_peaks.narrowPeak", sep="")
  ###############################################################################
  cmd <- paste("cat ", peakFile, " | wc -l", sep="")
  orig_peaks$default_na[i] <- as.numeric(system(cmd, intern=TRUE))
}

write.table(orig_peaks, outFile, col.names=T, row.names=F, quote=F, sep="\t")

################################
## Make data table template for novel/lost peak counts
res <- NULL
for (i in 1:dim(sampleTab)[1]){
  paper <- sampleTab$paper[i]
  sample_alias <- sampleTab$sample_alias[i]
  tmp.res <- data.frame(paper=rep(paper, 9*10), sample_alias=rep(sample_alias, 9*10), frac=rep(seq(10,90, by=10), each=10), repl=rep(1:10, 9))
  res <- rbind(res, tmp.res)
}

################################ 
## Number of Peaks & number of Novel Peaks to original peak files for DEFAULT BROAD and NARROW peaks
outFile <- paste(path, "NumPeaks_NovelLost_AllReps.txt", sep="")
res$numPeaks_def_br <- 0
res$numNovelPeaks_def_br <- 0
res$numLostPeaks_def_br <- 0

res$numPeaks_def_na <- 0
res$numNovelPeaks_def_na <- 0
res$numLostPeaks_def_na <- 0

for (i in 1:dim(sampleTab)[1]){
  paper <- sampleTab$paper[i]
  sample_alias <- sampleTab$sample_alias[i]
  ####### Original Broad and Narrow Peak Files need to be specified here ########
  br_peakFile <- paste("./", sample_alias, ".default_broad.macs2_peaks.broadPeak", sep="")
  na_peakFile <- paste("./", sample_alias, ".default_narrow.macs2_peaks.narrowPeak", sep="")
  ###############################################################################

  ### With Broad Peaks
  o_peakFile <- br_peakFile
  for (f in seq(10,90,by=10)){
    cat(i, "\tbroad\t", sample_alias, "\t", f, "\n")
    for (r in 1:10){
      ####### Sampled Broad Peak Files need to be specified here ########
      s_peakFile <- paste("./", sample_alias, ".f", f, "_r", r, ".DefaultBroad.macs2_broadPeak.gz", sep="")
      ###############################################################################
      cmd <- paste("zcat ", s_peakFile, " | wc -l", sep="")
      res$numPeaks_def_br[which(res$paper==paper & res$sample_alias==sample_alias & res$frac==f & res$repl==r)] <- as.numeric(system(cmd, intern=TRUE))

      cmd <- paste("bedtools intersect -a ", s_peakFile, " -b ", o_peakFile, " -v | wc -l", sep="")
      res$numNovelPeaks_def_br[which(res$paper==paper & res$sample_alias==sample_alias & res$frac==f & res$repl==r)]  <- as.numeric(system(cmd, intern=TRUE))

      cmd <- paste("bedtools intersect -a ", o_peakFile, " -b ", s_peakFile, " -v | wc -l", sep="")
      res$numLostPeaks_def_br[which(res$paper==paper & res$sample_alias==sample_alias & res$frac==f & res$repl==r)]  <- as.numeric(system(cmd, intern=TRUE))
    }
  }

  ### With Narrow Peaks
  o_peakFile <- na_peakFile
  for (f in seq(10,90,by=10)){
    cat(i, "\tnarrow\t", sample_alias, "\t", f, "\n")
    for (r in 1:10){
      ####### Sampled Narrow Peak Files need to be specified here ########
      s_peakFile <- paste(path, paper, "/", sample_alias, "/", sample_alias, ".f", f, "_r", r, ".DefaultNarrow.macs2_narrowPeak.gz", sep="")
      ###############################################################################
      cmd <- paste("zcat ", s_peakFile, " | wc -l", sep="")
      res$numPeaks_def_na[which(res$paper==paper & res$sample_alias==sample_alias & res$frac==f & res$repl==r)] <- as.numeric(system(cmd, intern=TRUE))
      
      cmd <- paste("bedtools intersect -a ", s_peakFile, " -b ", o_peakFile, " -v | wc -l", sep="")
      res$numNovelPeaks_def_na[which(res$paper==paper & res$sample_alias==sample_alias & res$frac==f & res$repl==r)]  <- as.numeric(system(cmd, intern=TRUE))

      cmd <- paste("bedtools intersect -a ", o_peakFile, " -b ", s_peakFile, " -v | wc -l", sep="")
      res$numLostPeaks_def_na[which(res$paper==paper & res$sample_alias==sample_alias & res$frac==f & res$repl==r)]  <- as.numeric(system(cmd, intern=TRUE))
    }
  }
}

write.table(res, outFile, col.names=T, row.names=F, quote=F, sep="\t")


################ Summary and Plot #####################
## Summary per Fraction
inFile <- paste(path, "NumPeaks_NovelLost_AllReps.txt.gz", sep="")
orig_file <- paste(path, "Original_number_peaks.txt", sep="")
orig_peaks <- read.table(orig_file, header=T, as.is=T, sep="\t")
res <- read.table(inFile, header=T, as.is=T, sep="\t")

## Novel Peaks
res$pctNovelPeaks_def_br <- res$numNovelPeaks_def_br/res$numPeaks_def_br * 100
res$pctNovelPeaks_def_na <- res$numNovelPeaks_def_na/res$numPeaks_def_na * 100

## summary over replicates
res2 <- aggregate(cbind(pctNovelPeaks_def_br, pctNovelPeaks_def_na) ~ paper + sample_alias + frac, data=res, FUN="mean")
res2 <- res2[order(res2$paper, res2$sample_alias, res2$frac),]

## summary over fraction (per study)
res3 <- aggregate(cbind(pctNovelPeaks_def_br, pctNovelPeaks_def_na) ~ paper + frac, data=res2, FUN="mean")
res3 <- res3[order(res3$paper, res3$frac),]

##### PLOTTING
max_pct_br <- max(res3$pctNovelPeaks_def_br)
min_pct_br <- min(res3$pctNovelPeaks_def_br)

max_pct_na <- max(res3$pctNovelPeaks_def_na)
min_pct_na <- min(res3$pctNovelPeaks_def_na)

colTab <- sampleTab[,c("paper", "color")]
colTab <- colTab[which(! duplicated(colTab)),]
colTab$fd_col <- sapply(colTab$col, "fadeColor", fade=fd)

res4 <- merge(res3, colTab, by="paper")
res4 <- res4[order(res4$paper, res4$frac),]

## spline
spl_res <- data.frame(frac=seq(10,90,by=10))
spl_res$all_avg.pctNovelPeaks_br <- -9
for (i in 1:9){
  f <- i*10
  spl_res$all_avg.pctNovelPeaks_br[i] <- mean(res4$pctNovelPeaks_def_br[which(res4$frac==f)])
  spl_res$all_avg.pctNovelPeaks_na[i] <- mean(res4$pctNovelPeaks_def_na[which(res4$frac==f)])
}

ss_br <- smooth.spline(spl_res$frac, spl_res$all_avg.pctNovelPeaks_br)
ss_na <- smooth.spline(spl_res$frac, spl_res$all_avg.pctNovelPeaks_na)

pdf("./pctNovelpeaks_vs_sampling.BrNa.pdf")
## Broad
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(1, type="n", xlab="Sampling %", ylab="Percentage of Novel Peaks", main="Avg Percentage of Novel Peaks (Broad)", xlim=c(10, 90), ylim=c(0, max_pct_br+2), xaxt="n")

for (i in unique(sampleTab$paper)){
  paper <- i
  tmp <- res4[which(res4$paper==paper),]
  points(tmp$frac, tmp$pctNovelPeaks_def_br, col=tmp$fd_col, type="b")
}
lines(ss_br, col="red", lwd=2)
axis(side=1, at=seq(10, 90, by=10), labels=seq(10, 90, by=10))
legend("topright", colTab$paper, col=colTab$fd_col, pch=1,  inset=c(-0.3,0))

## Narrow
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(1, type="n", xlab="Sampling %", ylab="Percentage of Novel Peaks", main="Avg Percentage of Novel Peaks (Narrow)", xlim=c(10, 90), ylim=c(0, max_pct_na+2), xaxt="n")

for (i in unique(sampleTab$paper)){
  paper <- i
  tmp <- res4[which(res4$paper==paper),]
  points(tmp$frac, tmp$pctNovelPeaks_def_na, col=tmp$fd_col, type="b")
}
lines(ss_na, col="red", lwd=2)
axis(side=1, at=seq(10, 90, by=10), labels=seq(10, 90, by=10))
legend("topright", colTab$paper, col=colTab$fd_col, pch=1,  inset=c(-0.3,0))

dev.off()


