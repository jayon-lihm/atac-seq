## 2018.02.15
## Calculate percentage of novel and lost peaks
## Percentage of novel peaks: (# new peaks) / (# sampled peaks)
## Percentage of lost peaks: (# peaks disappeared) / (# of original peaks)
## Process per fraction

args <- commandArgs(TRUE)
frac <- args[1]
sampleFile <- "./SRA_sample_info.color_20170913.txt"
geneFile <- "./Merged_mm9_refSeq_TSS_1kb.txt"

sampleTab <- read.table(sampleFile, header=T, as.is=T, sep="\t")
nsample <- dim(sampleTab)[1]

geneTab <- read.table(geneFile, header=T, as.is=T, sep="\t") ## Selected TSS based on strand, merged set (within same name2 field)

origNumPeakFile <- "./Original_number_peaks.20161028.txt"
origNumPeaksTab <- read.table(origNumPeakFile, header=T, as.is=T, sep="\t")

## Merged number of original peaks and sample information to be consistent in ordering of paper & sample_alias
a <- merge(sampleTab, origNumPeaksTab[,c("sample_alias", "numPeaks")], by="sample_alias")

## For all sampling fractions and all replicates
res <- data.frame(paper=rep(a$paper, each=10), sample_alias=rep(a$sample_alias, each=10), fractions=rep(frac, nsample*10), replicates=rep(1:10, dim(a)[1]),
                  orig_numPeaks=rep(a$numPeaks, each=10), numPeaks=-9, novelPeaks=-9, lostPeaks=-9, stringsAsFactors=FALSE)

for (i in 1:dim(sampleTab)[1]){
  cat("f", frac, i, ",")
  paper <- sampleTab$paper[i]
  sample_alias <- sampleTab$sample_alias[i]
  ###### PEAK FILE PATH needs to be specified here #######
  o_peakFile <- paste("./", sample_alias, ".merged_genome.broad.macs2_peaks.broadPeak", sep="")
  ########################################################
  orig_numPeaks <- a$numPeaks[which(a$sample_alias==sample_alias)]
  
  for (r in 1:10){
    ####### SAMPLED PEAK FILE PATH needs to be specified here ########
    sep_peakFile <- paste("./", sample_alias, ".f", frac, "_r", r, ".broad.macs2_peaks.broadPeak.gz", sep="")
    ##################################################################
    
    ## a. number of peaks
    cmd <- paste("zcat ", sep_peakFile, " | wc -l", sep="")
    numPeaks <- as.numeric(system(cmd, intern=TRUE))

    ## b. novel peaks
    cmd <- paste("bedtools intersect -a ", sep_peakFile, " -b ", o_peakFile, " -v | wc -l", sep="")
    novelPeaks <- as.numeric(system(cmd, intern=TRUE))

    ## c. lost peaks
    cmd <- paste("bedtools intersect -a ", o_peakFile, " -b ", sep_peakFile, " -v | wc -l", sep="")
    lostPeaks <- as.numeric(system(cmd, intern=TRUE))

    ## d. For double-checking numbers, count overlapping peaks
    cmd <- paste("bedtools intersect -a ", o_peakFile, " -b ", sep_peakFile, " -u | wc -l", sep="")
    ovlpPeaks_o <- as.numeric(system(cmd, intern=TRUE))
      
    cmd <- paste("bedtools intersect -a ", sep_peakFile, " -b ", o_peakFile, " -u | wc -l", sep="")
    ovlpPeaks_s <- as.numeric(system(cmd, intern=TRUE))
    
    if ( ovlpPeaks_o + lostPeaks != orig_numPeaks ){ cat("ERROR 1 :: Original number of peaks should be the sum of overlapping peaks and lost peaks\n") }
    if ( ovlpPeaks_s + novelPeaks != numPeaks ){ cat("ERROR 2 :: Number of ", frac, "% sampling peaks should be the sum of overlapping peaks and novel peaks\n") }
    
    res$numPeaks[which(res$paper==paper & res$sample_alias==sample_alias & res$replicates==r)] <- numPeaks
    res$novelPeaks[which(res$paper==paper & res$sample_alias==sample_alias & res$replicates==r)] <- novelPeaks
    res$lostPeaks[which(res$paper==paper & res$sample_alias==sample_alias & res$replicates==r)] <- lostPeaks
  }
}

outfile <- paste("./Sampling_f", frac, "_withRep_NumPeaks_AllReps.txt", sep="")
write.table(res, outfile, col.names=T, row.names=F, quote=F, sep="\t")
  

