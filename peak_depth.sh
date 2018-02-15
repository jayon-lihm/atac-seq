## 2018.02.08
## Figure 2(b)
bamfile="/path/X.bam"
peak_file="/path/sampleID.broadPeak"

## 1. Get depth for a given bam file
samtools depth ${bamfile} | gzip -c > ${bamfile}.depth.txt.gz

## 2. Extract depth within peak regions
depth_file="${bamfile}.depth.txt.gz" # input
depth_peak="${bamfile}.depth.sample_peaks.txt.gz" # output

zcat ${depth_file} | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $2, $3}' | bedtools intersect -a stdin -b ${peak_file} -wa | gzip -c > ${depth_peak}

## 3. FINAL STEP: Compute the sum of depth per peak
out_file="${peak_file}.depth.txt"

bedtools intersect -a ${peak_file} -b ${depth_peak} -wa -wb | python ./depth_per_peak.py | gzip -c > ${out_file}.gz

