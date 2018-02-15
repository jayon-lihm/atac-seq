#!/usr/bin/env python
import sys
import numpy

## Figure 2(b)
## Used within "peak_depth.sh"

##bedtools intersect -a /seq/jlihm/ATAC-seq/paper/data_processed/LaraAstiaso/GSM1463172/GSM1463172.merged_genome.broad.macs2_peaks.broadPeak -b /seq/jlihm/ATAC-seq/paper/data_processed/LaraAstiaso/GSM1463172/GSM1463172.depth.sample_peaks.txt.gz -wa -wb |  python ~/ATAC-seq/paper/test.py

old_name="."
dp_vec=[]
print "chrom" + "\t" + "start" + "\t" + "end" + "\t" + "peak_name" + "\t" + "sum_depth" + "\t" + "mean_depth" + "\t" + "std_depth"
for line in sys.stdin : 
    fields=line.strip().split()
    peak_chrom=fields[0]
    peak_start=fields[1]
    peak_end=fields[2]
    peak_name=fields[3]
    peak_x1=fields[4]
    peak_x2=fields[5]
    peak_x3=fields[6]
    peak_x4=fields[7]
    peak_x5=fields[8]
    depth_chrom=fields[9]
    depth_start=fields[10]
    depth_end=fields[11]
    depth=int(fields[12])
     
    if old_name==".":
        old_name=peak_name
        old_chrom=peak_chrom
        old_start=peak_start
        old_end=peak_end
        dp_vec.append(depth)
    else:
        if (old_name==peak_name):
            dp_vec.append(depth)
        else:
            sum_DP=sum(dp_vec)
            sd_DP=numpy.round(numpy.std(dp_vec), decimals=2)
            mean_DP=numpy.round(numpy.mean(dp_vec), decimals=2)
            print old_chrom + "\t" + old_start + "\t" + old_end + "\t" + old_name + "\t" + str(sum_DP) + "\t" + str(mean_DP) + "\t" + str(sd_DP)

            old_name=peak_name
            old_chrom=peak_chrom
            old_start=peak_start
            old_end=peak_end
            dp_vec=[depth]




sum_DP=sum(dp_vec)
sd_DP=numpy.around(numpy.std(dp_vec), decimals=2)
mean_DP=numpy.around(numpy.mean(dp_vec), decimals=2)
print old_chrom + "\t" + old_start + "\t" + old_end + "\t" + old_name + "\t" + str(sum_DP) + "\t" + str(mean_DP) + "\t" + str(sd_DP)

