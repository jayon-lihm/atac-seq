## 2018.02.08
## Non-overlapping 10% sampling (10 reps)
bamfile="/path/X.bam"
sampleID="X"
tmpdir="/path/tmp/"
out_dir="/output/path/"

echo "START: $(date +"%D"), $(date +"%T")"
original_reads=$(samtools view ${bamfile} | wc -l)

echo "1. Split into 10 data sets"
samtools view -H ${bamfile} > ${tmpdir}original_header.txt
samtools view ${bamfile} | gawk 'BEGIN {srand()}{rdn=rand();
                                 if (rdn<=0.1) print $0 > "'"${tmpdir}${sampleID}.random.1"'";
                                 else if (rdn>0.1 && rdn<=0.2) print $0 > "'"${tmpdir}${sampleID}.random.2"'";
                                 else if (rdn>0.2 && rdn<=0.3) print $0 > "'"${tmpdir}${sampleID}.random.3"'";
                                 else if (rdn>0.3 && rdn<=0.4) print $0 > "'"${tmpdir}${sampleID}.random.4"'";
                                 else if (rdn>0.4 && rdn<=0.5) print $0 > "'"${tmpdir}${sampleID}.random.5"'";
                                 else if (rdn>0.5 && rdn<=0.6) print $0 > "'"${tmpdir}${sampleID}.random.6"'";
                                 else if (rdn>0.6 && rdn<=0.7) print $0 > "'"${tmpdir}${sampleID}.random.7"'";
                                 else if (rdn>0.7 && rdn<=0.8) print $0 > "'"${tmpdir}${sampleID}.random.8"'";
                                 else if (rdn>0.8 && rdn<=0.9) print $0 > "'"${tmpdir}${sampleID}.random.9"'";
                                 else if (rdn>0.9 && rdn<=1) print $0 > "'"${tmpdir}${sampleID}.random.10"'"; }'

echo "2. Peak calling for each data-subset"
for spr in {1..10} ## 10 seperate files
do
    cat ${tmpdir}original_header.txt ${tmpdir}${sampleID}.random.${spr} | samtools view -S - -bh >  ${tmpdir}${sampleID}.random.${spr}.bam
    rm ${tmpdir}${sampleID}.random.${spr} ## Remove samfile
    numRead=$(samtools view ${tmpdir}${sampleID}.random.${spr}.bam | wc -l)
    echo "${sampleID}.random.${spr}.bam :: $numRead"
    macs2 callpeak --nomodel --broad --keep-dup all -g mm -f BAM -t ${tmpdir}${sampleID}.random.${spr}.bam -n ${sampleID}.random_${spr}.broad.macs2 --outdir ${out_dir}
    rm ${tmpdir}${sampleID}.random.${spr}.bam
    gzip ${out_dir}${sampleID}.random_${spr}.broad.macs2_peaks.broadPeak
    gzip ${out_dir}${sampleID}.random_${spr}.broad.macs2_peaks.gappedPeak
    echo "...random split file ${spr} done:: $(date +"%D"), $(date +"%T")"
done

