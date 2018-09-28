## 2018.02.15
## Downsampling to 10%, 20%, ... 90% (with replacement)
## 10 replicates at each sampling
bamfile="/path/X.bam"
out_prefix= "/dir/sampleID" ## prefix for bams from sampled reads
out_dir="./" ## MACS2 peaks output directory


echo "START: $(date +"%D"), $(date +"%T")"
for j in {1..9} ## fraction 10%~90%
do
    (( frac=$j*10 ))
    echo "****fraction $frac %****"
    for i in {1..10} ## repeats 1-10
    do
	seed_num="${RANDOM}"
	echo "replicate: $i"
	echo "seed number: $seed_num"
	samtools view -s ${seed_num}.${frac} $bamfile -b > ${out_prefix}.f${frac}_r${i}.bam
	numRead=$(samtools view ${out_prefix}.f${frac}_r${i}.bam | wc -l)
	echo "number of reads: $numRead"

	## peak calling
	echo "default, broad, $(date +"%D"), $(date +"%T")"
        macs2 callpeak --nomodel --broad --keep-dup all -g mm -f BAM -t ${out_prefix}.f${frac}_r${i}.bam -n ${sample_alias}.f${frac}_r${i}.DefaultBroad.macs2 --outdir ${out_dir}
        echo "default, narrow, $(date +"%D"), $(date +"%T")"
        macs2 callpeak --nomodel --keep-dup all -g mm -f BAM -t ${out_prefix}.f${frac}_r${i}.bam -n ${sample_alias}.f${frac}_r${i}.DefaultNarrow.macs2 --outdir ${out_dir}
	rm ${out_prefix}.f${frac}_r${i}.bam
	echo "...r$i done:: $(date +"%D"), $(date +"%T")"
    done
    echo ""
done

