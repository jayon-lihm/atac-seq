## 2018.02.08
##1. Qu et al script: Offset "+" strand +4bp, "-" strand -5bp 
##2. MACS2: Peak calling on GENOME bam file

sampleID="X"
out_dir="/dir/sampleID/"
tmpdir="/path/tmpdir/tmp"
picard_prog_path="/path/picard-tools-1.98/"
file_prefix="${out_dir}${sampleID}.merged"
input_bam=${file_prefix}_genome.sorted.bam 


## 1. Offset positions
## + strand: + 4 bp
## - strand: - 5 bp
## Use Qu et al's script
## input: sam file
## output: sam file
outfile=${file_prefix}_genome.sorted.offset
samtools view -h ${input_bam} | perl ./shift_sam_bases_JL.pl | samtools view -bhS - > ${outfile}.bam

## 2. Peak Calling (by MACS2) :: Parameters set as in Denny et al
echo "3. Peak Calling by MACS2"
bam_in=${outfile}.bam
macs2 callpeak --nomodel --broad --keep-dup all -g mm -f BAM -t ${bam_in} -n ${sample_alias}.merged_genome.broad.macs2 --outdir ${out_dir}

