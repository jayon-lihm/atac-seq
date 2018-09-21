## 2018.02.08                                                                                                                                                                                           
##1. Qu et al script: Offset "+" strand +4bp, "-" strand -5bp                                                                                                                                           
##2. MACS2: Peak calling on GENOME bam file                                                                                                                                                             

sampleID="X"
out_dir="/dir/sampleID/"
tmpdir="/path/tmpdir/tmp"
input_bam=${out_dir}${sampleID}_genome.sorted.bam
out_bam=${out_dir}${sampleID}_genome.sorted.offset.bam

## 1. Offset positions                                                                                                                                                                                  
## + strand: + 4 bp                                                                                                                                                                                     
## - strand: - 5 bp                                                                                                                                                                                     
## Use Qu et al's script                                                                                                                                                                                
echo "Shift read positions"
samtools view -h ${input_bam} | perl ./shift_sam_bases_JL.pl | samtools view -bhS - > ${out_bam}

## 2. Peak Calling (by MACS2)                                                                                                                                                                           
echo "Peak Calling by MACS2"
bam_in=${out_bam}

## Default Parameters                                                                                                                                                                                   
macs2 callpeak --nomodel --broad --keep-dup all -g mm -f BAM -t ${bam_in} -n ${sampleID}.default_broad.macs2 --outdir ${out_dir}
macs2 callpeak --nomodel --keep-dup all -g mm -f BAM -t ${bam_in} -n ${sampleID}.default_narrow.macs2 --outdir ${out_dir}

## ENCODE parameters                                                                                                                                                                                    
macs2 callpeak --nomodel --shift -37 --extsize 73 -p 0.01 --broad --keep-dup all -g mm -f BAM -t ${bam_in} -n ${sampleID}.encode_broad.macs2 --outdir ${dir}
macs2 callpeak --nomodel --shift -37 --extsize 73 -p 0.01 --keep-dup all -g mm -f BAM -t ${bam_in} -n ${sampleID}.encode_narrow.macs2 --outdir ${dir}

## BAMPE parameters                                                                                                                                                                                     
macs2 callpeak --broad --keep-dup all -g mm -f BAMPE -t ${bam_in} -n ${sampleID}.bampe_broad.macs2  --outdir ${dir}
macs2 callpeak --keep-dup all -g mm -f BAMPE -t ${bam_in} -n ${sampleID}.bampe_narrow.macs2  --outdir ${dir}
