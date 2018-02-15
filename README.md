# atac-seq
Analysis scripts for ATAC-seq data

Software Requirements
- Trimmomatics v0.32
- Picard v1.88
- Bowtie2

1. pipeline_paired.sh & pipeline_single.sh
Trim adapters, mapping, remove duplicates, generate .bam and .bai for genome (chr1-19,X,Y) and chrM.  
NOTE: FASTQ files need to be gzipped in advance  

a) "./pipeline_paired.sh ${file_prefix}" will process ${file_prefix}_1.fastq.gz & ${file_prefix}_2.fastq.gz.  
b) "./pipeline_single.sh ${file_prefix}" will process ${file_prefix}.fastq.gz.  

Output: ${file_prefix}.MQ30.genome.sorted.bam, ${file_prefix}.MQ30.chrM.sorted.bam  
