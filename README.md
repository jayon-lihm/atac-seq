# atac-seq
Analysis scripts for ATAC-seq data

1. **pipeline_paired.sh & pipeline_single.sh**  
  Trim adapters, mapping, remove duplicates, generate .bam and .bai for genome (chr1-19,X,Y) and chrM.  
  Output: ${file_prefix}.MQ30.genome.sorted.bam, ${file_prefix}.MQ30.chrM.sorted.bam  
  \[Software Requirements\]
    - Trimmomatics v0.32
    - Picard v1.88
    - Bowtie2

  NOTE: Software path needs to be set within the scripts.  
  NOTE: FASTQ files need to be gzipped in advance  

2. **peak_calling.sh**  
  Offset read start and end positions and call peaks by MACS2.  
  Input: merged bam file (NOTE: Technical replicates were merged before running this script)  
  Output: adjusted ("+"strand +4bp, "-" strand -5bp) bam file, MACS2 broad peaks
  \[Software Requirements\]
    - samtools
    - MACS2
  
3. **peak_depth.sh**  
  Compute the sum of depth per peak
  \[Software Requirements\]
    - samtools
  
  
