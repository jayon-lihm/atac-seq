# atac-seq
Analysis scripts for ATAC-seq data

1. **pipeline_paired.sh & pipeline_single.sh**  
  Trim adapters, mapping, remove duplicates, generate .bam and .bai for genome (chr1-19,X,Y) and chrM.  
  Output: ${file_prefix}.MQ30.genome.sorted.bam, ${file_prefix}.MQ30.chrM.sorted.bam  
  \[Software Requirements\]
    - Trimmomatics v0.32
    - Picard v1.88
    - Bowtie2
    - Samtools

    NOTE: Software path needs to be set within the scripts.  
    NOTE: FASTQ files need to be gzipped in advance  

2. **peak_calling.sh**  
  Offset read start and end positions and call peaks by MACS2.  
  Input: merged bam file (NOTE: Technical replicates were merged before running this script)  
  Output: adjusted ("+"strand +4bp, "-" strand -5bp) bam file, MACS2 broad peaks  
  \[Software Requirements\]
    - Samtools
    - MACS2
  
3. **peak_depth.sh**  
  Compute the sum of depth per peak. Used for Figure 2(B).  
  \[Software Requirements\]
    - Samtools
    - Bedtools
    - Python
 4. **count_TSS_accessible_genes.R**  
  Used for Figure 2(C). Generate three tables:  
    1) Numbe of accessible TSS, Number of TSS-accessible genes (genes with peaks within TSS regions) - output: Num_TSS_accessible_genes.txt,  
    2) Binary table for each TSS region, row: TSS, column: samples, 1 if the TSS is accessible, otherwise 0 - output: TSS_numSamples.txt,  
    3) Binary table for each gene, row: gene, column: samples - output: Gene_numSamples.txt   
  
    \[Software Requirements\]  
      - Bedtools
      
    NOTE: "Merged_mm9_refSeq_TSS_1kb.txt" is a table of TSS regions. Overlapping TSS regions were merged within the same gene.  
    NOTE: "SRA_sample_info.color_20170913.txt" is a table of 193 sample information used in the meta analysis.  
  
  5. **num_genes_stacked_histogram.R**  
    Generates Figure 2(C); Stacked histogram of number of TSS-accessible genes. Color represents each study.  
    Input: "SRA_sample_info.color_20170913.txt" for sample information and color codes and "Num_TSS_accessible_genes.txt" for the number of TSS-accessible genes  
    Output: Number_of_Genes_Distribution_histogram.pdf  
    
