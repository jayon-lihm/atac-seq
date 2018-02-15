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
    
  6. **downsampling.sh**  
    Downsampling reads to 10%, 20%, ..., 90%, generating ten replicates are generated at each sampling fraciton.  
    Peaks are called right after the sampling and the sampled bam files are removed.  
    Output: ${sample_alias}.f${frac}\_r${i}.broad.macs2 where sample_alias is sample ID, frac is sampling fraction, i is a replicate number.  
    \[Software Requirements\]
      - MACS2
      - Samtools  
      
      NOTE: Input bam file and output prefix needs to be set within the script.  
   
   7. **Novel_Lost_Peak_Counts.R**  
      Used for Figure 3(B). Calculate % of novel and lost peaks. Novel peak is a newly appeared peak after sampling and lost peak is a disappeared peak after sampling. Peaks from the original bam file before sampling are called "original peaks". This script runs for each fraction for parallelization.       
      Percentage of novel peaks: (# novel peaks) / (Total # peaks after sampling)  
      Percentage of lost peaks: (# lost peaks) / (Total # peaks before sampling)  
      Output: Sampling_f\<frac\>\_withRep_NumPeaks_AllReps.txt - a table with sample ID, study name, sampling fraction, replicate number, original number of peaks, number of peaks at this sampling, number of novel peaks, number of lost peaks.  
      NOTE: Paths for original peak file and sample peak file need to be specified within the script.  
      
      ```
      Rscript ./Novel_Lost_Peak_Counts.R 30
      ```
      \[Software Requirements\]
        - Bedtools
        
   8. **Summary_Gene_withPeaks_frac1090.R**  
      Used for Figure 3(C). Generate table per sample with TSS/Gene info + original + f10_r1 + ... + f10_r10 + f90_r1 + f90_r10.  
      1 if the TSS or gene has overlapping peaks, 0 otherwise.  
      It takes the index number of sample in the sample info table and temporary output directory ("tmpdir") as arguments.  
      ```
      Rscript ./Summary_Gene_withPeaks_frac1090.R 10 "./"
      ```
      NOTE: Again, the location of peak files for 10% sampling and 90% sampling needs to be specificed within scripts.  
      
      Input files:  
        - ./Merged_mm9_refSeq_TSS_1kb.txt
        - ./mm9_genename.txt: unique list of gene names and its chromosomal locations
        - ./Called_TSS_original_peaks.txt.gz: Table with the list of TSS and all 193 samples for whether a sample has peaks or not.  
        
      
