# atac-seq
Analysis scripts for paper "How low can you go? Calling robust ATAC-seq peaks through read down-sampling"

Authors:  Jayon Lihm1, Sandra Ahrens1, Sara Ballouz1, Hayan Lee2,3, Megan Crow1, Jessica Tollkuhn1, Shane McCarthy1, Bo Li1, W.R. McCombie1, Jesse Gillis1*


Affiliations:  
1 The Stanley Institute for Cognitive Genomics, Cold Spring Harbor Laboratory, Cold Spring Harbor, NY, 11724, USA  
2 Department of Energy Joint Genome Institute, Walnut Creek, CA 94598, USA  
3 Department of Genetics, Stanford University School of Medicine, Stanford University, Stanford, CA 94304, USA  


For questions, please contact Jayon Lihm (jlihm@cshl.edu) or Jesse Gillis (JGillis@cshl.edu).  

===============================================================================

Input and output files (tables, graphs, texts) are located under "files" folder. 
List and pubmed link of 12 studies in our meta analysis is in "List_of_12_studies_in_meta_analysis.txt".  

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
    1) Numbe of accessible TSS, Number of TSS-accessible genes (genes with peaks within TSS regions) - output: Num_TSS_accessible_genes_BrNa.txt,  
    2) Binary table for each TSS region, row: TSS, column: samples, 1 if the TSS is accessible, otherwise 0 - output: TSS_AllSamples_Br.txt, TSS_AllSamples_Na.txt,  
    3) Binary table for each gene, row: gene, column: samples - output: Gene_AllSamples_Br.txt, Gene_AllSamples_Na.txt     
  
    \[Software Requirements\]  
      - Bedtools
      
    NOTE: "Merged_mm9_refSeq_TSS_1kb.txt" is a table of TSS regions. Overlapping TSS regions were merged within the same gene.  
    NOTE: "SRA_sample_info.color_20170913.txt" is a table of 193 sample information used in the meta analysis.  
  
  5. **num_genes_stacked_histogram.R**  
    Generates Figure 2(C); Stacked histogram of number of TSS-accessible genes. Color represents each study.  
    Input files:
      - SRA_sample_info.color_20170913.txt: sample information and color codes
      - Num_TSS_accessible_genes_BrNa.txt: the number of TSS-accessible genes with broad and narrow peaks  
      
    Output file: Number_of_Genes_Distribution_histogram_BrNa.pdf  
    
  6. **downsampling.sh**  
    Downsampling reads to 10%, 20%, ..., 90%, generating ten replicates are generated at each sampling fraciton.  
    Peaks are called right after the sampling and the sampled bam files are removed.  
    Output: ${sample_alias}.f${frac}_r${i}.DefaultBroad.macs2 & ${sample_alias}.f${frac}_r${i}.DefaultNarrow.macs2 where sample_alias is sample ID, frac is sampling fraction, i is a replicate number.  
    
    [Software Requirements]
      - MACS2
      - Samtools  
      
      NOTE: Input bam file and output prefix needs to be set within the script.  
   
   7. **Novel_Lost_Peak_Counts.R**  
      Used for Figure 3(A). Calculate % of novel and lost peaks. Novel peak is a newly appeared peak after sampling and lost peak is a disappeared peak after sampling. Peaks from the original bam file before sampling are called "original peaks". This script runs for each fraction for parallelization.       
      Percentage of novel peaks: (# novel peaks) / (Total # peaks after sampling)  
      Percentage of lost peaks: (# lost peaks) / (Total # peaks before sampling)  
      Output  
      - NumPeaks_NovelLost_AllReps.txt - A table with sample ID, study name, sampling fraction, replicate number, original number of peaks, number of peaks at this sampling, number of novel peaks, number of lost peaks.  
      - pctNovelpeaks_vs_sampling.BrNa.pdf - A summary plot of % of novel peaks per study
      NOTE: Paths for original peak file and sample peak file need to be specified within the script.  
      
      ```
      Rscript ./Novel_Lost_Peak_Counts.R
      ```
      \[Software Requirements\]
        - Bedtools
        
   8. **novel_and_lost_peaks_genes_f90.R**  
      Used for Figure 3(B&C). Generate table per sample with TSS/Gene info + original + f90_r1 + f90_r10.  
      1 if the TSS or gene has overlapping peaks, 0 otherwise.  
      NOTE: Again, the location of peak files from original reads and 90% sampling needs to be specificed within scripts. (multiple times)  
      
      Input files:  
        - SRA_sample_info.color_20170913.txt
        - Merged_mm9_refSeq_TSS_1kb.txt
        - mm9_genename.txt: unique list of gene names and its chromosomal locations
        - Num_TSS_accessible_genes_BrNa.txt, TSS_AllSamples_Br.txt.gz, TSS_AllSamples_Na.txt.gz, 
        
      Output files:
        - INTERMEDIATE FILES
        - Br_outfile_tss <- paste("./", sample_alias, ".TSS_overlaps_Broad.byTSS.txt", sep="")
        - Br_outfile_gene <- paste("./", sample_alias, ".TSS_overlaps_Broad.byGene.txt", sep="")
        - Na_outfile_tss <- paste("./", sample_alias, ".TSS_overlaps_Narrow.byTSS.txt", sep="")
        - Na_outfile_gene <- paste("./", sample_alias, ".TSS_overlaps_Narrow.byGene.txt", sep="")
        - Novel_lost_genes_f90_AllReps.txt: count of novel and lost genes from 90% sampled reads, all 10 reps
        - ...
        - SUMMARY OUTPUTS
        - Average_novel_lost_genes_f90.txt: averaged over 10 replicates
        - histogram_frac_lost_peaks_genes_f90.pdf: plot of novel and lost genes at 90% sampling
        
   9. **Gene_Frequency_Distribution.R**  
      Used for Figure 4(A). Generate a reverse cumulative distribution P(X>=x) graph.  
      NOTE: Peak file needs to be specified within the script.  
      Input files:  
        - SRA_sample_info.color_20170913.txt
        - Merged_mm9_refSeq_TSS.txt, Merged_mm9_refSeq_TSS_1kb.txt, Merged_mm9_refSeq_Gene_1kb.txt, mm9_genename.txt  
        
      Output files:  
        - RevCumulative_Histogram_gene_TSS_peaks.pdf: Reverse cumulative distribution graph
        - peaks_overlap_with_Genes_perGeneSummary.txt: List of genes with frequency of samples at TSS, TSS1kb, Gene1kb
        - CommonGenes_451.txt: List of commonly accessible genes (TSS+/- 1kb regions)
   10. **common_genes_analysis.R**  
      Used for Figure 4(B). Generates a png graph plotting the rank of mean expression vs. the ranked SD, a list of 451 commonly accessible genes, and the results of GO enrichment analysis based on MannWhitney test.  
      GO enrichment analysis: In this script, we perform Mann Whitney test for broader trend and Hypergeometric test for enrichment.  
      Input files:
        - SRA_sample_info.color_20170913.txt
        - Merged_mm9_refSeq_TSS_1kb.txt
        - Gene_numSamples.txt.gz
        - mouse.ranked_mean_sd_expression.txt.gz: Mouse expression data
        - EGAD sources (https://www.bioconductor.org/packages/release/bioc/html/EGAD.html): EGADlite.RData, biogrid.RData, GO.mouse.RData, GO.voc.RData, gene_set_enrichment.r  
       
      Output files:
        - CommonAccessibleGenes_451_RankedGeneExpression.txt
        - All_and_CommonGenes_mean_expression_vs_sd.png
        - MannWhitney_rankTest_CommonGenes_GO_enrichment.txt
   
   11. **non_overlapping_10pct_sampling.sh**  
      Used for our mouse amygdala and cortex samples processing. This script will generate non-overlapping 10 splits of the input bam files, followed by MACS2 peak calling. The paths to bamfile, sampleID, tmp dir, and output dir need to be specified within the script. After the peak calling, the bam files for the sampled reads will be removed to save space.  
      \[Software Requirements\]  
        - MACS2
        - Samtools

   12. **MouseBrain_DiffAcc_ttest.R**  
      Perform Student t-tests for selected genes that are called robustly at least in one sample. Compare fear-conditioning vs control, ErbB4 knock-out vs wildtype, and amygdala vs cortex.  
      Output Figures include Figure 5C, 5D, and Figure S7.  
      
      Input files:
        - MouseBrain_SampleInfo.txt: Sample info
        - Merged_mm9_refSeq_TSS_1kb.txt
        - mm9_genename.txt
        - Peak calling results: MouseBrain_Broad_Original_TSS_Hits.txt, MouseBrain_Broad_f10_TSS_Hits.txt, MouseBrain_Narrow_Original_TSS_Hits.txt, MouseBrain_Narrow_f10_TSS_Hits.txt
        - TSS ATAC-seq signal: MouseBrain_TSS_1k.depth_summary.txt.gz
        
      Output files:
        - T-test results: MouseBrain_DA_Test_DefBroad_original.txt, MouseBrain_DA_Test_DefBroad_f10.txt, MouseBrain_DA_Test_DefNarrow_original.txt, MouseBrain_DA_Test_DefNarrow_f10.txt
        - Figures: MouseBrain_nzd_reads_amyg_vs_cortex.pdf, MouseBrain_nzd_reads_vs_pvalues.pdf
  
