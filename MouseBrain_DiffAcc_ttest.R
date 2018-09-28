### Differential Accessibility T-test for mouse brain samples
## Figure 5C,D
s_info_file <- paste("./MouseBrain_SampleInfo.txt", sep="")
sampleTab <- read.table(s_info_file, header=T, as.is=T, sep="\t")

geneFile <- "./Merged_mm9_refSeq_TSS_1kb.txt"
geneNameFile <- "./mm9_genename.txt" ## 24,361 gene names
geneTab <- read.table(geneFile, header=T, as.is=T, sep="\t") ## Selected TSS based on strand, merged set (within same name2 field)
gene_list <- unique(geneTab$name2)
genename <- read.table(geneNameFile, header=T, as.is=T, sep="\t")

### Manually typed ###
FC_case_id_list <- c("AMYGDALA_FEAR_1", "AMYGDALA_FEAR_2", "CORTEX_FEAR_1", "CORTEX_FEAR_2")
FC_ctrl_id_list <- c("AMYGDALA_CTRL_1", "AMYGDALA_CTRL_2", "CORTEX_CTRL_1", "CORTEX_CTRL_2")
KO_case_id_list <- c("AMYGDALA_ERRB4_KO_1", "AMYGDALA_ERRB4_KO_2", "AMYGDALA_ERRB4_KO_3", "CORTEX_ERRB4_KO_1", "CORTEX_ERRB4_KO_2", "CORTEX_ERRB4_KO_3")
KO_ctrl_id_list <- c("AMYGDALA_ERRB4_WT_1", "AMYGDALA_ERRB4_WT_2", "CORTEX_ERRB4_WT_1", "CORTEX_ERRB4_WT_2")
AM_case_id_list <- c("AMYGDALA_CTRL_1", "AMYGDALA_CTRL_2", "AMYGDALA_FEAR_1", "AMYGDALA_FEAR_2" , "AMYGDALA_ERRB4_KO_1", "AMYGDALA_ERRB4_KO_2", "AMYGDALA_ERRB4_KO_3",
                                          "AMYGDALA_ERRB4_WT_1", "AMYGDALA_ERRB4_WT_2")
AM_ctrl_id_list <- c("CORTEX_FEAR_1", "CORTEX_FEAR_2", "CORTEX_CTRL_1", "CORTEX_CTRL_2", "CORTEX_ERRB4_KO_1", "CORTEX_ERRB4_KO_2", "CORTEX_ERRB4_KO_3", "CORTEX_ERRB4_WT_1", "CORTEX_ERRB4_WT_2")


filter_samples <- c("AMYGDALA_FEAR_2", "AMYGDALA_ERRB4_WT_2")
eligible_samples <- setdiff(sampleTab$merged_ID, filter_samples)

## Read peak calling results
br_org_res <- read.table("MouseBrain_Broad_Original_TSS_Hits.txt", header=T, as.is=T, sep="\t")
br_f10_res <- read.table("MouseBrain_Broad_f10_TSS_Hits.txt", header=T, as.is=T, sep="\t")

na_org_res <- read.table("MouseBrain_Narrow_Original_TSS_Hits.txt", header=T, as.is=T, sep="\t")
na_f10_res <- read.table("MouseBrain_Narrow_f10_TSS_Hits.txt", header=T, as.is=T, sep="\t")

## Read summed depth info within TSS region
TSS_1k_depth <- read.table("MouseBrain_TSS_1k.depth_summary.txt.gz", header=T, as.is=T, sep="\t") ## 26,728 TSS regions

## Normalize TSS_1k_depth
tot_tss1k_depth <- colSums(TSS_1k_depth[,sampleTab$merged_ID])
nzd_tss1k_depth <- TSS_1k_depth[,c("gene.chrom", "gene.start", "gene.end", "gene.name")]

for (i in sampleTab$merged_ID){
    nzd_tss1k_depth[,i] <- TSS_1k_depth[,i]/tot_tss1k_depth[i]
  }

## Aggregate at gene level
nzd_gene_depth <- aggregate( . ~ gene.name, nzd_tss1k_depth[,c("gene.name", sampleTab$merged_ID)], FUN="sum")

## Filter genes on multiple chromosomes
filter_genes <- genename[grep(",", genename$chrom), "name2"] ## 19 genes
genename2 <- genename[which(! genename$name2 %in% filter_genes),] ## 24,342 genes ==> on single chromosome

nzd_gene_depth2 <- nzd_gene_depth[which(! nzd_gene_depth$gene.name %in% filter_genes),]
nzd_gene_depth2 <- merge(genename2, nzd_gene_depth2, by.x="name2", by.y="gene.name")

######## Differential T-test with Original and sampled peaks ################
## Select regions that are called in at least one sample
br_org_test <- br_org_res[which(rowSums(br_org_res[,eligible_samples])>0),] ## 16,462
br_f10_test <- br_f10_res[which(rowSums(br_f10_res[,paste("f10.", eligible_samples, sep="")])>0),] ## 11,521 (all overlapping with original broad peaks)

na_org_test <- na_org_res[which(rowSums(na_org_res[,eligible_samples])>0),] ## 15,715
na_f10_test <- na_f10_res[which(rowSums(na_f10_res[,paste("f10.", eligible_samples, sep="")])>0),] ## 10,439 (all overlapping with original narrow peaks)

br_org_test2 <- br_org_test[which(! br_org_test$name2 %in% filter_genes), ] ## 16,451
br_f10_test2 <- br_f10_test[which(! br_f10_test$name2 %in% filter_genes), ] ## 11,515

na_org_test2 <- na_org_test[which(! na_org_test$name2 %in% filter_genes), ] ## 15,706
na_f10_test2 <- na_f10_test[which(! na_f10_test$name2 %in% filter_genes), ] ## 10,434

######################################################################################
### FUNCTION ###
diff_atac_test <- function(test_genes, atac_signal, groupA, groupB, groupA.name, groupB.name){
    ### test_genes: gene name and chromosome
    ### atac_signal: normalized read depth of ATAC-seq data. gene name, chrom, normalized depth per sample
    ### groupA: sample ID list of one group
    ### groupB: sample ID list of the other group

    tmp <- sapply(test_genes[,1], function(x, test_data=atac_signal, A=groupA, B=groupB){
          data_x <- test_data[which(test_data[,1]==x), A];
              data_y <- test_data[which(test_data[,1]==x), B];
              res <- t.test(data_x, data_y);
              return(c(rowMeans(data_x), rowMeans(data_y), res$statistic, res$p.value))
        })

      test_genes[, c( paste("mean_", c(groupA.name, groupB.name), sep=""), "t.stat", "p.value")] <- t(tmp)
      test_genes$adj.p.value <- p.adjust(test_genes$p.value, method="BH")
      test_genes$fold_change_AB <- test_genes[,paste("mean_", groupA.name, sep="")] / test_genes[,paste("mean_", groupB.name, sep="")]
      return(test_genes)
  }

#### TEST Three conditions: Amygdala vs cortex, fear vs control, and ERBB4 KO vs WT
#### A. AMYGDALA vs CORTEX
AM_case_id_list <- c("AMYGDALA_CTRL_1", "AMYGDALA_CTRL_2", "AMYGDALA_FEAR_1", "AMYGDALA_FEAR_2",
                                          "AMYGDALA_ERRB4_KO_1", "AMYGDALA_ERRB4_KO_2", "AMYGDALA_ERRB4_KO_3", "AMYGDALA_ERRB4_WT_1", "AMYGDALA_ERRB4_WT_2")
AM_ctrl_id_list <- c("CORTEX_FEAR_1", "CORTEX_FEAR_2", "CORTEX_CTRL_1", "CORTEX_CTRL_2", "CORTEX_ERRB4_KO_1", "CORTEX_ERRB4_KO_2", "CORTEX_ERRB4_KO_3", "CORTEX_ERRB4_WT_1", "CORTEX_ERRB4_WT_2")

amygdala <- setdiff(AM_case_id_list, filter_samples) ## 7
cortex <- setdiff(AM_ctrl_id_list, filter_samples) ## 9

br_orig_tissue <- diff_atac_test(br_org_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=amygdala, groupB=cortex, groupA.name="amygdala", groupB.name="cortex")
br_f10_tissue <- diff_atac_test(br_f10_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=amygdala, groupB=cortex, groupA.name="amygdala", groupB.name="cortex")

na_orig_tissue <- diff_atac_test(na_org_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=amygdala, groupB=cortex, groupA.name="amygdala", groupB.name="cortex")
na_f10_tissue <- diff_atac_test(na_f10_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=amygdala, groupB=cortex, groupA.name="amygdala", groupB.name="cortex")

length(which(br_orig_tissue$adj.p.value<0.05)) ## 727
length(which(br_f10_tissue$adj.p.value<0.05)) ## 402
length(which(na_orig_tissue$adj.p.value<0.05)) ## 723
length(which(na_f10_tissue$adj.p.value<0.05)) ## 274

## Output Results
write.table(br_orig_tissue, "MouseBrain_DA_Test_DefBroad_original.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(br_f10_tissue, "MouseBrain_DA_Test_DefBroad_f10.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(na_orig_tissue, "MouseBrain_DA_Test_DefNarrow_original.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(na_f10_tissue, "MouseBrain_DA_Test_DefNarrow_f10.txt", col.names=T, row.names=F, quote=F, sep="\t")

## pvalue vs normalized reads
orig_sig_genes <- na_orig_tissue[which(na_orig_tissue$adj.p.value<0.05),]
f10_sig_genes <- na_f10_tissue[which(na_f10_tissue$adj.p.value<0.05),]

orig_tmp <- merge(orig_sig_genes, nzd_gene_depth2[which(nzd_gene_depth2$name2 %in% orig_sig_genes$name2),], by=c("name2", "chrom"))
f10_tmp <- merge(f10_sig_genes, nzd_gene_depth2[which(nzd_gene_depth2$name2 %in% f10_sig_genes$name2),], by=c("name2", "chrom"))

orig_amyg <- rowMeans(orig_tmp[,amygdala])
orig_cort <- rowMeans(orig_tmp[,cortex])
f10_amyg <- rowMeans(f10_tmp[,amygdala])
f10_cort <- rowMeans(f10_tmp[,cortex])

### Figure 5C
pdf("MouseBrain_nzd_reads_amyg_vs_cortex.pdf")
plot(log(orig_amyg), log(orig_cort), col="black", xlim=log(c(ymin, ymax)), pch=19, main="log(normalized reads) of cortex vs amygdala")
points(log(f10_amyg), log(f10_cort), col="red", xlim=log(c(ymin, ymax)), pch=19)
abline(a=0, b=1, col="blue", lwd=2, lty="dashed")
legend("topleft", pch=19, col=c("black", "red"), legend=c("original", "sampled"))
dev.off()

## Supplementary Figure
pdf("MouseBrain_nzd_reads_vs_pvalues.pdf")
ymax <- max(orig_amyg, orig_cort, f10_amyg, f10_cort)
ymin <- min(orig_amyg, orig_cort, f10_amyg, f10_cort)

xmax <- max(orig_tmp$p.value, f10_tmp$p.value)
xmin <- min(orig_tmp$p.value, f10_tmp$p.value)

plot(-log(orig_tmp$p.value), log(orig_amyg), col="black", xlim=-log(c(xmax, xmin)), ylim=log(c(ymin, ymax)), pch=19,
          main="avg of normalized reads in amygdala", xlab="-log(pvalue)", ylab="log(nzd reads)")
points(-log(f10_tmp$p.value), log(f10_amyg), col="red", pch=19)
legend("topright", pch=19, col=c("black", "red"), legend=c("original amygdala", "sampled amygdala"))

plot(-log(orig_tmp$p.value), log(orig_cort), col="black", xlim=-log(c(xmax, xmin)), ylim=log(c(ymin, ymax)), pch=19,
          main="avg of normalized reads in cortex", xlab="-log(pvalue)", ylab="log(nzd reads)")
points(-log(f10_tmp$p.value), log(f10_cort), col="blue", pch=19)
legend("topright", pch=19, col=c("black", "blue"), legend=c("original cortex", "sampled cortex"))
dev.off()

#### B. Fear conditioning vs control
FC_case_id_list <- c("AMYGDALA_FEAR_1", "AMYGDALA_FEAR_2", "CORTEX_FEAR_1", "CORTEX_FEAR_2")
FC_ctrl_id_list <- c("AMYGDALA_CTRL_1", "AMYGDALA_CTRL_2", "CORTEX_CTRL_1", "CORTEX_CTRL_2")

fear <- setdiff(FC_case_id_list, filter_samples) ## 3
control <- setdiff(FC_ctrl_id_list, filter_samples) ## 4

br_orig_fear <- diff_atac_test(br_org_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=fear, groupB=control, groupA.name="fear", groupB.name="control")
br_f10_fear <- diff_atac_test(br_f10_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=fear, groupB=control, groupA.name="fear", groupB.name="control")

na_orig_fear <- diff_atac_test(na_org_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=fear, groupB=control, groupA.name="fear", groupB.name="control")
na_f10_fear <- diff_atac_test(na_f10_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=fear, groupB=control, groupA.name="fear", groupB.name="control")

length(which(br_orig_fear$adj.p.value<0.05)) ## 0
length(which(br_f10_fear$adj.p.value<0.05)) ## 0
length(which(na_orig_fear$adj.p.value<0.05)) ## 0
length(which(na_f10_fear$adj.p.value<0.05)) ## 0

#### C. KO vs WT
KO_case_id_list <- c("AMYGDALA_ERRB4_KO_1", "AMYGDALA_ERRB4_KO_2", "AMYGDALA_ERRB4_KO_3", "CORTEX_ERRB4_KO_1", "CORTEX_ERRB4_KO_2", "CORTEX_ERRB4_KO_3")
KO_ctrl_id_list <- c("AMYGDALA_ERRB4_WT_1", "AMYGDALA_ERRB4_WT_2", "CORTEX_ERRB4_WT_1", "CORTEX_ERRB4_WT_2")

ko <- setdiff(KO_case_id_list, filter_samples) ## 6
wt <- setdiff(KO_ctrl_id_list, filter_samples) ## 3

br_orig_ko <- diff_atac_test(br_org_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=ko, groupB=wt, groupA.name="ko", groupB.name="wt")
br_f10_ko <- diff_atac_test(br_f10_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=ko, groupB=wt, groupA.name="ko", groupB.name="wt")

na_orig_ko <- diff_atac_test(na_org_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=ko, groupB=wt, groupA.name="ko", groupB.name="wt")
na_f10_ko <- diff_atac_test(na_f10_test2[,c("name2", "chrom")], atac_signal=nzd_gene_depth2, groupA=ko, groupB=wt, groupA.name="ko", groupB.name="wt")

length(which(br_orig_ko$adj.p.value<0.05)) ## 0
length(which(br_f10_ko$adj.p.value<0.05)) ## 0
length(which(na_orig_ko$adj.p.value<0.05)) ## 0
length(which(na_f10_ko$adj.p.value<0.05)) ## 0


