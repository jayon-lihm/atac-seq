## 2018.02.18
## f10 : called 1 if >=9 replicates have peaks

## 1. FC vs CTRL
## 2. AMYG vs CORTEX - FC design & ERBB4 KO design
## 3. ERBB4 KO vs WT

s_info_file <- "MouseBrain_SampleInfo.txt"
geneFile <- "Merged_mm9_refSeq_TSS_1kb.txt"  ## Selected TSS based on strand, merged set (within same name2 field)


###
sampleTab <- read.table(s_info_file, header=T, as.is=T, sep="\t")

ctrl_id_list <- c("AMYGDALA_CTRL_1", "AMYGDALA_CTRL_2", "CORTEX_CTRL_1", "CORTEX_CTRL_2", "AMYGDALA_ERRB4_WT_1", "AMYGDALA_ERRB4_WT_2", "CORTEX_ERRB4_WT_1", "CORTEX_ERRB4_WT_2")
FC_id_list <- c("AMYGDALA_CTRL_1", "AMYGDALA_CTRL_2", "AMYGDALA_FEAR_1", "AMYGDALA_FEAR_2", "CORTEX_CTRL_1", "CORTEX_CTRL_2", "CORTEX_FEAR_1", "CORTEX_FEAR_2")
KO_id_list <- c("AMYGDALA_ERRB4_KO_1", "AMYGDALA_ERRB4_KO_2", "AMYGDALA_ERRB4_KO_3", "AMYGDALA_ERRB4_WT_1", "AMYGDALA_ERRB4_WT_2",
                "CORTEX_ERRB4_KO_1", "CORTEX_ERRB4_KO_2", "CORTEX_ERRB4_KO_3", "CORTEX_ERRB4_WT_1", "CORTEX_ERRB4_WT_2")

## Gene calling table with 10% non-overlapping sampling
## 1 if >=9 replicates have peaks, 0 if < 9 replicates have peaks
f10_all.res <- read.table("MouseBrain_frac10_ALL.TSS_1kb_genes.txt.gz", header=T, as.is=T, sep="\t")

source("./hyp.geo.test.R")
filter_samples <- c("AMYGDALA_FEAR_2", "AMYGDALA_ERRB4_WT_2")

## Remove two low quality samples
a <- apply(as.matrix(f10_all.res[,setdiff(sampleTab$merged_ID, filter_samples)]), 1, "sum")
## table(a)

pdf("MouseBrain_numGenes_distributions_numSamples16.pdf")
b <- barplot(table(a), xlab="number of samples", ylab="number of genes", xaxt="n")
axis(1, at=b, c(0:16), cex=0.5)
dev.off()

## 1. FC vs CTRL (7 samples)
FC_id_list2 <- setdiff(FC_id_list, filter_samples)

num_FC <- length(grep("FEAR", FC_id_list2)) ## 3
num_CTRL <- length(grep("CTRL", FC_id_list2)) ## 4
FC_colnames <- FC_id_list2[grep("FEAR", FC_id_list2)]
CTRL_colnames <- FC_id_list2[grep("CTRL", FC_id_list2)]

FC.f10_HypGeo.res <- hyp_geo_test_selected(df=f10_all.res, x="FEAR", y="CTRL", colname_list=FC_id_list2, lowB=min(num_FC, num_CTRL), upB=max(num_FC, num_CTRL))
length(which(FC.f10_HypGeo.res$xy_p.value_adj<0.05))
length(which(FC.f10_HypGeo.res$yx_p.value_adj<0.05))

## 2. AMYGDALA vs. CORTEX (16 samples)
TIS_id_list <- unique(c(FC_id_list, KO_id_list))
TIS_id_list2 <- setdiff(TIS_id_list, filter_samples)
num_am <- length(grep("AMYGDALA", TIS_id_list2)) ## 7
num_ct <- length(grep("CORTEX", TIS_id_list2)) ## 9
am_colnames <- TIS_id_list2[grep("AMYGDALA", TIS_id_list2)]
ct_colnames <- TIS_id_list2[grep("CORTEX", TIS_id_list2)]

TIS.f10_HypGeo.res <- hyp_geo_test_selected(df=f10_all.res, x="AMYGDALA", y="CORTEX", colname_list=TIS_id_list2, lowB=min(num_am, num_ct), upB=max(num_am, num_ct))
length(which(TIS.f10_HypGeo.res$xy_p.value_adj<0.05))
length(which(TIS.f10_HypGeo.res$yx_p.value_adj<0.05))

#write.table(TIS.f10_HypGeo.res, "./hyp_geo_test.TISSUE.txt", col.names=T, row.names=F, quote=F, sep="\t")

## 3. ERBB4 KO vs ERBB4 WT (9 samples)
KO_id_list2  <- setdiff(KO_id_list, filter_samples)
num_ko <- length(grep("KO", KO_id_list2)) ## 4
num_wt <- length(grep("WT", KO_id_list2)) ## 5
ko_colnames <- KO_id_list2[grep("KO", KO_id_list2)]
wt_colnames <- KO_id_list2[grep("WT", KO_id_list2)]

KO.f10_HypGeo.res <- hyp_geo_test_selected(df=f10_all.res, x="KO", y="WT", colname_list=TIS_id_list2, lowB=min(num_ko, num_wt), upB=max(num_ko, num_wt))
length(which(KO.f10_HypGeo.res$xy_p.value_adj<0.05))
length(which(KO.f10_HypGeo.res$yx_p.value_adj<0.05))

###################################################################
#### GO enrichment test
load("EGADlite.RData")
load("biogrid.RData")
load("GO.mouse.RData")
load("GO.voc.RData")
source("gene_set_enrichment.r")

geneTab <- read.table(geneFile, header=T, as.is=T, sep="\t")
gene_list <- unique(geneTab$name2)

## remove "IEA"
filt <- GO.mouse[,4] != "IEA"

## intersection of gene_list & GO.mouse
temp <- GO.mouse[filt,]
filt2 <- temp$name %in% gene_list
temp2 <- temp[filt2,]

## make table
go.terms <- unique(temp2[,3])
gene.symb <- unique(temp2[,1])
annotations <- make_annotations(temp2[ ,c(1,3)], gene.symb, go.terms ) ## row: gene, column: GO

## annotations: 17,045 genes & 17,724 GO terms
labels.counts <- colSums(annotations)
filt.GO_names <- colnames(annotations)[which(labels.counts>=10 & labels.counts<=300)]## 5,579 GO terms
annotations2 <- annotations[,filt.GO_names]

####### TEST
sig_genes <-  TIS.f10_HypGeo.res$name2[which( TIS.f10_HypGeo.res$yx_p.value_adj<0.05 )]

enrichment_results <- gene_set_enrichment(sig_genes, annotations2, GO.voc)
res <- enrichment_results[order(enrichment_results[,"pvals.adj[v.f]"]),]
sigGO <- res[,"pvals.adj[v.f]"]<0.05
sig.res <- res[sigGO,]

filt.GO_names <- colnames(annotations)[which(labels.counts>=5 & labels.counts<=100)]## 7,388 GO terms
annotations2 <- annotations[,filt.GO_names]
enrichment_results <- gene_set_enrichment(sig_genes, annotations2, GO.voc)
res <- enrichment_results[order(enrichment_results[,"pvals.adj[v.f]"]),]
sigGO <- res[,"pvals.adj[v.f]"]<0.05
sig.res <- res[sigGO,]

filt.GO_names <- colnames(annotations)[which(labels.counts<=500)]## 17,291 GO terms
annotations2 <- annotations[,filt.GO_names]
enrichment_results <- gene_set_enrichment(sig_genes, annotations2, GO.voc)
res <- enrichment_results[order(enrichment_results[,"pvals.adj[v.f]"]),]
sigGO <- res[,"pvals.adj[v.f]"]<0.05
sig.res <- res[sigGO,]

