#Loading libraries
library(Rsubread)
library(dplyr)
library(tidyr)
library(tibble)
library(DESeq2)
library(plyr)
library(purrr)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(RColorBrewer)
library(grid)

#Creating output directories
dir.create("Results")
dir.create("Results/EWAT")
dir.create("Results/IWAT")
dir.create("Results/RWAT")

# Load Metadata

Metadata <- read.csv(file = "data/Metadata.csv")

Metadata <- Metadata %>%  #mutate(Tissue_And_Diet = revalue(as.factor(paste(tissue, Diet, sep = "_")), 
                                                           # c("epididymal white adipose tissue_45% high-fat diet" = "EWAT_45",
                                                           #   "epididymal white adipose tissue_standard rodent diet" = "EWAT_STANDARD",
                                                           #   "inguinal white adipose tissue_45% high-fat diet" = "IWAT_45",
                                                           #   "inguinal white adipose tissue_standard rodent diet" = "IWAT_STANDARD",
                                                           #   "retroperitoneal white adipose tissue_45% high-fat diet" = "RWAT_45",
                                                           #   "retroperitoneal white adipose tissue_standard rodent diet" = "RWAT_STANDARD"))) %>% 
  # filter(tissue == "epididymal white adipose tissue") %>% 
  mutate(Diet = revalue(as.factor(Diet),c(
    "45% high-fat diet" = "45_High_Fat_Diet",
    "standard rodent diet" = "Standard"
    )
    )) %>%
  mutate(tissue = revalue(as.factor(tissue), c(
    "epididymal white adipose tissue" = "EWAT",
    "inguinal white adipose tissue" = "IWAT",
    "retroperitoneal white adipose tissue" = "RWAT"
  )))

# Counting features in all alignments performed

f <- c("data/FilteredReads/SRR6984609_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984610_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984611_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984612_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984613_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984614_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984615_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984616_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984617_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984618_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984619_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984620_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984621_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984622_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984623_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984624_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984625_trimmed_sorted_unmapped_reads_filtered.bam",
       "data/FilteredReads/SRR6984626_trimmed_sorted_unmapped_reads_filtered.bam")

DF <- featureCounts(files=f,
                      annot.ext="data/genome/grcm38/Mus_musculus.GRCm38.102.gtf",
                      isGTFAnnotationFile=TRUE,
                      GTF.featureType="exon",
                      GTF.attrType="gene_id")


# Saving counts in Rda file

save(DF, file = "data/Features.Rda")

# Creating DSEq2 object
design_formula <- ~ tissue + Diet + tissue:Diet

dds <- DESeqDataSetFromMatrix(countData = DF$counts, colData = Metadata, design = design_formula)
# Removing genes with low counts
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]
# Creating DesqObject
dds <- DESeq(dds)

save(dds, file = "data/DDS.Rda")

################################################################################
# Exploratory analysis
################################################################################

# Stabilising RNA-Seq variance for visualisation
vsd <- vst(dds, blind=FALSE)

# PCA components 1 and 2
PCA <- plotPCA(vsd, intgroup=c("Diet", "tissue")) + scale_color_viridis_d()

ggsave(PCA, file = "Results/PCA.pdf", width = 10, height = 10)
ggsave(PCA, file = "Results/PCA.png", width = 10, height = 10)
save(PCA, file = "data/PCA.Rda")

# Heatmap

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

colnames(vsd) <- stringr::str_remove(colnames(vsd), "_trimmed_sorted_unmapped_reads_filtered.bam")

pdf("Results/Pheatmap.pdf", width = 13, height = 7)


p <- pheatmap(assay(vsd)[select,], 
         cluster_rows=FALSE, 
         show_rownames=TRUE,
         cluster_cols=FALSE,
         cellheight = 12,
         cellwidth = 20)

grid.text("Samples", x=0.4, y=0.16, gp=gpar(fontsize=12))
grid.text("Top 20 Genes with more counts",y=0.56, x= 0.2, rot=90, gp=gpar(fontsize=12))

p

dev.off()

################################################################################
# EWAT
################################################################################

res <- results(dds, name = "Diet_Standard_vs_45_High_Fat_Diet",alpha = 0.05)

# Getting differentialy expressed genes

sum(res$padj < 0.05, na.rm=TRUE)
resOrdered <- res[order(res$padj),]

# PLotting top 10 differentially expressed genes

top10_Genes <- rownames(resOrdered)[1:10]

for (i in top10_Genes){
  pdf(file.path("Results/EWAT", paste0(i, ".pdf")))
    plotCounts(dds, gene=i, intgroup="Diet")
  dev.off()
}

print(paste("Gene with lowest pvalue:", top10_Genes[1]))

# Volcano plot
pdf("Results/EWAT/VolcanoPlot.pdf", width = 10, 8)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano plot of DESeq2 EWAT',
                subtitle = 'Differential expression analysis',
                caption = 'Log2 fold change vs. p-value', 
                drawConnectors = TRUE,
                widthConnectors = 1, boxedLabels = TRUE,
                selectLab = top10_Genes)

dev.off()

# Geting significant differentially expressed genes

res_filtered <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

# GO Enrichment
go_results <- enrichGO(gene = rownames(res_filtered), OrgDb = org.Mm.eg.db,
                       ont = "BP", pvalueCutoff = 0.05, readable = TRUE,
                       universe = rownames(dds),
                       qvalueCutoff = 0.2, keyType = "ENSEMBL")

# GO Group

go_group_BP <- groupGO(gene = rownames(res_filtered), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                    ont = "BP", readable = TRUE)

# Barplot enrichment analysis

pdf("Results/EWAT/EnrichmentGO.pdf", width = 10, height = 10)

Bar_Ewat <- barplot(go_results, 
        drop = TRUE, 
        showCategory = 40,
        font.size = 6)

Bar_Ewat

dev.off()

save(Bar_Ewat, file = "data/EWAT_GO_Barplot.Rda")

# dotplot enrichment analysis

pdf("Results/EWAT/EnrichmentGO_Dotplot.pdf", width = 10, height = 10)

dotplot(go_results, orderBy = "p.adjust", showCategory = 40)

dev.off()



# CSV GO enrichment
write.csv(as.data.frame(go_results), file = "Results/EWAT/Top10_differentially_expressed_genes_GO.csv")

save(go_results, file = "data/EWAT_GO.Rda")

# Barplot and csv for GO groups

pdf("Results/EWAT/GO_BP_Group.pdf")
barplot(go_group_BP, drop = TRUE)
dev.off()

write.csv(as.data.frame(go_group_BP), file = "Results/EWAT/Group_BP.csv")

################################################################################
# IWAT
################################################################################

res <- results(dds, contrast = list(c("Diet_Standard_vs_45_High_Fat_Diet", 
                                      "tissueIWAT.DietStandard")),
               alpha = 0.05)

# Getting differentialy expressed genes

sum(res$padj < 0.05, na.rm=TRUE)
resOrdered <- res[order(res$padj),]

# Plotting top 10 differentially expressed genes

top10_Genes <- rownames(resOrdered)[1:10]

for (i in top10_Genes){
  pdf(file.path("Results/IWAT", paste0(i, ".pdf")))
  plotCounts(dds, gene=i, intgroup="Diet")
  dev.off()
}

print(paste("Gene with lowest pvalue:", top10_Genes[1]))

# Volcano plot
pdf("Results/IWAT/VolcanoPlot.pdf", width = 10, 8)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano plot of DESeq2 EWAT',
                subtitle = 'Differential expression analysis',
                caption = 'Log2 fold change vs. p-value', 
                drawConnectors = TRUE,
                widthConnectors = 1, boxedLabels = TRUE,
                selectLab = top10_Genes)

dev.off()

# Geting significant differentially expressed genes

res_filtered <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

# GO Analysis
go_results <- enrichGO(gene = rownames(res_filtered), OrgDb = org.Mm.eg.db,
                       ont = "BP", pvalueCutoff = 0.05, readable = TRUE,
                       universe = rownames(dds),
                       qvalueCutoff = 0.2, keyType = "ENSEMBL")

go_group_BP <- groupGO(gene = rownames(res_filtered), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                       ont = "BP", readable = TRUE)

# Barplot enrichment analysis

pdf("Results/IWAT/EnrichmentGO.pdf", width = 10, height = 10)

Bar_Iwat <- barplot(go_results, 
        drop = TRUE, 
        showCategory = 40,
        font.size = 6)

Bar_Iwat

dev.off()

save(Bar_Iwat, file = "data/IWAT_GO_Barplot.Rda")

# CSV GO enrichment
write.csv(as.data.frame(go_results), file = "Results/IWAT/Top10_differentially_expressed_genes_GO.csv")

save(go_results, file = "data/IWAT_GO.Rda")

# Barplot and csv for GO groups

pdf("Results/IWAT/GO_BP_Group.pdf")
barplot(go_group_BP, drop = TRUE)
dev.off()

write.csv(as.data.frame(go_group_BP), file = "Results/IWAT/Group_BP.csv")

################################################################################
# RWAT
################################################################################

res <- results(dds, contrast = list(c("Diet_Standard_vs_45_High_Fat_Diet", 
                                      "tissueRWAT.DietStandard")),
               alpha = 0.05)

# Getting differentialy expressed genes

sum(res$padj < 0.05, na.rm=TRUE)
resOrdered <- res[order(res$padj),]

# PLotting top 10 differentially expressed genes

top10_Genes <- rownames(resOrdered)[1:10]

for (i in top10_Genes){
  pdf(file.path("Results/RWAT", paste0(i, ".pdf")))
  plotCounts(dds, gene=i, intgroup="Diet")
  dev.off()
}

print(paste("Gene with lowest pvalue:", top10_Genes[1]))

# Volcano plot
pdf("Results/RWAT/VolcanoPlot.pdf", width = 10, 8)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano plot of DESeq2 EWAT',
                subtitle = 'Differential expression analysis',
                caption = 'Log2 fold change vs. p-value', 
                drawConnectors = TRUE,
                widthConnectors = 1, boxedLabels = TRUE,
                selectLab = top10_Genes)

dev.off()

# Geting significant differentially expressed genes

res_filtered <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

# GO Analysis
go_results <- enrichGO(gene = rownames(res_filtered), OrgDb = org.Mm.eg.db,
                       ont = "BP", pvalueCutoff = 0.05, readable = TRUE,
                       universe = rownames(dds),
                       qvalueCutoff = 0.2, keyType = "ENSEMBL")

go_group_BP <- groupGO(gene = rownames(res_filtered), OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                       ont = "BP", readable = TRUE)

# Barplot enrichment analysis

pdf("Results/RWAT/EnrichmentGO.pdf", width = 10, height = 10)

Bar_RWAT <- barplot(go_results, 
        drop = TRUE, 
        showCategory = 40,
        font.size = 6)

Bar_RWAT

dev.off()

save(Bar_RWAT, file = "data/RWAT_GO_Barplot.Rda")

# CSV GO enrichment
write.csv(as.data.frame(go_results), file = "Results/RWAT/Top10_differentially_expressed_genes_GO.csv")

save(go_results, file = "data/RWAT_GO.Rda")

# Barplot and csv for GO groups

pdf("Results/RWAT/GO_BP_Group.pdf")
barplot(go_group_BP, drop = TRUE)
dev.off()

write.csv(as.data.frame(go_group_BP), file = "Results/RWAT/Group_BP.csv")



