# DESCRIPTION ----------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: March 25th, 2022

#Purpose: SEC EV COVID-19 Project differential expression

#Imports:
# deduplicated raw counts (countTable_dedup.tsv)
# sample data

#Exports: 
# DEG analysis

# import packages -------------------------------------------------------------------------------
library(tidyverse) #42; base for data/table manipulation
library(ggpubr) #pretty graphs
library(DESeq2)
library(EnhancedVolcano)
library(devtools)
library(fst)
source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

groupPalette <- getPalette(5, "complementary") #this is the color palette for your groupings
dev.off()

dir.create("analysis/byHospitalization/DESeq2/")

# volcano plots ----------------------
createVolcano <- function(df, plotTitle, pval, logfc) {
  EnhancedVolcano(df,
                  lab = df$gene_name,
                  #selectLab = c("ENSG00000100985.7"),
                  selectLab = c(""),
                  x = 'log2FoldChange',
                  y = 'padj',
                  drawConnectors = TRUE,
                  title = plotTitle,
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  legendLabels = c('NS', bquote(~Log[2]~ 'fold change'), 'Adjusted p-value',
                                   bquote('Adjusted p-value and ' ~Log[2]~ 'fold change')),
                  pCutoff = pval,
                  FCcutoff = logfc,
                  pointSize = 2.5,
                  labSize = 5,
                  legendPosition = "bottom") +
    #theme_bw(base_size = 18) +
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_blank())
}

#function for running DE ------
genDEResults <- function(contrastName, param, contrast1, contrast2, dds, logfc, pval, outDir){
  res <- results(dds, contrast = c(param, contrast1, contrast2))
  res <- res[order(res$padj),]
  summary(res)
  res.df <- rownames_to_column(as.data.frame(res), var = "gene_id")
  res.df <- right_join(genes.anno, res.df)
  write_csv(res.df, paste(outDir, contrastName,"_allGenes.csv", sep = ""))
  
  sig <- subset(res, baseMean > 10)
  sig <- subset(res, abs(log2FoldChange) > logfc)
  sig <- subset(sig, padj < pval)
  summary(sig)
  

  #output cutoff table
  sig.df <- rownames_to_column(as.data.frame(sig), var = "gene_id")
  sig.df <- right_join(genes.anno, sig.df)
  sig.df.withMeans <- left_join(sig.df, meanNormCounts)
  print(nrow(sig.df))
  write_csv(sig.df.withMeans, paste(outDir, contrastName,"_differentiallyExpressedGenes.csv", sep = ""))
  
  #output prettier volcano plot
  volPlot <- createVolcano(res.df, contrastName, pval, logfc)
  ggsave(paste(outDir, contrastName, "_enhancedVolcano.png", sep = ""), width = 10, height = 8)
  

}


#load data -------------------------------------------------------------------------------
sampleData <- read_tsv("data/sampleData.txt")
sampleData$group <- factor(paste(sampleData$GroupCode, sampleData$Visit, sep = "_"))
sampleData$Group2 <- NULL
sampleData$Group2[sampleData$GroupCode == "icu"] <- "hospitalized"
sampleData$Group2[sampleData$GroupCode == "non_icu"] <- "hospitalized"
sampleData$Group2[sampleData$GroupCode == "out_pt_pos"] <- "not_hospitalized"
sampleData$Group2[sampleData$GroupCode == "out_pt_neg"] <- "not_hospitalized"

countTable.anno <- read_tsv("data/countTable_dedup_anno.txt")

countData <- countTable.anno %>% select(-gene_type, -gene_name, -biotype_class)
countData <- countData %>% column_to_rownames(var = "gene_id")

countData.mRNA <- countTable.anno %>% filter(biotype_class == "protein_coding") %>%
  select(-gene_type, -gene_name, -biotype_class) %>%
  column_to_rownames(var = "gene_id")

countData.lncRNA <- countTable.anno %>% filter(biotype_class == "lncRNA") %>%
  select(-gene_type, -gene_name, -biotype_class) %>%
  column_to_rownames(var = "gene_id")

countsSub <- countData

#order to match
i <- match(sampleData$Sample, names(countsSub), nomatch = 0)
countsSub <- countsSub[,i]

identical(names(countsSub), sampleData$Sample) #double checking

#differential expression -------------

dds <- DESeqDataSetFromMatrix(countData = countsSub, colData = sampleData, design = ~1 + Group2)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

resultsNames(dds)

normCounts <- counts(dds, normalized = TRUE)

vsd <- vst(dds)
plotPCA(vsd, intgroup = c("Group2")) +
  scale_color_manual(values = groupPalette[c(5,2)]) +
  geom_text_repel(aes(label=name)) +
  theme_bw(base_size = 18)
ggsave(paste("analysis/byHospitalization/DESeq2/PCA_DE_byGroup.png", sep = ""), width = 8, height = 8)

#normalized counts per group -----------
hospGroup <- sampleData$Sample[sampleData$Group2 %in% c("hospitalized")]
nonHospGroup <- sampleData$Sample[sampleData$Group2 %in% c("not_hospitalized")]

hospCounts <- normCounts[,hospGroup]
nonHopsCounts <- normCounts[,nonHospGroup]


meanNormCounts <- data.frame("gene_id" = row.names(normCounts),
                             "AvgExpr_Hosp" = rowMeans(hospCounts),
                             "AvgExpr_nonHosp" = rowMeans(nonHopsCounts)
                             )




#output differential expression -----------
outdir <- "analysis/byHospitalization/DESeq2_padj0.05/"
dir.create(outdir)
genDEResults(contrastName = "not_hospitalized_vs_hospitalized", param = "Group2", contrast1 = "hospitalized", contrast2 = "not_hospitalized", dds = dds, logfc = log2(1.5), pval = 0.05, outDir = outdir)

outdir <- "analysis/DESeq2_padj0.1/"
dir.create(outdir)
genDEResults(contrastName = "not_hospitalized_vs_hospitalized", param = "Group2", contrast1 = "hospitalized", contrast2 = "not_hospitalized", dds = dds, logfc = log2(1.5), pval = 0.1, outDir = outdir)
