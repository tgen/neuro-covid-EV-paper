# DESCRIPTION ----------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: April 27th, 2022

#Purpose: TEMPLATE SCRIPT FOR COMPLEX HEATMAP

#Imports:
# count table, sample data, genes of interest

#Exports: 
# mapping summaries, diversity, PCA, etc.

# import packages -------------------------------------------------------------------------------

library(tidyverse) #42; basis for data manipulation
library(devtools) #load scripts from github
library(fst) #load .fst files
library(ComplexHeatmap) #heatmaps
library(dendsort) #sort dendrograms
library(DESeq2) #for count table normalization
library(circlize) #for heatmap colorRamp2
source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

groupPalette <- getPalette(5, "complementary") #this is the color palette for your groupings
dev.off()

#set heatmap colors
col_fun = colorRamp2(c(-2, 0, 2), c("#0571b0", "white", "#ca0020"))

#calculate Z score function  -----------------------------------------------------------------------
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

#import count data, sample data, and gene annotations -----------------------------------------------
countTable <- read_tsv("data/countTable.txt")
sampleData <- read_tsv("data/sampleData.txt")
genes.anno <- read_fst("~/Dropbox/repos/longRNA-ping/sourceFiles/GRCh38_GENCODE29LNCipedia_geneAnnotations.fst")

#SET YOUR GENES OF INTEREST HERE -----------------------------------------------------
myGenes <- c("ACAT1",
             "ACE2",
             "CTSL",
             "FURIN",
             "NRP1",
             "SCARB1",
             "SLC6A19",
             "TLR7",
             "TMPRSS2",
             "TMPRSS4")

#set colors for groups in sampleData -----------------------------------------------------
myColors <- groupPalette
names(myColors) <- c("OutPatient Covid(-) A",
                     "OutPatient Covid(+) A",
                     "Hosp.NON-ICU A",
                     "Hosp. ICU B",
                     "Hosp. ICU A")

#import into DESeq2 for vsd  -----------------------------------------------------------------------
countData <- column_to_rownames(countTable, var = "Geneid")

#order to match
countsSub <- countData[,names(countData) %in% c(sort(sampleData$Sample))]
i <- match(sampleData$Sample, names(countsSub), nomatch = 0)
countsSub <- countsSub[,i]
identical(names(countsSub), as.character(sampleData$Sample)) #double checking

#normalized counts and vst
dds <- DESeqDataSetFromMatrix(countData = countsSub, colData = sampleData, design = ~GroupCode)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

vsd <- varianceStabilizingTransformation(dds)


vst.myGenes.anno <- rownames_to_column(as.data.frame(assay(vsd)), var = "gene_id") %>% left_join(genes.anno) %>% filter(gene_name %in% myGenes) %>%
  mutate(id_name = paste(gene_id, gene_name, sep = " : "))
vst.myGenes <- vst.myGenes.anno %>%
  dplyr::select(-gene_id, -gene_type, -biotype_class, -id_name) %>%
  column_to_rownames(var = "gene_name")

vst.myGenes_norm <- t(apply(vst.myGenes, 1, cal_z_score))
vst.myGenes_norm[is.na(vst.myGenes_norm)] <- 0
vst.myGenes_norm <- vst.myGenes_norm[rowSums(vst.myGenes_norm) != 0,]

# my_hclust_sample <- hclust(dist(t(vst.myGenes)), method = "mcquitty")
# #my_hclust_gene <- hclust(dist(vst.myGenes), method = "mcquitty")
# mat_cluster_cols <- sort_hclust(my_hclust_sample)

#SET UP YOUR ANNOTATION HERE  -----------------------------------------------------------------------
ann <- data.frame(Sample = names(as.data.frame(vst.myGenes_norm))) %>%
  left_join(sampleData[,c("Sample", "Group")])
ha = HeatmapAnnotation(df = column_to_rownames(ann, var = "Sample"),
                       col = list(Group = myColors),
                       which = "col",
                       gp = gpar(col = "black"),
                       simple_anno_size = unit(0.5, "cm"),
                       annotation_legend_param = list(Group = list(title = "",
                                                                   #title_gp = gpar(fontsize = 16),
                                                                   labels_gp = gpar(fontsize = 16),
                                                                   legend_width = unit(5, "cm")
                       ))
)



#create dendograms and sort with dendsort  -----------------------------------------------------------------------
row_dend = dendsort(hclust(dist(vst.myGenes_norm)))
col_dend = dendsort(hclust(dist(t(vst.myGenes_norm))))

#create heatmap and draw
ht2 <- Heatmap(vst.myGenes_norm,
               col = col_fun,
               cluster_rows = row_dend,
               cluster_columns = col_dend,
               show_row_names = TRUE,
               show_column_names = FALSE,
               heatmap_legend_param = list(
                 #legend_direction = "horizontal", 
                 legend_width = unit(5, "cm"),
                 title = "Z-score",
                 title_gp = gpar(fontsize = 16),
                 labels_gp = gpar(fontsize = 16)),
               top_annotation = ha)

png("analysis/WP_SARSCOV2_AND_COVID19_PATHWAY_Heatmap.png", width = 600, height = 400)
draw(ht2, merge_legend = TRUE, adjust_annotation_extension = TRUE)
dev.off()

#heatmaps of DEG  -----------------------------------------------------------------------
ICUA_vs_nonICU <- read_csv("analysis/limma_voom_qualityWeights_padj0.05/ICUA_vs_nonICU_differentiallyExpressedGenes.csv")
ICUA_vs_OUTNEG <- read_csv("analysis/limma_voom_qualityWeights_padj0.05/ICUA_vs_OUTNEG_differentiallyExpressedGenes.csv")
ICUA_vs_OUTPOS <- read_csv("analysis/limma_voom_qualityWeights_padj0.05/ICUA_vs_OUTPOS_differentiallyExpressedGenes.csv")
nonICU_vs_OUTNEG <- read_csv("analysis/limma_voom_qualityWeights_padj0.05/nonICU_vs_OUTNEG_differentiallyExpressedGenes.csv")
hosp_nonHosp <- read_csv("analysis/byHospitalization/limma_voom_qualityWeights_padj0.05/limma_hosp_vs_nonHosp_differentiallyExpressedGenes.csv")

myInputList <- list(ICUA_vs_nonICU, ICUA_vs_OUTNEG, ICUA_vs_OUTPOS, nonICU_vs_OUTNEG,
                    hosp_nonHosp)
myCompList <- list("ICUA_vs_nonICU", "ICUA_vs_OUTNEG", "ICUA_vs_OUTPOS", "nonICU_vs_OUTNEG",
                   "hosp_nonHosp")


myHeatmap <- function(myGenes, myComp){
  vst.myGenes.anno <- rownames_to_column(as.data.frame(assay(vsd)), var = "gene_id") %>% left_join(genes.anno) %>% filter(gene_name %in% myGenes) %>%
    mutate(id_name = paste(gene_id, gene_name, sep = " : "))
  vst.myGenes <- vst.myGenes.anno %>%
    dplyr::select(-gene_id, -gene_type, -biotype_class, -gene_name) %>%
    column_to_rownames(var = "id_name")
  
  vst.myGenes_norm <- t(apply(vst.myGenes, 1, cal_z_score))
  vst.myGenes_norm[is.na(vst.myGenes_norm)] <- 0
  vst.myGenes_norm <- vst.myGenes_norm[rowSums(vst.myGenes_norm) != 0,]
  
  row_dend = dendsort(hclust(dist(vst.myGenes_norm)))
  col_dend = dendsort(hclust(dist(t(vst.myGenes_norm))))
  
  ht2 <- Heatmap(vst.myGenes_norm,
                 col = col_fun,
                 cluster_rows = row_dend,
                 cluster_columns = col_dend,
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 heatmap_legend_param = list(
                   #legend_direction = "horizontal", 
                   legend_width = unit(5, "cm"),
                   title = "Z-score",
                   title_gp = gpar(fontsize = 16),
                   labels_gp = gpar(fontsize = 16)),
                 top_annotation = ha)
  png(paste("analysis/limma_voom_qualityWeights_padj0.05/", myComp, "_heatmap.png", sep = ""), width = 600, height = 600)
  draw(ht2, merge_legend = TRUE, adjust_annotation_extension = TRUE)
  dev.off()
}

x <- 1
while(x <= length(myCompList)){
  myGenes <- myInputList[[x]]$gene_name
  myHeatmap(myGenes, myCompList[x])
  
  x <- x + 1
}
