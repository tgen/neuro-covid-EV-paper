# DESCRIPTION ----------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: June 3rd, 2022

#Purpose: cluster profiler pathway analysis

#Imports:
# DE tables

#Exports: 
# 

# import packages -------------------------------------------------------------------------------
library(tidyverse) #42; basis for data manipulation
library(ggpubr) #pretty graphs
library(gprofiler2) #pathway analysis tool that the proteomics group uses
library(devtools) #load R scripts from github
library(clusterProfiler) #supports over-representation tests and GSEA of GO
library(DOSE) #suppose disease ontology semantic and enrichment analysis
library(enrichplot)
library(org.Hs.eg.db) #load homo sapiens database
library(ReactomePA) #reactome pathway analysis


dir.create("analysis/clusterProfiler")

#function to output results and generate plots  -------------------------------------
genResults <- function(myGeneList, enrichObject, analysisName, compName){
  
  #create output directory for this comparison if it doesn't exist
  outDir <- paste(outdir, "/", compName, "/", sep = "")
  
  if(!dir.exists(outDir)){
    dir.create(outDir)
  }
  
  write_tsv(as.data.frame(enrichObject), paste(outDir, compName, "_", analysisName, ".tsv", sep = ""))
  
  dotplot(enrichObject, showCategory = 20)
  ggsave(paste(outDir, compName, "_", analysisName, "_dotplot.png", sep = ""), width = 8, height = 12)
  
  cnetplot(enrichObject, categorySize="pvalue", foldChange=myGeneList) + ggtitle(compName)
  ggsave(paste(outDir, compName, "_", analysisName, "_cnetplot.png", sep = ""), width = 12, height = 12)
  
  # heatplot(enrichObject, foldChange=myGeneList) + ggtitle(compName)
  # ggsave(paste(outDir, compName, "_", analysisName, "_heatplot.png", sep = ""), width = 20, height = 6)
  
  upsetplot(enrichObject) + ggtitle(compName)
  ggsave(paste(outDir, compName, "_", analysisName, "_upsetplot.png", sep = ""), width = 12, height = 8)
  
  x2 <- pairwise_termsim(enrichObject)
  emapplot(x2, pie_scale=1.5, pie = "count") + ggtitle(compName)
  ggsave(paste(outDir, compName, "_", analysisName, "_emapplot.png", sep = ""), width = 12, height = 12)
  
  
  treeplot(x2) + ggtitle(compName) #+ theme(legend.position = "bottom",
  #plot.margin = margin(2,10,0,0, "cm"))
  ggsave(paste(outDir, compName, "_", analysisName, "_treeplot.png", sep = ""), width = 12, height = 8)
  
}

#function to run cluster profiler on GO and KEGG enrichment  -------------------------------------
runClusterProfiler <- function(myGeneList, backgroundGenes, myGeneList.entrez, backgroundGenes.entrez, compName, pcutoff){
  
  # GO, KEGG, and DO: sig categories ----------------------------------------------------
  #run enrichGO
  ego.BP <- enrichGO(gene          = names(myGeneList),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL',
                     universe      = backgroundGenes,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = pcutoff,
                     qvalueCutoff  = pcutoff)
  
  
  if(nrow(ego.BP) != 0){
    head(as.data.frame(ego.BP))
    genResults(myGeneList, ego.BP, "GOBP_enrich_sig", compName)
  }
  
  ego.MF <- enrichGO(gene          = names(myGeneList),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL',
                     universe      = backgroundGenes,
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = pcutoff,
                     qvalueCutoff  = pcutoff)
  
  
  if(nrow(ego.MF) != 0){
    head(as.data.frame(ego.MF))
    genResults(myGeneList, ego.MF, "GOMF_enrich_sig", compName)
  }
  
  ego.CC <- enrichGO(gene          = names(myGeneList),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL',
                     universe      = backgroundGenes,
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = pcutoff,
                     qvalueCutoff  = pcutoff)
  
  
  if(nrow(ego.CC) != 0){
    head(as.data.frame(ego.CC))
    genResults(myGeneList, ego.CC, "GOCC_enrich_sig", compName)
  }
  
  
  enDO <- enrichDO(gene          = names(myGeneList.entrez),
                   ont           = "DO",
                   pvalueCutoff  = pcutoff,
                   pAdjustMethod = "BH",
                   universe      = backgroundGenes.entrez$ENTREZID,
                   minGSSize     = 5,
                   maxGSSize     = 500,
                   qvalueCutoff  = pcutoff)
  if(nrow(enDO) != 0){
    head(as.data.frame(enDO))
    genResults(myGeneList, enDO, "DO_enrich_sig", compName)
  }
  
  
  kk <- enrichKEGG(gene         = names(myGeneList.entrez),
                   organism     = 'hsa',
                   universe     = backgroundGenes.entrez$ENTREZID,
                   pvalueCutoff = pcutoff)
  if(nrow(kk) != 0){
    head(as.data.frame(kk))
    genResults(myGeneList, kk, "KEGG_enrich_sig", compName)
  }
  
  reac <- enrichPathway(gene         = names(myGeneList.entrez),
                        organism     = 'human',
                        universe     = backgroundGenes.entrez$ENTREZID,
                        pvalueCutoff = pcutoff)
  if(nrow(reac) != 0){
    head(as.data.frame(reac))
    genResults(myGeneList, reac, "REACTOME_enrich", compName)
    
  }
  
}


#clusterProfiler ----------
outdir <- "analysis/clusterProfiler"

#IMPORT GENE LIST OF INTEREST AND SET HERE, IE
#ICUA_vs_ICUB <- read_csv("analysis/limma_voom_qualityWeights_padj0.1/ICUA_vs_ICUB_differentiallyExpressedGenes.csv")
ICUA_vs_nonICU <- read_csv("analysis/limma_voom_qualityWeights_padj0.1/ICUA_vs_nonICU_differentiallyExpressedGenes.csv")
ICUB_vs_nonICU <- read_csv("analysis/limma_voom_qualityWeights_padj0.1/ICUB_vs_nonICU_differentiallyExpressedGenes.csv")
ICUA_vs_OUTNEG <- read_csv("analysis/limma_voom_qualityWeights_padj0.1/ICUA_vs_OUTNEG_differentiallyExpressedGenes.csv")
ICUA_vs_OUTPOS <- read_csv("analysis/limma_voom_qualityWeights_padj0.1/ICUA_vs_OUTPOS_differentiallyExpressedGenes.csv")
#OUTPOS_vs_OUTNEG <- read_csv("analysis/limma_voom_qualityWeights_padj0.1/OUTPOS_vs_OUTNEG_differentiallyExpressedGenes.csv")
nonICU_vs_OUTNEG <- read_csv("analysis/limma_voom_qualityWeights_padj0.1/nonICU_vs_OUTNEG_differentiallyExpressedGenes.csv")
#nonICU_vs_OUTPOS <- read_csv("analysis/limma_voom_qualityWeights_padj0.1/nonICU_vs_OUTPOS_differentiallyExpressedGenes.csv")

hosp_nonHosp <- read_csv("analysis/byHospitalization/limma_voom_qualityWeights_padj0.1/hosp_vs_nonHosp_differentiallyExpressedGenes.csv")

backgroundGenesA <- read_csv("analysis/limma_voom_qualityWeights_padj0.1/ICUA_vs_ICUB_allGenes.csv")
backgroundGenesB <- read_csv("analysis/byHospitalization/limma_voom_qualityWeights_padj0.1/hosp_vs_nonHosp_allGenes.csv")

#if DESeq:
#myGeneList <- myGenes$log2FoldChange

myInputList <- list(ICUA_vs_ICUB, ICUA_vs_nonICU, ICUA_vs_OUTNEG, ICUA_vs_OUTPOS, OUTPOS_vs_OUTNEG, nonICU_vs_OUTNEG, nonICU_vs_OUTPOS,
                    hosp_nonHosp)
myCompList <- list("ICUA_vs_nonICU", "ICUA_vs_OUTNEG", "ICUA_vs_OUTPOS", "nonICU_vs_OUTNEG",
                   "hosp_nonHosp")
myBackgroundList <- list(backgroundGenesA, backgroundGenesA, backgroundGenesA, backgroundGenesA,
                         backgroundGenesB)


myInputList <- list(ICUB_vs_nonICU)
myCompList <- list("ICUB_vs_nonICU")
myBackgroundList <- list(backgroundGenesA)

x <- 1
while(x <= length(myInputList)){
  myGenes <- myInputList[x][[1]]
  myCompName <- myCompList[x][[1]]
  myBackground <- myBackgroundList[x][[1]]
  
  myGeneList <- myGenes$logFC
  
  names(myGeneList) <- myGenes$gene_name
  
  backgroundGenes <- myBackground
  
  #convert to entrez ids
  backgroundGenes.entrez <- bitr(backgroundGenes$gene_name, fromType = "SYMBOL", toType = "ENTREZID", org.Hs.eg.db)
  length(backgroundGenes.entrez$ENTREZID)
  
  myGeneList.entrez.tmp <- bitr(names(myGeneList), fromType = "SYMBOL", toType = "ENTREZID", org.Hs.eg.db)
  myGeneList.entrez.tmp <- left_join(myGeneList.entrez.tmp, data.frame(SYMBOL = names(myGeneList), FC = myGeneList))
  myGeneList.entrez <- myGeneList.entrez.tmp$FC
  names(myGeneList.entrez) <- myGeneList.entrez.tmp$ENTREZID
  myGeneList.entrez <- sort(myGeneList.entrez, decreasing = TRUE)
  
  
  runClusterProfiler(myGeneList, as.character(backgroundGenes$gene_name), myGeneList.entrez, backgroundGenes.entrez, myCompName, pcutoff = 0.05)
  
  x <- x + 1
}



#look at a few comparisons by hand ---------
hosp_nonHosp.up <- hosp_nonHosp %>% filter(logFC > 0)
hosp_nonHosp.down <- hosp_nonHosp %>% filter(logFC < 0)

#hosp vs non hosp upregulated ---------
myGenes <- hosp_nonHosp.up
myCompName <- "hosp_nonHosp_upregulated"
myBackground <- backgroundGenesB

myGeneList <- myGenes$logFC

names(myGeneList) <- myGenes$gene_name

backgroundGenes <- myBackground

#convert to entrez ids
backgroundGenes.entrez <- bitr(backgroundGenes$gene_name, fromType = "SYMBOL", toType = "ENTREZID", org.Hs.eg.db)
length(backgroundGenes.entrez$ENTREZID)

myGeneList.entrez.tmp <- bitr(names(myGeneList), fromType = "SYMBOL", toType = "ENTREZID", org.Hs.eg.db)
myGeneList.entrez.tmp <- left_join(myGeneList.entrez.tmp, data.frame(SYMBOL = names(myGeneList), FC = myGeneList))
myGeneList.entrez <- myGeneList.entrez.tmp$FC
names(myGeneList.entrez) <- myGeneList.entrez.tmp$ENTREZID
myGeneList.entrez <- sort(myGeneList.entrez, decreasing = TRUE)


runClusterProfiler(myGeneList, as.character(backgroundGenes$gene_name), myGeneList.entrez, backgroundGenes.entrez, myCompName, pcutoff = 0.99)



#hosp vs non hosp downregulated  ---------
myGenes <- hosp_nonHosp.down
myCompName <- "hosp_nonHosp_downregulated"
myBackground <- backgroundGenesB

myGeneList <- myGenes$logFC

names(myGeneList) <- myGenes$gene_name

backgroundGenes <- myBackground

#convert to entrez ids
backgroundGenes.entrez <- bitr(backgroundGenes$gene_name, fromType = "SYMBOL", toType = "ENTREZID", org.Hs.eg.db)
length(backgroundGenes.entrez$ENTREZID)

myGeneList.entrez.tmp <- bitr(names(myGeneList), fromType = "SYMBOL", toType = "ENTREZID", org.Hs.eg.db)
myGeneList.entrez.tmp <- left_join(myGeneList.entrez.tmp, data.frame(SYMBOL = names(myGeneList), FC = myGeneList))
myGeneList.entrez <- myGeneList.entrez.tmp$FC
names(myGeneList.entrez) <- myGeneList.entrez.tmp$ENTREZID
myGeneList.entrez <- sort(myGeneList.entrez, decreasing = TRUE)


runClusterProfiler(myGeneList, as.character(backgroundGenes$gene_name), myGeneList.entrez, backgroundGenes.entrez, myCompName, pcutoff = 0.99)


#generate plots for hosp vs non hosp ---------
myGenes <- hosp_nonHosp
myCompName <- "hosp_nonHosp"
myBackground <- backgroundGenesB

myGeneList <- myGenes$logFC

names(myGeneList) <- myGenes$gene_name

backgroundGenes <- myBackground

#convert to entrez ids
backgroundGenes.entrez <- bitr(backgroundGenes$gene_name, fromType = "SYMBOL", toType = "ENTREZID", org.Hs.eg.db)
length(backgroundGenes.entrez$ENTREZID)

myGeneList.entrez.tmp <- bitr(names(myGeneList), fromType = "SYMBOL", toType = "ENTREZID", org.Hs.eg.db)
myGeneList.entrez.tmp <- left_join(myGeneList.entrez.tmp, data.frame(SYMBOL = names(myGeneList), FC = myGeneList))
myGeneList.entrez <- myGeneList.entrez.tmp$FC
names(myGeneList.entrez) <- myGeneList.entrez.tmp$ENTREZID
myGeneList.entrez <- sort(myGeneList.entrez, decreasing = TRUE)
