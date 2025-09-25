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
library(edgeR) #expression normalization
library(limma) #more expression normalization, differential expression
library(ggpubr) #paper level visuals
library(devtools) #load scripts from github
library(fst) #import compressed fst txt files
library(EnhancedVolcano) #prettier volcano plots

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

groupPalette <- getPalette(5, "complementary") #this is the color palette for your groupings
dev.off()

dir.create("analysis/byHospitalization/")

# volcano plots ----------------------
createVolcano <- function(df, plotTitle, pval, logfc) {
  EnhancedVolcano(df,
                  lab = df$gene_name,
                  #selectLab = c("ENSG00000100985.7"),
                  selectLab = c(""),
                  x = 'logFC',
                  y = 'adj.P.Val',
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

#resultsOut function ---------------------------------------------------------
resultsOut <- function(limmaFit, contr.matrix, outDir, logfc, pval){
  
  x <- 1
  while(x <= length(colnames(contr.matrix))){
    contrName <- colnames(contr.matrix)[x]
    print(contrName)
    
    dt <- decideTests(limmaFit)
    
    #output full table
    df <- rownames_to_column(topTable(efit, coef = x, n = Inf), var = "gene_id")
    df <- right_join(genes.anno, df)
    df.withMeans <- left_join(df, meanCPM)
    write_csv(df.withMeans, paste(outDir, contrName,"_allGenes.csv", sep = ""))
    
    #output volcano plot
    png(paste(outDir, contrName, "_volcano.png", sep = ""))
    with(df, plot(logFC, -log10(P.Value), pch = 20, main = paste("Volcano Plot: ", contrName, sep = "")))
    with(subset(df, adj.P.Val < pval & abs(logFC) > logfc), points(logFC, -log10(P.Value), pch=20, col="red"))
    dev.off()
    
    #output prettier volcano plot
    volPlot <- createVolcano(df, contrName, pval, logfc)
    ggsave(paste(outDir, contrName, "_enhancedVolcano.png", sep = ""), width = 10, height = 8)
    
    #output MD plot
    png(paste(outDir, contrName, "_MDPlot.png", sep = ""))
    plotMD(limmaFit, column = x, status = dt[,x], main = contrName,
           xlim=c(-10,15))
    dev.off()
    
    #output cutoff table
    df.cut <- rownames_to_column(topTable(efit, coef = x, n = Inf, lfc = logfc, p.value = pval), var = "gene_id")
    df.cut <- right_join(genes.anno, df.cut)
    df.cut.withMeans <- left_join(df.cut, meanCPM)
    dim(df.cut.withMeans)
    print(nrow(df.cut))
    write_csv(df.cut.withMeans, paste(outDir, contrName,"_differentiallyExpressedGenes.csv", sep = ""))
    
    
    x <- x + 1
  }
  
}


#load data -------------------------------------------------------------------------------
genes.anno <- read_fst("~/Dropbox/repos/longRNA-ping/sourceFiles/GRCh38_GENCODE29LNCipedia_geneAnnotations.fst")

sampleData <- read_tsv("data/sampleData.txt")
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

## sample design --------------------------------------------------------------------------
participant <- factor(as.character(sampleData$Participant))
group <- factor(sampleData$Group2)

design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))

contr.matrix <- makeContrasts(
  hosp_vs_nonHosp = hospitalized - not_hospitalized,
  levels = colnames(design)
)

## normalization and filtering -------------------------------------------------------------------

keep <- filterByExpr(countsSub, group = group)
sum(keep) #10393 genes

dge <- DGEList(countsSub[keep,])
dge  <- calcNormFactors(dge)

logCPM <- cpm(dge, log = TRUE, prior.count = 3)

## MDS plot ----------------------------------------------------------------------------------

col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

plotMDS(logCPM, labels=group, col=col.group)

## limma/voom ----------------------------------------------------------------------------------
v <- voomWithQualityWeights(dge, design, plot=TRUE) #check voom plot for curve to see if we need to do more filtering
#v <- voom(dge, design, plot=TRUE) #check voom plot for curve to see if we need to do more filtering
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)
summary(decideTests(efit))

#make table of mean expression for each group
hospGroup <- sampleData$Sample[sampleData$Group2 %in% c("hospitalized")]
nonHospGroup <- sampleData$Sample[sampleData$Group2 %in% c("not_hospitalized")]

hospCPM <- logCPM[,hospGroup]
nonHospCPM <- logCPM[,nonHospGroup]


meanCPM <- data.frame("gene_id" = row.names(logCPM),
                      "AvgExpr_Hosp" = rowMeans(hospCPM),
                      "AvgExpr_nonHosp" = rowMeans(nonHospCPM)
                      )


#output results
outdir <- "analysis/byHospitalization/limma_voom_qualityWeights_padj0.1/"
dir.create(outdir)
resultsOut(efit, contr.matrix, outdir, logfc = log2(1.5), pval = 0.1)

outdir <- "analysis/byHospitalization/limma_voom_qualityWeights_padj0.05/"
dir.create(outdir)
resultsOut(efit, contr.matrix, outdir, logfc = log2(1.5), pval = 0.05)
