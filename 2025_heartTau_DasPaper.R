# import packages -------------------------------------------------------------------------------
library(tidyverse) #42; basis for data manipulation
library(ggpubr) #pretty graphs
library(DESeq2) #data normalization, pca
library(ggrepel) #connectors to data points
library(devtools) #load scripts from github
library(fst) #import compressed fst txt files
library(UpSetR) #upset plots

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes
source_gist("ab773041ac27834c767f5eff1ddf6dc8") #biotype plotter

groupPalette <- getPalette(5, "complementary")[c(5,3)] #this is the color palette for your groupings
dev.off()

countTable <- read_tsv("data/countTable.txt") #counts (gene-level)

genes.anno <- read_fst("/data/NGD/resources/longRNA-ping/ref/grch38_hg38/gencode40_lncipedia5.2/GRCh38_GENCODE40LNCipedia_geneAnnotations.fst") %>%
  dplyr::select(gene_id, gene_type, gene_name, biotype_class) %>%
  unique()
genes.anno.ens <- read_fst("/data/NGD/resources/longRNA-ping/ref/grch38_hg38/gencode40_lncipedia5.2/GRCh38_GENCODE40LNCipedia_geneAnnotations.fst") %>%
  dplyr::select(gene_id, gene_type, gene_name, biotype_class, ensembl_id) %>%
  mutate(ens_base = sapply(str_split(ensembl_id, "\\."), '[', 1)) %>%
  select(-ensembl_id) %>%
  unique()

names(countTable) <- gsub("_CMB", "", names(countTable))

batch3IDs <- c(8543, 8559, 8642, 8999)
exclude <- c("33123_2_PD_Plasma_combined_2_down75M", "33111_2_PD_Plasma_combined_2_down75M",
             "33151_1_Control_Urine_down75M", "33151_2_Control_Urine_down75M",
             "33112_1_PD_Urine_down75M", "33112_2_PD_Urine_down75M",
             "33142_1_AD_Urine_down75M", "33142_2_AD_Urine_down75M",
             "33124_1_PD_Urine_down75M", "33124_2_PD_Urine_down75M"
)

samples <- names(countTable)[-1]
sampleData <- data.frame(Sample = samples,
                         Participant = sapply(str_split(samples, "_"), '[', 1),
                         Visit = sapply(str_split(samples, "_"), '[', 2),
                         Diag = sapply(str_split(samples, "_"), '[', 3),
                         Biofluid = sapply(str_split(samples, "_"), '[', 4),
                         Group = "Participant",
                         Type = "EVs"
                         
) %>%
  filter(!grepl("70EtOH", Sample)) %>%
  filter(!Sample %in% exclude) %>%
  filter(!Participant %in% batch3IDs)
sampleData[sampleData$Participant == "33131" & sampleData$Visit == 1,]$Visit <- 0
sampleData[sampleData$Participant == "33131" & sampleData$Visit == 2,]$Visit <- 1
sampleData[sampleData$Participant == "33131" & sampleData$Visit == 3,]$Visit <- 2

sampleData <- sampleData %>% filter(!Visit == 0)

#sampleData$Sample <- gsub("_100EtOH", "_100EtOH", sampleData$Sample)


sampleData$Group <- factor(sampleData$Group, levels = c("Participant", "Pool"))
sampleData$Group[grep("PC", sampleData$Participant)] <- "Pool"
sampleData$Biofluid[grep("PC", sampleData$Participant)] <- "Plasma"
sampleData$Group[grep("UC", sampleData$Participant)] <- "Pool"
sampleData$Biofluid[grep("UC", sampleData$Participant)] <- "Urine"
sampleData$Visit[grep("Pool", sampleData$Group)] <- "Pool"
sampleData$Diag[grep("Pool", sampleData$Group)] <- "Pool"

demo <- read_csv("data/demographics.csv")
demo$Participant <- as.character(demo$Participant)
head(demo)

sampleData <- left_join(sampleData, demo)

sampleData_forPlot <- sampleData %>%
  filter(!Group == "Pool") %>%
  filter(Biofluid == "Plasma") %>%
  filter(Visit == 1)
write_csv(sampleData_forPlot %>% dplyr::select(-Sample, -Visit), "data/sampleData_forPlot.csv")
table(sampleData_forPlot$Diag, sampleData_forPlot$sex)
summary(sampleData_forPlot$age_at_exam)
summary(sampleData_forPlot$days_between_visits)

#DESeq2 normalized counts  -------------------------------------------------------------------------------
#order to match
countData <- countTable %>% column_to_rownames(var = "Geneid")

countsSub <- countData[,names(countData) %in% c(sort(sampleData$Sample))]
i <- match(sampleData$Sample, names(countsSub), nomatch = 0)
countsSub <- countsSub[,i]
identical(names(countsSub), as.character(sampleData$Sample)) #double checking

#normalized counts and vst
dds <- DESeqDataSetFromMatrix(countData = countsSub, colData = sampleData, design = ~Biofluid + Diag)
dds <- estimateSizeFactors(dds)

#vsd <- vst(dds)
normCounts <- counts(dds, normalized = TRUE)
# filter <- rowSums(normCounts >= 20) >= 10
# normCounts.filt <- normCounts[filter,]
normCounts.anno <- right_join(genes.anno, rownames_to_column(as.data.frame(normCounts), var = "gene_id"))
normCounts.anno.ens <- right_join(genes.anno.ens, rownames_to_column(as.data.frame(normCounts), var = "gene_id"))

nc.longer <- normCounts %>% as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "Sample", values_to = "nc") %>%
  left_join(sampleData)


#import tau and subset to heart specific genes --------
x <- 0.6
tau.spec <- read_csv("/data/NGD/elizabeth/from_tgen/projects/TissueSpecificity/tau_tissue_specificity_gtex.csv")
i <- grep(paste("Heart", "$", sep = ""), names(tau.spec))
myGenes <- dplyr::filter(tau.spec, tau.spec[,i] >= x)
myGenes <- myGenes[,c("gene_name", "gene_id", "gene_biotype")] %>%
  mutate(ens_base = sapply(str_split(gene_id, "\\."), '[', 1))

tau.spec.ens <- tau.spec %>%
  mutate(ens_base = sapply(str_split(gene_id, "\\."), '[', 1))
names(tau.spec.ens) <- gsub("Heart", "heart_tau", names(tau.spec.ens))
#scatterplots  -------------------------------------------------------------------------------
mean.nc.long <- nc.longer %>%
  group_by(gene_id, Biofluid, Visit) %>%
  summarise(mean = mean(nc)) %>%
  replace(is.na(.), 0) %>%
  left_join(genes.anno.ens[,c("gene_id", "ens_base")])

mean.nc.long.gr10 <- nc.longer %>%
  filter(nc > 10) %>%
  group_by(gene_id, Biofluid, Visit) %>%
  summarise(mean = mean(nc)) %>%
  replace(is.na(.), 0) %>%
  left_join(genes.anno.ens[,c("gene_id", "ens_base")])

sp.input <- mean.nc.long.gr10 %>% filter(ens_base %in% myGenes$ens_base) %>%
  filter(Biofluid == "Plasma")

sp <- ggscatter(sp.input %>%
                  pivot_wider(names_from = "Visit", values_from = "mean", names_prefix = "Visit_"),
                x = "Visit_1", y = "Visit_2", size = 3, shape = 21,
                alpha = 0.5,
                #color = "Biofluid",
                fill = "Biofluid",
                add = "reg.line", conf.int = TRUE,
                palette = groupPalette) +
  #expand_limits(x = 1, y = 1) +
  stat_cor(aes(color = Biofluid, label = ..rr.label..), label.y = 4, method = "spearman", size = 16) + 
  facet_grid(. ~ Biofluid) +
  theme_pubr(base_size = 22)

ggpar(sp,
      title = "Heart specific genes (tau > 0.6)",
      legend = "none",
      xlab = "Visit 1 normalized counts \n (log10 scale)",
      ylab = "Visit 2 normalized counts \n (log10 scale)",
      xscale = c("log10"),
      yscale = c("log10")
) 
ggsave("analysis/heartTauGenes_scatterplot_visits_basePlot.png", width = 16, height = 8)


#fancier
sp2 <- ggscatter(sp.input %>%
                  pivot_wider(names_from = "Visit", values_from = "mean", names_prefix = "Visit_") %>%
                  left_join(tau.spec.ens[,c("ens_base", "heart_tau")]),
                x = "Visit_1", y = "Visit_2", size = 3, shape = 21,
                alpha = 0.7,
                color = "heart_tau",
                fill = "heart_tau",
                add = "reg.line", conf.int = TRUE,
                palette = groupPalette) +
  #expand_limits(x = 1, y = 1) +
  #aes(color = Heart) + 
  stat_cor(aes(label = ..rr.label..), label.y = 3, method = "spearman", size = 8) + 
  scale_color_viridis(option = "rocket", direction = -1) +
  scale_fill_viridis(option = "rocket", direction = -1) +
  #facet_grid(. ~ Biofluid) +
  scale_x_continuous(breaks = c(10, 100, 1000, 1000)) +
  scale_x_log10() +
  scale_y_log10() +
  #xlim(10, 3000) +
  #ylim(10, 3000) +
  xlab("Visit 1 normalized counts \n (log10 scale)") +
  ylab("Visit 2 normalized counts \n (log10 scale)") +
  ggtitle("Heart specific genes (tau > 0.6)") +
  theme_pubr(base_size = 22) +
  theme(legend.position = "right")
sp2

# ggpar(sp2,
#       title = "Heart specific genes (tau > 0.6)",
#       legend = "none",
#       xlab = "Visit 1 normalized counts \n (log10 scale)",
#       ylab = "Visit 2 normalized counts \n (log10 scale)",
#       xscale = c("log10"),
#       yscale = c("log10")
# ) 
ggsave("analysis/heartTauGenes_scatterplot_visits_tau0.6.png", width = 16, height = 8)


sp3 <- ggscatter(sp.input %>%
                   pivot_wider(names_from = "Visit", values_from = "mean", names_prefix = "Visit_") %>%
                   left_join(tau.spec.ens[,c("ens_base", "heart_tau")]) %>%
                   filter(heart_tau > 0.9),
                 x = "Visit_1", y = "Visit_2", size = 3, shape = 21,
                 alpha = 0.7,
                 color = "heart_tau",
                 fill = "heart_tau",
                 add = "reg.line", conf.int = TRUE,
                 palette = groupPalette) +
  #expand_limits(x = 1, y = 1) +
  #aes(color = Heart) + 
  stat_cor(aes(label = ..rr.label..), label.y = 2.5, method = "spearman", size = 8) + 
  scale_color_viridis(option = "rocket", direction = -1) +
  scale_fill_viridis(option = "rocket", direction = -1) +
  #facet_grid(. ~ Biofluid) +
  scale_x_continuous(breaks = c(10, 100, 1000, 1000)) +
  scale_x_log10() +
  scale_y_log10() +
  #xlim(10, 3000) +
  #ylim(10, 3000) +
  xlab("Visit 1 normalized counts \n (log10 scale)") +
  ylab("Visit 2 normalized counts \n (log10 scale)") +
  ggtitle("Heart specific genes (tau > 0.9)") +
  theme_pubr(base_size = 22) +
  theme(legend.position = "right")
sp3

# ggpar(sp2,
#       title = "Heart specific genes (tau > 0.)",
#       legend = "none",
#       xlab = "Visit 1 normalized counts \n (log10 scale)",
#       ylab = "Visit 2 normalized counts \n (log10 scale)",
#       xscale = c("log10"),
#       yscale = c("log10")
# ) 
ggsave("analysis/heartTauGenes_scatterplot_visits_tau0.9.png", width = 16, height = 8)
