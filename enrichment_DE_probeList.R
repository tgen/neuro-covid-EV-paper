
library(tidyverse)

limma.DE <- read_csv("analysis/byHospitalization/limma_voom_qualityWeights_padj0.05/limma_hosp_vs_nonHosp_differentiallyExpressedGenes.csv")
deseq.DE <- read_csv("analysis/byHospitalization/DESeq2_padj0.05/DESeq2_not_hospitalized_vs_hospitalized_differentiallyExpressedGenes.csv")

tauSpec <- read_csv("/Volumes/kjensen/elizabeth/projects/TissueSpecificity/tau_tissue_specificity_gtex.csv") %>%
  mutate(ens_base = sapply(str_split(gene_id, "\\."), '[', 1))

tauSpec2 <- read_csv("/Volumes/kjensen/elizabeth/projects/TissueSpecificity/tau_tissue_specificity_gtex_noReproductive.csv") %>%
  mutate(ens_base = sapply(str_split(gene_id, "\\."), '[', 1))

DE.intersect <- inner_join(limma.DE, deseq.DE, by = c("gene_id", "gene_type", "gene_name", "biotype_class")) %>%
  mutate(ens_base = sapply(str_split(gene_id, "\\."), '[', 1)) %>%
  left_join(tauSpec2, by = c("ens_base"))
write_csv(DE.intersect, "analysis/byHospitalization/limma_DESeq2_intersection_not_hospitalized_vs_hospitalized_differentiallyExpressedGenes.csv")

table(DE.intersect$tau)

DESeq2.withTau <- deseq.DE %>%
  mutate(ens_base = sapply(str_split(gene_id, "\\."), '[', 1)) %>%
  left_join(tauSpec, by = c("ens_base"))
write_csv(DESeq2.withTau, "analysis/byHospitalization/DESeq2_not_hospitalized_vs_hospitalized_differentiallyExpressedGenes_withTau.csv")
