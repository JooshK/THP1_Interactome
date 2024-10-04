library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(tidyverse)

load("../../GSE207487/THP1_DGE/out/DGEList_UNvsWT.Rdata")
load("../../GSE207487/THP1_DGE/out/exactTest_UNvsWT.Rdata")

perox_genes <- read.csv("perox_lipid_genes.csv") |> 
  column_to_rownames(var = "ENSEMBL")
genes <- rownames(perox_genes)
i <- which(rownames(et$table) %in% genes) # pull out gene logFC
perox_expression <- et$table[i,]

perox_expression <- merge(perox_expression, perox_genes, by=0, all.x = T) |> 
  column_to_rownames(var = "Row.names")

head(perox_expression)

# KEGG ------------------
# Enrichment
kk <- enrichKEGG(gene         = perox_expression$ENTREZ_ACC,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
View(kk@result)

# Visualization
significant_pathways <- kk@result[kk@result$p.adjust < 0.01, ]
enriched_pathways <- significant_pathways$ID

perox_expression_vector <- setNames(perox_expression$logFC, perox_expression$ENTREZ_ACC)
head(perox_expression_vector)

pathview(gene.data = perox_expression_vector,
         pathway.id = "hsa04146",
         species = "hsa",
         limit = list(gene = 2, cpd = 1), low = c("red", "red"), high = c("green", "green"), kegg.native = T)

