library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

load("~/THP1_Interactome/GSE207487/THP1_DGE/out/DGEList_UNvsWT.Rdata")
load("~/THP1_Interactome/GSE207487/THP1_DGE/out/exactTest_UNvsWT.Rdata")

perox_genes <- read.table("perox_genes.txt", sep = '\t', header = T)
genes <- perox_genes$GeneID
i <- which(rownames(et$table) %in% genes) # pull out FA gene logFC
perox_expression <- et$table[i,]

perox_expression$entrezid <- bitr(rownames(perox_expression), fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")$ENTREZID
head(perox_expression)

# KEGG ------------------
# Enrichment
ids <- bitr(genes, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
kk <- enrichKEGG(gene         = ids$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
View(kk@result)

# Visualization
significant_pathways <- kk@result[kk@result$p.adjust < 0.01, ]
enriched_pathways <- significant_pathways$ID

perox_expression_vector <- setNames(perox_expression$logFC, perox_expression$entrezid)
head(perox_expression_vector)

pathview(gene.data = perox_expression_vector,
         pathway.id = "hsa03320",
         species = "hsa",
         limit = list(gene = 2, cpd = 1), low = c("red", "red"), high = c("green", "green"), kegg.native = T)


