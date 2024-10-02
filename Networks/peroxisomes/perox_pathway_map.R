library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

load("~/THP1_Interactome/GSE207487/THP1_DGE/out/DGEList_UNvsWT.Rdata")
load("~/THP1_Interactome/GSE207487/THP1_DGE/out/exactTest_UNvsWT.Rdata")

fatty_acid_cat_genes <- read.table("fatty_acid_catabolisim_autoAnn.txt", sep = '\t', header = T)
genes <- fatty_acid_cat_genes$Gene
i <- which(rownames(et$table) %in% genes) # pull out FA gene logFC
fa_expression <- et$table[i,]
fa_expression$entrezid <- bitr(rownames(fa_expression), fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")$ENTREZID
head(fa_expression)

# KEGG ------------------
# Enrichment
ids <- bitr(genes, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
kk <- enrichKEGG(gene         = ids$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
kk

# Visualization
significant_pathways <- kk@result[kk@result$p.adjust < 0.01, ]
enriched_pathways <- significant_pathways$ID

fa_expression_vector <- setNames(fa_expression$logFC, fa_expression$entrezid)
head(fa_expression_vector)

for (pathway in enriched_pathways) {
  tryCatch(
    {
    pathview(gene.data = fa_expression_vector,
            pathway.id = pathway,
            species = "hsa",
            limit = list(gene = 2, cpd = 1), low = c("red", "red"), high = c("green", "green"), kegg.native = T)
    }, 
    error = function(cond) {
      next
    }
  )
}


