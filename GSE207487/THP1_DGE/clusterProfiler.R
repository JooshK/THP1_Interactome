library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)

# clusterProfiler based gene annotations, the numbers correspond to part/chapter in 
# the documentation

# GO Annotations - ORA --------------------
geneList <- read.csv("out/UN_vs_WT_DGE.csv")
genes_entrez <- bitr(geneList$ENSEMBL, fromType="ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
genes_entrez <- merge(geneList, genes_entrez, by.x=1, by.y="ENSEMBL")

sigGenes <- genes_entrez[genes_entrez$FDR < 0.05,]$ENTREZID

# GO enrichment analysis 
ego <- enrichGO(gene          = sigGenes,
                universe      = genes_entrez$ENTREZID, # specify the universe, different than g:Profiler
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)s

head(ego)

# Gene set testing 
gsea_rnks <- read.table("out/UNvsWT.rnk", sep = "\t", header = T)
rnks_entrez <- bitr(gsea_rnks$GeneName, fromType="SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
gsea_rnks <- merge(gsea_rnks, rnks_entrez, by.x=1, by.y = "SYMBOL")

ranks <- gsea_rnks[,2]
names(ranks) <- as.character(gsea_rnks[,3])
ranks <- sort(ranks, decreasing = T)
head(ranks)

ego3 <- gseGO(geneList     = ranks,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 1000,
              pvalueCutoff = 0.05,
              verbose      = FALSE)


dotplot(ego, showCategory=30)
ridgeplot(ego3)



