library(pheatmap)
library(RColorBrewer)
library(pathview)

load("~/Documents/THP1_Interactome/GSE207487/THP1_DGE/out/DGEList_UNvsWT.Rdata")
load("~/Documents/THP1_Interactome/GSE207487/THP1_DGE/out/exactTest_UNvsWT.Rdata")

genes_list <- read.table("fatty_acid_catabolisim_autoAnn.txt", sep = '\t', header = T)

lcpm_norm <- as.data.frame(cpm(x,log=TRUE))
i <- which(rownames(et$genes) %in% genes_list$Gene) # pull out the genes of interest
group <- as.factor(c(rep("UN", 3), rep("WT", 3)))

pheatmap(lcpm_norm[i,], 
         cluster_rows = T, cluster_cols = T, trace = "none",
         angle_col = 45, 
         fontsize_row = 6, fontsize_col = 8, show_rownames = T, 
         main = "Expression of Peroxisome Genes", labels_row = genes_list$Description, 
         labels_col = group)


perox_genes <- read.table("perox_genes.txt", sep = '\t', header = T)

i <- which(rownames(et$genes) %in% genes_list$GeneID) # pull out the genes of interest



