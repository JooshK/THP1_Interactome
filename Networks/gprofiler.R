library(gprofiler2)

# THP1 RNAseq data - q value < 0.05
sig_genes <- read.table("./data/THP1_RNAseq/qval005_UNvsWT.txt")

upload_GMT_file(gmtfile = "./data/gene_sets/GO_AllPathways_withPFOCR_noiea_KEGG_Sept_16.gmt") # returns a file ID 
gostres <- gost(sig_genes$V1, organism = 'gp__sELz_ES6x_jaM')

gem <- gost$result[,c("term_id", "term_name", "p_value", "intersection_size")]
colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem$Phenotype = "+1"
gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
head(gem, 3)

write.table(gem, file = "./data/gprofiler/BaderLab_combined_gem.txt", sep = "\t", quote = F, row.names = F)

gostres2 <- gost(query = sig_genes$V1, 
                 organism = "hsapiens", ordered_query = FALSE, 
                 multi_query = FALSE, significant = TRUE, exclude_iea = TRUE, 
                 measure_underrepresentation = FALSE, evcodes = TRUE, 
                 user_threshold = 0.05, correction_method = "g_SCS", 
                 domain_scope = "annotated", custom_bg = NULL, 
                 numeric_ns = "", sources = NULL, highlight = TRUE)
