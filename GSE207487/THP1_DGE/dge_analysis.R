library(edgeR)
library(Homo.sapiens)

# GO and KEGG analysis with edgeR ----------
load("out/DGEList_UNvsWT.Rdata")
load("out/exactTest_UNvsWT.Rdata")

entrez_ids <- select(org.Hs.eg.db, keys=rownames(et),
                 columns=c("ENSEMBL","ENTREZID"), keytype="ENSEMBL")

dge_table <- merge(et$table, entrez_ids, by.x = 0, by.y = "ENSEMBL")
dge_table
