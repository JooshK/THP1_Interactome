library(edgeR)
library(Homo.sapiens)
library(gage)
library(pathview)
library(SBGNview)

# GO and KEGG analysis with edgeR ----------
load("out/DGEList_UNvsWT.Rdata")
load("out/exactTest_UNvsWT.Rdata")

# convert to entrez ids 
entrez_ids <- select(org.Hs.eg.db, keys=rownames(et),
                 columns=c("ENSEMBL","ENTREZID"), keytype="ENSEMBL")

et$table <- merge(et$table, entrez_ids, by.x = 0, by.y = "ENSEMBL")
et$genes <- merge(et$genes, entrez_ids, by.x = 1, by.y = "ENSEMBL")

# GO enrichment
go <- goana(et, geneid = 'ENTREZID')
write.csv(go, file = "out/GO_UNvsWT.csv")
topGO(go, ontology = "BP", n = 20)

# KEGG enrichment
kk <- kegga(et, geneid = "ENTREZID")
topKEGG(kk, sort = "Up")
write.csv(kk, file = "out/KEGG_UNvsWT.csv")

# GAGE and Pathview ------------
data("kegg.gs") # load the kegg databases
data("hsa_pathwayCommons_ENSEMBL")

edger.fc=et$table$logFC
names(edger.fc)=et$table$ENTREZID

# Analysis for up-regulated pathways
fc.kegg.p <- gage(edger.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &
  !is.na(fc.kegg.p$greater[, "q.val"]) # select for significant pathways
path.ids <- rownames(fc.kegg.p$greater)[sel]

sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
   !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)

# visualize the top 3 pathways
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(
  gene.data = edger.fc, pathway.id = pid, 
  low = c("gene" = "red"), high = "green"))


# SBGNview -----------
load("out/exactTest_UNvsWT.Rdata") # reload exact test
data("sbgn.xmls") # loads the sbgn xml collection
data("pathways.info")

edger.fc <- et$table$logFC
names(edger.fc) <- rownames(et$table)

ensembl.pathway <- sbgn.gsets(id.type = "ENSEMBL",
                              mol.type = "gene",
                              output.pathway.name = TRUE
)

# re run gage with ENSEMBL 
degs <- gage(edger.fc, gsets = ensembl.pathway)
up.pathways <- row.names(degs$greater)[1:10]
up.pathways <- sapply(strsplit(up.pathways,"::"), "[", 1)
head(up.pathways)

#  SBGN view driver function
sbgnview.obj <- SBGNview(
  gene.data = et$table$logFC,
  gene.id.type = "ENSEMBL",
  input.sbgn = up.pathways[1:2],
  output.file = "UNvsWT",
  show.pathway.name = TRUE,
  max.gene.value = 2,
  min.gene.value = -2,
  mid.gene.value = 0,
  node.sum = "mean",
  output.format = c("pdf")
)






















