library(SBGNview)
library(gage)
library(tidyverse)
library(SummarizedExperiment)

load("out/DGEList_UNvsWT.Rdata")
load("out/exactTest_UNvsWT.Rdata")
data("pathways.info", "IFNg")
data(sbgn.xmls)
data('mapped.ids')

expression_data <- read.table("out/EM_expressionFile_RNAseq.txt", header = T) |> 
  select(c(1, 3:8)) |> 
  column_to_rownames(var="Name")

expression_data_log <- log(expression_data)

logFC <- read.table("out/UNvsWT_exactTest.txt", header=T,sep = '\t') |> 
  rownames_to_column(var="Name") |> 
  select(c(1:2))
head(logFC)
  
head(expression_data)

wt.cols <- which(x$samples$group == "WT")
un.cols <- which(x$samples$group == "UN")

ensembl.pathway <- sbgn.gsets(id.type = "ENSEMBL",
                              species = "hsa",
                              mol.type = "gene",
                              output.pathway.name = TRUE
)

degs <- gage(exprs = expression_data,
             gsets = ensembl.pathway,
             ref = un.cols,
             samp = wt.cols,
             compare = "paired" #"as.group"
)
head(degs$greater)[,3:5]
head(degs$less)[,3:5]
down.pathways <- row.names(degs$less)[1:10]
head(down.pathways)

down.pathways <- sapply(strsplit(down.pathways,"::"), "[", 1)
head(down.pathways)

# logCPM differences
ensembl.unVsWT <- expression_data_log[,wt.cols] - expression_data_log[,un.cols]
head(ensembl.unVsWT)

sbgnview.obj <- SBGNview(
  gene.data = ensembl.unVsWT,
  gene.id.type = "ENSEMBL",
  input.sbgn = "R-HSA-1799339",
  output.file = "./plots/sbgn/UNvsWT.sbgn_down",
  show.pathway.name = TRUE,
  max.gene.value = 2,
  min.gene.value = -2,
  mid.gene.value = 0,
  node.sum = "mean",
  output.format = c("pdf"),
  
  font.size = 2.3,
  org = "HSA",
  
  text.length.factor.complex = 3,
  if.scale.compartment.font.size = TRUE,
  node.width.adjust.factor.compartment = 0.04, 
  
)
sbgnview.obj









