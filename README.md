# The Legionella and THP1 Interactome 

## Getting the Counts Matrix and Quality Control

From the Short Read Archive (SRA), for Gene Expression Omnibus(GEO) Experiment GSE207487 the RNAseq raw gene level counts were obtained from the [GREIN](https://www.ilincs.org/apps/grein/?gse=GSE207487) server analysis. This uses the hg38 assembly from Ensembl release 91 (2017) and SALMON to quantify the reads.  The metadata and fastqc reports (as summarized by MultiQC) were downloaded as well. 

Another option for this is to do the alignment via HISAT2/Salmon manually. Benefits include the ability to use a more modern chromosome assembly and to have more control over the data. The drawbacks are the significant time/compute/storage requirement. 