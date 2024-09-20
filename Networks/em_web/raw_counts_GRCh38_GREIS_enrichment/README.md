EnrichmentMap::RNA-Seq - raw_counts_GRCh38_GREIS
------------------------------------------------

Network Permalink: https://enrichmentmap.org/document/9aad8240-46bb-47c8-8701-75f54f12a795

EnrichmentMap is a web-app that allows you to perform functional enrichment analysis on 
gene lists derived from RNA-seq experiments and visualise the results as a network.

This archive contains the following files:
* images/enrichment_map_large.png
* images/enrichment_map_medium.png
* images/enrichment_map_small.png
  * Network PNG images in various sizes.
* images/node_color_legend_NES.svg
  * An SVG image of the NES color legend used for the nodes in the network.
* data/enrichment_results.txt
  * Results of Gene Set Enrichment Analysis from the FGSEA R package.
* data/ranks.txt
  * Gene ranks.
* data/gene_sets.gmt
  * Gene sets (pathways) that correspond to nodes in the network.


How to cite EnrichmentMap
-------------------------
To cite this app in a paper, for now, please cite this Nature Protocols 
article (an article specific to this app will be published shortly):
https://doi.org/10.1038/s41596-018-0103-9
Reimand, J., Isserlin, R., ..., Bader, G.
Pathway enrichment analysis and visualization of omics data using g:Profiler, GSEA, Cytoscape and EnrichmentMap.
Nat Protoc 14, 482–517 (2019).


Importing data into the Cytoscape EnrichmentMap App
---------------------------------------------------
* Download and install Cytoscape
  * https://cytoscape.org/download.html
* Download and install the EnrichmentMap App
  * https://apps.cytoscape.org/apps/enrichmentmap
* (optional) If your network was created from an RNA-seq expression file you may copy 
  it to the 'data' folder.
* Start Cytoscape.
* Go to the main menu and select *Apps > EnrichmentMap*
* Click the button that says *Add* then select *Scan a folder for enrichment data*.
* Select the 'data' folder.
* The three files contained in this archive will show up in the *Enrichments*, *GMT* 
  and *Ranks* fields in the dialog.
  * Note, if you copied the RNA-seq expression file to the zip output 'data' 
    folder it should also appear in the *Expressions* field. If it does not then rename 
    the file to include the word “expression”, then try again.
* Click the *Build* button.
* Documentation for the EnrichmentMap Cytoscape App is available here...
  * https://enrichmentmap.readthedocs.io/


Gene-set filtering parameters
-----------------------------
* The Baderlab gene set database used for this network was:
  * Human_GOBP_AllPathways_noPFOCR_no_GO_iea_May_01_2024_symbol.gmt
* The following cutoff parameters were used to filter the results of enrichment analysis.
  * Gene-sets with q-value (padj) greater than 0.05 are removed from the network.
  * Edges represent similarity between gene-sets. 
    * Similarity is calculated using the JACCARD method and must have a 
      value of at least 0.25.
    * JACCARD coefficient is computed as the size of the intersection divided by the size of the union.