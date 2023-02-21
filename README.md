# williams_et_al
 Analysis, files, and figures for the following publication: "Genome-wide characterization of the Caenorhabditis elegans intestine GATA transcription factor ELT-2"

This repository serves to simplify and supersede the exploratory analysis performed for this publication. Code related to initial exploratory analysis can be found in the following repository: [ELT-2-ChIP-revision](https://github.com/meekrob/ELT-2-ChIP-revision)

Analysis for this publication was performed using the R language. To utilize and explore the code in this repository, clone the repository and open the `williams_et_al.Rproj` file using RStudio. Relevant package installation commands are listed where necessary. Output files and plots were generated using R version 4.1.0 and RStudio version 1.4.1106.

Below describes the `analysis` directory structure:

- `01_tissue_specific_genes`: Generate a list of "tissue specific genes" based on [WormBase Ontology Browser](https://wormbase.org/tools/ontology_browser).
- `02_emb_L1_L3_intestine_RNAseq`: RNA-seq analysis of FACS isolated intestine samples using the DESeq2 package.
- `03_elt2_RNAseq`: RNA-seq analysis of elt-2 (-) vs. wt samples ([Dineen and Nishimura et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29360433/)) using the DESeq2 package.
- `04_promoters`: Analysis of ELT-2 ChIP-seq data downloaded from the [modERN](https://epic.gs.washington.edu/modERN/) resource. Analysis in this directory was performed by David King.
- `05_elt2_target_analysis`: Quantification and visualization gene set categories corresponding to intestine enrichment and ELT-2 regulation. Additional analysis includes Gene Ontology (GO) analysis and quantification of microscopy images from translation reporter worm strains.
- `06_intestine_enriched_genes`: GO analysis of intestine enriched genes evaluated generated in the intestine FACS RNA-seq dataset.
- `07_expression_heatmaps`: generate heatmaps for transcript abundance in the intestine transcriptome and elt-2 (-) RNA-seq datasets.
- `08_elt2_promoter_regulation`: evaluate if there is a higher abundance of RNA-seq reads aligning to the elt-2 promoter in elt-2(-) compared to wildtype

Within each directory are three sub-directories:

- `01_input`: contains any input data necessary to perform the analysis
- `02_scripts`: contains the R markdown scripts that were utilized to perform the analysis. Corresponding HTML and PDF files are associated with each R markdown script to facilitate documentation of the analysis performed.
- `03_output`: Contains any output files or figures generated in the R markdown scripts.
