# Alternative splicing analysis of single cell RNA-Sequencing (scRNA-Seq) dataset

This repository contain resources and guides for the single-cell alternative splicing project.  

The aim of this project is to identify regulated alternative splicing events between cell clusters from scRNA-Seq experiments. Here are the key objectives:

1. Understand a typical scRNA-seq workflow (Smart-Seq vs 10X Genomics)
2. Review current approaches to analyze single-cell alternative splicing and its limitations
3. Process published scRNA-seq dataset to identify cell clusters
4. Merge sequencing reads from cell clusters and compare alternative splicing landscape between clusters
5. Identify enrichment in gene groups regulated by alternative splicing

### Resources
Below are several resources that are useful for this project.

- scRNA-seq datasets
	1. [Early emergence of cortical interneuron diversity in the mouse embryo, 2018; Mi et al](https://science.sciencemag.org/content/360/6384/81)
		- Data can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109796)
	2. [Shared and distinct transcriptomic cell types across neocortical areas, 2018; Tasic et al](https://www.nature.com/articles/s41586-018-0654-5)
		- Data can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115746) or [here](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq)
- Programs/packages for bioinformatics analysis
	1. [Seurat, for QC and clustering of single-cell data](https://satijalab.org/seurat/)
	2. [Whippet, for analyzing alternative splicing changes](https://github.com/timbitz/Whippet.jl)

### Tasks
Some of the current tasks that can be done:

- [ ] Literature review of scRNA-seq workflows (Smart-Seq vs 10X Genomics)
- [ ] Literature review of single-cell alternative splicing analysis
- [x] Download count matrix from scRNA-seq dataset
- [x] Familiarize with Seurat package and R programming
- [x] Import scRNA-seq matrix into R
- [x] Carry QC on dataset and normalise+scale data
- [x] Perform dimensional reduction (preferably UMAP/tSNE) and create clusters 
- [ ] Identify cluster biomarkers and infer its cell type (if possible)
- [ ] Label tsne/UMAP plot with cell types annotated from main paper
- [ ] Plan pipeline for whippet analysis on clustered scRNA-seq transcriptome
