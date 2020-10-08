# Step 1: import dataset into R


# Step 2: re-annotate MT genes
# I will assist with this 

# Step 2.1: download mouse reference genome
# firstly, we need to install a package that contain genome databases 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("AnnotationHub")

# load mouse genome
library(AnnotationHub) 
ah <- AnnotationHub()

# retrieve latest GRCm38 annotation from ensembl
mouse_genome <- ah[['AH60127']]

# Step 2.2: prepare dataframe containing a list of gene-if, gene_name and chromosome origin
mouse_genes <- distinct(data.frame(gene_id = mouse_genome$gene_id, gene_name = mouse_genome$gene_name, chromosome = as.character(seqnames(mouse_genome)))) 
#to be continued

# Step 3: Create Seurat object 


# Step 4: Perform QC on dataset and subset data further 


# Step 5: Normalize, identify highly variable features and scale data 


# Step 6: carry out linear dimensional reduction


# Step 7: Cluster cells


# Step 8: Perform non-linear dimensional reduction 

