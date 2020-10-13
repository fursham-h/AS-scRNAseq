# Step 1: import dataset into R as matrix
# hint: a matrix should contain only count integers.
## The first row of the interneuron data has column names, 
## and the first column of the data has feature names.
## We need to import the data in such a way that the first row will be converted to column names
## and the first column converted into row names
> library(Seurat)
> library(patchwork)
> library(dplyr)

> setwd("~/Desktop/Seurat")

input.matrix <- read.table("GSE109796_Oscar.GEO.singleCell.gene.count.txt", header = TRUE, sep = "", row.names = 1)
> rownames(input.matrix) %>% head(n=5)
[1] "ENSMUSG00000000001|GNAI3"
[2] "ENSMUSG00000000003|PBSN" 
[3] "ENSMUSG00000000028|CDC45"
[4] "ENSMUSG00000000031|H19"  
[5] "ENSMUSG00000000037|SCML2"
> colnames(input.matrix) %>% head(n=5)
[1] "C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all"
[2] "C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all" 
[3] "C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all" 
[4] "C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all" 
[5] "C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all" 

# Step 2: Create Seurat object 
> mydata <- CreateSeuratObject(counts = input.matrix, min.cells = 3, min.genes = 200, project = "interneuron")
Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
Warning: Feature names cannot have pipe characters ('|'), replacing with dashes ('-')
> mydata
An object of class Seurat 
25616 features across 2669 samples within 1 assay 
Active assay: RNA (25616 features, 0 variable features)

# Step 3: Perform QC on dataset and subset data further 
## I just realize that MT genes have been annotated in the list of features.
## However, it is named like this "ENSMUSG00000064336-MT-TF"
## Therefore, the pattern input for PercentageFeatureSet function will be different



# Step 4: Normalize, identify highly variable features and scale data 




# Step 5: carry out linear dimensional reduction




# Step 6: Cluster cells




# Step 7: Perform non-linear dimensional reduction 



