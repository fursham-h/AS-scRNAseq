# Step 1: import dataset into R as matrix
# hint: a matrix should contain only count integers.
## The first row of the interneuron data has column names, 
## and the first column of the data has feature names.
## We need to import the data in such a way that the first row will be converted to column names
## and the first column converted into row names
library(Seurat)
library(patchwork)
library(dplyr)

setwd("~/Desktop/Seurat")

input.matrix <- read.table("GSE109796_Oscar.GEO.singleCell.gene.count.txt", header = TRUE, sep = "", row.names = 1) #F: Good
rownames(input.matrix) %>% head(n=5)
# [1] "ENSMUSG00000000001|GNAI3"
# [2] "ENSMUSG00000000003|PBSN" 
# [3] "ENSMUSG00000000028|CDC45"
# [4] "ENSMUSG00000000031|H19"  
# [5] "ENSMUSG00000000037|SCML2"
colnames(input.matrix) %>% head(n=5)
# [1] "C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all"
# [2] "C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all" 
# [3] "C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all" 
# [4] "C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all" 
# [5] "C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all" 

## F: The above code heads the names of columns and rows. 
# If you want to preview the first 5 column and rows, do this:
input.matrix[1:5,1:5]
# C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all
# ENSMUSG00000000001|GNAI3                                        0                                       0                                       0
# ENSMUSG00000000003|PBSN                                         0                                       0                                       0
# ENSMUSG00000000028|CDC45                                        0                                       0                                    1439
# ENSMUSG00000000031|H19                                        189                                       1                                      84
# ENSMUSG00000000037|SCML2                                        0                                       0                                       0
# C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all
# ENSMUSG00000000001|GNAI3                                       1                                       0
# ENSMUSG00000000003|PBSN                                        0                                       0
# ENSMUSG00000000028|CDC45                                      31                                       0
# ENSMUSG00000000031|H19                                       151                                    2231
# ENSMUSG00000000037|SCML2                                       0                                       0

# Step 2: Create Seurat object 
mydata <- CreateSeuratObject(counts = input.matrix, min.cells = 3, min.genes = 200, project = "interneuron") #F: Good
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Warning: Feature names cannot have pipe characters ('|'), replacing with dashes ('-')

mydata
# An object of class Seurat 
# 25616 features across 2669 samples within 1 assay 
# Active assay: RNA (25616 features, 0 variable features)

# Step 3: Perform QC on dataset and subset data further 
## I just realize that MT genes have been annotated in the list of features.
## However, it is named like this "ENSMUSG00000064336-MT-TF"
## Therefore, the pattern input for PercentageFeatureSet function will be different
mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "^ENSMUSG00000064336-MT-TF")
head(mydata@meta.data, 5)
plot1 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#F: The pattern above is not correct. 
# Remember, in the Seurat vignette, the authors used "^MT-" because the MT genes start with the letters MT- (for example MT-TF)
# The ^ character means that you are finding a pattern of test that BEGINS with MT-
# In this interneuron dataset, the MT genes do not start with MT-, but this pattern -MT- is found in the middle of its name.
# Try correcting the pattern one more time
mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "CHANGE THIS PATTERN")
VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #F: Plot this violin plot out
head(mydata@meta.data, 5)
plot1 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

mydata <- subset(mydata, subset = nFeature_RNA > ? & nFeature_RNA < ? & percent.mt < ?) #F: filter the data up a little after viewing the violin plot

#F: One more thing, can you add the plots generated in the folder "Plots" that I made in the parent directory
# You can export the plots using the export button above the graphs and drop into the folder


##########################F: You can stop here, and we can discuss on the plots before moving further

# Step 4: Normalize, identify highly variable features and scale data 
mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 10000)
mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mydata), 10)
plot1 <- VariableFeaturePlot(mydata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# When using repel, set xnudge and ynudge to 0 for optimal results
# Warning message:
# Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
# Please use `as_label()` or `as_name()` instead.
# This warning is displayed once per session. 
plot1 + plot2
# Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : 
#   Viewport has zero dimension(s)

all.genes <- rownames(mydata)
# There were 13 warnings (use warnings() to see them)
mydata <-ScaleData(mydata, features = all.genes)

#F: Try Scaling the data with only the highly-variant features too
                   
# Step 5: carry out linear dimensional reduction




# Step 6: Cluster cells




# Step 7: Perform non-linear dimensional reduction 



