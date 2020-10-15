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
input.matrix[1:5,1:5]
                         C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all
ENSMUSG00000000001|GNAI3                                        0
ENSMUSG00000000003|PBSN                                         0
ENSMUSG00000000028|CDC45                                        0
ENSMUSG00000000031|H19                                        189
ENSMUSG00000000037|SCML2                                        0
                         C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all
ENSMUSG00000000001|GNAI3                                       0
ENSMUSG00000000003|PBSN                                        0
ENSMUSG00000000028|CDC45                                       0
ENSMUSG00000000031|H19                                         1
ENSMUSG00000000037|SCML2                                       0
                         C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all
ENSMUSG00000000001|GNAI3                                       0
ENSMUSG00000000003|PBSN                                        0
ENSMUSG00000000028|CDC45                                    1439
ENSMUSG00000000031|H19                                        84
ENSMUSG00000000037|SCML2                                       0
                         C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all
ENSMUSG00000000001|GNAI3                                       1
ENSMUSG00000000003|PBSN                                        0
ENSMUSG00000000028|CDC45                                      31
ENSMUSG00000000031|H19                                       151
ENSMUSG00000000037|SCML2                                       0
                         C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all
ENSMUSG00000000001|GNAI3                                       0
ENSMUSG00000000003|PBSN                                        0
ENSMUSG00000000028|CDC45                                       0
ENSMUSG00000000031|H19                                      2231
ENSMUSG00000000037|SCML2                                       0

## F: The above code heads the names of columns and rows. 
# If you want to preview the first 5 column and rows, do this:

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
mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "-MT")
head(mydata@meta.data, 5)
             orig.ident
C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all interneuron
C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all  interneuron
C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all  interneuron
C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all  interneuron
C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all  interneuron
                                         nCount_RNA
C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all     651271
C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all     1027481
C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all     1462016
C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all     1220828
C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all      553652
                                         nFeature_RNA
C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all         3198
C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all          2847
C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all          3444
C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all          3495
C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all          3371
                                         percent.mt
C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all  0.5816319
C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all   2.4223319
C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all   2.6754153
C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all   3.2504988
C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all   1.1888334
plot1 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mydata <- subset(mydata, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 10)
#F: The pattern above is not correct. 
# Remember, in the Seurat vignette, the authors used "^MT-" because the MT genes start with the letters MT- (for example MT-TF)
# The ^ character means that you are finding a pattern of test that BEGINS with MT-
# In this interneuron dataset, the MT genes do not start with MT-, but this pattern -MT- is found in the middle of its name.
# Try correcting the pattern one more time

mouse.gtf <- rtracklayer::import("Mus_musculus.GRCm38.101.gtf")
gtf_df <- as.data.frame(mouse.gtf1)
geneid_df <- dplyr::select(gtf_df,c(gene_name,gene_id,gene_biotype))

ribo.genes <- subset(geneid_df, subset = gene_biotype == "rRNA")

ribo.genes[1:5,]
     gene_name            gene_id gene_biotype
9476  n-R5s209 ENSMUSG00000095256         rRNA
9477  n-R5s209 ENSMUSG00000095256         rRNA
9478  n-R5s209 ENSMUSG00000095256         rRNA
9866   Gm27362 ENSMUSG00000098529         rRNA
9867   Gm27362 ENSMUSG00000098529         rRNA

miRNA <- subset(geneid_df, subset = gene_biotype == "miRNA")
> miRNA[1:5,]
     gene_name            gene_id gene_biotype
889    Gm22463 ENSMUSG00000093015        miRNA
890    Gm22463 ENSMUSG00000093015        miRNA
891    Gm22463 ENSMUSG00000093015        miRNA
1988   Gm23358 ENSMUSG00000093970        miRNA
1989   Gm23358 ENSMUSG00000093970        miRNA

#C: Hmm how do I make it so that it can subset with several conditions? As in gene_biotype == "miRNA" OR "rRNA" OR "processed_pseudogene" etc? I tried |



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







#F: Try Scaling the data with only the highly-variant features too
                   
# Step 5: carry out linear dimensional reduction




# Step 6: Cluster cells




# Step 7: Perform non-linear dimensional reduction 



