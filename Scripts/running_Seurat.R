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
#                          C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all
# ENSMUSG00000000001|GNAI3                                        0
# ENSMUSG00000000003|PBSN                                         0
# ENSMUSG00000000028|CDC45                                        0
# ENSMUSG00000000031|H19                                        189
# ENSMUSG00000000037|SCML2                                        0
#                          C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all
# ENSMUSG00000000001|GNAI3                                       0
# ENSMUSG00000000003|PBSN                                        0
# ENSMUSG00000000028|CDC45                                       0
# ENSMUSG00000000031|H19                                         1
# ENSMUSG00000000037|SCML2                                       0
#                          C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all
# ENSMUSG00000000001|GNAI3                                       0
# ENSMUSG00000000003|PBSN                                        0
# ENSMUSG00000000028|CDC45                                    1439
# ENSMUSG00000000031|H19                                        84
# ENSMUSG00000000037|SCML2                                       0
#                          C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all
# ENSMUSG00000000001|GNAI3                                       1
# ENSMUSG00000000003|PBSN                                        0
# ENSMUSG00000000028|CDC45                                      31
# ENSMUSG00000000031|H19                                       151
# ENSMUSG00000000037|SCML2                                       0
#                          C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all
# ENSMUSG00000000001|GNAI3                                       0
# ENSMUSG00000000003|PBSN                                        0
# ENSMUSG00000000028|CDC45                                       0
# ENSMUSG00000000031|H19                                      2231
# ENSMUSG00000000037|SCML2                                       0

# Step 2: Clean up input dataset
# We will apply several filters to select for high quality cells and contributing features for this analysis

## Filter 1: Select coding genes
### We need to import a mouse GRCm38 annotation and parse data for gene_ids of coding genes
mouse.gtf <- rtracklayer::import("Mus_musculus.GRCm38.101.gtf")

#F: first things first, you can preview the type of genes in the annotation
unique(mouse.gtf$gene_biotype)
# [1] "TEC"                                "snRNA"                              "protein_coding"                    
# [4] "processed_pseudogene"               "antisense"                          "sense_intronic"                    
# [7] "lincRNA"                            "processed_transcript"               "miRNA"                             
# [10] "snoRNA"                             "misc_RNA"                           "transcribed_unprocessed_pseudogene"
# [13] "unprocessed_pseudogene"             "sense_overlapping"                  "rRNA"                              
# [16] "transcribed_processed_pseudogene"   "ribozyme"                           "unitary_pseudogene"                
# [19] "scaRNA"                             "pseudogene"                         "polymorphic_pseudogene"            
# [22] "bidirectional_promoter_lncRNA"      "transcribed_unitary_pseudogene"     "macro_lncRNA"                      
# [25] "3prime_overlapping_ncRNA"           "translated_unprocessed_pseudogene"  "TR_V_gene"                         
# [28] "TR_V_pseudogene"                    "TR_D_gene"                          "TR_J_gene"                         
# [31] "TR_C_gene"                          "TR_J_pseudogene"                    "IG_LV_gene"                        
# [34] "IG_V_gene"                          "IG_V_pseudogene"                    "IG_J_gene"                         
# [37] "IG_C_gene"                          "sRNA"                               "scRNA"                             
# [40] "IG_C_pseudogene"                    "IG_D_gene"                          "IG_D_pseudogene"                   
# [43] "IG_pseudogene"                      "Mt_tRNA"                            "Mt_rRNA"    


#F: As you can see, there are alot of non-coding genes in the annotation. 
# it sort of make more sense to select only protein_coding genes since it makes up the bulk of the annotation

#F: Here, I am going to generate a table on the number of genes in each gene_biotype category
# You don't need to understand it, but maybe you could try to find out what does "%>%" does
mouse.gtf %>% 
  as.data.frame() %>% 
  filter(type == "gene") %>% 
  group_by(gene_biotype) %>% 
  tally() %>% 
  arrange(desc(n))
# # A tibble: 45 x 2
# gene_biotype               n
# <chr>                  <int>
#   1 protein_coding         21936
# 2 processed_pseudogene   10003
# 3 lincRNA                 5629
# 4 TEC                     3238
# 5 antisense               2991
# 6 unprocessed_pseudogene  2723
# 7 miRNA                   2207
# 8 snoRNA                  1507
# 9 snRNA                   1385
# 10 processed_transcript     779
# # … with 35 more rows

#F: Make a dataframe containing gene_name and gene_id of protein coding genes
coding.genes <- mouse.gtf %>% 
  as.data.frame() %>% 
  dplyr::select(gene_name, gene_id, gene_biotype) %>% 
  filter(gene_biotype == "protein_coding")

#F: Now we will filter input.matrix to only select genes containing protein_coding gene_ids
# you might need tidyverse for this
install.packages("tidyverse")

# 1) extract the feature names from input.matrix into a dataframe
features.df <- data.frame(names = rownames(input.matrix))
head(features.df)
# names
# 1 ENSMUSG00000000001|GNAI3
# 2  ENSMUSG00000000003|PBSN
# 3 ENSMUSG00000000028|CDC45
# 4   ENSMUSG00000000031|H19
# 5 ENSMUSG00000000037|SCML2
# 6  ENSMUSG00000000049|APOH

# 2) I want to split the gene_id and gene_name from each other.
features.df <- features.df %>% 
  separate(names, c("gene_id","gene_name"), sep = "\\|")
head(features.df)
# gene_id gene_name
# 1 ENSMUSG00000000001     GNAI3
# 2 ENSMUSG00000000003      PBSN
# 3 ENSMUSG00000000028     CDC45
# 4 ENSMUSG00000000031       H19
# 5 ENSMUSG00000000037     SCML2
# 6 ENSMUSG00000000049      APOH

# 3) add a new column to annotate if the gene_id corresponds to a coding gene
features.df <- features.df %>% 
  mutate(coding = ifelse(gene_id %in% coding.genes$gene_id,T,F))
head(features.df)
# gene_id gene_name coding
# 1 ENSMUSG00000000001     GNAI3   TRUE
# 2 ENSMUSG00000000003      PBSN   TRUE
# 3 ENSMUSG00000000028     CDC45   TRUE
# 4 ENSMUSG00000000031       H19  FALSE
# 5 ENSMUSG00000000037     SCML2   TRUE
# 6 ENSMUSG00000000049      APOH   TRUE

#F: you can calculate how many features are retained
sum(features.df$coding)
# [1] 20087

#4) subset input.matrix to keep only coding genes
input.matrix <- input.matrix[features.df$coding,]
nrow(input.matrix)
# [1] 20087  #F: See it works ;P

#F: NEW, I want to rename the rows to just use the gene_names
input.matrix <- input.matrix %>% 
  rownames_to_column("ID") %>% 
  separate(names, c("gene_id","ID"), sep = "\\|") %>%
  dplyr::select(-gene_id) %>%
  group_by(ID) %>% 
  mutate(across(everything(), sum)) %>% 
  ungroup() %>% 
  distinct(ID, .keep_all = T) %>% 
  column_to_rownames("ID") %>% 
  as.matrix()

#Filter 2: This next part is a abit tedious to understand, but bear with me. We can go through this
# you need matrixStats package
install.packages("matrixStats")

#F: we will filter for cells that contain at least 50000 reads mapping to coding features(which we have already subsetted)
input.matrix <- input.matrix[,colSums2(as.matrix(input.matrix)) > 50000]
ncol(input.matrix)
# [1] 2658


#Filter 3: Next, We will try to remove features that are expressed in less than 10 cells with less than 5CPM
# 5CPM means that instead of counts, each feature in each gene have been normalized.
# So we will use Seurat NormalizeData function to get this normalized value, but we will later re-normalize after filtering is done
testdata <- CreateSeuratObject(counts = input.matrix, project = "filter")
testdata <- NormalizeData(testdata, normalization.method = "LogNormalize", scale.factor = 1000000)

# we will select features that have 10 cells with at least than 5CPM expression
features.to.keep <-  apply(testdata[["RNA"]]@counts, 1, function(x){sum(x>=5) >=10})

# and then filter the input.matrix
input.matrix <- input.matrix[features.to.keep,]
input.matrix[1:5,1:5]
# C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all
# ENSMUSG00000000001|GNAI3                                        0                                       0
# ENSMUSG00000000028|CDC45                                        0                                       0
# ENSMUSG00000000037|SCML2                                        0                                       0
# ENSMUSG00000000056|NARF                                       193                                       0
# ENSMUSG00000000078|KLF6                                        13                                       0
# C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all
# ENSMUSG00000000001|GNAI3                                       0                                       1
# ENSMUSG00000000028|CDC45                                    1439                                      31
# ENSMUSG00000000037|SCML2                                       0                                       0
# ENSMUSG00000000056|NARF                                        0                                       1
# ENSMUSG00000000078|KLF6                                        0                                       0
# C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all
# ENSMUSG00000000001|GNAI3                                       0
# ENSMUSG00000000028|CDC45                                       0
# ENSMUSG00000000037|SCML2                                       0
# ENSMUSG00000000056|NARF                                        0
# ENSMUSG00000000078|KLF6                                        0


# Step 3: Create Seurat object 
mydata <- CreateSeuratObject(counts = input.matrix, min.cells = 3, min.genes = 200, project = "interneuron") #F: Good
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Warning: Feature names cannot have pipe characters ('|'), replacing with dashes ('-')

mydata
# An object of class Seurat 
# 25616 features across 2669 samples within 1 assay 
# Active assay: RNA (25616 features, 0 variable features)

# Step 4: Perform QC on dataset and subset data further 
## I just realize that MT genes have been annotated in the list of features.
## However, it is named like this "ENSMUSG00000064336-MT-TF"
## Therefore, the pattern input for PercentageFeatureSet function will be different
mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "^MT-") #F: I changed the pattern, since the naming is changed
head(mydata@meta.data, 5)
#              orig.ident
# C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all interneuron
# C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all  interneuron
# C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all  interneuron
# C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all  interneuron
# C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all  interneuron
#                                          nCount_RNA
# C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all     651271
# C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all     1027481
# C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all     1462016
# C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all     1220828
# C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all      553652
#                                          nFeature_RNA
# C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all         3198
# C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all          2847
# C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all          3444
# C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all          3495
# C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all          3371
#                                          percent.mt
# C1.101.A10_CGAGGCTG.GCGTAAGA_L008_R1_all  0.5816319
# C1.101.A1_TAAGGCGA.GCGTAAGA_L008_R1_all   2.4223319
# C1.101.A4_TCCTGAGC.GCGTAAGA_L008_R1_all   2.6754153
# C1.101.A5_GGACTCCT.GCGTAAGA_L008_R1_all   3.2504988
# C1.101.A6_TAGGCATG.GCGTAAGA_L008_R1_all   1.1888334
plot1 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mydata <- subset(mydata, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 5)


# Step 5: Normalize, identify highly variable features and scale data 
mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 1000000)


#F: THe authors used a different method for selection.method. Another free gift from me ;P
mydata <- FindVariableFeatures(mydata, selection.method = "mvp", mean.cutoff = c(0.5,8), dispersion.cutoff = c(0.5,5))
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

#F: I scaled the data based on highly variable features
mydata <- ScaleData(mydata)



# Step 6: Regress out cell cycle determinants
## Collect s and g2m genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Score Cell cycle genes
mydata <- CellCycleScoring(mydata, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# RUn a PCA plot to see if cells cluster based on mitotic phase
mydata <- RunPCA(mydata, features = c(s.genes, g2m.genes))
DimPlot(mydata, reduction = "pca")

# Regress out cell cycle factors
mydata <- ScaleData(mydata, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(mydata))

# Re-check PCA plot
mydata <- RunPCA(mydata, features = c(s.genes, g2m.genes))
DimPlot(mydata, reduction = "pca")



# Step 5: carry out linear dimensional reduction




# Step 6: Cluster cells




# Step 7: Perform non-linear dimensional reduction 



