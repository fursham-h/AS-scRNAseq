# Step 1: import dataset into R as matrix
# hint: a matrix should contain only count integers.
## The first row of the interneuron data has column names, 
## and the first column of the data has feature names.
## We need to import the data in such a way that the first row will be converted to column names
## and the first column converted into row names
input.matrix <- # complete this code



# Step 2: Create Seurat object 


  
  
# Step 3: Perform QC on dataset and subset data further 
## I just realize that MT genes have been annotated in the list of features.
## However, it is named like this "ENSMUSG00000064336-MT-TF"
## Therefore, the pattern input for PercentageFeatureSet function will be different



# Step 4: Normalize, identify highly variable features and scale data 




# Step 5: carry out linear dimensional reduction




# Step 6: Cluster cells




# Step 7: Perform non-linear dimensional reduction 



