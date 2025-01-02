# Exploratory Data Analysis Utilities
Custom functions that can be leveraged for conducting exploratry data analysis for RNA-Seq data and MGX / 16S. 

# Table of Contents
- [Dendrogram](#dendrogram)
- [Distances](#distances)
- [PCA](#pca)
- [Eigencorrelations](#eigencorrelations)

# Dendrogram
Function to plot dendrogram of samples to visualize clustering

Arguments:

- VSD: a variance stabilized or rlog DESeq2 transformation
- META: a metadata table
- FEATURE: the column in your metadata to color by within ""
- COLORS: a color palette
- OUTPUT: the output folder to use

# Distances
Function to plot Euclidean distances to visualize feature specific clustering

Arguments: 

- VSD: a variance stabilized or rlog DESeq2 transformation
- META: a metadata table
- FEATURES: a comma separated vector of META column names to include as annotations
- palettes: a list of color palettes for each feature you are looking at
- OUTPUT: the output folder name encased with ""

# PCA
Functions to determine number of variablel features to include and calculation of principal components/PCA plot

**determineVarFeatures**

arguments:
- VSD: a variance stabilized or rlog DESeq2 transformation

**explorePCA**

arguments:

- VSD: a variance stabilized or rlog DESeq2 transformation
- NFEATURES: the number of features to include


# Eigencorrelations
Function to calculate eigencorrelations (correlations between prinicipal components and metadata features

Arguments: 

- VSD: a variance stabilized transformed object
- Metadata: a metadata table
- num_pcs: the number of pcs to look at
