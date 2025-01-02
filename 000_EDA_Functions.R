#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Functions for exploratory data analysis (bulk RNASeq)
#
# Lab: Dartmouth Genomic Data Science Core
# Author: Mike Martinez
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#-----Required libraries
# library(tidyverse)
# library(dplyr)
# library(DESeq2)
# library(RColorBrewer)
# library(pheatmap)
# library(dendextend)
# library(ggplot2)
# library(circlize)
# library(viridis)
# library(scales)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# DENDROGRAM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Exploratory Data Analysis Function: Dendrogram
# VSD is a variance stabilized or rlog DESeq2 transformation
# META is a metadata table
# FEATURE is the column in your metadata to color by within ""
# COLORS is a color palette
# OUTPUT is the output folder to use
dendrogram <- function(VSD, META, FEATURE, COLORS, OUTPUT) {
  
  # Calculate sample distances and perform hierarchical clustering
  sampleDists <- dist(t(assay(VSD)))
  hclustRes <- hclust(sampleDists)
  
  # Create a dendrogram based on the hierarchical clustering
  dend <- as.dendrogram(hclustRes)
  
  # Define a color palette
  color_palette <- COLORS
  
  # Map the metadata to colors
  colors <- color_palette[META[labels(dend), FEATURE]]
  
  # Color the branches
  dend2 <- color_branches(dend, col = colors)
  
  # Output folder
  opFolder <- paste("Outputs", OUTPUT, sep = "/")
  fileName <- paste("HClust", FEATURE, "tiff", sep = ".")
  
  # Initialize a plot
  tiff(paste(opFolder, fileName, sep = "/"), 
       width = 10, height = 10, units = "in", res = 200)

  # Adjust the line thickness
  par(lwd = 5)
  
  # Plot the dendrogram with colored branches
  plot(dend2,
       main = paste("HClust", gsub("_", " ", FEATURE), sep = " "),
       ylab = "Height")
  
  # Add a legend
  legend("topright",  # You can change the position as needed
         legend = names(color_palette),  # Group names
         col = color_palette,  # Corresponding colors
         pch = 15,  # Type of point (15 is a filled square)
         bty = "n")  # No box around the legend
  
  # Close the plot
  dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SAMPLE DISTANCE HEATMAP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Exploratory Data Analysis Function: sample distances
# VSD is a variance stabilized or rlog DESeq2 transformation
# META is a metadata table
# FEATURES is a comma separated vector of META column names to include as 
# annotations
# palettes is a list of color palettes for each feature you are looking at
# OUTPUT is the output folder name encased with ""
# Function to create distance heatmap with annotated colors
distanceHM <- function(VSD, META, FEATURES = c(...), palettes, OUTPUT) {
  
  # Calculate the sample distances
  sampleDists <- dist(t(assay(VSD)))
  
  # Convert distances into matrix
  distMat <- as.matrix(sampleDists)
  
  # Set the diagonal to zero (self:self comparisons)
  diag(distMat) <- 0
  
  # Create a data frame for annotations based on features
  groupdf <- META[, colnames(META) %in% FEATURES]
  groupdf <- as.data.frame(groupdf)
  rownames(groupdf) <- rownames(META)
  
  # Create a color mapping data frame for annotations
  annotation_colors <- list()  # Initialize list to store colors for each feature
  
  for (feature in FEATURES) {
    if (feature %in% names(palettes)) {
      # Get the color mapping for the feature
      colors <- palettes[[feature]]
      # Ensure the feature is treated as a factor
      groupdf[[feature]] <- as.factor(groupdf[[feature]])
      # Create a named vector of colors based on the factor levels
      annotation_colors[[feature]] <- colors[levels(groupdf[[feature]])]
    }
  }
  
  # Set colors for the heatmap
  heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  # Output folder
  opFolder <- paste("Outputs/", OUTPUT, sep = "")
  fileName <- c("Sample_Distance_HM.tiff")
  path <- paste(opFolder, fileName, sep = "/")
  
  print(path)
  
 
  tiff(path, width = 8, height = 8, units = "in", res = 150)
  
  # Visualize the heatmap with colored annotations
  H <- pheatmap(distMat,
                 #clustering_distance_rows = sampleDists,
                 annotation_col = groupdf,
                 annotation_colors = annotation_colors,  # Use the colors for annotations
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 col = heatmap_colors,
                 fontsize = 8,             
                 fontsize_row = 8,        
                 fontsize_col = 8,        
                 legend = TRUE)
  draw(H)

  dev.off()
}

# Example usage
# distanceHM(VSD, META, FEATURES = c("Age", "Sex", "Subtype", "Subtype_Class"), palettes = palette1, "241008")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PCA: DETERMINE VARIABLE FEATURES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Exploratory Data Analysis Function: Principal component analysis
#----- finding suitable number of variable features
# VSD is a variance stabilized or rlog DESeq2 transformation
# Function to plot variances to determine a suitable number of features
determineVarFeatures <- function(VSD) {
  
  # Extract the vsd assay as a matrix
  expMat <- assay(VSD)
  
  # Calculate variance
  var <- apply(expMat, 1, var)
  
  # Sort the variance in decreasing order
  var <- var[order(var, decreasing = TRUE)]
  
  # reset plotting window to 1 row vs 1 columns
  par(mfrow=c(1,1))
  
  # Initialize a plot
  tiff("Variable_Features.tiff", width = 8, height = 8, units ="in", res = 200)
  
  # plot variance for genes accross samples
  varPlot <- plot(var,
                  las = 1,
                  main="Sample gene expression variance",
                  xlab = "Gene", 
                  ylab = "Variance")
  
  # add vertical lines at specific gene number indexes
  abline(v=1000, col="red")
  abline(v=500, col="green")
  abline(v=250, col="blue")

  # Close the graphics device
  dev.off()
  
  # Return the variances
  return(var)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PCA: GENERATE PCS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Exploratory Data Analysis Function: Calculate principal components
# VSD is a variance stabilized or rlog DESeq2 transformation
# NFEATURES is the number of features to include
explorePCA <- function(VSD, VARS, NFEATURES) {
  
  # Set a variable for the number of genes (features) to be used for PCA and clustering
  var_features_n <- NFEATURES
  
  # Order variance and select the rows (genes) with the most variance
  select <- order(VARS, decreasing = TRUE)[1:var_features_n]
  
  # Subset vsd values for genes by top variance ranks
  vsd_sub <- assay(VSD)[select,]
  
  # Transpose the matrix
  vsd_sub <- t(vsd_sub)
  
  # Run principal component analysis
  pca <- prcomp(vsd_sub)
  
  # extract the variance explained by each PC
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  # subset for first 5 elemets
  percentVar <- percentVar[1:5]
  
  # add names to the percentVar vector
  cat("Percent Variations")
  cat("\n")
  percentVar <- paste(round(percentVar*100,2), "%", sep = " ")
  names(percentVar) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
  print(percentVar)
  
  # Construct a data frame with PC loadings and add sample labels
  pca_df <- as.data.frame(pca$x)
  
  return(pca_df)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PCA: EIGENCORRELATIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Function to generate Eigencorrelations and P-values
# VSD is a variance stabilized transformed object
# Metadata is metadata
# num_pcs are the number of pcs to look at
getEigenCorr <- function(vsd, metadata, num_pcs = 10) {
  
  # Convert all metadata categories to numeric
  for (i in colnames(metadata)) {
    if (!is.numeric(metadata[[i]])) {
      metadata[[i]] <- as.numeric(as.factor(metadata[[i]]))
    }
  }
  
  # Perform PCA on the VSD data
  pca_res <- prcomp(t(assay(vsd)))
  
  # Get the principal components (PCs)
  pcs <- as.data.frame(pca_res$x[, 1:num_pcs])  # Use only the top PCs
  
  # Ensure metadata rownames match those of the vsd object
  if (!all(rownames(metadata) == colnames(vsd))) {
    stop("Rownames of metadata must match the colnames of the VSD object.")
  }
  
  # Create empty matrices to store correlations and p-values
  cor_matrix <- matrix(NA, ncol = ncol(pcs), nrow = ncol(metadata))
  pval_matrix <- matrix(NA, ncol = ncol(pcs), nrow = ncol(metadata))
  rownames(cor_matrix) <- rownames(pval_matrix) <- colnames(metadata)
  colnames(cor_matrix) <- colnames(pval_matrix) <- colnames(pcs)
  
  # Loop through metadata variables and calculate correlation for each PC
  for (meta_var in colnames(metadata)) {
    for (pc in colnames(pcs)) {
      if (is.numeric(metadata[[meta_var]])) {
        # For continuous variables, use Pearson correlation and calculate p-values
        test <- cor.test(pcs[[pc]], metadata[[meta_var]], method = "pearson", use = "complete.obs")
        cor_matrix[meta_var, pc] <- test$estimate
        pval_matrix[meta_var, pc] <- test$p.value
      } else if (is.factor(metadata[[meta_var]])) {
        # For categorical variables, use ANOVA to calculate F-statistic and p-value
        aov_res <- aov(pcs[[pc]] ~ metadata[[meta_var]])
        summary_aov <- summary(aov_res)
        cor_matrix[meta_var, pc] <- summary_aov[[1]][["F value"]][1]  # F-statistic
        pval_matrix[meta_var, pc] <- summary_aov[[1]][["Pr(>F)"]][1]  # p-value
      }
    }
  }
  
  # Remove NA
  cor_matrix <- na.omit(cor_matrix)
  pval_matrix <- na.omit(pval_matrix)
  
  # Create significance stars based on p-value thresholds
  stars <- ifelse(pval_matrix < 0.001, "***", 
                  ifelse(pval_matrix < 0.01, "**", 
                         ifelse(pval_matrix < 0.05, "*", "")))
  
  # Set colors for the heatmap
  heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
  
  
  # Return list
  return_list <- list(cor_matrix = cor_matrix,
                      pval_matrix = pval_matrix,
                      stars = stars
  )
  
  return(return_list)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# EIGENCORRELATION HEATMAP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Function to plot Eigencorrelation heatmap
# corrs is the "cor_matrix" slot in the return list from getEigenCorr()
# stars is the "stars" slot in the return list from getEigenCorr()
# Output is the output folder name to save heatmap encased in ""
plotEigenCorrs <- function(corrs, stars, output) {
  
  # Generate output file path
  fileName <- c("EigenCorrelations.tiff")
  opPath <- paste("Outputs", output, sep = "/")
  path <- paste(opPath, fileName, sep = "/")
  
  # Initialize image
  tiff(path, width = 10, height = 10, units = "in", res = 200)
  
  # Plot heatmap
  H <- Heatmap(corrs,
               cluster_rows = TRUE,
               cluster_columns = FALSE,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%.2f", corrs[i, j]), x, y, gp = gpar(fontsize = 10))
                 if (stars[i, j] != "") {
                   grid.text(stars[i, j], x, y + unit(5, "mm"), gp = gpar(col = "black", fontsize = 12))
                 }
               },
               name = "Correlation",
               show_row_names = TRUE,
               show_column_names = TRUE,
               column_title = "Principal Components",
               row_title = "Metadata Variables")
  draw(H)
  dev.off()
  return(H)
  
}





