### ----- WGCNA Analysis -----
# Reference: https://www.jianshu.com/p/e9cc3f43441d
# Install required packages (commented for sharing)
# BiocManager::install("WGCNA")
# BiocManager::install('impute')
# BiocManager::install('preprocessCore')

library(WGCNA)
library(stringr)
library(reshape2)

# Set network type and correlation method
type <- "signed"
corType <- "bicor"  # Recommended method
corFnc <- ifelse(corType == "pearson", cor, bicor)
maxPOutliers <- ifelse(corType == "pearson", 1, 0.05)
robustY <- ifelse(corType == "pearson", TRUE, FALSE)

# Load expression data (genes × samples)
dataExpr <- read.csv('bulk.csv', header = TRUE, row.names = 1)

# Filter genes by standard deviation (threshold = 0.5)
threshold <- 0.5
gene_sds <- apply(dataExpr, 1, sd)
filtered_indices <- which(gene_sds >= threshold)
dataExprVar <- dataExpr[filtered_indices, ]

# Transpose to samples × genes matrix
dataExpr <- as.data.frame(t(dataExprVar))
nSamples <- nrow(dataExpr)

# Sample clustering to detect outliers
sampleTree <- hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub = "", xlab = "")

# Determine soft threshold power
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
sft <- pickSoftThreshold(dataExpr, powerVector = powers, 
                         networkType = type, verbose = 5)

# Plot scale independence and mean connectivity
par(mfrow = c(1, 2))
cex1 <- 0.9

plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R²",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red")  # Cutoff line

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     labels = powers, cex = cex1, col = "red")

power <- sft$powerEstimate
power

# Use empirical power if automatic selection fails
if (is.na(power)) {
  power <- ifelse(nSamples < 20, ifelse(type == "unsigned", 9, 18),
                  ifelse(nSamples < 30, ifelse(type == "unsigned", 8, 16),
                         ifelse(nSamples < 40, ifelse(type == "unsigned", 7, 14),
                                ifelse(type == "unsigned", 6, 12))))
}

# One-step network construction and module detection
exprMat <- "myExpressionData"
net <- blockwiseModules(dataExpr, 
                        power = 10,  # Adjust based on soft threshold analysis
                        maxBlockSize = ncol(dataExpr),
                        TOMType = "signed",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        corType = corType,
                        maxPOutliers = maxPOutliers,
                        loadTOMs = TRUE,
                        saveTOMFileBase = paste0(exprMat, ".tom"),
                        verbose = 3)

# Module summary
table(net$colors)

# Convert labels to colors
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

# Plot dendrogram with module colors
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

# Module eigengenes
MEs <- net$MEs
MEs_col <- MEs
colnames(MEs_col) <- paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs), "ME", ""))))
MEs_col <- orderMEs(MEs_col)

# Eigengene network heatmap
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3, 3, 2, 4),
                      marHeatmap = c(3, 4, 2, 2),
                      plotDendrograms = TRUE,
                      xLabelsAngle = 90)

# Load TOM matrix
load(net$TOMFiles[1], verbose = TRUE)
TOM <- as.matrix(TOM)

# Prepare trait data
traitData <- read.csv('batch.csv', header = TRUE, row.names = 1)
traitData <- data.frame(lactation = traitData[, "lactation"], 
                        row.names = rownames(traitData))
sampleName <- rownames(dataExpr)
traitData <- traitData[match(sampleName, rownames(traitData)), , drop = FALSE]

# Convert to binary variables
traitData$lactation_bin <- ifelse(traitData$lactation == "lactating", 1, 0)
traitData$dry_bin <- ifelse(traitData$lactation == "dry", 1, 0)
traitData <- traitData[, c("lactation_bin", "dry_bin")]

# Calculate module-trait correlations
if (corType == "pearson") {
  modTraitCor <- cor(MEs_col, traitData, use = "p")
  modTraitP <- corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP <- bicorAndPvalue(MEs_col, traitData, robustY = robustY)
  modTraitCor <- modTraitCorP$bicor
  modTraitP <- modTraitCorP$p
}

# Create labeled heatmap
textMatrix <- paste(signif(modTraitCor, 2), "\n(", 
                    signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) <- dim(modTraitCor)
dark_blue_to_red <- colorRampPalette(
  c("#0099CC", "#99CCFF", "#FFFFFF", "#FFCCCC", "#CC0000"))(50)

labeledHeatmap(Matrix = modTraitCor,
               xLabels = colnames(traitData),
               yLabels = colnames(MEs_col),
               cex.lab = 0.5,
               ySymbols = colnames(MEs_col),
               colorLabels = FALSE,
               colors = dark_blue_to_red,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-trait relationships")

# Calculate gene-module membership
if (corType == "pearson") {
  geneModuleMembership <- as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA <- bicorAndPvalue(dataExpr, MEs_col, robustY = robustY)
  geneModuleMembership <- geneModuleMembershipA$bicor
  MMPvalue <- geneModuleMembershipA$p
}

# Calculate gene-trait significance
if (corType == "pearson") {
  geneTraitCor <- as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP <- as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA <- bicorAndPvalue(dataExpr, traitData, robustY = robustY)
  geneTraitCor <- as.data.frame(geneTraitCorA$bicor)
  geneTraitP <- as.data.frame(geneTraitCorA$p)
}

# Analyze specific modules
probes <- colnames(dataExpr)
modNames <- substring(colnames(MEs_col), 3)

# Function to extract hub genes from a module
extractHubGenes <- function(module, pheno, traitData, geneModuleMembership, 
                            geneTraitCor, moduleColors, probes, MMs_col) {
  module_column <- match(module, modNames)
  pheno_column <- match(pheno, colnames(traitData))
  moduleGenes <- moduleColors == module
  modProbes <- probes[moduleGenes]
  
  # Scatter plot
  sizeGrWindow(7, 7)
  verboseScatterplot(
    abs(geneModuleMembership[moduleGenes, module_column]),
    abs(geneTraitCor[moduleGenes, pheno_column]),
    xlab = paste("Module Membership in", module, "module"),
    ylab = paste("Gene significance for", pheno),
    main = paste("Module membership vs. gene significance\n"),
    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module
  )
  
  # Create result dataframe
  MM <- abs(geneModuleMembership[moduleGenes, module_column])
  GS <- abs(geneTraitCor[moduleGenes, pheno_column])
  key_MMGS <- data.frame(MM = MM, GS = GS, row.names = modProbes)
  key_MMGS$hub <- abs(key_MMGS$MM) > 0.8 & abs(key_MMGS$GS) > 0.2
  key_MMGS$module <- module
  key_MMGS$pheno <- pheno
  
  return(key_MMGS)
}

# Analyze modules of interest
turquoise_result <- extractHubGenes("turquoise", "lactation_bin", 
                                    traitData, geneModuleMembership, 
                                    geneTraitCor, moduleColors, probes, MEs_col)

tan_result <- extractHubGenes("tan", "lactation_bin", 
                              traitData, geneModuleMembership, 
                              geneTraitCor, moduleColors, probes, MEs_col)

black_result <- extractHubGenes("black", "dry_bin", 
                                traitData, geneModuleMembership, 
                                geneTraitCor, moduleColors, probes, MEs_col)

# Combine results
all_results <- rbind(turquoise_result, tan_result, black_result)
hub_genes <- subset(all_results, hub == TRUE)

# Save results
write.csv(all_results, "wgcna_keymodule.csv")
write.csv(hub_genes, "wgcna_hub_genes.csv")