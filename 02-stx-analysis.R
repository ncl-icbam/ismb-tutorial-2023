## ----02_loadPackages, message=FALSE-----------------------------------------------------------------
## Load packages {-}
library(SpatialExperiment)
library(STexampleData)
library(ggspavis)
library(ggplot2)
library(scater)
library(scran)
library(igraph)
library(pheatmap)
library(ggExtra)


## ----02_reloadData----------------------------------------------------------------------------------
## Reload the example dataset
spe <- Visium_humanDLPFC()


## ----02_plot-maps-gTruth, fig.show = 'hold', out.width="50%", fig.height=5, fig.width=4-------------
## Plot spatial coordinates without annotations
plotSpots(spe)

## Plot spatial coordinates with annotations
plotSpots(spe,
          annotate = "ground_truth")


## ----02_keep_on-tissue------------------------------------------------------------------------------
## Dataset dimensions before the filtering
dim(spe)

## Subset to keep only on-tissue spots
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)


## ----02_find-mitoGenes------------------------------------------------------------------------------
## Classify genes as "mitochondrial" (is_mito == TRUE) 
## or not (is_mito == FALSE)
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
rowData(spe)$gene_name[is_mito]


## ----02_perCellQCs----------------------------------------------------------------------------------
## Calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
head(colData(spe))


## ----02_plot-libSize-histo, fig.height=3.5, warning=FALSE, message=FALSE----------------------------
## Density and histogram of library sizes
ggplot(data = as.data.frame(colData(spe)),
       aes(x = sum)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Library size") + 
  ylab("Density") + 
  theme_classic()


## ----02_plot-libSizeVScelNo, fig.width=6, fig.height=5, warning=FALSE, message=FALSE----------------
## Scatter plot, library size against number of cells per spot
plotQC(spe, type = "scatter", 
       metric_x = "cell_count", metric_y = "sum", 
       threshold_y = 700)


## ----02_ggplot-libSizeVScelNo, fig.width=6, fig.height=5, warning=FALSE, message=FALSE--------------
p = ggplot(as.data.frame(colData(spe)), aes(x=cell_count, y=sum)) +
  geom_point(size=0.5) + 
  geom_smooth(se=FALSE) +
  geom_hline(yintercept = 700, colour='red') + 
  theme_minimal()
ggMarginal(p, type='histogram', margins = 'both')


## ----02_libSize-thresh, fig.height=4----------------------------------------------------------------
## Select library size threshold
qc_lib_size <- colData(spe)$sum < 700
## Check how many spots are filtered out
table(qc_lib_size)
## Add threshold in colData
colData(spe)$qc_lib_size <- qc_lib_size

## Check putative spatial patterns of removed spots
plotQC(spe, type = "spots", 
       discard = "qc_lib_size")


## ----02_exercise01, fig.height=4, eval=FALSE--------------------------------------------------------
## ## Select library size threshold
## code...
## ## Check how many spots are filtered out
## code...
## ## Add threshold in colData
## code...
## 
## ## Check putative spatial patterns of removed spots
## plotQC(...)


## ----02_plot-genesInSpot-histo, fig.height=4, warning=FALSE, message=FALSE--------------------------
## Density and histogram of expressed genes
ggplot(data = as.data.frame(colData(spe)),
       aes(x = detected)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Genes expressed in each spot") + 
  ylab("Density") + 
  theme_classic()


## ----02_genesInSpot-scatter, fig.width=6, fig.height=5, warning=FALSE, message=FALSE----------------
# plot number of expressed genes vs. number of cells per spot
p = ggplot(as.data.frame(colData(spe)), aes(x=cell_count, y=detected)) +
  geom_point(size=0.5) + 
  geom_smooth(se=FALSE) +
  geom_hline(yintercept = 500, colour='red') + 
  theme_minimal()
ggMarginal(p, type='histogram', margins = 'both')


## ----02_genesInSpot-thresh, fig.height=4------------------------------------------------------------
## Select expressed genes threshold
qc_detected <- colData(spe)$detected < 500
## Check how many spots are filtered out
table(qc_detected)
## Add threshold in colData
colData(spe)$qc_detected <- qc_detected

## Check for putative spatial pattern of removed spots
plotQC(spe, type = "spots", 
       discard = "qc_detected")


## ----02_exercise02, fig.height=4, eval=FALSE--------------------------------------------------------
## ## Select library size threshold
## code...
## ## Check how many spots are filtered out
## code...
## ## Add threshold in colData
## code...
## 
## ## Check putative spatial patterns of removed spots
## plotQC(...)


## ----02_plot-mitoPercent-histo, fig.height=4, warning=FALSE, message=FALSE--------------------------
## Density and histogram of percentage of mitochondrial expression
ggplot(data = as.data.frame(colData(spe)),
       aes(x = subsets_mito_percent)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Percentage of mitochondrial expression") + 
  ylab("Density") + 
  theme_classic()


## ----02_mitoPercent-scatter, fig.width=6, fig.height=5, warning=FALSE, message=FALSE----------------
# plot mitochondrial read proportion vs. number of cells per spot
p = ggplot(as.data.frame(colData(spe)), aes(x=cell_count, y=subsets_mito_percent)) +
  geom_point(size=0.5) + 
  geom_smooth(se=FALSE) +
  geom_hline(yintercept = 28, colour='red') + 
  theme_minimal()
ggMarginal(p, type='histogram')



## ----02_mitoPercent-thresh, fig.height=4------------------------------------------------------------
## Select expressed genes threshold
qc_mito <- colData(spe)$subsets_mito_percent > 28
## Check how many spots are filtered out
table(qc_mito)
## Add threshold in colData
colData(spe)$qc_mito <- qc_mito

## Check for putative spatial pattern of removed spots
plotQC(spe, type = "spots", 
       discard = "qc_mito")


## ----02_exercise03, fig.height=4, eval=FALSE--------------------------------------------------------
## ## Select library size threshold
## code...
## ## Check how many spots are filtered out
## code...
## ## Add threshold in colData
## code...
## 
## ## Check putative spatial patterns of removed spots
## plotQC(...)


## ----02_plot-cellsPerSpot-histo, fig.height=4, warning=FALSE, message=FALSE-------------------------
## Density and histogram of the number of cells in each spot
ggplot(data = as.data.frame(colData(spe)),
       aes(x = cell_count)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 1,
                 colour = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.5,
               adjust = 1.5,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Number of cells per spot") + 
  ylab("Density") + 
  theme_classic()

## Have a look at the values
table(colData(spe)$cell_count)


## ----02_cellsPerSpot-scatter, fig.width=6, fig.height=5, warning=FALSE, message=FALSE---------------
# plot number of expressed genes vs. number of cells per spot
p = ggplot(as.data.frame(colData(spe)), aes(x=cell_count, y=detected)) +
  geom_point() + 
  geom_smooth(se=FALSE) +
  geom_vline(xintercept = 10, colour='red') + 
  theme_minimal()
ggMarginal(p, type='histogram')


## ----02_cellsPerSpot-thresh, fig.height=4-----------------------------------------------------------
## Select expressed genes threshold
qc_cell_count <- colData(spe)$cell_count > 10
## Check how many spots are filtered out
table(qc_cell_count)
## Add threshold in colData
colData(spe)$qc_cell_count <- qc_cell_count

## Check for putative spatial pattern of removed spots
plotQC(spe, type = "spots", 
       discard = "qc_cell_count")


## ----02_checkQC-thresh, fig.height=4----------------------------------------------------------------
## Check the number of discarded spots for each metric
apply(cbind(qc_lib_size, qc_detected, qc_mito, qc_cell_count), 2, sum)
## Combine together the set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
## Store the set in the object
colData(spe)$discard <- discard

## Check the spatial pattern of combined set of discarded spots
plotQC(spe, type = "spots", 
       discard = "discard")


## ----02_notAnnotSpots, fig.height=4-----------------------------------------------------------------
## Select locations without annotation
qc_NA_spots <- is.na(colData(spe)$ground_truth)
## Combine together the set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count | qc_NA_spots
## Store the set in the object
colData(spe)$discard <- discard

## Check the spatial pattern of combined set of discarded spots
plotQC(spe, type = "spots", 
       discard = "discard")


## ----02_applyFilter---------------------------------------------------------------------------------
## remove combined set of low-quality spots
spe <- spe[, !colData(spe)$discard]


## ----02_libraryFactors, message=FALSE---------------------------------------------------------------
## Calculate library size factors
spe <- computeLibraryFactors(spe)
## Have a look at the size factors
summary(sizeFactors(spe))


## ----02_plot-labfact-histo, fig.height=4, warning=FALSE, message=FALSE------------------------------
## Density and histogram of library sizes
ggplot(data = data.frame(sFact = sizeFactors(spe)), 
       aes(x = sFact)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey") +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Library size") + 
  ylab("Density") + 
  theme_classic()


## ----02_logNormCounts-------------------------------------------------------------------------------
## Calculate logcounts and store in the spe object
spe <- logNormCounts(spe)

## Check that a new assay has been added
assayNames(spe)


## ----02_features_remvMito, message=FALSE------------------------------------------------------------
## Remove mitochondrial genes
spe <- spe[!is_mito, ]


## ----02_features_FitModel, message=FALSE, fig.height=5----------------------------------------------
## Fit mean-variance relationship
dec <- modelGeneVar(spe)
## Visualize mean-variance relationship
fit <- metadata(dec)
fit_df <- data.frame(mean = fit$mean,
                     var = fit$var,
                     trend = fit$trend(fit$mean))

ggplot(data = fit_df, 
       aes(x = mean, y = var)) + 
  geom_point() + 
  geom_line(aes(y = trend), colour = "dodgerblue", linewidth = 1.5) + 
  labs(x = "mean of log-expression",
       y = "variance of log-expression") + 
  theme_classic()



## ----02_features_selectHVGs-------------------------------------------------------------------------
## Select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)

## How many HVGs?
length(top_hvgs)

## Plot the HVGs on the mean:variance trend
fit_df <- fit_df %>%
    tibble::rownames_to_column(var = "row.names") %>%
    dplyr::mutate(topHVGs = ifelse(row.names %in% top_hvgs, TRUE, FALSE)) %>%
    tibble::column_to_rownames("row.names")

ggplot(data = fit_df,
       aes(x = mean, y = var, colour = topHVGs)) + 
    geom_point() +
    geom_line(aes(y = trend), colour = "dodgerblue", linewidth = 1.5) + 
    scale_colour_manual(values = c("black", "red")) + 
    labs(x = "Mean of log-expression",
         y = "Variance of log-expression",
         colour = "Top HVGs") + 
    theme_classic()


## ----02_dimRedct_PCA, message=FALSE, warning=FALSE--------------------------------------------------
## Set seed
set.seed(987)
## Compute PCA
spe <- runPCA(spe, subset_row = top_hvgs)
## Check correctness - names
reducedDimNames(spe)
## Check correctness - dimensions
dim(reducedDim(spe, "PCA"))


## ----02_dimRed_UMAP---------------------------------------------------------------------------------
## Set seed
set.seed(987)
## Compute UMAP on top 50 PCs
spe <- runUMAP(spe, dimred = "PCA")
## Check correctness - names
reducedDimNames(spe)
## Check correctness - dimensions
dim(reducedDim(spe, "UMAP"))
## Update column names for easier plotting
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)


## ----02_dimRed_UMAP-vis, fig.width=6, fig.height=5, fig.show='hold'---------------------------------
## Plot top 2 PCA dimensions
# plotDimRed(spe, type = "PCA")

ggplot(data = as.data.frame(spe@int_colData@listData$reducedDims$PCA),
       aes(x = PC1, y = PC2, colour = spe@colData$ground_truth)) + 
  geom_point(size = 0.5) + 
  scale_colour_brewer(type = "qual") + 
  labs(title = "Reduced dimensions: PCA",
       x = "PC1",
       y = "PC2",
       colour = "Layers") +
  theme_classic()

## Plot top 2 UMAP dimensions
# plotDimRed(spe, type = "UMAP")

ggplot(data = as.data.frame(spe@int_colData@listData$reducedDims$UMAP),
       aes(x = UMAP1, y = UMAP2, colour = spe@colData$ground_truth)) + 
  geom_point(size = 0.5) + 
  scale_colour_brewer(type = "qual") + 
  labs(title = "Reduced dimensions: UMAP",
       x = "UMAP1",
       y = "UMAP2",
       colour = "Layers") +
  theme_classic()


## ----02_clustering----------------------------------------------------------------------------------
## Set seed
set.seed(987)
## Set number of Nearest-Neighbours (NNs)
k <- 10
## Build the k-NN graph
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
## Run walktrap clustering
g_walk <- igraph::cluster_walktrap(g)
## Get the cluster labels
clus <- g_walk$membership
## Check how many
table(clus)
## Store cluster labels in column 'label' in colData
colLabels(spe) <- factor(clus)


## ----02_clust_vis-map, fig.show = 'hold', out.width="50%", fig.height=5, fig.width=4----------------
## Plot in tissue map
plotSpots(spe, annotate = "label", 
          palette = "libd_layer_colors")

## Plot ground truth in tissue map
plotSpots(spe, annotate = "ground_truth", 
          palette = "libd_layer_colors")


## ----02_clust_vis-DimRed, message=FALSE, fig.show = 'hold', out.width="50%", fig.height=4, fig.width=4----
## Plot clusters in PCA space
plotDimRed(spe, type = "PCA", 
           annotate = "label", palette = "libd_layer_colors")

## Plot clusters in UMAP space
plotDimRed(spe, type = "UMAP", 
           annotate = "label", palette = "libd_layer_colors")


## ----02_dges----------------------------------------------------------------------------------------
## Set gene names as row names ease of plotting
rownames(spe) <- rowData(spe)$gene_name
## Test for DGEs
markers <- findMarkers(spe, test = "binom", direction = "up")
## Check output
markers


## ----02_dges_vis-clst1, message=FALSE, fig.show = 'hold', out.width="50%", fig.height=4, fig.width=4----
## Select cluster 1 genes
interesting <- markers[[1]]
## Get the top genes
best_set <- interesting[interesting$Top <= 5, ]
## Calculate the effect
logFCs <- getMarkerEffects(best_set)
## Plot a heat map
pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101))


## ----02_dges_vis-2, message=FALSE, fig.width=7, fig.height=7----------------------------------------
## Select genes
top_genes <- head(rownames(interesting))
## Plot expression
plotExpression(spe, x = "label", features = top_genes)


## ----all_together, eval=FALSE-----------------------------------------------------------------------
## # clear workspace from previous chapters
## rm(list = ls(all = TRUE))
## 
## # LOAD DATA
## 
## library(SpatialExperiment)
## library(STexampleData)
## spe <- Visium_humanDLPFC()
## 
## # QUALITY CONTROL (QC)
## 
## library(scater)
## # subset to keep only spots over tissue
## spe <- spe[, colData(spe)$in_tissue == 1]
## # identify mitochondrial genes
## is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
## # calculate per-spot QC metrics
## spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
## # select QC thresholds
## qc_lib_size <- colData(spe)$sum < 600
## qc_detected <- colData(spe)$detected < 400
## qc_mito <- colData(spe)$subsets_mito_percent > 28
## qc_cell_count <- colData(spe)$cell_count > 10
## # combined set of discarded spots
## discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
## colData(spe)$discard <- discard
## # filter low-quality spots
## spe <- spe[, !colData(spe)$discard]
## 
## # NORMALIZATION
## 
## library(scran)
## # calculate logcounts using library size factors
## spe <- logNormCounts(spe)
## 
## # FEATURE SELECTION
## 
## # remove mitochondrial genes
## spe <- spe[!is_mito, ]
## # fit mean-variance relationship
## dec <- modelGeneVar(spe)
## # select top HVGs
## top_hvgs <- getTopHVGs(dec, prop = 0.1)
## 
## # DIMENSIONALITY REDUCTION
## 
## # compute PCA
## set.seed(123)
## spe <- runPCA(spe, subset_row = top_hvgs)
## # compute UMAP on top 50 PCs
## set.seed(123)
## spe <- runUMAP(spe, dimred = "PCA")
## # update column names
## colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)
## 
## # CLUSTERING
## 
## # graph-based clustering
## set.seed(123)
## k <- 10
## g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
## g_walk <- igraph::cluster_walktrap(g)
## clus <- g_walk$membership
## colLabels(spe) <- factor(clus)
## 
## # MARKER GENES
## # test for marker genes
## rownames(spe) <- rowData(spe)$gene_name
## markers <- findMarkers(spe, test = "binom", direction = "up")

