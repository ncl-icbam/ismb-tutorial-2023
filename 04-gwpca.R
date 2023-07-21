## ----04_loadPackages, include=FALSE-----------------------------------------------------------------
library(SpatialFeatureExperiment)
library(tidyverse)
library(scran)
library(scater)
library(sf)
library(spdep)
library(GWmodel)
library(future)
library(doFuture)
library(foreach)
library(progressr)
library(parallel)
library(cols4all)
library(pheatmap)
library(RColorBrewer)
library(STExplorerDev)


## ----04_load_data-----------------------------------------------------------------------------------
sfe <- readRDS(file = "./data/to_load/practical03_sfe.rds")
top_hvgs <- readRDS(file = "./data/to_load/practical03_topHVGs.rds")


## ----04_set_parameters, eval=FALSE------------------------------------------------------------------
## ## Get the gene names that are going to be evaluated
## vars = top_hvgs
## ## Set a fixed bandwidth
## bw = 6*sfe@metadata[["spotDiameter"]][["JBO019"]][["spot_diameter_fullres"]]
## ## Set the number of components to be retained
## k = 20
## ## Set the kernel to be used
## kernel = "gaussian"
## ## Set the Minkowski distance power: p = 2 --> Euclidean
## p = 2
## ## Is the bandwidth adaptive?: No because spots are fixed
## adaptive = FALSE
## ## Cross-Validate GWPCA?
## cv = TRUE
## ## Calculate PCA scores?
## scores = FALSE
## ## Run a robust GWPCA?
## robust = FALSE
## ## Make a cluster for parallel computing (otherwise GWPCA is slow!)
## my.cl <- parallel::makeCluster(parallelly::availableCores() - 1, type = 'FORK')


## ----04_run_gwpca1, eval=FALSE----------------------------------------------------------------------
## # DO NOT RUN THIS CHUNK
## 
## pcagw <- gwpcaSTE(sfe = sfe,
##                   assay = "logcounts",
##                   vars = vars,
##                   p = p,
##                   k = k,
##                   bw = bw,
##                   kernel = kernel,
##                   adaptive = adaptive,
##                   scores = scores,
##                   robust = robust,
##                   cv = cv,
##                   future = FALSE,
##                   strategy = "cluster",
##                   workers = my.cl,
##                   verbose = FALSE)
## saveRDS(pcagw, file = "./data/to_load/practical04_pcagw.rds")


## ----04_run_gwpca2, eval=TRUE-----------------------------------------------------------------------
pcagw <- readRDS(file = "./data/to_load/practical04_pcagw.rds")


## ----04_scree_plot, eval=TRUE, fig.height=3, fig.width=8--------------------------------------------
plotGWPCA_global(gwpca = pcagw,
                 comps = 1:10,
                 type = "scree",
                 point_args = list(size = 3, colour = "red"),
                 line_args = list(linewidth = 1, colour = "dodgerblue"))


## ----leading_genes1, eval=TRUE----------------------------------------------------------------------
## Extract leading genes
pcagw <- gwpca_LeadingGene(gwpca = pcagw, 
                           sfe = sfe, 
                           pc_nos = 1:4, 
                           type = "single", 
                           names = "gene_names")

pcagw <- gwpca_LeadingGene(gwpca = pcagw, 
                           sfe = sfe, 
                           pc_nos = 1:4, 
                           genes_n = 4, 
                           type = "multi", 
                           method = "membership", 
                           names = "gene_names")


## ----leading_genes2, eval=TRUE, fig.show = 'hold', out.width='.49\\linewidth', fig.asp=1, fig.ncol = 1----
## Plot leading genes
plotGWPCA_leadingG(gwpca = pcagw,
                   comps = 1:2,
                   type = "single",
                   arrange = FALSE)

plotGWPCA_leadingG(gwpca = pcagw,
                   comps = 1,
                   type = "multi",
                   arrange = FALSE)


## ----leading_genes3, eval=TRUE, fig.show = 'hold', out.width='.49\\linewidth', fig.asp=1, fig.ncol = 1----
### Plot multi type (alternative)
## The data
leadingGsMulti <- pcagw$leadingGeneMulti
## The Legend labels
spot_labels <- data.frame(table(leadingGsMulti[1])) %>%
    dplyr::rename(LeadingGs = colnames(leadingGsMulti)[1], 
                  count = Freq) %>%
    dplyr::arrange(desc(count)) %>%
    mutate(show = ifelse(count > 12, TRUE, FALSE))
    
## The legend breaks:
spot_breaks <- spot_labels %>%
    dplyr::filter(show == TRUE) %>% 
    dplyr::arrange(LeadingGs) %>%
    dplyr::select(LeadingGs) %>% 
    .[["LeadingGs"]] %>% 
    as.vector()
    
## The colours:
col_No <- sum(spot_labels$show)
colour_values <- getColours(col_No)
names(colour_values) <- spot_labels$LeadingGs[spot_labels$show]
pc <- "PC1"
    
## The Plot:
ggplot() + 
    geom_sf(data = leadingGsMulti, 
            aes(geometry = geometry$geometry,
                fill = .data[[pc]]),
            colour = "grey30", 
            show.legend = TRUE) + 
    scale_fill_manual(values = colour_values,
                      breaks = spot_breaks,
                      na.value = "gray95") +
    labs(title = NULL,
         fill = "Group of\nLeading\nGenes") + 
    theme_void() +
    theme(legend.position = "bottom", legend.text = element_text(size=6)) +
    guides(fill = guide_legend(ncol = 3, byrow = TRUE))


## ----04_ptv, eval=TRUE, fig.show='hold'-------------------------------------------------------------
## Calculate the PTV for multiple Components
pcagw <- gwpca_PropVar(gwpca = pcagw, n_comp = 2:10, sfe = sfe)

## Plot PTV
plotGWPCA_ptv(gwpca = pcagw,
              comps = 1:10,
              type = "violin")

## Map PTV
plotGWPCA_ptv(gwpca = pcagw,
              comps = 1:6,
              type = "map")


## ----04_discrep1, eval=TRUE, fig.height=3, fig.width=8----------------------------------------------
## Plot the discrepancies as boxplot
plotGWPCA_discr(pcagw, type = "box")


## ----04_discrep2, eval=TRUE-------------------------------------------------------------------------
## Plot the discrepancies map
plotGWPCA_discr(pcagw, type = "map")


## ----04_discrep3, eval=TRUE-------------------------------------------------------------------------
## Get location data for the discrepancies
discrepancy_loc_dt <- getDiscrepancyLocData(sfe = sfe, 
                                            gwpca = pcagw, 
                                            sample_id = "JBO019")


## ----04_discrep4, eval=TRUE, message=FALSE, fig.show='hold', fig.height=15, fig.width=9-------------
head(discrepancy_loc_dt)
focus <- discrepancy_loc_dt$barcodes[1:2]
bw = 3*sfe@metadata[["spotDiameter"]][["JBO019"]][["spot_diameter_fullres"]]

# Plot the heatmap to visualise the genes that make this location an outlier
plotGWPCA_discrHeatmap(sfe = sfe,
                       assay = "logcounts",
                       vars = NULL,
                       focus = focus,
                       dMetric = "euclidean", 
                       sample_id = "JBO019",
                       bw = bw, 
                       mean.diff = 1, 
                       show.vars = "top", 
                       scale = "row", 
                       gene.names = TRUE,
                       color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(1000)),
                       fontsize_row = 3)


## ----04_discrep5, message=FALSE---------------------------------------------------------------------
discrepancy_gene_dt <- getDiscrepancyGeneData(sfe = sfe,
                                              assay = "logcounts",
                                              vars = NULL,
                                              focus = focus[2],
                                              dMetric = "euclidean", 
                                              sample_id = "JBO019",
                                              bw = bw, 
                                              mean.diff = 1, 
                                              show.vars = "top",
                                              exportExpression = TRUE)
head(discrepancy_gene_dt)

