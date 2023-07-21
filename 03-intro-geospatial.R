## ----03_loadPackages, include=FALSE, eval=TRUE------------------------------------------------------
library(SpatialFeatureExperiment)
library(tidyverse)
library(scran)
library(scater)
library(ggspavis)
library(sf)
library(spdep)
library(GWmodel)
library(tidyterra)
# STExplorer is our own, under-development package for STx analysis
library(STExplorerDev)


## ----03_sf_LoadTest---------------------------------------------------------------------------------
nc <- st_read(system.file("shape/nc.shp", package = "sf"))


## ----03_sf_TestClass--------------------------------------------------------------------------------
class(nc)


## ----03_sf_Test_sf_column---------------------------------------------------------------------------
attr(nc, "sf_column")


## ----03_sf_Test_print, echo=TRUE, eval=FALSE--------------------------------------------------------
## print(nc[9:15], n = 3)


## ----Sf-overview, echo=FALSE, out.width = "100%", fig.align="center", fig.cap="Overview of the `sf` object."----
knitr::include_graphics("images/sf_xfig.png")


## ----03_sf_no.sf------------------------------------------------------------------------------------
nc.no_sf <- as.data.frame(nc)
class(nc.no_sf)


## ----03_sfc_Test1-----------------------------------------------------------------------------------
(nc_geom <- st_geometry(nc))


## ----03_sfc_Test2-----------------------------------------------------------------------------------
nc_geom[[1]]


## ----03_sfc_Test3, fig.height=3---------------------------------------------------------------------
par(mar = c(0,0,1,0))
plot(nc[1], reset = FALSE) # reset = FALSE: we want to add to a plot with a legend
plot(nc[1,1], col = 'grey', add = TRUE)


## ----03_sfc_Test4-----------------------------------------------------------------------------------
nc_geom[[4]][[2]][[1]][1:3,]


## ----03_sfc_Test5-----------------------------------------------------------------------------------
class(nc_geom)


## ----03_sf_Test6, echo=FALSE, eval=TRUE, message=FALSE----------------------------------------------
p <- rbind(c(3.2,4), c(3,4.6), c(3.8,4.4), c(3.5,3.8), c(3.4,3.6), c(3.9,4.5))
(mp <- st_multipoint(p))
s1 <- rbind(c(0,3),c(0,4),c(1,5),c(2,5))
(ls <- st_linestring(s1))
s2 <- rbind(c(0.2,3), c(0.2,4), c(1,4.8), c(2,4.8))
s3 <- rbind(c(0,4.4), c(0.6,5))
(mls <- st_multilinestring(list(s1,s2,s3)))
p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
p2 <- rbind(c(1,1), c(1,2), c(2,2), c(1,1))
pol <-st_polygon(list(p1,p2))
p3 <- rbind(c(3,0), c(4,0), c(4,1), c(3,1), c(3,0))
p4 <- rbind(c(3.3,0.3), c(3.8,0.3), c(3.8,0.8), c(3.3,0.8), c(3.3,0.3))[5:1,]
p5 <- rbind(c(3,3), c(4,2), c(4,3), c(3,3))
(mpol <- st_multipolygon(list(list(p1,p2), list(p3,p4), list(p5))))
(gc <- st_geometrycollection(list(mp, mpol, ls)))


## ----03_sf_Test7, echo=FALSE, eval=TRUE-------------------------------------------------------------
par(mar = c(0.1, 0.1, 1.3, 0.1), mfrow = c(2, 3))
plot(mp, col = 'red')
box()
title("MULTIPOINT")
plot(ls, col = 'red')
box()
title("LINESTRING")
plot(mls, col = 'red')
box()
title("MULTILINESTRING")
plot(pol, border = 'red', col = 'grey', xlim = c(0,4))
box()
title("POLYGON")
plot(mpol, border = 'red', col = 'grey')
box()
title("MULTIPOLYGON")
plot(gc, border = 'grey', col = 'grey')
box()
title("GEOMETRYCOLLECTION")
par(mfrow = c(1, 1))


## ----03_sf_Test8, collapse=TRUE---------------------------------------------------------------------
(x <- st_geometrycollection())
length(x)


## ----03_load_sfe, warning=FALSE, message=FALSE------------------------------------------------------
sampleDir <- "./data/spaceranger_outs/Human_Liver_Steatotic/JBO019_Results"
sampleNames <- "JBO019"
sfe <- read10xVisiumSFE(samples = sampleDir, 
                        sample_id = sampleNames, 
                        type = "sparse", 
                        data = "filtered", 
                        images = "lowres", 
                        style = "W", 
                        zero.policy = TRUE)

ground_truth <- read_table("./data/to_load/spotzonationGroup.txt")



## ----03_QC_sfe1, message=FALSE, warning=FALSE-------------------------------------------------------
is_mito <- grepl("(^MT-)|(^mt-)", rowData(sfe)$symbol)
sfe <- addPerLocQC(sfe, gTruth = ground_truth, assay = "counts", 2, subsets = list(mito = is_mito))
sfe <- addGeometries(sfe, samples = sampleDir, sample_id = sampleNames, res = "fullres")
sfe <- addPerGeneQC(sfe, assay = "counts", version = NULL, mirror = NULL)

colData(sfe)
rowData(sfe)
colGeometries(sfe)


## ----03_QC_sfe2, message=FALSE, warning=FALSE-------------------------------------------------------
ggplot() + 
  geom_sf(aes(geometry = colGeometries(sfe)$spotHex$geometry, fill = colData(sfe)$annotation)) + 
  theme_void() + 
  theme(legend.position = "right") + 
  labs(fill = "Annotation")


## ----03_QC_sfe3-------------------------------------------------------------------------------------
# ----------------------------------------------- #
## Density and histogram of library sizes
ggplot(data = as.data.frame(colData(sfe)),
       aes(x = sum)) +
    geom_histogram(aes(y = after_stat(density)), 
                   colour = "black", 
                   fill = "grey",
                   bins = 50) +
    geom_density(alpha = 0.5,
                 adjust = 0.5,
                 fill = "#A0CBE8",
                 colour = "#4E79A7") +
    geom_vline(xintercept = c(1000, NA),
               colour = "red", 
               linetype = "dashed") + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    xlab("Library size") + 
    ylab("Density") + 
    theme_classic()
## Select library size threshold
qc_lib_size <- colData(sfe)$sum < 1000 #| colData(sfe)$sum > 45000
## Check how many spots are filtered out
table(qc_lib_size)
## Add threshold in colData
colData(sfe)$qc_lib_size <- qc_lib_size
## Check putative spatial patterns of removed spots
ggplot() + 
    geom_sf(data = colGeometry(sfe, "spotHex"),
            aes(geometry = geometry)) + 
    geom_sf(data = colGeometry(sfe, "spotHex"),
            aes(geometry = geometry, fill = colData(sfe)$qc_lib_size)) +
    scale_fill_manual(values = c("grey95", "red")) + 
    labs(fill = "Discarded") + 
    theme_bw()


## ----03_QC_sfe4-------------------------------------------------------------------------------------
# ----------------------------------------------- #
## Density and histogram of expressed genes
ggplot(data = as.data.frame(colData(sfe)),
       aes(x = detected)) +
    geom_histogram(aes(y = after_stat(density)), 
                   colour = "black", 
                   fill = "grey",
                   bins = 50) +
    geom_density(alpha = 0.5,
                 adjust = 0.5,
                 fill = "#A0CBE8",
                 colour = "#4E79A7") + 
    geom_vline(xintercept = c(550, NA),
               colour = "red", 
               linetype = "dashed") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    xlab("Genes expressed in each spot") + 
    ylab("Density") + 
    theme_classic()
## Select expressed genes threshold
qc_detected <- colData(sfe)$detected < 550 #| colData(sfe)$detected > 6000
## Check how many spots are filtered out
table(qc_detected)
## Add threshold in colData
colData(sfe)$qc_detected <- qc_detected
## Check for putative spatial pattern of removed spots
ggplot() + 
    geom_sf(data = colGeometry(sfe, "spotHex"),
            aes(geometry = geometry)) + 
    geom_sf(data = colGeometry(sfe, "spotHex"),
            aes(geometry = geometry, fill = colData(sfe)$qc_detected)) +
    scale_fill_manual(values = c("grey95", "red")) + 
    labs(fill = "Discarded") + 
    theme_bw()


## ----03_QC_sfe5-------------------------------------------------------------------------------------
# ----------------------------------------------- #
## Density and histogram of percentage of mitochondrial expression
ggplot(data = as.data.frame(colData(sfe)),
       aes(x = subsets_mito_percent)) +
    geom_histogram(aes(y = after_stat(density)), 
                   colour = "black", 
                   fill = "grey",
                   bins = 50) +
    geom_density(alpha = 0.5,
                 adjust = 0.5,
                 fill = "#A0CBE8",
                 colour = "#4E79A7") + 
    geom_vline(xintercept = c(22, NA),
               colour = "red", 
               linetype = "dashed") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    xlab("Percentage of mitochondrial expression") + 
    ylab("Density") + 
    theme_classic()
## Select mitochondrial percentage threshold
qc_mito <- colData(sfe)$subsets_mito_percent > 22
## Check how many spots are filtered out
table(qc_mito)
## Add threshold in colData
colData(sfe)$qc_mito <- qc_mito
## Check for putative spatial pattern of removed spots
ggplot() + 
    geom_sf(data = colGeometry(sfe, "spotHex"),
            aes(geometry = geometry)) + 
    geom_sf(data = colGeometry(sfe, "spotHex"),
            aes(geometry = geometry, fill = colData(sfe)$qc_mito)) +
    scale_fill_manual(values = c("grey95", "red")) + 
    labs(fill = "Discarded") + 
    theme_bw()


## ----03_QC_sfe6-------------------------------------------------------------------------------------
# ----------------------------------------------- #
## Check the number of discarded spots for each metric
apply(cbind(qc_lib_size, qc_detected, qc_mito), 2, sum)
## Combine together the set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito
table(discard)
## Store the set in the object
colData(sfe)$discard <- discard
## Check for putative spatial pattern of removed spots
ggplot() + 
    geom_sf(data = colGeometry(sfe, "spotHex"),
            aes(geometry = geometry)) + 
    geom_sf(data = colGeometry(sfe, "spotHex"),
            aes(geometry = geometry, fill = colData(sfe)$discard)) +
    scale_fill_manual(values = c("grey95", "red")) + 
    labs(fill = "Discarded") + 
    theme_bw()

# ----------------------------------------------- #
## remove combined set of low-quality spots
sfe <- sfe[, !colData(sfe)$discard]


## ----03_LogNorm_sfe---------------------------------------------------------------------------------
## Calculate library size factors
sfe <- computeLibraryFactors(sfe)
## Have a look at the size factors
summary(sizeFactors(sfe))
## Density and histogram of library sizes
ggplot(data = data.frame(sFact = sizeFactors(sfe)), 
       aes(x = sFact)) +
    geom_histogram(aes(y = after_stat(density)), 
                   colour = "black", 
                   fill = "grey",
                   bins = 40) +
    geom_density(alpha = 0.5,
                 adjust = 0.5,
                 fill = "#A0CBE8",
                 colour = "#4E79A7") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    xlab("Library size") + 
    ylab("Density") + 
    theme_classic()

# calculate logcounts using library size factors
sfe <- logNormCounts(sfe)


## ----03_GeneQC_sfe1---------------------------------------------------------------------------------
rowData(sfe)[["JBO019.s_logMean"]] <- rowSums(assay(sfe, "logcounts")) / rowData(sfe)[["JBO019.nLocations"]]


## ----03_GeneQC_sfe2---------------------------------------------------------------------------------
is_zero <- rowData(sfe)$total == 0
is_logLow <- rowData(sfe)[["JBO019.s_logMean"]] <= 1
discard_gs <- is_zero | is_mito | is_logLow
table(discard_gs)

rowData(sfe)$discard <- discard_gs

## FEATURE SELECTION
## remove mitochondrial and other genes
sfe <- sfe[!rowData(sfe)$discard, ]


## ----03_HVGs_sfe------------------------------------------------------------------------------------
## Fit mean-variance relationship
dec <- modelGeneVar(sfe,
                    assay.type = "logcounts")

## Visualize mean-variance relationship
fit <- metadata(dec)
fit_df <- data.frame(mean = fit$mean,
                     var = fit$var,
                     trend = fit$trend(fit$mean))

## Select top HVGs
top_hvgs <- getTopHVGs(dec, 
                       var.field = "bio", 
                       prop = 0.5,
                       var.threshold = 0,
                       fdr.threshold = 0.1)

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



## ----shortcut_code, eval=FALSE----------------------------------------------------------------------
## ## Import data
## sampleDir <- "./data/spaceranger_outs/Human_Liver_Steatotic/JBO019_Results"
## sampleNames <- "JBO019"
## sfe <- read10xVisiumSFE(samples = sampleDir,
##                         sample_id = sampleNames,
##                         type = "sparse",
##                         data = "filtered",
##                         images = "lowres",
##                         style = "W",
##                         zero.policy = TRUE)
## # ----------------------------------------------- #
## ground_truth <- read_table("./data/to_load/spotzonationGroup.txt")
## ## Add QC metrics
## is_mito <- grepl("(^MT-)|(^mt-)", rowData(sfe)$symbol)
## sfe <- addPerLocQC(sfe, gTruth = ground_truth, assay = "counts", 2, subsets = list(mito = is_mito))
## sfe <- addGeometries(sfe, samples = sampleDir, sample_id = sampleNames, res = "fullres")
## sfe <- addPerGeneQC(sfe, assay = "counts", version = NULL, mirror = NULL)
## # ----------------------------------------------- #
## ## SPOT SELECTION
## ## Select library size threshold
## qc_lib_size <- colData(sfe)$sum < 1000
## ## Add threshold in colData
## colData(sfe)$qc_lib_size <- qc_lib_size
## ## Select expressed genes threshold
## qc_detected <- colData(sfe)$detected < 550
## ## Add threshold in colData
## colData(sfe)$qc_detected <- qc_detected
## ## Select mitochondrial percentage threshold
## qc_mito <- colData(sfe)$subsets_mito_percent > 22
## ## Add threshold in colData
## colData(sfe)$qc_mito <- qc_mito
## ## Combine together the set of discarded spots
## discard <- qc_lib_size | qc_detected | qc_mito
## ## Store the set in the object
## colData(sfe)$discard <- discard
## ## Remove combined set of low-quality spots
## sfe <- sfe[, !colData(sfe)$discard]
## # ----------------------------------------------- #
## ## FEATURE SELECTION
## ## Calculate library size factors
## sfe <- computeLibraryFactors(sfe)
## ## Calculate logcounts using library size factors
## sfe <- logNormCounts(sfe)
## ## Calculate log-counts sample mean
## rowData(sfe)[["JBO019.s_logMean"]] <- rowSums(assay(sfe, "logcounts")) / rowData(sfe)[["JBO019.nLocations"]]
## ## Set and apply filters
## is_zero <- rowData(sfe)$total == 0
## is_logLow <- rowData(sfe)[["JBO019.s_logMean"]] <= 1
## discard_gs <- is_zero | is_mito | is_logLow
## rowData(sfe)$discard <- discard_gs
## ## Remove mitochondrial and other genes
## sfe <- sfe[!rowData(sfe)$discard, ]
## ## Fit mean-variance relationship
## dec <- modelGeneVar(sfe,
##                     assay.type = "logcounts")
## ## Select top HVGs
## top_hvgs <- getTopHVGs(dec,
##                        var.field = "bio",
##                        prop = 0.5,
##                        var.threshold = 0,
##                        fdr.threshold = 0.05)


## ----GWmodelFig1, echo=FALSE, out.width = "100%", fig.align="center", fig.cap="The math equations that define the kernels."----
knitr::include_graphics("images/gwmodel_kernel_math.png")


## ----GWmodelFig2, echo=FALSE, out.width = "100%", fig.align="center", fig.cap="Examples from using each kernel."----
knitr::include_graphics("images/gwmodel_kernel_graphs.png")


## ----03_spatial_weights_to_sfe----------------------------------------------------------------------
## add a neighbour graph using a weighted distance matrix
sfe <- addSpatialNeighGraphs(sfe, "JBO019", type = "knearneigh", style = "W", distMod = "raw", k = 6)

colGraphs(sfe)

## Calculate a simple distance matrix
sfe <- addDistMat(sfe, p = 2)



## ----03_visualise_neighbours------------------------------------------------------------------------
## Retrieve the tissue image
sfei <- getImg(sfe, image_id = "lowres")
## Extract the spot locations
spot_coords <- spatialCoords(sfe) %>% as.data.frame()

## Set limits
xlim <- c(min(spot_coords$pxl_col_in_fullres) - 100, 
          max(spot_coords$pxl_col_in_fullres) + 100)
ylim <- c(min(spot_coords$pxl_row_in_fullres) - 100, 
          max(spot_coords$pxl_row_in_fullres) + 100)
nbs <- colGraph(sfe)
ggplot() + 
    geom_spatraster_rgb(data = imgRaster(sfei)) + 
    geom_sf(data = as(nb2lines(nbs$neighbours, coords = spatialCoords(sfe)), "sf")) + 
    lims(x = xlim, y = ylim) +
    coord_sf() + 
    theme_void()


## ----all_in_one, message=FALSE, warning=FALSE, eval=FALSE-------------------------------------------
## ## Import data
## sampleDir <- "./data/spaceranger_outs/Human_Liver_Steatotic/JBO019_Results"
## sampleNames <- "JBO019"
## sfe <- read10xVisiumSFE(samples = sampleDir,
##                         sample_id = sampleNames,
##                         type = "sparse",
##                         data = "filtered",
##                         images = "lowres",
##                         style = "W",
##                         zero.policy = TRUE)
## # ----------------------------------------------- #
## ground_truth <- read_table("./data/to_load/spotzonationGroup.txt")
## ## Add QC metrics
## is_mito <- grepl("(^MT-)|(^mt-)", rowData(sfe)$symbol)
## sfe <- addPerLocQC(sfe, gTruth = ground_truth, assay = "counts", 2, subsets = list(mito = is_mito))
## sfe <- addGeometries(sfe, samples = sampleDir, sample_id = sampleNames, res = "fullres")
## sfe <- addPerGeneQC(sfe, assay = "counts", version = NULL, mirror = NULL)
## # ----------------------------------------------- #
## ## SPOT SELECTION
## ## Select library size threshold
## qc_lib_size <- colData(sfe)$sum < 1000
## ## Add threshold in colData
## colData(sfe)$qc_lib_size <- qc_lib_size
## ## Select expressed genes threshold
## qc_detected <- colData(sfe)$detected < 550
## ## Add threshold in colData
## colData(sfe)$qc_detected <- qc_detected
## ## Select mitochondrial percentage threshold
## qc_mito <- colData(sfe)$subsets_mito_percent > 22
## ## Add threshold in colData
## colData(sfe)$qc_mito <- qc_mito
## ## Combine together the set of discarded spots
## discard <- qc_lib_size | qc_detected | qc_mito
## ## Store the set in the object
## colData(sfe)$discard <- discard
## ## Remove combined set of low-quality spots
## sfe <- sfe[, !colData(sfe)$discard]
## # ----------------------------------------------- #
## ## FEATURE SELECTION
## ## Calculate library size factors
## sfe <- computeLibraryFactors(sfe)
## ## Calculate logcounts using library size factors
## sfe <- logNormCounts(sfe)
## ## Calculate log-counts sample mean
## rowData(sfe)[["JBO019.s_logMean"]] <- rowSums(assay(sfe, "logcounts")) / rowData(sfe)[["JBO019.nLocations"]]
## ## Set and apply filters
## is_zero <- rowData(sfe)$total == 0
## is_logLow <- rowData(sfe)[["JBO019.s_logMean"]] <= 1
## discard_gs <- is_zero | is_mito | is_logLow
## rowData(sfe)$discard <- discard_gs
## ## Remove mitochondrial and other genes
## sfe <- sfe[!rowData(sfe)$discard, ]
## ## Fit mean-variance relationship
## dec <- modelGeneVar(sfe,
##                     assay.type = "logcounts")
## ## Select top HVGs
## top_hvgs <- getTopHVGs(dec,
##                        var.field = "bio",
##                        prop = 0.5,
##                        var.threshold = 0,
##                        fdr.threshold = 0.05)
## # ----------------------------------------------- #
## ## ADD GEOGRAPHY
## ## Add a neighbour graph using a weighted distance matrix
## sfe <- addSpatialNeighGraphs(sfe, "JBO019", type = "knearneigh", style = "W", distMod = "raw", k = 6)
## ## Calculate a simple distance matrix
## sfe <- addDistMat(sfe, p = 2)

