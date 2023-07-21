## ----01_load-data, message=FALSE--------------------------------------------------------------------
library(SpatialExperiment)
library(STexampleData)
library(ggplot2)
library(ggspavis)

# Load the object
spe <- Visium_humanDLPFC()


## ----SpExp-overview, echo=FALSE, out.width = "100%", fig.align="center", fig.cap="Overview of the `SpatialExperiment` object class structure."----
knitr::include_graphics("images/SpatialExperiment.png")


## ----01_data-inspect, message=FALSE-----------------------------------------------------------------
## Check the object's structure
spe

## Check number of features/genes (rows) and spots (columns)
dim(spe)

## Check names of 'assay' tables
assayNames(spe)


## ----01_counts-inspect, message=FALSE---------------------------------------------------------------
## Have a look at the counts table
assay(spe)[1:6,1:4]


## ----01_counts-chunks, message=FALSE----------------------------------------------------------------
assay(spe)[20:40, 2000:2010]

assay(spe)[33488:33508, 2000:2010]


## ----01_gene-metaData, message=FALSE----------------------------------------------------------------
## Have a look at the genes metadata
head(rowData(spe))


## ----01_coordinates-inspect, message=FALSE----------------------------------------------------------
## Check the spatial coordinates
head(spatialCoords(spe))

## spot-level metadata
head(colData(spe))



## ----01_image-inspect, message=FALSE----------------------------------------------------------------
## Have a look at the image metadata
imgData(spe)


## ----01_plot-image, message=FALSE, fig.height=8, fig.width=8----------------------------------------
## retrieve the image
spi <- getImg(spe)
## "plot" the image
plot(imgRaster(spi))


## ----01_plot-spots, message=FALSE, fig.height=8, fig.width=8----------------------------------------
## "Plot" the image
plot(imgRaster(spi))
## Extract the spot locations
spot_coords <- spatialCoords(spe) %>% as.data.frame
## Scale by low-res factor
lowres_scale <- imgData(spe)[imgData(spe)$image_id == 'lowres', 'scaleFactor']
spot_coords$x_axis <- spot_coords$pxl_col_in_fullres * lowres_scale
spot_coords$y_axis <- spot_coords$pxl_row_in_fullres * lowres_scale
## lowres image is 600x600 pixels
dim(imgRaster(spi))
## flip the Y axis
spot_coords$y_axis <- abs(spot_coords$y_axis - (ncol(imgRaster(spi)) + 1))
points(x=spot_coords$x_axis, y=spot_coords$y_axis)


## ----01_ggplot-spots, message=FALSE, fig.height=8, fig.width=8--------------------------------------
ggplot(mapping = aes(1:600, 1:600)) +
  annotation_raster(imgRaster(spi), xmin = 1, xmax = 600, ymin = 1, ymax = 600) +
  geom_point(data=spot_coords, aes(x=x_axis, y=y_axis), alpha=0.2) + xlim(1, 600) + ylim(1, 600) +
  coord_fixed() + 
  theme_void()


## ----01_ggplot_ontissue, message=FALSE, fig.height=8, fig.width=8-----------------------------------
## Add the annotation to the coordinate data frame
spot_coords$on_tissue <- as.logical(colData(spe)$in_tissue)

ggplot(mapping = aes(1:600, 1:600)) +
  annotation_raster(imgRaster(spi), xmin = 1, xmax = 600, ymin = 1, ymax = 600) +
  geom_point(data=spot_coords, aes(x=x_axis, y=y_axis, colour=on_tissue), alpha=0.2) + xlim(1, 600) + ylim(1, 600) +
  coord_fixed() + 
  theme_void()


## ----01_ggspavis-ontissue, message=FALSE, fig.height=8, fig.width=8---------------------------------
plotSpots(spe, in_tissue = NULL, annotate='in_tissue', size=0.5)

