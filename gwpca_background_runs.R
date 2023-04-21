library(spdep)
library(sf)
library(GWmodel)
library(ggplot2)
library(tidyverse)
load("/Users/b9047753/Downloads/ismbTutorial.RData")
countsNormCenter <- readRDS(file = "./data/countsNormScaled.rds")
rowDATA <- readRDS(file = "./data/rowDATA.rds")
colDATA <- readRDS(file = "./data/colDATA.rds")
## Prepare for Geographically Weighted PCA (GWPCA)
countsNormCenter <- countsNormCenter %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "rowname") %>%
  arrange("rowname") %>%
  column_to_rownames("rowname")

## Get the coordinates
coords <- colDATA[, c("Barcode", "pixel_x", "pixel_y")] %>%
  arrange(Barcode) %>%
  column_to_rownames(var = "Barcode")

## Get the data into a SpatialPointsDataFrame object
inputPCAgw <- SpatialPointsDataFrame(coords, 
                                     countsNormCenter, 
                                     match.ID = TRUE)
## Get the gene names that are going to be evaluated
vars <- colnames(inputPCAgw@data)

## Set the number of components to be retained
k <- 20

## Set the kernel to be used
kernel = "gaussian"

## Set a bandwidth -in pixels- for neighbourhood 
dist.Mat <- gw.dist(dp.locat = st_coordinates(colDATA$geom_cntd), p = 2)
# bw.choice <- readRDS(file = "./data/bw.choice.rds")
bw.choice <- 19

pcaGW <- gwpca(inputPCAgw,
               vars = vars,
               bw = bw.choice,
               k = k,
               dMat = dist.Mat,
               adaptive = TRUE,
               kernel = "gaussian")

saveRDS(pcaGW, file = "./data/pcaGW.rds")

