library(rstac)
library(terra)
library(sf)
library(jsonlite)
library(lubridate)
library(uuid)
library(dplyr)

BASE_DIR <- "/Volumes/T7SSD/"
DEM_PATH <- paste0(BASE_DIR, "/lynker-spatial/dem/vrt/USGS_1/ned_USGS_1.vrt")
CATALOG_DIR <- paste0(BASE_DIR, "/hydrofabric-stac")
